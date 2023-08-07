# Functions to call variants from Sanger sequencing data

#' Title
#'
#' @param sampleKey_file File name to a table specifying the samples.
#' This table must have the following columns:
#' * SampleID: a unique ID for each bacterial sample
#' * Species
#' * Strain
#' * Gene
#' * Primer: the name of the primer used for the sequencing
#' * FileName: the name of the fasta file (without file extension)
#' * Success: whether or not the sequencing was successful (TRUE or FALSE)
#' * Reference: the name of the corresponding reference sequence
#' @param fasta_path The path where the fasta files are located.
#' @param fasta_extension File extension of the fasta files, usually either ".fa" (default) or ".fasta".
#' @param reference_files The path and file name(s) of the fasta file(s) with the reference sequences
#' @param hardTrim A vector specifying how many nucleotides on each side of a sequence should be hard-trimmed.
#'
#' @return A list containing two tables: "summary" list the number of variants of different types detected for each sample,
#' whereas "details" provides detailed information for each variant.
#' @export
#'
callSangerVariants_fasta <- function(sampleKey_file,
                                     fasta_path,
                                     fasta_extension = ".fa",
                                     reference_files,
                                     hardTrim = c(0, 0)) {

  # read in, check and sort sampleKey:
  sampleKey <- readr::read_csv(sampleKey_file, show_col_types = FALSE)
  if (!all(names(sampleKey)[1:8] == c("SampleID",
                                      "Species",
                                      "Strain",
                                      "Gene",
                                      "Primer",
                                      "FileName",
                                      "Success",
                                      "Reference"))) {
    stop("Columns in sampleKey file don't confirm to what is expected - refer to documentation.")
  }
  sampleKey <- sampleKey |>
    dplyr::arrange(Gene, SampleID)

  # read in and prepare fasta files for samples:
  fastaFileNames <- paste0(fasta_path,
                           sampleKey$FileName[sampleKey$Success],
                           fasta_extension)
  sampleSeqs <- Biostrings::readDNAStringSet(fastaFileNames) |>
    DECIPHER::RemoveGaps() |>
    XVector::subseq(hardTrim[1] + 1, -hardTrim[2] - 1) # trim base pairs on both sides

  # read in and prepare reference sequences:
  referenceSeqs <- Biostrings::readDNAStringSet(reference_files)

  # calculate coordinates - hard-coded for now to only work for E. coli MG1655 ref:
  genes <- unique(sampleKey$Gene)
  coordinates <- vector(mode = "list", length = length(genes))
  names(coordinates) <- genes
  for(i in 1:length(genes)) {
    refs <- referenceSeqs[unique(sampleKey$Reference[sampleKey$Gene == sampleKey$Gene[i]])]
    ref_Ecoli <- paste0(genes[i], "_Escherichia_coli_MG1655")
    coordinates[[i]] <- getAllCoordinates(refs, ref_Ecoli)
  }

  # prepare variants data frames:
  variantsSummary <- sampleKey |>
    dplyr::mutate(ContigLength = NA,
                  MapStart = NA,
                  MapEnd = NA,
                  nVariants = NA,
                  nMismatches = NA,
                  nInsertions = NA,
                  nDeletions = NA)

  variantsDetails <- data.frame(SampleID = rep(NA, 10000),  # make long data frame for more efficiency
                         Species = NA,
                         Strain = NA,
                         Gene = NA,
                         Primer = NA,
                         FileName = NA,
                         Success = NA,
                         Reference = NA,
                         ContigLength = NA,
                         MapStart = NA,
                         MapEnd = NA,
                         VariantType = NA,
                         Nt_pos = NA,
                         Nt_pos_Ecoli = NA,
                         Nt_original = NA,
                         Nt_mutation = NA,
                         Nt_OtoM = NA,
                         Nt_mut_name = NA,
                         Nt_mut_name_Ecoli = NA,
                         Codon_original = NA,
                         Codon_mutation = NA,
                         Codon_OtoM = NA,
                         AA_pos = NA,
                         AA_pos_Ecoli = NA,
                         AA_original = NA,
                         AA_mutation = NA,
                         AA_OtoM = NA,
                         AA_mut_name = NA,
                         AA_mut_name_Ecoli = NA)

  j <- 1 # tracker for row in variantsDetails table
  for(i in (1:nrow(sampleKey))[sampleKey$Success]) {
    svMisc::progress(i, progress.bar = TRUE, max.value = nrow(sampleKey))

    # mapping of contig to reference:
    contig <- sampleSeqs[[sampleKey$FileName[i]]]
    ref <- referenceSeqs[[sampleKey$Reference[i]]]
    alig <- Biostrings::pairwiseAlignment(contig, ref, type = "local")
    mismatches <- Biostrings::mismatchTable(alig)
    indels <- Biostrings::indel(alig)

    # fill in variantsSummary table:
    variantsSummary$ContigLength[i] <- length(contig)
    variantsSummary$MapStart[i] <- alig@subject@range@start
    variantsSummary$MapEnd[i] <- alig@subject@range@end
    variantsSummary$nVariants[i] <- nrow(mismatches) +
      length(indels@insertion[[1]]) + length(indels@deletion[[1]])
    variantsSummary$nMismatches[i] <- nrow(mismatches)
    variantsSummary$nInsertions[i] <- length(indels@insertion[[1]])
    variantsSummary$nDeletions[i] <- length(indels@deletion[[1]])

    # fill in variantsDetails table:
    if (variantsSummary$nMismatches[i] > 0) { # any mismatches?

      for(k in 1:nrow(mismatches)) {
        variantsDetails[j, 1:11] <- variantsSummary[i, 1:11] # same first columns
        variantsDetails$VariantType[j] <- "substitution"
        variantsDetails$Nt_pos[j] <- mismatches$SubjectStart[k]
        variantsDetails$Nt_original[j] <- mismatches$SubjectSubstring[k]
        variantsDetails$Nt_mutation[j] <- mismatches$PatternSubstring[k]
        variantsDetails <- fillMutationsTableRow_DNA(variantsDetails,
                                                     j,
                                                     sampleKey$Reference[i],
                                                     referenceSeqs,
                                                     coordinates[[sampleKey$Gene[i]]])
        variantsDetails <- fillMutationsTableRow_Protein(variantsDetails,
                                                         j,
                                                         sampleKey$Reference[i],
                                                         referenceSeqs,
                                                         coordinates[[sampleKey$Gene[i]]])
        j <- j + 1
      }
    }
  }
  variantsDetails <- variantsDetails[1:(j - 1),] |>
    select(-Success)
  return(list(summary = variantsSummary,
              details = variantsDetails))
}
