# Functions to call variants from Sanger sequencing data

#' Calling variants in Sanger sequencing data
#'
#' @param sampleKey_file File name to a table specifying the samples.
#' This table must have the following columns:
#' * Sample_ID: a unique ID for each bacterial sample
#' * Species
#' * Strain
#' * Gene
#' * Primer: the name of the primer used for the sequencing
#' * File_name: the name of the fasta file (without file extension)
#' * Success: whether or not the sequencing was successful (TRUE or FALSE)
#' * Reference: the name of the corresponding reference sequence
#' @param fasta_path The path where the fasta files are located.
#' @param reference_files The path and file name(s) of the fasta file(s) with the reference sequences
#' @param fasta_extension File extension of the fasta files, usually either ".fa" (default) or ".fasta".
#' @param hard_trim A vector specifying how many nucleotides on each side of a sequence should be hard-trimmed.
#'
#' @return A list containing two tables: "summary" list the number of variants of different types detected for each sample,
#' whereas "details" provides detailed information for each variant.
#' @export
#'
callSangerVariants_fasta <- function(sampleKey_file,
                                     fasta_path,
                                     reference_files,
                                     fasta_extension = ".fa",
                                     hard_trim = c(0, 0)) {

  # read in, check and sort sampleKey:
  sampleKey <- readr::read_csv(sampleKey_file, show_col_types = FALSE)
  if (!("Success" %in% names(sampleKey))) { # add success column if not there
    sampleKey |> mutate(Success = TRUE)
  }
  if (!("Direction" %in% names(sampleKey))) { # add success column if not there
    sampleKey |> mutate(Direction = "F")
  }
  requiredSampleKeyColumns <- c("Sample_ID",
                                "Species",
                                "Strain",
                                "Gene",
                                "Primer",
                                "Direction",
                                "File_name",
                                "Success",
                                "Reference")
  if (!all(requiredSampleKeyColumns %in% names(sampleKey))) {
    stop("Columns in sampleKey file don't conform to what is expected - refer to documentation.")
  }
  if (!all(sampleKey$Direction %in% c("F", "R"))) {
    stop("Values in Direction column in sampleKey file need to be either F or R.")
  }

  sampleKey <- sampleKey |>
    dplyr::select(all_of(requiredSampleKeyColumns),
                  tidyselect::everything()) |>
    dplyr::arrange(Gene, Primer, Sample_ID)

  # read in and prepare fasta files for samples:
  fastaFileNames <- paste0(fasta_path,
                           sampleKey$File_name[sampleKey$Success],
                           fasta_extension)
  sampleSeqs <- Biostrings::readDNAStringSet(fastaFileNames) |>
    DECIPHER::RemoveGaps() |>
    XVector::subseq(hard_trim[1] + 1, -hard_trim[2] - 1) # trim base pairs on both sides

  # read in and prepare reference sequences:
  referenceSeqs <- Biostrings::readDNAStringSet(reference_files)

  # calculate coordinates - hard-coded for now to only work for E. coli MG1655 ref:
  genes <- unique(sampleKey$Gene)
  coordinates <- vector(mode = "list", length = length(genes))
  names(coordinates) <- genes
  for (i in 1:length(genes)) {
    ref_Ecoli <- paste0(genes[i], "_Escherichia_coli_MG1655")
    refs <- referenceSeqs[unique(c(
      ref_Ecoli,
      sampleKey$Reference[sampleKey$Gene == sampleKey$Gene[i]]))]
    coordinates[[i]] <- getAllCoordinates(refs, ref_Ecoli)
  }

  # prepare variants data frames:
  variantsSummary <- sampleKey |>
    dplyr::mutate(Contig_length = NA,
                  Map_start = NA,
                  Map_end = NA,
                  n_variants = NA,
                  n_mismatches = NA,
                  n_insertions = NA,
                  n_deletions = NA)

  variantsDetails <- data.frame(Contig_length = rep(NA, 10000), # make long for efficiency
                         Map_start = NA,
                         Map_end = NA,
                         Variant_type = NA,
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
                         AA_mut_name_Ecoli = NA,
                         RefSeq_ID = NA)
  variantsDetails[, names(sampleKey)] <- NA
  variantsDetails <- variantsDetails |>
    dplyr::select(all_of(names(sampleKey)),
                  tidyselect::everything())

  cat("\nProcessing samples.\n")
  j <- 1 # tracker for row in variantsDetails table
  for(i in (1:nrow(sampleKey))[sampleKey$Success]) {
    svMisc::progress(i, progress.bar = TRUE, max.value = nrow(sampleKey))

    # mapping of contig to reference:
    contig <- sampleSeqs[[sampleKey$File_name[i]]]
    if (sampleKey$Direction[i] == "R")
      contig <- Biostrings::reverseComplement(contig)
    ref <- referenceSeqs[[sampleKey$Reference[i]]]
    alig <- Biostrings::pairwiseAlignment(contig, ref, type = "local")
    mismatches <- Biostrings::mismatchTable(alig)
    indels <- Biostrings::indel(alig)

    # fill in variantsSummary table:
    variantsSummary$Contig_length[i] <- length(contig)
    variantsSummary$Map_start[i] <- alig@subject@range@start
    variantsSummary$Map_end[i] <- alig@subject@range@start + alig@subject@range@width - 1
    variantsSummary$n_variants[i] <- nrow(mismatches) +
      length(indels@insertion[[1]]) + length(indels@deletion[[1]])
    variantsSummary$n_mismatches[i] <- nrow(mismatches)
    variantsSummary$n_insertions[i] <- length(indels@insertion[[1]])
    variantsSummary$n_deletions[i] <- length(indels@deletion[[1]])

    # fill in variantsDetails table:
    if (variantsSummary$n_mismatches[i] > 0) { # any mismatches?
      for(k in 1:nrow(mismatches)) {
        variantsDetails[j, 1:11] <- variantsSummary[i, 1:11] # same first columns
        variantsDetails$Variant_type[j] <- "substitution"
        variantsDetails$Nt_pos[j] <- mismatches$SubjectStart[k]
        variantsDetails$Nt_original[j] <- mismatches$SubjectSubstring[k]
        variantsDetails$Nt_mutation[j] <- mismatches$PatternSubstring[k]
        variantsDetails$RefSeq_ID[j] <- sampleKey$Reference[i]

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
    if (variantsSummary$n_insertions[i] > 0) { # any insertions?
      warning(paste0("Insertions detected for sample ", sampleKey$Sample_ID[i], ", but processing these is currently not supported."))
      # for(k in 1:length(indels@insertion[[1]])) {
      #   variantsDetails[j, 1:11] <- variantsSummary[i, 1:11] # same first columns
      #   variantsDetails$Variant_type[j] <- "insertion"
      #   variantsDetails$Nt_pos[j] <- indels@insertion[[1]]@start[k]
      #   variantsDetails$Nt_original[j] <- mismatches$SubjectSubstring[k]
      #   variantsDetails$Nt_mutation[j] <- mismatches$PatternSubstring[k]
      #   variantsDetails <- fillMutationsTableRow_DNA(variantsDetails,
      #                                                j,
      #                                                sampleKey$Reference[i],
      #                                                referenceSeqs,
      #                                                coordinates[[sampleKey$Gene[i]]])
      #   j <- j + 1
      # }
    }
  }
  variantsDetails <- variantsDetails[1:(j - 1),] |>
    dplyr::select(-RefSeq_ID)
  variantsSummary <- variantsDetails |>
    dplyr::group_by(dplyr::across(tidyselect::all_of(names(sampleKey)))) |>
    dplyr::summarise(Nt_variants = paste(Nt_mut_name, collapse = ", ")) |>
    dplyr::right_join(variantsSummary) |>
    dplyr::ungroup() |>
    dplyr::relocate(Nt_variants, .after = last_col()) |>
    as.data.frame()
  return(list(summary = variantsSummary,
              details = variantsDetails))
}
