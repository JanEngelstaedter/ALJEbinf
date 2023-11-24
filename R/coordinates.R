# Functions relating to coordinates, i.e. tables relating positions in one focal sequence
# to positions in a homologous reference sequence (usually in E. coli).


#' Calculate a table of coordinates
#'
#' This function takes two DNA sequences, aligns the translated amino acid sequences,
#' and calculates a table of corresponding positions in this alignment ("coordinates").
#'
#' @param dnaFocal Focal DNA sequence.
#' @param dnaRef Reference DNA sequence, e.g. Escherichia coli.
#' @param aligOutput If FALSE (default), only the actual coordinates are returned.
#' If TRUE then the alignment is also returned.
#'
#' @return By default, a table with two columns: the position in the focal sequence (posFocal),
#' and the position in the reference sequence (posRef).
#' If aligOutput==TRUE, a list containing the coordinates and the alignment.
#'
#' @export
#'
getCoordinates <- function(dnaFocal,
                           dnaRef,
                           aligOutput = FALSE) {
  AAFocal <- Biostrings::translate(dnaFocal,
                                   genetic.code = Biostrings::getGeneticCode("Bacterial",
                                                                             full.search = TRUE))
  AARef <- Biostrings::translate(dnaRef,
                                 genetic.code = Biostrings::getGeneticCode("Bacterial",
                                                                           full.search = TRUE))
  alig <- Biostrings::pairwiseAlignment(AAFocal, AARef)

  # extract the aligned sequences from the alignment
  aligFocal <- alig |>
    Biostrings::alignedPattern() |>
    toString() |>
    stringr::str_replace_all("[a-zA-Z*]", "NNN") |>
    stringr::str_replace_all("[-]", "---") |>
    Biostrings::DNAString()

  aligRef <- alig |>
    Biostrings::alignedSubject() |>
    toString() |>
    stringr::str_replace_all("[a-zA-Z*]", "NNN") |>
    stringr::str_replace_all("[-]", "---") |>
    Biostrings::DNAString()

  posFocal <- 0L
  posRef <- 0L
  coordinates <- data.frame(posFocal = rep(NA_integer_, length(aligFocal)),
                            posRef   = NA_integer_)
  for (i in 1:length(aligFocal)) {
    letterFocal <- Biostrings::extractAt(aligFocal, IRanges::IRanges(i)) |>
      as.character()
    letterRef <- Biostrings::extractAt(aligRef, IRanges::IRanges(i)) |>
      as.character()
    if (letterFocal != "-") {
      posFocal <- posFocal + 1L
      coordinates$posFocal[i] <- posFocal
    }
    if (letterRef != "-") {
      posRef <- posRef + 1L
      coordinates$posRef[i] <- posRef
    }
  }
  if (aligOutput) {
    return(list(coordinates = coordinates,
                alignment = alig))
  } else {
    return(coordinates)
  }
}

#' Obtain tables of coordinates.
#'
#' Takes a list of homologous sequences of which one is designated a reference,
#' aligns them all to the reference (at amino-acid level), and calculates for each sequence
#' a table of corresponding coordinates.
#'
#' @param seqs A DNAStringSet object containing a set of homologous sequences.
#' @param refSeq_ID The name of the sequence designated as the reference.
#'
#' @return A list of tables of coordinates. Each of these tables has two columns:
#' the position in the focal sequence (posFocal), and the position in the reference sequence (posRef).
#' @export
#'
getAllCoordinates <- function(seqs, refSeq_ID) {
  coordinates <- vector(mode = "list", length = length(seqs))
  names(coordinates) <- names(seqs)
  for(i in 1:length(seqs)) {
    cat(paste0("Determining coordinates for ", names(seqs)[i], ".\n"))
    coordinates[[i]] <- getCoordinates(seqs[[i]], seqs[[refSeq_ID]])
  }
  return(coordinates)
}

#' Translates between positions in two sequences.
#'
#' Using coordinates calculated using the `getCoordinates` or `getAllCoordinates` function,
#' this function translates between a position in a focal sequence to a position in
#' the reference sequence, or vice versa.
#'
#' @param pos Position within the sequence.
#' @param coordinates Table of coordinates.
#' @param direction Either "FocalToRef" (the default) or "RefToFocal".
#' @param AAinput When `FALSE` (the default), the `pos` argument will be assumed to be a
#' nucleotide position. Otherwise, `pos` will be treated as an amino acid position.
#' @param AAoutput When `FALSE` (the default), the returned translated position will be a nucleotide
#' position, otherwise an amino acid position.
#'
#' @return The translated position.
#' @export
#'
translateCoordinate <- function(pos,
                                coordinates,
                                direction = "FocalToRef",
                                AAinput = FALSE,
                                AAoutput = FALSE) {
  if (AAinput) {
    pos <- pos * 3L - 2L
  }
  if(length(pos) > 1L) {
    pos <- pos[1]
    warning("More than one position in pos vector, only the first one will be used.")
  }
  if (direction == "FocalToRef") {
    coordinates$posRef[duplicated(coordinates$posRef) |
                           (coordinates$posRef == 0L)] <- NA
    trans <- coordinates |>
      dplyr::filter(posFocal == pos) |>
      dplyr::pull(posRef) |>
      dplyr::first()
  } else if (direction == "RefToFocal") {
    coordinates$posFocal[duplicated(coordinates$posFocal) |
                           (coordinates$posFocal == 0L)] <- NA
    trans <- coordinates |>
      dplyr::filter(posRef == pos) |>
      dplyr::pull(posFocal) |>
      dplyr::first()
  } else {
    stop("Argument direction must be either FocalToRef or RefToFocal")
  }
  if (AAoutput == TRUE)
    trans <- ((trans - 1L) %/% 3L) + 1L
  return(trans)
}
