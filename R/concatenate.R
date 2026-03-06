## **************************************************************************
##
##    (c) 2025 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    ** Sequence concatenation **
##
##    This file is part of RPEM
##
##    RPEM is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    RPEM is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with RPEM. If not, see <https://www.gnu.org/licenses/>.
##
##    R source code file
##
## **************************************************************************
##
#' @name concatenate
#' 
#' @title Sequence Concatenation
#' 
#' @description Converts rows of a \code{\link{raw}} matrix—used to store DNA or
#' RNA sequences—into a character vector of strings. Optionally removes
#' specified characters (e.g., gaps).
#' 
#' @param x A \code{\link{raw}} matrix where each row represents a sequence of
#' nucleotides (or other molecular states) encoded as raw bytes.
#' @param discard A character vector specifying which symbols to remove from the
#' sequences before conversion (e.g., \code{"-"} for gaps). Default is
#' \code{NULL}, meaning no characters are removed.
#' 
#' @return A character vector of length \code{nrow(x)}, where each element is
#' the string representation of the corresponding row in \code{x}, optionally
#' filtered by \code{discard}.
#' 
#' @details This function applies \code{\link{rawToChar}} to each row of the
#' input matrix \code{x}, effectively converting raw-encoded sequences into
#' readable character strings. If \code{discard} is provided, all occurrences of
#' the specified characters (after conversion to raw) are removed prior to
#' string construction.
#'
#' This is particularly useful after simulating DNA sequences using
#' \code{molEvolSim()}, where sequences are stored efficiently as \code{raw}
#' matrices, and need to be exported or visualized as standard strings.
#'
#' @author \packageAuthor{RPEM}
#' 
#' @import magrittr
#' 
#' @examples
#' # Define example sequences
#' sequences <- c(
#'   Sequence_1 = "ATCG-TTTCG--CCCCA--TTA--TTAA-GTAA-GTAATCTTTCA",
#'   Sequence_2 = "TTGGCTTCC--TC--CAT-TTC--TTCA-GT-ACG-ATTCTTTTA",
#'   Sequence_3 = "TTCG-TACC-T-T---A-ATAA--T-AA-GTCTTGTAATCGTTCA"
#' )
#'
#' # Convert to a raw matrix (each row is a sequence)
#' raw_matrix <- sapply(sequences, charToRaw, USE.NAMES = TRUE) %>%
#'   t()
#'
#' # Display raw matrix
#' print(raw_matrix)
#'
#' # Convert to character strings
#' concatenate(raw_matrix)
#'
#' # Convert to character strings with gaps removed
#' concatenate(raw_matrix, discard = "-")
#'
#' # Clean up
#' rm(sequences, raw_matrix)
#'
#' @export
concatenate <- function(x, discard = NULL) {
  apply(
    x,
    1L,
    function(x, discard, linebreak, sep) {
      if(!is.null(discard)) x <- x[!(x %in% charToRaw(discard))]
      rawToChar(x)
    },
    discard = discard
  )
}
