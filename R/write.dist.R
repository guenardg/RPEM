## **************************************************************************
##
##    (c) 2024-2025 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    ** Write distances to a text file in PHYLIP format **
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
#' @name write.dist
#' 
#' @title Write Distances to a PHYLIP Text File
#' 
#' @description A function to write a \code{\link{dist}}-class object into a
#' text file.
#' 
#' @param x A distance matrix in the form of a \code{\link{dist}}-class object.
#' @param file A character string describing the connection (see the description
#' of argument \code{description} in function \code{\link{file}} help page).
#' @param ... Further arguments to be internally passed to function
#' \code{\link{file}}.
#' 
#' @returns \code{NULL} (invisibly).
#' 
#' @author \packageAuthor{RPEM}
#' 
#' @export
write.dist <- function(x, file, ...) {
  
  if(!inherits(x,"dist"))
    stop("Argument 'x' must be of class 'dist'")
  
  n <- attr(x,"Size")
  
  unlist(
    lapply(
      strsplit(attr(x,"Labels")," "),
      paste,
      collapse = "_"
    )
  ) -> attr(x,"Labels")
  
  out <- file(description=file, open="wt", ...)
  
  cat(sprintf("%d\n", attr(x,"Size")), file=out)
  
  for(i in 1L:n) {
    cat(
      sprintf(
        "%s    %s\n",
        attr(x,"Labels")[i],
        paste(sprintf("%0.12f", x[i,]), collapse=" ")
      ),
      file = out
    )
  }
  
  close(out)
  
  invisible(NULL)
}
