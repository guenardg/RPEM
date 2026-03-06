## **************************************************************************
##
##    (c) 2025 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    ** Pairwise distance matrix extraction operator **
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
#' @name dist-operator
#'
#' @title Subsetting Operator for Distance Objects
#'
#' @description Provides subsetting functionality for objects of class
#' \code{"dist"} (as created by \code{\link[stats]{dist}}), enabling extraction
#' of subsets of distances by index, name, or logical condition.
#'
#' @param x A distance object of class \code{"dist"}, typically returned by
#' \code{\link[stats]{dist}}.
#' @param i A numeric, logical, or character vector specifying the rows (or
#' elements) to extract. If \code{j} is missing, the result is a new
#' \code{"dist"} object containing only the pairwise distances among the
#' selected objects.
#' @param j A numeric, logical, or character vector specifying the columns when
#' extracting a rectangular subset (i.e., distances from set \code{i} to set
#' \code{j}). If omitted, the operation is symmetric and returns a \code{"dist"}
#' object.
#' @param drop A logical value indicating whether dimensions should be dropped
#' if the result is a single row or column. Default is \code{TRUE}.
#' @param ... Further arguments passed to or from other methods (currently
#' unused).
#'
#' @return
#' \itemize{
#'   \item If only \code{i} is provided: a \code{\link[stats]{dist}} object
#'   containing distances among the selected elements.
#'   \item If both \code{i} and \code{j} are provided: a numeric vector (if
#'   \code{drop = TRUE} and one index has length 1) or a numeric matrix
#'   (otherwise), with distances from elements in \code{i} to those in \code{j}.
#'   \item If no indices are provided (\code{x[]}): the full distance vector,
#'   coerced to numeric if \code{drop = TRUE}.
#' }
#'
#' @details
#' The base \code{stats} package lacks a proper subsetting method for
#' \code{"dist"} objects. This method fills that gap, allowing intuitive
#' indexing of distance matrices in their compact form.
#'
#' The function supports:
#' \itemize{
#'   \item Indexing by position (\code{1:5}),
#'   \item Logical indexing (\code{c(TRUE, TRUE, FALSE)}),
#'   \item Label-based indexing (\code{c("A", "B", "C")}) if labels are present,
#'   \item Negative indexing (e.g., \code{-c(2,3)}),
#'   \item Rectangular extraction (asymmetric subsets via \code{i} and
#'   \code{j}).
#' }
#'
#' When \code{i} and \code{j} are both given, the result is no longer symmetric
#' and is returned as a matrix or vector. This is useful for extracting specific
#' blocks of distances (e.g., between two clades).
#'
#' @author \packageAuthor{RPEM}
#'
#' @examples
#' # Create a distance object
#' x <- dist(c(A = 1, B = 2, C = 5, D = 9, E = 4, F = 2, G = 6, H = 4, I = 0, J = 3))
#' x
#'
#' # Subset by index
#' x[c(1, 2, 3, 8, 10, 6)]
#' x[1:4]
#' x[-(2:4)]
#'
#' # Subset by logical vector
#' x[c(TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE)]
#'
#' # Subset by labels
#' x[c("A", "B", "F", "J")]
#'
#' # Rectangular extraction: distances from first 3 to next 7
#' x[1:3, 4:10]
#'
#' # Extract single row as vector
#' x[1, 2:10]
#'
#' # Keep matrix form
#' x[1:3, 2:4, drop = FALSE]
#'
#' # Full extraction
#' x[]
#'
#' @method [ dist
#'
#' @export
`[.dist` <- function(x, i, j, ..., drop=TRUE) {
  
  cl <- match.call()
  sx <- attr(x, "Size")
  na <- nargs() - !missing(drop)
  
  ## Processing indices i, if any:
  if(missing(i)) {
    i <- NULL
  } else
    if(is.character(i)) {
      if(is.null(attr(x, "Labels")))
        stop("labels are provided as references, but the object has no labels")
      if(!all(i %in% attr(x, "Labels")))
        stop(
          "\nunknown references: ",
          paste(i[!(i %in% attr(x, "Labels"))], collapse=", ")
        )
      i <- match(i, attr(x, "Labels"))
    } else if(is.logical(i)) {
      if(length(i) != sx)
        stop(length(i), " logical references for ", sx, " objects")
      i <- which(i)
    } else if(is.numeric(i)) {
      storage.mode(i) <- "integer"
      if(any(i < 0) && any(i > 0))
        stop("negative indices can only be mixed with 0s")
      if(all(i > 0)) {
        if(any(i > sx))
          stop("\nunknown indices: ", paste(i[i > sx], collapse=", "))
      } else {
        if(any(-i > sx))
          stop("\nunknown indices: ", paste(-i[-i > sx], collapse=", "))
        i <- (1L:sx)[i]
      }
    } else
      stop("\nunsuitable storage mode for indices: ", storage.mode(i))
  
  ## Processing indices j, if any:
  if(missing(j)) {
    j <- NULL
  } else {
    if(is.character(j)) {
      if(is.null(attr(x, "Labels")))
        stop("labels are provided as references, but the object has no labels")
      if(!all(j %in% attr(x, "Labels")))
        stop(
          "\nunknown references: ",
          paste(j[!(j %in% attr(x, "Labels"))], collapse=", ")
        )
      j <- match(j, attr(x, "Labels"))
    } else if(is.logical(j)) {
      if(length(j) != sx)
        stop(length(j), " logical references for ", sx, " objects")
      j <- which(j)
    } else if(is.numeric(j)) {
      storage.mode(j) <- "integer"
      if(any(j < 0) && any(j > 0))
        stop("negative indices can only be mixed with 0s")
      if(all(j > 0)) {
        if(any(j > sx))
          stop("\nunknown indices: ", paste(j[j > sx], collapse=", "))
      } else {
        if(any(-j > sx))
          stop("\nunknown indices: ", paste(-j[-j > sx], collapse=", "))
        j <- (1L:sx)[j]
      }
    } else
      stop("\nunsuitable storage mode for indices: ", storage.mode(j))
  }
  
  if(is.null(i) && is.null(j))
    return(if(drop) as.numeric(x) else x)
  
  if((na < 3L) && !is.null(i)) {
    
    sy <- length(i)
    storage.mode(x) <- "double"
    y <- double(sy*(sy - 1)/2)
    
    return(
      structure(
        .C("extractDistC", sx, x, sy, i, y, NAOK=TRUE, PACKAGE="RPEM")[[5L]],
        Size = sy,
        Labels = if(!is.null(attr(x, "Labels"))) attr(x, "Labels")[i] else NULL,
        Diag = attr(x, "Diag"),
        Upper = attr(x, "Upper"),
        method = attr(x, "method"),
        call = c(attr(x, "call"), cl),
        class = "dist"
      )
    )
  }
  
  if(!(na < 3L)) {
    
    if(is.null(i)) i <- 1L:sx
    if(is.null(j)) j <- 1L:sx
    
    ni <- length(i)
    nj <- length(j)
    
    .C("extractDistC2", sx, x, ni, i, nj, j, double(ni*nj), NAOK=TRUE,
       PACKAGE="RPEM")[[7L]] -> out
    
    if(drop && ((ni == 1L) || (nj == 1L))) {
      return(out)
    } else
      return(
        matrix(
          data = out,
          nrow = ni,
          ncol = nj,
          byrow = FALSE,
          dimnames = list(attr(x,"Labels")[i], attr(x,"Labels")[j])
        )
      )
  }
  
  invisible(NULL)
}

