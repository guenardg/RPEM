## **************************************************************************
##
##    (c) 2025 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    ** Method Definition **
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
#' @name RPEM-generics
#' 
#' @aliases generics
#'
#' @title Generic Functions and Methods for RPEM
#'
#' @description Defines a set of S3 generic functions used internally and
#' externally by the \pkg{RPEM} package for graph manipulation, phylogenetic
#' modeling, and trait evolution analysis. These methods support the package's
#' object-oriented design, particularly around the \code{\link{graph-class}} and
#' \code{PEM-class} objects.
#' 
#' Users typically interact with these generics through their methods on
#' specific classes (e.g., \code{"graph"}, \code{"PEM"}), rather than calling
#' the generics directly.
#'
#' @param object A \code{\link{graph-class}} object or any object that can be
#' converted to one.
#' @param x A \code{\link{graph-class}} object, or a numeric vector/matrix of
#' auxiliary traits used in evolutionary model estimation. If used as a trait,
#' default is \code{NULL} (no auxiliary traits).
#' @param value For assignment methods: a two-column \code{data.frame} of
#' coordinates (e.g., x, y) to assign to a graph; or new value appropriate to
#' the method.
#' @param v A character string (vertex name) or integer (index) specifying the
#' origin (root) vertex.
#' @param link A data frame containing information to add cross-links (e.g.,
#' reticulation edges) to a graph.
#' @param target A numeric or character vector specifying one or more target
#' species (vertices) in graph \code{x}, typically for placement or prediction.
#' @param y A numeric vector or matrix of response trait(s) used to estimate
#' parameters in evolutionary models (e.g., in \code{evolution.model}).
#' @param ... Additional arguments passed to or from methods, depending on the
#' class and implementation.
#'
#' @details
#' This file declares S3 generics and associated method dispatchers for key
#' operations in \pkg{RPEM}, including:
#' \itemize{
#'   \item Graph structure access and modification (e.g., \code{nedge},
#'   \code{edge}, \code{edgenames}),
#'   \item Geometry and spatial layout handling (e.g., \code{geometry},
#'   \code{geometry<-}),
#'   \item Graph conversion and labeling (e.g., \code{as.graph},
#'   \code{labels<-}),
#'   \item Phylogenetic operations (e.g., \code{locate}, \code{set.origin},
#'   \code{crosslink}),
#'   \item Evolutionary model estimation (e.g., \code{evolution.model},
#'   \code{vcv}).
#' }
#'
#' @return The return value depends on the specific method and class. For
#' example:
#' \itemize{
#'   \item \code{nedge()}: integer scalar (number of edges),
#'   \item \code{edge()}: list or data structure representing edges,
#'   \item \code{as.graph()}: an object of class \code{"graph"},
#'   \item \code{evolution.model()}: a list with estimated parameters and fit
#'   statistics.
#' }
#'
#' @author \packageAuthor{RPEM}
#'
#' @keywords internal methods
#'
NULL
#'
#'
#' @describeIn RPEM-generics
#'
#' @title Number of Edges
#'
#' @description Access the number of edges in a graph-like object.
#'
#' @export
nedge <- function(x) UseMethod("nedge")
#'
#'
#' @describeIn RPEM-generics
#'
#' @title Edge Extraction
#'
#' @description Extract the edges of a graph-like object as a list or data
#' structure.
#'
#' @export
edge <- function(x) UseMethod("edge")
#'
#'
#' @describeIn RPEM-generics
#'
#' @title Edge Assignment
#'
#' @description Replace the edges of an object with a new set of edges.
#'
#' @export
`edge<-` <- function(x, value)  UseMethod("edge<-")
#'
#'
#' @describeIn RPEM-generics
#' 
#' Edge Names Extraction
#'
#' A method to extract edge names from an object.
#'
#' @export
edgenames <- function(x) UseMethod("edgenames")
#'
#'
#' @describeIn RPEM-generics
#'
#' @title Edge Names Assignment
#'
#' @description Assign names (labels) to the edges of a graph-like object.
#'
#' @export
`edgenames<-` <- function(x, value)  UseMethod("edgenames<-")
#'
#'
#' @describeIn RPEM-generics
#'
#' @title Geometry Extraction
#'
#' @description Obtain the geometric coordinates (e.g., x, y) associated with an
#' object.
#'
#' @export
geometry <- function(x)  UseMethod("geometry")
#'
#'
#' @describeIn RPEM-generics
#'
#' @title Geometry Assignment
#'
#' @description Assign geometric coordinates to an object, such as spatial
#' layout for plotting.
#' 
#' @export
`geometry<-` <- function(x, value)  UseMethod("geometry<-")
#'
#'
#' @describeIn RPEM-generics
#'
#' @title Labels Assignment
#'
#' @description Assign labels (e.g., species names) to the vertices of an
#' object.
#'
#' @export
`labels<-` <- function(object, value)  UseMethod("labels<-")
#'
#'
#' @describeIn RPEM-generics
#'
#' @title Graph-class Conversion
#'
#' @description Convert an object (e.g., phylo, network) into a
#' \code{\link{graph-class}} object.
#'
#' @export
as.graph <- function(x, ...)  UseMethod("as.graph")
#'
#'
#' @describeIn RPEM-generics
#'
#' @title Set Origin Vertex
#'
#' @description Designate a specific vertex as the origin (root) of the graph.
#'
#' @export
set.origin <- function(x, v, ...) UseMethod("set.origin")
#'
#'
#' @describeIn RPEM-generics
#'
#' @title Add Cross Links
#'
#' @description Add reticulation edges (e.g., hybridization, lateral transfer)
#' to a graph.
#'
#' @export
crosslink <- function(x, link, ...) UseMethod("crosslink")
#'
#'
#' @describeIn RPEM-generics
#'
#' @title Locate Target Species
#'
#' @description Determine the phylogenetic placement (LCA and distance) of
#' target species in a graph.
#'
#' @export
locate <- function(x, target, ...) UseMethod("locate")
#'
#'
#' @describeIn RPEM-generics
#'
#' @title Estimate Evolution Model Parameters
#'
#' @description Estimate the parameters of a trait evolution model (e.g.,
#' steepness, rate) using response and auxiliary traits on a PEM-class object.
#'
#' @export
evolution.model <- function(object, y, ..., x = NULL)
  UseMethod("evolution.model")
#'
#'
#' @describeIn RPEM-generics
#'
#' @title Estimate Variance-Covariance Matrix
#'
#' @description Compute the variance-covariance matrix of an object, typically
#' reflecting phylogenetic or evolutionary structure.
#'
#' @export
vcv <- function(x, ...)
  UseMethod("vcv")
#'
