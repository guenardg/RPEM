## **************************************************************************
##
##    (c) 2025 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    ** Package RPEM description **
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
#' \packageTitle{RPEM}
#' 
#' @name RPEM-package
#' 
#' @description \packageDescription{RPEM}
#' 
#' @details Phylogenetic eignevector maps (PEM) is a method for using phylogeny
#' to model features of organism, most notably quantitative traits. It consists
#' in calculating sets of explanatory variables (eigenvectors) that are meant to
#' represent different patterns in trait values that are likely to have been
#' inducted by evolution. These patterns are used to model the data, using a
#' linear model for instance.
#' 
#' If one in interested in a \sQuote{target} species (i.e. a species for which
#' the trait value is unknown), and provided that we know the phylogenetic
#' relationships between that species and those of the model, the method allows
#' us to obtain the scores of that new species on the phylogenetic
#' eigenfunctions underlying a PEM. These scores are used to make empirical
#' predictions of trait values for the target species on the basis of those
#' observed for the species used in the model.
#' 
#' Function \code{\link{PEM}}, and its methods \code{\link{update.PEM}} and
#' \code{\link{evolution.model.PEM}} allow one to build a PEM, update it, and
#' estimate its evolution model (selection steepness and local evolution rates),
#' respectively.
#' 
#' For making predictions, \code{\link{graph-class}} method
#' \code{\link{locate.graph}} calculate the residual graph and graph locations
#' for a set of target vertices, whereas \code{\link{PEM-class}} method
#' \code{\link{predict.PEM}} enables one to obtain the values of the PEM
#' eigenfunctions (ie., scores) at any arbitrary graph location. RPEM also
#' provides a linear modeling interface called \code{\link{pemlm}}. It enables
#' one to fit a base linear model, with or without auxiliary traits, and then
#' run a forward selection of the most relevant phylogenetic eigenfunctions on
#' the basis of the (corrected) Akaike Information Criterion. A
#' \code{\link{predict}} method is also available for \code{\link{pemlm-class}}
#' objects (\code{\link{predict.pemlm}}), which enables one to make predictions
#' directly for a set of target graph locations using an auxiliary trait data
#' frame and a \code{\link{PEM-class}} object.
#' 
#' The package provides low-level utility functions for performing operations on
#' graphs (see \link{graph-functions}), calculate influence matrix
#' (\code{\link{InflMat}}), and so on.
#' 
#' A phylogenetic modelling tutorial using \code{RPEM} is available as a
#' package vignette. See example below.
#' 
#' The DESCRIPTION file:
#' \packageDESCRIPTION{RPEM}
#' \packageIndices{RPEM}
#' 
#' @author \packageAuthor{RPEM}
#' Maintainer: \packageMaintainer{RPEM}
#' 
#' @references
#' Guénard, G., Legendre, P., and Peres-Neto, P. 2013. Phylogenetic eigenvector
#' maps: a framework to model and predict species traits. Methods in Ecology 
#' and Evolution 4: 1120-1131
#' 
#' @seealso
#' Makarenkov, V., Legendre, P. & Desdevise, Y. 2004. Modelling phylogenetic
#' relationships using reticulated networks. Zoologica Scripta 33: 89-96
#' 
#' Blanchet, F. G., Legendre, P. & Borcard, D. 2008. Modelling directional
#' spatial processes in ecological data. Ecological Modelling 215: 325-336
#' 
#' @examples
#' ## To view RPEM tutorial
#' vignette("RPEM_tutorial", package="RPEM")
#' vignette("Simulating_phylogenetic_network", package="RPEM")
#'
NULL
##
