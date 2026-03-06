## **************************************************************************
##
##    (c) 2025 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    ** Global variable definitions and residual code **
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
#' 
#' @useDynLib RPEM, .registration = TRUE
#' 
## Global variable definition:
globalVariables(c(
  ".",                                  ## Used by package magrittr.
  "candidate", "included", "pemModel"   ## Used by pemlm.
))
## 
PEMvar <- function(d, a = 0, psi = 1) {
  
  nd <- length(d)
  a <- rep(a, length.out = nd)
  psi <- rep(psi, length.out = nd)
  
  .C(
    "PEMvarC",
    as.double(d),
    as.integer(nd),
    as.double(a),
    as.double(psi),
    res = double(nd),
    PACKAGE = "RPEM"
  )$res
}
## 
generateVertexNames <- function(x, prefix)
  paste(
    prefix,
    sum(substr(rownames(x),1L,nchar(prefix)) == prefix) + 1L,
    sep=""
  )
## 
generateEdgeNames <- function(x, prefix)
  paste(
    prefix,
    sum(substr(edgenames(x),1L,nchar(prefix)) == prefix) + 1L,
    sep=""
  )
## 
joinEdge <- function(x, va_edge, va_loc, vb_edge, vb_loc, data,
                     vnames = list()) {
  
  ## Processing the edge origin:
  if(va_loc == 0) {
    
    ## Location: at the beginning of the origin edge.
    va <- edge(x)[[1L]][va_edge]
    
  } else if(va_loc == edge(x)$distance[va_edge]) {
    
    ## Location: at the end of the origin edge.
    va <- edge(x)[[2L]][va_edge]
    
  } else {
    
    ## Location: somewhere in the middle of the origin edge.
    
    ## Create an intermediary vertex to anchor the origin edge.
    ifelse(
      is.null(vnames$va),
      generateVertexNames(x, "IVA"),
      vnames$va
    ) -> vnm
    add.vertex(
      x,
      data.frame(
        row.names = vnm,
        species = FALSE
      )
    ) -> x
    va <- which(rownames(x) == vnm)
    
    ## Create the edge connecting the intermediary vertex to the vertex at the
    ## end of the origin edge:
    ifelse(
      is.null(vnames$ea),
      generateEdgeNames(x, "IEA"),
      vnames$ea
    ) -> enm
    add.edge(
      x,
      from = va,
      to = edge(x)[[2L]][va_edge],
      data = data.frame(
        row.names = enm,
        distance = edge(x)$distance[va_edge] - va_loc
      )
    ) -> x
    
    ## Reconnecting the upstream edge:
    edge(x)[[2L]][va_edge] <- va
    edge(x)$distance[va_edge]  <- va_loc
  }
  
  ## Processing the edge destination
  if(vb_loc == 0) {
    
    ## Location: at the beginning of the destination edge.
    vb <- edge(x)[[1L]][vb_edge]
    
  } else if(vb_loc == edge(x)$distance[vb_edge]) {
    
    ## Location: at the end of the destination edge.
    vb <- edge(x)[[2L]][vb_edge]
    
  } else {
    
    ## Location: somewhere in the middle of the destination edge.
    
    ## Create an intermediary vertex to anchor the destination edge.
    ifelse(
      is.null(vnames$vb),
      generateVertexNames(x, "IVB"),
      vnames$vb
    ) -> vnm
    add.vertex(
      x,
      data.frame(
        row.names = vnm,
        species = FALSE
      )
    ) -> x
    vb <- which(rownames(x) == vnm)
    
    ## Create the edge connecting the intermediary vertex to the vertex at the
    ## end of the destination edge:
    ifelse(
      is.null(vnames$eb),
      generateEdgeNames(x, "IEB"),
      vnames$eb
    ) -> enm
    add.edge(
      x,
      from = vb,
      to = edge(x)[[2L]][vb_edge],
      data = data.frame(
        row.names = enm,
        distance = edge(x)$distance[vb_edge] - vb_loc
      )
    ) -> x
    
    ## Reconnecting the upstream edge:
    edge(x)[[2L]][vb_edge] <- vb
    edge(x)$distance[vb_edge]  <- vb_loc
  }
  
  ## Adding the connecting edge:
  add.edge(
    x,
    from = va,
    to = vb,
    data = data
  )
}
## 
NULL
## 
