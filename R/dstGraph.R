## **************************************************************************
##
##    (c) 2025 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    ** Distance-based directed graph **
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
#' @name dstGraph
#'
#' @title Distance-Based Directed Graph
#'
#' @description Constructs a directed graph from a dissimilarity matrix by
#' connecting vertices within a user-defined threshold distance, starting from a
#' specified origin vertex. Optionally allows threshold stretching to connect
#' otherwise isolated vertices.
#'
#' @param d A dissimilarity matrix, such as those produced by
#' \code{\link[stats]{dist}} or \code{\link[ape]{dist.dna}}.
#' @param th Numeric. Threshold value for dissimilarity: vertices are connected
#' if their pairwise distance is \eqn{\leq} \code{th}.
#' @param origin Integer. Index of the origin (root) vertex from which the graph
#' grows via directed edges.
#' @param stretch Numeric (optional). If provided, this value is used to stretch
#' the threshold for vertices that are initially unreachable, allowing the
#' shortest edge to them to be included if its distance is \eqn{\leq}
#' \code{stretch}.
#'
#' @return A \code{\link{graph-class}} object representing the directed
#' phylogenetic-like graph. The graph includes:
#' \itemize{
#'   \item Directed edges from ancestors to descendants,
#'   \item Edge attribute \code{distance},
#'   \item Vertex attribute \code{order}, indicating the sequence in which
#'   vertices were added,
#'   \item Vertex attribute \code{species} (logical), set to \code{TRUE} for all
#'   vertices by default.
#' }
#' Unconnected vertices (excluding the origin) are reported via a message.
#'
#' @details
#' The algorithm builds the graph incrementally, starting from the specified
#' \code{origin} vertex. At each step, it connects all unvisited vertices within
#' the dissimilarity threshold \code{th} to any already-connected vertex. This
#' process repeats until no more vertices can be reached.
#'
#' If \code{stretch} is provided, isolated vertices are re-evaluated using the
#' stretched threshold to ensure full connectivity, which can be useful in
#' fragmented or noisy data.
#'
#' This method is particularly useful for:
#' \itemize{
#'   \item Constructing phylogenetic hypotheses from genetic distances,
#'   \item Building reticulate networks,
#'   \item Creating graphs for use in Phylogenetic Eigenvector Maps (PEM).
#' }
#'
#' @author \packageAuthor{RPEM}
#'
#' @references
#' Guénard, G., Legendre, P., and Peres-Neto, P. (2013).
#' Phylogenetic eigenvector maps: a framework to model and predict species
#' traits.
#' \emph{Methods in Ecology and Evolution}, 4: 1120–1131.
#'
#' Makarenkov, V., Legendre, L., & Desdevise, Y. (2004).
#' Modelling phylogenetic relationships using reticulated networks.
#' \emph{Zoologica Scripta}, 33: 89–96.
#'
#' @import magrittr
#'
#' @examples
#' # Set seed for reproducibility
#' set.seed(7653401)
#'
#' # Simulate spatial coordinates (e.g., representing genetic variation)
#' N <- 100
#' coords <- cbind(x = runif(N, -1, 1), y = runif(N, -1, 1))
#' rownames(coords) <- sprintf("N%d", 1:N)
#'
#' # Compute Euclidean dissimilarity matrix
#' dst <- dist(coords)
#'
#' # Build directed graph with threshold = 0.25, origin at vertex 15
#' gr <- dstGraph(d = dst, th = 0.25, origin = 15)
#'
#' # Plot the graph
#' plot(coords, type = "n", asp = 1, main = "Directed Graph from Dissimilarity")
#' col <- head(rainbow(max(gr$order)), max(gr$order))
#' arrows(
#'   x0 = coords[edge(gr)[[1]], 1],
#'   y0 = coords[edge(gr)[[1]], 2],
#'   x1 = coords[edge(gr)[[2]], 1],
#'   y1 = coords[edge(gr)[[2]], 2],
#'   length = 0.05,
#'   col = col[gr$order[edge(gr)[[2]]]]
#' )
#' points(coords, pch = 21, bg = "black", cex = 0.25)
#'
#' # Try higher threshold
#' gr <- dstGraph(d = dst, th = 0.28, origin = 15)
#'
#' # Now allow stretching to connect isolated vertices
#' gr <- dstGraph(d = dst, th = 0.28, origin = 15, stretch = 0.5)
#'
#' # Verify all vertices are connected
#' if (!any((gr$order == 0)[-15])) cat("All vertices are now connected.\n")
#'
#' # Compute and visualize influence matrix
#' infl_matrix <- InflMat(gr)
#' image(t(infl_matrix[nrow(infl_matrix):1, ]), col = gray(c(1, 0)), asp = 1,
#'       main = "Influence Matrix (Transposed)")
#'
#' @export
dstGraph <- function(d, th, origin, stretch) {
  
  n <- attr(d, "Size")
  d <- as.numeric(d)
  
  if(!missing(stretch)) {
    dc <- d
    d[d > th] <- NA
    
    for(i in 1L:n) {
      idx <- dst_idx(n,i)
      if(all(is.na(d[idx]))) {
        wh <- which(dc[idx] <= stretch)
        d[idx[wh]] <- dc[idx[wh]]
      }
    }
    
  } else {
    d[d > th] <- NA
  }
  
  graph(
    data.frame(
      order = integer(n),
      row.names = attr(d,"Label")
    )
  ) -> gr
  
  ord <- 1L
  oo <- origin
  
  while(length(oo)) {
    for(i in oo) {
      idx <- dst_idx(n,i)
      dd <- d[idx]
      wh <- which(!is.na(dd))
      dd <- dd[wh]
      vv <- (1L:n)[-i][wh]
      gr$order[vv] <- ord
      nvv <- length(vv)
      if(nvv)
        add.edge(
          gr,
          from = rep(i,nvv),
          to = vv,
          data = data.frame(
            distance = dd,
            row.names = sprintf("E%d", nedge(gr) + (1L:nvv))
          )
        ) -> gr
      d[idx[wh]] <- NA
    }
    oo <- which(gr$order == ord)
    ord <- ord + 1L
  }
  
  unconnected <- gr$order == 0
  unconnected[origin] <- FALSE
  
  attr(gr,"processOrder") <- getProcessOrder(gr)
  
  if(any(unconnected)) {
    message(
      sum(unconnected),
      " vertices are not reachable from the specified origin:\n",
      sprintf("%s\t",rownames(gr)[unconnected])
    )
  }
  
  ## By default, all the vertices have observations.
  gr$species <- TRUE
  
  gr
}
