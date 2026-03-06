## **************************************************************************
##
##    (c) 2025 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    ** Class and method for directed graphs **
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
#' @name graph-class
#'
#' @title Class and Methods for Directed Phylogenetic Graphs
#'
#' @description
#' The \code{graph-class} is the core data structure in \pkg{RPEM} for
#' representing phylogenetic and reticulated evolutionary graphs. It supports
#' both tree-like and network-based topologies, enabling modeling of vertical
#' descent and lateral evolutionary processes (e.g., hybridization, horizontal
#' gene transfer).
#'
#' This class extends \code{data.frame} and includes an \code{edge} attribute (a
#' data frame) that defines the directed connections between vertices. It is
#' designed for integration with phylogenetic modeling, trait evolution
#' simulation, and PEM computation.
#'
#' @docType class
#'
#' @keywords classes
#'
#' @section Structure:
#' A \code{graph-class} object is a \code{data.frame} where:
#' \itemize{
#'   \item Each row represents a vertex (tip or internal node),
#'   \item Row names are vertex labels,
#'   \item Columns represent vertex properties (e.g., \code{species},
#'   \code{type}),
#'   \item An attribute \code{edge} (a data frame) stores edge information:
#'   \code{from}, \code{to} (vertex indices), and edge properties (e.g.,
#'   \code{distance}).
#' }
#'
#' @section Required Properties:
#' For full functionality (e.g., in \code{PEM} or \code{locate}):
#' \itemize{
#'   \item A logical vertex property \code{species} indicating which vertices
#'   have trait data,
#'   \item A numeric edge property \code{distance} representing evolutionary
#'   time or divergence.
#' }
#'
#' @section Methods:
#' The class provides methods for:
#' \itemize{
#'   \item Conversion: \code{as.graph()}, \code{as.phylo()},
#'   \item Manipulation: \code{add.edge()}, \code{crosslink()},
#'   \item Querying: \code{nedge()}, \code{edge()}, \code{edgenames()},
#'   \item Geometry: \code{geometry()}, \code{geometry<-},
#'   \item Plotting: \code{plot.graph()},
#'   \item Topology: \code{set.origin()}, \code{locate()}.
#' }
#'
#' @section Cross-Linking:
#' The \code{crosslink()} method adds reticulation edges (e.g., hybridization
#' events) between locations on existing edges. The \code{link} data frame must
#' include:
#' \itemize{
#'   \item \code{$from}, \code{$loc_from}: origin edge and position,
#'   \item \code{$to}, \code{$loc_to}: destination edge and position,
#'   \item \code{$length}: evolutionary distance of the new edge,
#'   \item Optional: \code{$va}, \code{$ea}, \code{$vb}, \code{$eb} for custom
#'   vertex/edge names.
#' }
#' If a location is internal to an edge, a new bisecting vertex is inserted.
#'
#' @section Plotting:
#' The \code{plot} method supports two modes:
#' \itemize{
#'   \item \strong{With geometry}: Uses spatial coordinates (via \code{sf}) to
#'   plot the graph.
#'   \item \strong{Without geometry}: Coerces the graph to a tree using
#'   \code{as.phylo()}, plots the backbone, and overlays reticulation edges in
#'   red.
#' }
#' Vertex markers can be scaled and colored by trait values (\code{y}) using
#' \code{bg} and \code{cex.min}/\code{cex.max}.
#'
#' @format
#' A \code{graph-class} object inherits from \code{data.frame} and includes:
#' \itemize{
#'   \item Vertex data: rows = vertices, columns = vertex properties,
#'   \item Edge data: stored in \code{attr(x, "edge")} as a data frame with
#'   \code{from} and \code{to} (vertex indices) and optional edge properties.
#' }
#'
#' @param x,object A \code{graph-class} object.
#' @param value A vector or \code{data.frame} of values to assign (e.g.,
#' geometry, edge data).
#' @param y Numeric vector of trait values to visualize at vertices (optional).
#' @param ylim Range for scaling \code{y} values in plotting (default: auto).
#' @param use.geometry Logical; use embedded geometry if available (default:
#' \code{TRUE}).
#' @param pch Plot character for vertices (default: 21, filled circle).
#' @param bg Background color(s) for vertices; single color or vector for color
#' scale.
#' @param cex.min Minimum marker size (default: 2).
#' @param cex.max Maximum marker size (default: \code{cex.min}).
#' @param cex.lab Label size for vertex names (default: \code{cex.min/3}).
#' @param axes Whether to show plot axes (default: \code{FALSE}).
#' @param xlab,ylab Axis labels (default: empty).
#' @param edge.color List of two colors: backbone edges and reticulation edges
#' (default: \code{list("black", "red")}).
#' @param length Arrowhead size in inches (default: 0.05).
#' @param code Arrow style (see \code{\link{arrows}}; default: 2).
#' @param show.vertex.labels Whether to display vertex labels (default:
#' \code{TRUE}).
#' @param direction Direction of tree layout when no geometry (default:
#' \code{"rightwards"}).
#' @param name,v,link,target,... See individual method documentation.
#' 
#' @author \packageAuthor{RPEM}
#' 
#' @references
#' Guénard, G., Legendre, P., and Peres-Neto, P. (2013).
#' Phylogenetic eigenvector maps: a framework to model and predict species traits.
#' \emph{Methods in Ecology and Evolution}, 4: 1120–1131.
#'
#' @seealso
#' \code{\link{graph-functions}}, \code{\link[ape]{plot.phylo}},
#' \code{\link[sf]{st_as_sf}}.
#' 
#' @importFrom sf st_as_sf
#' @importFrom utils head tail
#' @importFrom ape as.phylo
#' @importFrom graphics arrows segments points text
#' 
#' @examples
#' # Create a simple graph
#' library(magrittr)
#' data.frame(
#'   species = rep(TRUE, 13),
#'   type = c(2,2,3,1,2,2,2,2,2,2,3,3,3),
#'   x = c(1,3,4,0,1.67,4,1,1.33,2.33,3.33,4.33,4,5),
#'   y = c(1,1,1,0,0.5,0,-1,0,0,-0.5,-1,-0.5,-0.5),
#'   row.names = sprintf("V%d", 1:13)
#' ) %>%
#'   st_as_sf(coords = c("x", "y"), crs = NA) %>%
#'   graph %>%
#'   add.edge(
#'     from = c(1,2,1,5,4,4,5,9,4,8,9,4,7,7,6,6,9,10,10),
#'     to = c(2,3,5,2,1,5,9,2,8,9,6,7,8,9,3,13,10,12,11),
#'     data = data.frame(
#'       distance = c(4.2,4.7,3.9,3.0,3.6,2.7,4.4,3.4,3.6,3.3,4.8,3.2,3.5,4.4,2.5,3.4,4.3,3.1,2.2),
#'       row.names = sprintf("E%d", 1:19)
#'     )
#'   ) -> x
#'
#' # Inspect
#' x
#' nedge(x)
#' plot(x)
#' plot(x, use.geometry = FALSE)
#'
#' # Convert to phylo
#' phy <- as.phylo(x)
#' plot(phy, show.node.label = TRUE)
#'
NULL
#' 
#' 
#' @describeIn graph-class
#'
#' @title Print Graph Object
#'
#' @description Prints a summary of the graph, including number of vertices and
#' edges, labels, and available properties.
#' 
#' @method print graph
#' 
#' @export
print.graph <- function (x, ...) {
  
  nv <- nrow(x)
  ne <- nedge(x)
  edge <- edge(x)
  
  cat("\nA graph with", nv, "vertices and", ne, "edges\n")
  cat("----------------------------------\n")
  
  vl <- rownames(x)
  
  if(length(vl)) {
    
    cat("Vertex labels: ")
    
    if(length(vl) > 10L) {
      cat(paste(c(head(vl,5L), paste("... +", length(vl) - 10L, "more ..."),
                  tail(vl,3L)), collapse=", "))
    } else cat(paste(vl, collapse=", "))
    
    cat("\n")
  }
  
  el <- rownames(edge)
  
  if(length(el)) {
    
    cat("Edge labels: ")
    
    if(length(el) > 10L) {
      cat(paste(c(head(el, 5L), paste("... +", length(el) - 10L, "more ..."),
                  tail(el, 3L)), collapse=", "))
    } else cat(paste(el, collapse=", "))
    
    cat("\n")
  }
  
  if(ncol(x) > 0L) {
    cat("Vertex information: ", paste(colnames(x), collapse = ", "), "\n")
  } else cat("No available vertex information\n")
  
  if(ncol(edge) > 2L) {
    cat("Available edge information: ",
        paste(colnames(edge)[-(1L:2L)], collapse = ", "), "\n")
  } else
    cat("No edge information available\n")
  
  cat("\n")
  
  invisible(NULL)
}
#' 
#' @describeIn graph-class
#'
#' @title Number of Edges
#'
#' @description Returns the number of edges in a graph.
#'
#' @return An integer.
#'
#' @method nedge graph
#'
#' @export
nedge.graph <- function(x) nrow(attr(x, "edge"))
#'
#' @describeIn graph-class
#'
#' @title Extract Edges
#'
#' @description Extracts the edge data frame from a graph.
#'
#' @return A data frame with columns \code{from}, \code{to}, and optional edge
#' properties.
#'
#' @method edge graph
#'
#' @export
edge.graph <- function(x) attr(x, "edge")
#'
#' @describeIn graph-class
#'
#' @title Assign Edges
#'
#' @description Replaces the edges of a graph with a new edge data frame.
#'
#' @return The updated graph.
#'
#' @method edge<- graph
#'
#' @export
`edge<-.graph` <- function(x, value) {
  
  if(is.null(value)) {
    
    attr(x, "edge") <- data.frame(from=integer(0L), to=integer(0L))
    
  } else if(is.data.frame(value) && !is.null(value[[1L]]) &&
            !is.null(value[[2L]])) {
    
      attr(x, "edge") <- value
      
    } else
      stop("The 'value' must be a data frame with at least two columns")
  
  x
}
#'
#' @describeIn graph-class
#'
#' @title Extract Edge Names
#'
#' @description Returns the names of the edges in the graph.
#'
#' @return A character vector of edge names.
#'
#' @method edgenames graph
#'
#' @export
edgenames.graph <- function(x) rownames(attr(x, "edge"))
#'
#' @describeIn graph-class
#'
#' @title Assign Edge Names
#'
#' @description Sets the names of the edges in the graph.
#'
#' @return The updated graph.
#'
#' @method edgenames<- graph
#'
#' @export
`edgenames<-.graph` <- function(x, value) {
  rownames(attr(x, "edge")) <- value
  x
}
#'
#' @describeIn graph-class
#'
#' @title Extract Geometry
#'
#' @description Retrieves the spatial geometry (e.g., coordinates) associated
#' with the graph.
#'
#' @return A geometry object (e.g., \code{sf} points), or \code{NULL} if none.
#'
#' @method geometry graph
#'
#' @export
geometry.graph <- function(x)
  if(is.null(attr(x,"sf_column"))) NULL else x[[attr(x,"sf_column")]]
#'
#' @describeIn graph-class
#'
#' @title Assign Geometry
#'
#' @description Assigns spatial coordinates to the graph vertices.
#'
#' @return The updated graph with geometry.
#'
#' @method geometry<- graph
#'
#' @export
`geometry<-.graph` <- function(x, value) {
  
  if(is.null(value)) {
    
    if(inherits(x, "sf")) {
      x[[attr(x, "sf_column")]] <- NULL
      attr(x, "sf_column") <- NULL
      class(x) <- class(x)[class(x) != "sf"]
    }
    
  } else {
    
    if(nrow(x) != length(value))
      stop("The number of coordinates provided (", length(value), ") does not ",
           "equal the number of vertices (", nrow(x),")")
    
    if(inherits(x, "sf")) {
    
      x[[attr(x, "sf_column")]] <- value
    
    } else {
    
      attr(x, "sf_column") <- "geometry"
      x[[attr(x, "sf_column")]] <- value
      cl <- class(x)
      class(x) <- c(cl[1L:which(cl == "graph")], "sf",
                    cl[(which(cl == "graph") + 1L):length(cl)])
      
    }
    
  }
  
  x
}
#'
#' @describeIn graph-class
#'
#' @title Insert or Replace Vertex Property
#'
#' @description Adds or replaces a vertex property (column) in the graph.
#'
#' @return The updated graph with in inserted property
#'
#' @method $<- graph
#'
#' @export
`$<-.graph` <- function(x, name, value)
  `[[<-.data.frame`(x, name, value=value)
#'
#' @describeIn graph-class
#'
#' @title Convert Graph to Phylo
#'
#' @description Coerces a \code{graph-class} object into a \code{phylo} object.
#' Reticulation edges are removed and stored in the \code{discarded} attribute.
#'
#' @method as.phylo graph
#'
#' @return An object of class \code{"phylo"} with attributes \code{discarded}
#' and \code{subst}.
#'
#' @export
as.phylo.graph <- function(x, ...) {
  
  nv <- nrow(x)
  
  po <- attr(x, "processOrder")
  if(is.null(po))
    po <- getProcessOrder(x)
  
  x <- reorderGraph(x, po)
  
  o <- getOrigin(x)
  
  if(length(o) > 1L)
    stop("The graph has to have a single origin to be transformable ",
         "into a tree.")
  
  ## The graph is coerced into as tree, if necessary.
  if(!isTree(x)) {
    
    discarded <- NULL
    
    d <- attr(x,"dist")
    if(is.null(d))
      d <- graphDist(x)
    
    do <- numeric(nv)
    do[-o] <- as.numeric(d)[dst_idx(nv, o)]
    
    for(i in 1L:nv) {
      wh <- which(edge(x)[[2L]] == i)
      if(length(wh) > 1L) {
        
        wh <- wh[-which.min(do[edge(x)[[1L]][wh]])]
        
        as.data.frame(
          lapply(edge(x), function(x) x[wh]),
          row.names = edgenames(x)[wh],
        ) -> df
        
        discarded <- rbind(discarded, df)
        
        x <- rm.edge(x, wh)
      }
    }
    
  } else
    discarded <- NULL
  
  sw <- integer(nv)
  tip <- getTerminal(x)
  
  sw[tip] <- 1L:length(tip)
  sw[-tip] <- (length(tip) + 1L):nv
  
  if(!is.null(discarded)) {
    discarded[[1L]] <- sw[discarded[[1L]]]
    discarded[[2L]] <- sw[discarded[[2L]]]
  }
  
  subst <- integer(nv)
  subst[sw] <- po
  
  structure(
    list(
      edge = cbind(sw[edge(x)[[1L]]], sw[edge(x)[[2L]]]),
      tip.label = rownames(x)[tip],
      node.label = rownames(x)[-tip],
      Nnode = nv - length(tip),
      edge.length = edge(x)$distance
    ),
    discarded = discarded,
    subst = subst,
    class = "phylo"
  )
}
#'
#' @describeIn graph-class
#'
#' @title Convert Phylo to Graph
#'
#' @description Converts a \code{phylo} object into a \code{graph-class} object.
#'
#' @method as.graph phylo
#'
#' @return A \code{graph-class} object.
#'
#' @export
as.graph.phylo <- function(x, ...) {
  
  args <- list(...)
  
  includeRoot <- if(is.null(args$includeRoot)) FALSE else args$includeRoot
  rootKnown <- if(is.null(args$rootKnown)) TRUE else args$rootKnown
  
  root <- includeRoot && !is.null(x$root.edge)
  
  if(is.null(x$node.label))
    x$node.label <- paste("N", 1L:x$Nnode, sep="")
  
  graph(
    data.frame(
      species = c(rep(TRUE, length(x$tip.label)),
                  rep(FALSE, x$Nnode), if(root) rootKnown),
      row.names = c(x$tip.label, x$node.label, if(root) "ROOT")
    )
  ) -> out
  
  add.edge(
    out,
    from = c(x$edge[,1L], if(root) nrow(out)),
    to = c(x$edge[,2L], if(root) x$edge[1L,1L]),
    data = data.frame(
      distance = c(x$edge.length, if(root) x$root.edge),
      row.names = sprintf("E%d", 1L:(nrow(x$edge) + as.integer(root)))
    )
  )
}
#'
#' @describeIn graph-class
#'
#' @title Convert Data Frame to Graph
#'
#' @description Converts a data frame into a \code{graph-class} object with no
#' edges.
#'
#' @return A \code{graph-class} object.
#'
#' @method as.graph data.frame
#'
#' @export
as.graph.data.frame <- function(x, ...) graph(data = x)
#'
#' @describeIn graph-class
#'
#' @title Summarize Graph
#'
#' @description Returns a structured summary of the graph.
#'
#' @return A list with \code{sumv} (vertex data) and \code{sume} (edge data).
#'
#' @method summary graph
#'
#' @export
summary.graph <- function(object, ...)
  structure(
    list(
      sumv = as.data.frame(object),
      sume = data.frame(
        row.names = edgenames(object),
        from = rownames(object)[edge(object)$from],
        to = rownames(object)[edge(object)$to],
        edge(object)[,-(1L:2L),drop=FALSE]
      )
    ),
    class = "summary.graph"
  )
#'
#' @describeIn graph-class
#'
#' Not to be documented.
#' @export
print.summary.graph <- function(x, ...) {
  cat("Graph summary\n------------------\n")
  cat(sprintf("Vertices (%d):\n",nrow(x$sumv)))
  cat(paste(rownames(x$sumv), collapse=", "), "\n")
  cat("Vertex descriptor(s):\n")
  print(summary(x$sumv))
  cat(sprintf("Edges (%d):\n",nrow(x$sume)))
  print(x$sume)
  invisible(NULL)
}
#'
#' @describeIn graph-class
#'
#' @title Plot Graph
#'
#' @description Plots the graph using either embedded geometry or a tree layout
#' with reticulation edges.
#'
#' @method plot graph
#'
#' @export
plot.graph <- function(x, y, ylim, use.geometry = TRUE, pch = 21L, bg = "white",
                       cex.min = 2, cex.max = cex.min, cex.lab = cex.min/3,
                       axes = FALSE, xlab = "", ylab = "",
                       edge.color = list("black","red"), length = 0.05,
                       code = 2L, show.vertex.labels = TRUE,
                       direction = c("rightwards", "leftwards", "upwards",
                                     "downwards"),
                       ...) {
  
  .midArrows <- function(x0, y0, x1, y1, length, code, ...) {
    arrows(x0=x0, y0=y0, x1=0.5*(x0 + x1), y1=0.5*(y0 + y1), length=length,
           code=code, ...)
    segments(x0=0.5*(x0 + x1), y0=0.5*(y0 + y1), x1=x1, y1=y1, ...)
  }
  
  if(use.geometry && !is.null(geometry(x))) {
    
    xy <- st_coordinates(geometry(x))
    
    plot(xy, asp=1, type="n", xlab=xlab, ylab=ylab, axes=FALSE, ...)
    
    .midArrows(
      x0 = xy[edge(x)[[1L]],1L],
      y0 = xy[edge(x)[[1L]],2L],
      x1 = xy[edge(x)[[2L]],1L],
      y1 = xy[edge(x)[[2L]],2L],
      col = edge.color[[1L]],
      length = length,
      code = code,
      ...
    )
    
    if(missing(y)) {
      
      points(x=xy[,1L], y=xy[,2L], pch=pch, cex=cex.min, bg=bg[1L], xpd=TRUE,
             ...)
      
    } else {
      
      if(missing(ylim))
        ylim <- range(y)
      
      cc <- (y - ylim[1L])/(ylim[2L] - ylim[1L])
      cc[cc < 0] <- 0
      cc[cc > 1] <- 1
      
      points(
        x = xy[,1L],
        y = xy[,2L],
        pch = pch,
        cex = cex.min + (cex.max - cex.min)*cc,
        bg = bg[1L + floor(cc*(length(bg) - 1L))],
        xpd = TRUE,
        ...
      )
      
    }
    
    if(show.vertex.labels)
      text(
        x = xy[,1L],
        y = xy[,2L],
        labels = rownames(x),
        cex = cex.lab,
        ...
      )
    
    invisible(x)
    
  } else {
    
    direction <- match.arg(direction)
    
    tre <- as.phylo(x)
    
    xy <- getPhyloXY(tre, direction, FALSE)
    attr(tre,"xy") <- xy
    
    plot(xy, type="n", xlab=xlab, ylab=ylab, axes=FALSE, ...)
    
    .midArrows(
      x0 = xy[tre$edge[,1L],1L],
      y0 = xy[tre$edge[,1L],2L],
      x1 = xy[tre$edge[,2L],1L],
      y1 = xy[tre$edge[,2L],2L],
      col = edge.color[[1L]],
      length = length,
      code = code,
      ...
    )
    
    if(!is.null(attr(tre,"discarded")))
      .midArrows(
        x0 = xy[attr(tre,"discarded")[[1L]],1L],
        y0 = xy[attr(tre,"discarded")[[1L]],2L],
        x1 = xy[attr(tre,"discarded")[[2L]],1L],
        y1 = xy[attr(tre,"discarded")[[2L]],2L],
        col = edge.color[[2L]],
        length = length,
        code = code,
        ...
      )
    
    if(missing(y)) {
      
      points(x=xy[,1L], y=xy[,2L], pch=pch, cex=cex.min, bg=bg[1L], xpd=TRUE,
             ...)
      
    } else {
      
      if(missing(ylim))
        ylim <- range(y)
      
      cc <- (y[attr(tre,"subst")] - ylim[1L])/(ylim[2L] - ylim[1L])
      cc[cc < 0] <- 0
      cc[cc > 1] <- 1
      
      points(
        x = xy[,1L],
        y = xy[,2L],
        pch = pch,
        cex = cex.min + (cex.max - cex.min)*cc,
        bg = bg[1L + floor(cc*(length(bg) - 1L))],
        xpd = TRUE,
        ...
      )
      
    }
    
    if(show.vertex.labels)
      text(
        x = xy[,1L],
        y = xy[,2L],
        labels = c(tre$tip.label,tre$node.label),
        cex = cex.lab,
        ...
      )
    
    invisible(tre)
    
  }
}
#' 
#' @describeIn graph-class
#'
#' @title Set Graph Origin
#'
#' @description Sets the origin (root) vertex of the graph by reversing edges as
#' needed.
#'
#' @return The updated graph with specified origin.
#'
#' @method set.origin graph
#'
#' @export
set.origin.graph <- function(x, v, ...) {
  
  if(length(v) > 1L) {
    warning("Only the first value will be taken as the origin.")
    v <- v[1L]
  }
  
  if(is.character(v)) {
    if(!(v %in% rownames(x)))
      stop("Unknown vertex: ", v)
    v <- which(v == rownames(x))
  } else
    if((v < 1L) || (v > nrow(x)))
      stop("Unknown vertex: ", v)
  
  erm <- logical(nedge(x))
  
  repeat{
    nextv <- c()
    for(i in v) {
      s <- !erm & (edge(x)$to == i)
      if(any(s)) {
        for(j in which(s)) {
          vf <- edge(x)$from[j]
          edge(x)$from[j] <- edge(x)$to[j]
          edge(x)$to[j] <- vf
          nextv <- c(nextv, vf)
        }
        erm[s] <- TRUE
      }
    }
    if(is.null(nextv)) break else v <- nextv
  }
  
  x
}
#'
#' @describeIn graph-class
#'
#' @title Add Cross Links (Reticulation)
#'
#' @description Adds reticulation edges (e.g., hybridization) between internal
#' points on edges.
#'
#' @return The updated graph with new reticulation edges.
#'
#' @method crosslink graph
#' 
#' @export
crosslink.graph <- function(x, link, ...) {
  
  ## Check whether 'x' is a graph:
  if(!inherits(x,"graph"))
    stop("argument 'x' must be a graph-class object")
  
  ## Graph 'x' must have an edge property called 'distance'
  if(is.null(edge(x)$distance))
    stop(
      "the graph provided as argument 'x' must have an edge property called ",
      "'distance'"
    )
  
  ## Graph 'x' must have a vertex property called 'species'
  if(is.null(x$species))
    stop(
      "the graph provided as argument 'x' must have a vertex property called ",
      "'species'"
    )
  
  c(
    is.data.frame(link),
    is.numeric(link$length),
    is.character(link$from),
    is.numeric(link$loc_from),
    is.character(link$to),
    is.numeric(link$loc_to),
    ifelse(is.null(link$va), TRUE, is.character(link$va)),
    ifelse(is.null(link$ea), TRUE, is.character(link$ea)),
    ifelse(is.null(link$vb), TRUE, is.character(link$vb)),
    ifelse(is.null(link$eb), TRUE, is.character(link$eb))
  ) -> itst
  
  if(!all(itst)) {
    cat("Storage modes:\n")
    print(
      data.frame(
        row.names = c("link","$length","$from","$loc_from","$to","$loc_to",
                      "$va","$ea","$vb","$eb"),
        Observed = c(
          storage.mode(link),
          storage.mode(link$length),
          storage.mode(link$from),
          storage.mode(link$loc_from),
          storage.mode(link$to),
          storage.mode(link$loc_to),
          ifelse(is.null(link$va), "null", storage.mode(link$va)),
          ifelse(is.null(link$ea), "null", storage.mode(link$ea)),
          ifelse(is.null(link$vb), "null", storage.mode(link$vb)),
          ifelse(is.null(link$eb), "null", storage.mode(link$eb))
        ),
        Expected = c("list","integer or double","character","integer or double",
                     "character","integer or double",
                     rep("character or null",4L)),
        OK = ifelse(itst, "yes", "no")
      )
    )
    stop("storage mode discrepancies observed")
  }
  
  
  if(!(is.character(link$from) && is.character(link$to)))
    stop(
      "edge coordinates provided through arguments 'link$from' and ",
      "'link$to' must both be character strings"
    )
  
  opt <- as.data.frame(lapply(list(...), rep, length.out=nrow(link)))
  
  ## i=1L
  for(i in 1L:nrow(link)) {
    
    if(rownames(link)[i] %in% edgenames(x))
      stop("attempting to create an edge with a name that already exists at ",
           "step ", i, " (", rownames(link)[i],")")
    
    if(!(link$from[i] %in% edgenames(x)))
      stop("unknown origin edge at step ", i, " (", link$from[i], ")")
    
    va_edge <- which(link$from[i] == edgenames(x))
    va_loc <- link$loc_from[i]
    
    if(va_loc > edge(x)$distance[va_edge])
      stop("location beyond the origin edge end at step ", i, " (",va_loc," > ",
           edge(x)$distance[va_edge],")")
    
    if(!(link$to[i] %in% edgenames(x)))
      stop("unknown destination edge at step ", i, " (", link$to[i], ")")
    
    vb_edge <- which(link$to[i] == edgenames(x))
    vb_loc <- link$loc_to[i]
    
    if(vb_loc > edge(x)$distance[vb_edge])
      stop("location beyond the destination edge end at step ", i, " (", vb_loc,
           " > ", edge(x)$distance[vb_edge], ")")
    
    data.frame(
      row.names = rownames(link)[i],
      distance = link$length[i],
      opt[i,,drop=FALSE]
    ) -> data
    
    list(
      va = link[["va"]][i],
      ea = link[["ea"]][i],
      vb = link[["vb"]][i],
      eb = link[["eb"]][i]
    ) -> vnames
    
    joinEdge(
      x = x,
      va_edge = va_edge,
      va_loc = va_loc,
      vb_edge = vb_edge,
      vb_loc = vb_loc,
      data = data,
      vnames = vnames
    ) -> x
  }
  
  x
}
#'
#' @describeIn graph-class
#'
#' @title Locate Target Species
#'
#' @description Removes target species and returns their placement (LCA and
#' distance) in the residual graph.
#' 
#' @return A list with \code{$x} (residual graph) and \code{$location}
#' (placement data).
#' 
#' @method locate graph
#' 
#' @export
locate.graph <- function(x, target, ...) {
  
  if(is.character(target)) {
    
    idx <- match(target, rownames(x))
    if(any(is.na(idx)))
      stop("Unknown target(s): ", paste(target[is.na(idx)], collapse=", "))
    
    data.frame(
      row.names = target,
      ref = idx,
      dist = as.double(rep(NA, length(target))),
      ladist = double(length(target))
    ) -> ttab
    
    target <- idx
    
  } else {
    
    outrange <- which((target < 0) | (target > nrow(x)))
    if(length(outrange))
      stop("Unknown target(s): ", paste(target[outrange], collapse=", "))
    
    data.frame(
      row.names = names(target),
      ref = target,
      dist = as.double(rep(NA, length(target))),
      ladist = double(length(target))
    ) -> ttab
  }
  
  if(is.null(x$species))
    stop("'x' has no vertex property called 'species'")
  
  if(!all(x$species[target]))
    stop("Non-species Target(s): ",
         paste(target[!x$species[target]], collapse=", "))
  
  ## Sanity check: is there any terminal vertex not marked as species?
  if(!all(x$species[getTerminal(x)]))
    stop("Sanity check failed: the graph has one or more terminal vertices ",
         "not marked as species.\nFunction purge.terminal() can be used to ",
         "discard them.")
  
  ## Sanity check: is there any median vertex not marked as a species?
  if(!all(x$species[getMedian(x)]))
    stop("Sanity check failed: the graph has one or more median vertices not ",
         "marked as species.\nFunction purge.median() can be used to discard ",
         "them automatically.")
  
  edge <- edge(x)
  
  if(is.null(edge$distance))
    stop("'x' has no edge property called 'distance'")
  
  vrm <- logical(nrow(x))
  erm <- logical(nrow(edge))
  trm <- logical(length(target))
  
  ## Step 1: removing any target that is terminal vertex.
  end <- FALSE
  while(!end) {
    end <- TRUE
    
    ## i=2L
    for(i in which(!trm))
      if(!any(!erm & (edge[[1L]] == target[i]))) {
        
        up <- which(!erm & (edge[[2L]] == target[i]))
        
        if(length(up) == 1L) {
          
          erm[up] <- TRUE
          vrm[target[i]] <- TRUE
          trm[i] <- TRUE
          
          s <- which(is.na(ttab$dist) & ttab$ref == edge[up,2L])
          ttab$ref[s] <- edge[up,1L]
          ttab$ladist[s] <- ttab$ladist[s] + edge$distance[up]
          
          end <- FALSE
          break
        }
      }
  }
  
  ## Step 2: removing any target that is a median vertex.
  end <- FALSE
  while(!end) {
    end <- TRUE
    
    ## i=2L
    for(i in which(!trm)) {
      
      down <- which(!erm & (edge[[1L]] == target[i]))
      
      if(length(down) == 1L) {
        
        up <- which(!erm & (edge[[2L]] == target[i]))
        
        if(length(up) == 1L) {
          
          if(!any(!erm & (edge[,1L] == edge[up,1L]) &
                  (edge[,2L] == edge[down,2L]))) {
            
            vrm[target[i]] <- TRUE
            
            edge[up,2L] <- edge[down,2L]
            dup <- edge$distance[up]
            edge$distance[up] <- sum(edge$distance[c(up,down)])
            
            erm[down] <- TRUE
            trm[i] <- TRUE
            
            ## Problem here...
            ## s <- which(ttab$ref == edge[down,1L])
            ## ttab$ref[s] <- up
            ## ttab$dist[s] <- ifelse(is.na(ttab$dist[s]), dup, ttab$dist[s] + dup)
            ## This is the old code.
            
            ## Treat vertices differently from edges
            ## Vertices:
            s <- which(is.na(ttab$dist) & (ttab$ref == edge[down,1L]))
            ttab$ref[s] <- up
            ttab$dist[s] <- dup
            
            ## Edges (if any):
            s <- which(!is.na(ttab$dist) & (ttab$ref == down))
            ttab$ref[s] <- up
            ttab$ladist[s] <- ttab$ladist[s] + ttab$dist[s]
            ttab$dist[s] <- dup
            
            end <- FALSE
            break
          }
        }
      }
    }
  }
  
  ## Step 3: Purging any new terminal vertex (ie., a previously non-terminal
  ## vertex that have become a terminal vertex following the removal of terminal
  ## target species) that is not marked as a species.
  end <- FALSE
  while(!end) {
    end <- TRUE
    
    ## i=10L
    for(i in which(!(vrm | x$species)))
      if(!any(!erm & (edge[[1L]] == i))) {
        
        up <- which(!erm & (edge[[2L]] == i))
        
        if(length(up) == 1L) {
          
          s <- which(is.na(ttab$dist) & ttab$ref == edge[up,2L])
          ttab$ref[s] <- edge[up,1L]
          ttab$ladist[s] <- ttab$ladist[s] + edge$distance[up]
          
          erm[up] <- TRUE
          vrm[edge[up,2L]] <- TRUE
          
          end <- FALSE
          break
        }
      }
  }
  
  ## Step 4: purging any new median vertex (ie., a previously non-median vertex
  ## that have become a median vertex following the previous removal of terminal
  ## vertice; either target species at step 1 or vertices not marked as a
  ## species at step 2-3) that is not-marked as a species vertex.
  end <- FALSE
  while(!end) {
    end <- TRUE
    
    ## i=3L
    for(i in which(!(vrm | x$species))) {
      
      down <- which(!erm & (edge[[1L]] == i))
      
      if(length(down) == 1L) {
        
        up <- which(!erm & (edge[[2L]] == i))
        
        if(length(up) == 1L) {
          
          if(!any(!erm & (edge[,1L] == edge[up,1L]) &
                  (edge[,2L] == edge[down,2L]))) {
            
            vrm[i] <- TRUE
            edge[up,2L] <- edge[down,2L]
            dup <- edge$distance[up]
            edge$distance[up] <- sum(edge$distance[c(up,down)])
            erm[down] <- TRUE
            
            ## Treat vertices differently from edges
            ## Vertices:
            s <- which(is.na(ttab$dist) & (ttab$ref == edge[down,1L]))
            ttab$ref[s] <- up
            ttab$dist[s] <- dup
            
            ## Edges (if any):
            ## Problem here: !is.na(dist), but ttab$ref is a vertex!
            ## s <- which(!is.na(ttab$dist) & (ttab$ref == edge[down,1L]))
            ## Old code above...
            s <- which(!is.na(ttab$dist) & (ttab$ref == down))
            ttab$ref[s] <- up
            ttab$ladist[s] <- ttab$ladist[s] + ttab$dist[s]
            ttab$dist[s] <- dup
            
            end <- FALSE
            break
          }
        }
      }
    }
  }
  
  ## Vertices that cannot be removed simply lose their species status.
  x$species[target[!trm]] <- FALSE
  
  ## Recalculating the new edge indices:
  mask <- rep(NA, nrow(edge))
  mask[!erm] <- 1L:sum(!erm)
  
  ## Changing the edge indices from the target table:
  ttab$ref[!is.na(ttab$dist)] <- mask[ttab$ref[!is.na(ttab$dist)]]
  
  ## Removing the edges marked for removal:
  edge <- edge[!erm,]
  
  ## Recalculating the new vertex indices:
  mask <- rep(NA, nrow(x))
  mask[!vrm] <- 1L:sum(!vrm)
  
  ## Changing the vertex indices from the target table:
  ttab$ref[is.na(ttab$dist)] <- mask[ttab$ref[is.na(ttab$dist)]]
  
  ## Reindexing the vertices:
  edge[[1L]] <- mask[edge[[1L]]]
  edge[[2L]] <- mask[edge[[2L]]]
  
  ## Removing the vertices marked for removal and reassign the edges:
  out <- x[!vrm,,drop=FALSE]
  edge(out) <- edge
  class(out) <- c("graph",class(out)[-which(class(out) == "graph")])
  
  if(!is.null(attr(out,"processOrder")))
    attr(out,"processOrder") <- getProcessOrder(out)
  
  if(!is.null(attr(out,"dist")))
    attr(out,"dist") <- graphDist(out)
  
  attr(out,"removedVertex") <- vrm
  attr(out,"removedEdge") <- erm
  
  list(
    x = out,
    location = ttab
  )
}
#'
