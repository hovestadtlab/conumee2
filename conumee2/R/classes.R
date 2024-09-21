##### CLASSES #####

#' CNV.anno class
#' @description Annotations required for CNV analysis are stored in this class.
#' @return \code{CNV.anno} class.
#' @details This class does not contain any sample data. Use \code{CNV.create_anno} to create.
#' @examples
#' # create object
#' anno <- CNV.create_anno()
#'
#' # general information
#' anno
#' show(anno)
#' @author Volker Hovestadt, Bjarne Daenekas \email{conumee@@hovestadt.bio}
#' @export
setClass("CNV.anno", representation(date = "character", args = "list",
                                    genome = "data.frame", gap = "GRanges", probes = "GRanges", exclude = "GRanges",
                                    detail = "GRanges", bins = "GRanges"))

#' @rdname CNV.anno-class
#' @importFrom methods show
#' @param object \code{CNV.anno} object
#' @export
setMethod("show", "CNV.anno", function(object) {
  cat("CNV annotation object\n")
  cat("   created  : ", object@date, "\n", sep = "")
  cat("  @genome   : ", nrow(object@genome), " chromosomes\n", sep = "")
  cat("  @gap      : ", length(object@gap), " regions\n", sep = "")
  cat("  @probes   : ", length(object@probes), " probes\n", sep = "")
  cat("  @exclude  : ", length(object@exclude), " regions (overlapping ",
      length(findOverlaps(object@probes, object@exclude)), " probes)\n",
      sep = "")
  cat("  @detail   : ", length(object@detail), " regions (overlapping ",
      length(findOverlaps(object@probes, object@detail)), " probes)\n",
      sep = "")
  cat("  @bins     : ", length(object@bins), " bins (min/avg/max size: ",
      object@args$bin_minsize/1000, "/", suppressWarnings(round(mean(width(object@bins))/1000,
                                                                1)), "/", object@args$bin_maxsize/1000, "kb, probes: ", object@args$bin_minprobes,
      "/", suppressWarnings(round(mean(values(object@bins)$probes), 1)),
      "/", max(values(object@bins)$probes), ")\n", sep = "")
})

#' CNV.data class
#' @description Intensities of one or multiple samples are stored in this class.
#' @return \code{CNV.data} class.
#' @details Use \code{CNV.load} to create.
#' @examples
#' # general information
#' d
#' show(d)
#'
#' # show or replace sample names
#' names(d)
#' names(d) <- toupper(names(d))
#'
#' # subset samples
#' d[1:2]
#' @author Volker Hovestadt, Bjarne Daenekas \email{conumee@@hovestadt.bio}
#' @export
setClass("CNV.data", representation(date = "character", intensity = "data.frame"))

#' @rdname CNV.data-class
#' @param object \code{CNV.data} object
setMethod("show", "CNV.data", function(object) {
  cat("CNV data object\n")
  cat("   created   : ", object@date, "\n", sep = "")
  if (length(object@intensity) == 0) {
    cat("  @intensity : unavailable, run CNV.load\n", sep = "")
  } else {
    cat("  @intensity : available (", ncol(object@intensity), " samples, ",
        nrow(object@intensity), " probes)\n", sep = "")
  }
})

#' @rdname CNV.data-class
#' @param x \code{CNV.data} object (defined by \code{Extract} generic).
#' @param i index. \code{logical}, \code{numeric} or \code{character}.
#' @export
setMethod("[", signature(x = "CNV.data"), function(x, i) {
  x@intensity <- x@intensity[, i, drop = FALSE]
  return(x)
})

#' @rdname CNV.data-class
setMethod("names", signature(x = "CNV.data"), function(x) {
  return(colnames(x@intensity))
})

#' @rdname CNV.data-class
#' @param value Replacement names.
setReplaceMethod("names", signature(x = "CNV.data"), function(x, value) {
  if (length(value) == ncol(x@intensity)) {
    colnames(x@intensity) <- value
  } else {
    stop("number of names does not fit number of samples.")
  }
  return(x)
})


#' CNV.analysis class
#' @description CNV analysis data of a single sample is stored in this class
#' @return \code{CNV.analysis} class.
#' @details Use \code{CNV.fit} to create. Modified by \code{CNV.bin}, \code{CNV.detail} and \code{CNV.segment}.
#' @examples
#' # general information
#' x
#' show(x)
#'
#' # coefficients of linear regression
#' coef(x)
#'
#' # show or replace sample name
#' names(x)
#' names(x) <- c('Sample.1', 'Sample.2')
#'
#' # subset samples
#'
#' x[1]
#' # output plots
#' CNV.genomeplot(x)
#' CNV.genomeplot(x, chr = 'chr6')
#' #CNV.detailplot(x, name = 'MYCN')
#' #CNV.detailplot_wrap(x)
#' CNV.write(x, what = 'segments')
#' @author Volker Hovestadt, Bjarne Daenekas \email{conumee@@hovestadt.bio}
#' @export
setClass("CNV.analysis", representation(name = "character", date = "character",anno = "CNV.anno", fit = "list", bin = "list", detail = "list", seg = "list"))

#' @rdname CNV.analysis-class
#' @param object \code{CNV.analysis} object
setMethod("show", "CNV.analysis", function(object) {
  cat("CNV analysis object\n")
  cat("   created   : ", object@date, "\n", sep = "")
  cat("  @name      :", colnames(object@fit$ratio), "\n", sep = " ")
  cat("  @anno      : ", nrow(object@anno@genome), " chromosomes, ",
      length(object@anno@probes), " probes, ", length(object@anno@bins),
      " bins\n", sep = "")
  if (length(object@fit) == 0) {
    cat("  @fit       : unavailable, run CNV.fit\n", sep = "")
  } else {
    cat("  @fit       : available (", ncol(object@fit$ratio), "sample(s), noise:", round(object@fit$noise,
                                                                                         2), ")\n", sep = " ")
  }
  if (length(object@bin) == 0) {
    cat("  @bin       : unavailable, run CNV.bin\n", sep = "")
  } else {
    cat("  @bin       : available (", ncol(object@fit$ratio), "sample(s), shift:", round(object@bin$shift,
                                                                                         2), ")\n", sep = " ")
  }
  if (length(object@detail) == 0) {
    cat("  @detail    : unavailable, run CNV.detail\n", sep = "")
  } else {
    cat("  @detail    : available (", length(object@detail$ratio[[1]]),
        " regions)\n", sep = "")
  }
  if (length(object@seg) == 0) {
    cat("  @seg       : unavailable, run CNV.segment\n", sep = "")
  } else {
    cat("  @seg       : available (", ncol(object@fit$ratio), " sample(s), ", sum(unlist(lapply(object@seg$summary, nrow))), " segments)\n",
        sep = "")
  }
})

#' @rdname CNV.data-class
setMethod("names", signature(x = "CNV.analysis"), function(x) {
  return(colnames(x@fit$ratio))
})

#' @rdname CNV.data-class
#' @param value Replacement names.
setReplaceMethod("names", signature(x = "CNV.analysis"), function(x, value) {
  if (length(value) == ncol(x@fit$ratio)) {
    x@name <- value
    colnames(x@fit$coef) <- value
    colnames(x@fit$ratio) <- value
    names(x@fit$noise) <- value
    if(!is.null(x@bin$ratio)){
      names(x@bin$ratio) <- value
      names(x@bin$variance) <- value
      names(x@bin$shift) <- value
    }
    if(!is.null(x@detail$ratio)){
      names(x@detail$ratio) <- value
    }
    if(!is.null(x@detail$amp.bins)){
      names(x@detail$amp.bins) <- value
      names(x@detail$del.bins) <- value
      names(x@detail$amp.detail.regions) <- value
      names(x@detail$del.detail.regions) <- value
      names(x@detail$amp.cancer.genes) <- value
      names(x@detail$del.cancer.genes) <- value
    }
    if(!is.null(x@seg$summary)){
      names(x@seg$summary) <- value
      names(x@seg$p) <- value
    }
  } else {
    stop("number of names must fit number of query samples.")
  }
  return(x)
})


#' @rdname CNV.analysis-class
#' @importFrom stats coef
#' @export
setMethod("coef", signature(object = "CNV.analysis"), function(object) {
  object@fit$coef
})

#' @rdname CNV.analysis-class
#' @param x \code{CNV.analysis} object (defined by \code{Extract} generic).
#' @param i index. \code{logical}, \code{numeric} or \code{character}.
#' @export
setMethod("[", signature(x = "CNV.analysis"), function(x, i) {
  x@fit$coef <- x@fit$coef[,i, drop = FALSE]
  x@fit$ratio <- x@fit$ratio[,i, drop = FALSE]
  x@fit$noise <- x@fit$noise[i, drop = FALSE]
  x@bin$ratio <- x@bin$ratio[i, drop = FALSE]
  x@bin$variance <- x@bin$variance[i, drop = FALSE]
  x@bin$shift <- x@bin$shift[i, drop = FALSE]
  x@detail$ratio <- x@detail$ratio[i, drop = FALSE]
  x@detail$amp.bins <- x@detail$amp.bins[i, drop = FALSE]
  x@detail$del.bins <- x@detail$del.bins[i, drop = FALSE]
  x@detail$amp.detail.regions <- x@detail$amp.detail.regions[i, drop = FALSE]
  x@detail$del.detail.regions <- x@detail$del.detail.regions[i, drop = FALSE]
  x@detail$amp.cancer.genes <- x@detail$amp.cancer.genes[i, drop = FALSE]
  x@detail$del.cancer.genes <- x@detail$del.cancer.genes[i, drop = FALSE]

  return(x)
})
