##### OUTPUT methods #####

.cumsum0 <- function(x, left = TRUE, right = FALSE, n = NULL) {
  xx <- c(0, cumsum(as.numeric(x)))
  if (!left)
    xx <- xx[-1]
  if (!right)
    xx <- head(xx, -1)
  names(xx) <- n
  xx
}

#'
#'
#' CNV.genomeplot
#' @description Create CNV plot for the whole genome or chromosomes. If the \code{CNV.analysis} object holds the information for multiple samples, the plots get either loaded individually in the graphical output or directly saved as .pdf or .png files.
#' @param object \code{CNV.analysis} object.
#' @param chr character vector. Which chromosomes to plot. Defaults to \code{'all'}.
#' @param centromere logical. Show dashed lines at centromeres? Defaults to \code{TRUE}.
#' @param detail logical. If available, include labels of detail regions? Defaults to \code{TRUE}.
#' @param bins_cex Default to \code{0.75}. Alternatively, the size of individual bin dots is inversely proportional to its variance of included probes' log2-ratios. Choose either \code{standardized} for fixed dot sizes (to make plots from different samples comparable) or \code{sample_level} (to scale the dot sizes for each sample individually)
#' @param sig_cgenes logical. Should the significant genes from the Cancer Gene Census be plotted that were identified with \code{CNV.focal}? Default to \code{FALSE}.
#' @param nsig_cgenes numeric. How many significant genes identified with \code{CNV.focal} should be plotted? Default to \code{3}. We do not recommend using values higher than 5 in order to avoid false positive results.
#' @param focal_thresholds logical. Should the dynamic thresholds be plotted? Default to \code{TRUE}. \code{sig_cgenes} must be \code{TRUE}.
#' @param main character vector. Title of the plot(s). Default to sample names. Please provide a vector of the same length as the number of samples.
#' @param ylim numeric vector. The y limits of the plot. Default to \code{c(-1.25, 1.25)}.
#' @param set_par logical. Use recommended graphical parameters for \code{oma} and \code{mar}? Defaults to \code{TRUE}. Original parameters are restored afterwards.
#' @param cols character vector. Colors to use for plotting intensity levels of bins. Centered around 0. Defaults to \code{c("darkblue","darkblue", "lightgrey", "#F16729", "#F16729")}.
#' @param directory character. Export directory for saving the files
#' @param output character. Choose between \code{output} (to display it in the graphical output) or \code{pdf} and \code{png} to save it. Defaults to \code{local}.
#' @param width numeric. Width in inches of the saved files. Defaults to \code{12}.
#' @param height numeric. Height in inches of the saved files. Defaults to \code{8}
#' @param res numeric. Resolution of the saved .png files. Defaults to \code{720}
#' @param ... Additional parameters (\code{CNV.detailplot} generic, currently not used).
#' @return \code{NULL}.
#' @details This method provides the functionality for generating CNV plots for the whole genome or defined chromosomes. Bins are shown as dots, segments are shown as lines. See parameters for more information.
#' @examples
#' # prepare
#' library(minfiData)
#' data(MsetEx)
#' d <- CNV.load(MsetEx)
#' data(detail_regions)
#' anno <- CNV.create_anno(detail_regions = detail_regions)
#'
#' # create/modify object
#' x <- CNV.segment(CNV.detail(CNV.bin(CNV.fit(query = d['GroupB_1'],
#'     ref = d[c('GroupA_1', 'GroupA_2', 'GroupA_3')], anno))))
#'
#' # output plots
#' CNV.genomeplot(x)
#' CNV.genomeplot(x, chr = 'chr6')
#' CNV.genomeplot(x, output = "pdf", directory = dir, sig_cgenes = TRUE)
#' CNV.detailplot(x, name = 'PTEN')
#' CNV.detailplot_wrap(x)
#' CNV.summaryplot(x)
#'
#' # output text files
#' CNV.write(x, what = 'segments')
#' CNV.write(x, what = 'detail')
#' CNV.write(x, what = 'bins')
#' CNV.write(x, what = 'probes')
#' @author Volker Hovestadt, Bjarne Daenekas \email{conumee@@hovestadt.bio}
#' @export
setGeneric("CNV.genomeplot", function(object, ...) {
  standardGeneric("CNV.genomeplot")
})

#' @rdname CNV.genomeplot
setMethod("CNV.genomeplot", signature(object = "CNV.analysis"), function(object, chr = "all", centromere = TRUE, detail = TRUE,
                                                                         main = NULL, sig_cgenes = FALSE, nsig_cgenes = 3, output = "local", directory = getwd(), ylim = c(-1.25, 1.25),
                                                                         bins_cex = 0.75, set_par = TRUE,
                                                                         width = 12, height = 6, res = 720, cols = c("darkblue","darkblue", "lightgrey", "#F16729", "#F16729")){

  # if(length(object@fit) == 0) stop('fit unavailable, run CNV.fit')
  if (length(object@bin) == 0)
    stop("bin unavailable, run CNV.bin")
  # if(length(object@detail) == 0) stop('bin unavailable, run CNV.detail')
  if (length(object@seg) == 0)
    stop("bin unavailable, run CNV.seg")
  if (nrow(object@fit$ratio) < 300000) {
    centromere = FALSE
  }

  if (set_par) {
    mfrow_original <- par()$mfrow
    mar_original <- par()$mar
    oma_original <- par()$oma
  }

  if (is.null(main)) {
    main = colnames(object@fit$ratio)
  }

  if (!is.null(main) & length(main) != ncol(object@fit$ratio)) {
    stop("please provide names for every sample")
  }

  for (i in 1:ncol(object@fit$ratio)) {

    message(main[i])

    if(output == "pdf"){
      p_names <- paste(directory,"/", main,"_genomeplot",".pdf",sep="")
      pdf(p_names[i], width = width, height = height)
      par(mfrow = c(1, 1), mar = c(4, 4, 4, 4), oma = c(0, 0, 0, 0))
    }

    if(output == "png"){
      p_names <- paste(directory,"/", main[i],"_genomeplot",".png",sep="")
      png(p_names[i], units = "in", width = width, height = height, res = res)
      par(mfrow = c(1, 1), mar = c(4, 4, 4, 4), oma = c(0, 0, 0, 0))
    }

    if (chr[1] == "all") {
      chr <- object@anno@genome$chr
    } else {
      chr <- intersect(chr, object@anno@genome$chr)
    }

    chr.cumsum0 <- .cumsum0(object@anno@genome[chr, "size"], n = chr)

    plot(NA, xlim = c(0, sum(as.numeric(object@anno@genome[chr, "size"])) -
                        0), ylim = ylim, xaxs = "i", xaxt = "n", yaxt = "n", xlab = NA,
         ylab = NA, main = main[i])
    abline(v = .cumsum0(object@anno@genome[chr, "size"], right = TRUE),
           col = "grey")
    if (centromere) {
      abline(v = .cumsum0(object@anno@genome[chr, "size"]) + object@anno@genome[chr,
                                                                                "pq"], col = "grey", lty = 2)
    }

    axis(1, at = .cumsum0(object@anno@genome[chr, "size"]) + object@anno@genome[chr,
                                                                                "size"]/2, labels = object@anno@genome[chr, "chr"], las = 2)
    if (all(ylim == c(-1.25, 1.25))) {
      axis(2, at = round(seq(-1.2, 1.2, 0.4), 1), las = 2)
    } else {
      axis(2, las = 2)
    }

    # ratio
    bin.ratio <- object@bin$ratio[[i]] - object@bin$shift[i]
    bin.ratio[bin.ratio < ylim[1]] <- ylim[1]
    bin.ratio[bin.ratio > ylim[2]] <- ylim[2]

    p_size <- 1/object@bin$variance[[i]][names(object@anno@bins)]

    if(bins_cex == "standardized") {
      p_size[p_size <15] <- 0.2
      p_size[p_size >= 15 & p_size <22.5] <- 0.3
      p_size[p_size >= 22.5 & p_size <30] <- 0.4
      p_size[p_size >= 30 & p_size <37.5] <- 0.5
      p_size[p_size >= 37.5 & p_size <45] <- 0.6
      p_size[p_size >= 45 & p_size <52.5] <- 0.7
      p_size[p_size >= 52.5 & p_size <60] <- 0.8
      p_size[p_size > 60] <- 0.9
    } else if(bins_cex == "sample_level") {
      b <- boxplot.stats(p_size)
      outliers <- names(b$out)
      p_size[outliers] <- as.numeric(b$stats[5])
      p_size <- round(0.7*((p_size - min(p_size))/(max(p_size) - min(p_size)))+ 0.2, digits = 2) #scaling from 0.1:0.8 for cex using predefined bins to enable comparability between plots
    } else {
      p_size <- bins_cex
    }

    bin.ratio.cols <- apply(colorRamp(cols)((bin.ratio + max(abs(ylim)))/(2 *max(abs(ylim)))),
                            1, function(x) rgb(x[1], x[2], x[3], maxColorValue = 255))

    lines(chr.cumsum0[as.vector(seqnames(object@anno@bins))] + values(object@anno@bins)$midpoint,
          bin.ratio, type = "p", pch = 16, cex = p_size, col = bin.ratio.cols)


    for (l in seq(length(object@seg$summary[[i]]$seg.median))) {
      lines(c(object@seg$summary[[i]]$loc.start[l] + chr.cumsum0[object@seg$summary[[i]]$chrom[l]],
              object@seg$summary[[i]]$loc.end[l] + chr.cumsum0[object@seg$summary[[i]]$chrom[l]]),
            rep(min(ylim[2], max(ylim[1], object@seg$summary[[i]]$seg.median[l])),
                2) - object@bin$shift[i], col = "darkblue", lwd = 2)
    }

    # detail

    if (detail & length(object@detail) > 0 & ncol(object@anno@genome) == 2) {        #mouse array

      detail.ratio <- object@detail$ratio[[i]] - object@bin$shift[i]
      detail.ratio[detail.ratio < ylim[1]] <- ylim[1]
      detail.ratio[detail.ratio > ylim[2]] <- ylim[2]
      detail.ratio.above <- (detail.ratio > 0 & detail.ratio < 0.85) |
        detail.ratio < -0.85

      lines(start(object@anno@detail) + (end(object@anno@detail) - start(object@anno@detail)) /2
            + chr.cumsum0[as.vector(seqnames(object@anno@detail))],
            detail.ratio, type = "p", pch = 16, col = "black")
      text(start(object@anno@detail) + (end(object@anno@detail) - start(object@anno@detail)) /2
           + chr.cumsum0[as.vector(seqnames(object@anno@detail))],
           ifelse(detail.ratio.above, detail.ratio, NA), labels = paste("  ", values(object@anno@detail)$name, sep = ""), adj = c(0,0.5),srt = 90, col = "black")
      text(start(object@anno@detail) + (end(object@anno@detail) - start(object@anno@detail)) /2
           + chr.cumsum0[as.vector(seqnames(object@anno@detail))],
           ifelse(detail.ratio.above, NA, detail.ratio), labels = paste(values(object@anno@detail)$name, "  ", sep = ""), adj = c(1, 0.5), srt = 90, col = "black")

    } else if(ncol(object@anno@genome) != 2 & detail & length(object@detail) > 0) { #human array

      detail.ratio <- object@detail$ratio[[i]] - object@bin$shift[i]
      detail.ratio[detail.ratio < ylim[1]] <- ylim[1]
      detail.ratio[detail.ratio > ylim[2]] <- ylim[2]
      detail.ratio.above <- (detail.ratio > 0 & detail.ratio < 0.85) |
        detail.ratio < -0.85

      lines(start(object@anno@detail) + (end(object@anno@detail) - start(object@anno@detail)) /2
            + chr.cumsum0[as.vector(seqnames(object@anno@detail))],
            detail.ratio, type = "p", pch = 16, col = "black")
      text(start(object@anno@detail) + (end(object@anno@detail) - start(object@anno@detail)) /2
           + chr.cumsum0[as.vector(seqnames(object@anno@detail))],
           ifelse(detail.ratio.above, detail.ratio, NA), labels = paste("  ", values(object@anno@detail)$name, sep = ""), adj = c(0,0.5),srt = 90, col = "black")
      text(start(object@anno@detail) + (end(object@anno@detail) - start(object@anno@detail)) /2
           + chr.cumsum0[as.vector(seqnames(object@anno@detail))],
           ifelse(detail.ratio.above, NA, detail.ratio), labels = paste(values(object@anno@detail)$name, "  ", sep = ""), adj = c(1, 0.5), srt = 90, col = "black")


      if(!length(object@detail$amp.detail.regions[[i]]) == 0 || !length(object@detail$del.detail.regions[[i]]) == 0){ #CNV focal was used

        sig.detail.regions.ratio <- c(object@detail$amp.detail.regions[[i]], object@detail$del.detail.regions[[i]])

        sig.detail.regions.ratio[sig.detail.regions.ratio < ylim[1]] <- ylim[1]
        sig.detail.regions.ratio[sig.detail.regions.ratio > ylim[2]] <- ylim[2]
        sig.detail.regions.ratio.above <- (sig.detail.regions.ratio > 0 & sig.detail.regions.ratio < 0.85) |
          sig.detail.regions.ratio < -0.85

        sig.detail.regions <- object@anno@detail[object@anno@detail$name %in% names(sig.detail.regions.ratio)]
        names(sig.detail.regions) <- sig.detail.regions$name
        sig.detail.regions.ratio <- sig.detail.regions.ratio[names(sig.detail.regions)]

        lines(start(sig.detail.regions) + (end(sig.detail.regions) - start(sig.detail.regions)) /2
              + chr.cumsum0[as.vector(seqnames(sig.detail.regions))],
              sig.detail.regions.ratio, type = "p", pch = 16, col = "red")
        text(start(sig.detail.regions) + (end(sig.detail.regions) - start(sig.detail.regions)) /2
             + chr.cumsum0[as.vector(seqnames(sig.detail.regions))],
             ifelse(sig.detail.regions.ratio.above, sig.detail.regions.ratio, NA), labels = paste("  ", names(sig.detail.regions), sep = ""), adj = c(0,0.5), srt = 90, col = "red")
        text(start(sig.detail.regions) + (end(sig.detail.regions) - start(sig.detail.regions)) /2
             + chr.cumsum0[as.vector(seqnames(sig.detail.regions))],
             ifelse(sig.detail.regions.ratio.above, NA, sig.detail.regions.ratio), labels = paste(names(sig.detail.regions), "  ", sep = ""), adj = c(1, 0.5), srt = 90, col = "red")

        if(sig_cgenes) {

          if(object@anno@args$genome == "hg38"){
            data("consensus_cancer_genes_hg38")
          } else {
            data("consensus_cancer_genes_hg19")
          }

          sig.cancer.genes.sorted <- names(sort(c(object@detail$amp.cancer.genes[[i]], abs(object@detail$del.cancer.genes[[i]])), decreasing = TRUE))

          if(nsig_cgenes > length(sig.cancer.genes.sorted)){
            nsig_cgenes <- length(sig.cancer.genes.sorted)
          }

          sig.cancer.genes.ratio <- c(object@detail$amp.cancer.genes[[i]], object@detail$del.cancer.genes[[i]])[sig.cancer.genes.sorted[1:nsig_cgenes]]

          sig.cancer.genes.ratio[sig.cancer.genes.ratio < ylim[1]] <- ylim[1]
          sig.cancer.genes.ratio[sig.cancer.genes.ratio > ylim[2]] <- ylim[2]
          sig.cancer.genes.ratio.above <- (sig.cancer.genes.ratio > 0 & sig.cancer.genes.ratio < 0.85) |
            sig.cancer.genes.ratio < -0.85

          sig.cancer.genes <- cancer_genes[names(sig.cancer.genes.ratio)]

          lines(start(sig.cancer.genes) + (end(sig.cancer.genes) - start(sig.cancer.genes)) /2
                + chr.cumsum0[as.vector(seqnames(sig.cancer.genes))],
                sig.cancer.genes.ratio, type = "p", pch = 16, col = "red")
          text(start(sig.cancer.genes) + (end(sig.cancer.genes) - start(sig.cancer.genes)) /2
               + chr.cumsum0[as.vector(seqnames(sig.cancer.genes))],
               ifelse(sig.cancer.genes.ratio.above, sig.cancer.genes.ratio, NA), labels = paste("  ", names(sig.cancer.genes), sep = ""), adj = c(0,0.5), srt = 90, col = "red")
          text(start(sig.cancer.genes) + (end(sig.cancer.genes) - start(sig.cancer.genes)) /2
               + chr.cumsum0[as.vector(seqnames(sig.cancer.genes))],
               ifelse(sig.cancer.genes.ratio.above, NA, sig.cancer.genes.ratio), labels = paste(names(sig.cancer.genes), "  ", sep = ""), adj = c(1, 0.5), srt = 90, col = "red")

        }}}

    if(is.element(output, c("pdf", "png"))){
      dev.off()
    }
  }

  if(is.element(output, c("pdf", "png"))){
    message(paste(ncol(object@fit$ratio)," files were created.", sep = ""))
  }

  if (set_par)
    par(mfrow = mfrow_original, mar = mar_original, oma = oma_original)

})



#' CNV.detailplot
#' @description Create CNV plot for detail region. If the \code{CNV.analysis} object holds the information for multiple samples, the plots get either loaded individually in the graphical output or directly saved as .pdf or .png files.
#' @param object \code{CNV.analysis} object.
#' @param name character. Name of detail region to plot.
#' @param yaxt character. Include y-axis? \code{'l'}: left, \code{'r'}: right, \code{'n'}: no. Defaults to \code{'l'}.
#' @param ylim numeric vector. The y limits of the plot. Defaults to \code{c(-1.25, 1.25)}.
#' @param directory character. Export directory for saving the files
#' @param output character. Choose between \code{pdf} and \code{png}. Defaults to \code{local}
#' @param width numeric. Width in inches of the saved files. Defaults to \code{12}.
#' @param height numeric. Height in inches of the saved files. Defaults to \code{8}
#' @param res numeric. Resolution of the saved .png files. Defaults to \code{720}
#' @param set_par logical. Use recommended graphical parameters for \code{oma} and \code{mar}? Defaults to \code{TRUE}. Original parameters are restored afterwards.
#' @param cols character vector. Colors to use for plotting intensity levels of bins. Centered around 0. Defaults to \code{c('red', 'red', 'lightgrey', 'green', 'green')}.
#' @param columns numeric. Needed for \code{detailplot_wrap}. Defaults to \code{NULL}. Do not manipulate.
#' @param ... Additional parameters (\code{CNV.detailplot} generic, currently not used).
#' @return \code{NULL}.
#' @details This method provides the functionality for generating detail regions CNV plots. Probes are shown as dots, bins are shown as lines. See parameters for more information.
#' @examples
#' # prepare
#' library(minfiData)
#' data(MsetEx)
#' d <- CNV.load(MsetEx)
#' data(detail_regions)
#' anno <- CNV.create_anno(detail_regions = detail_regions)
#'
#' # create/modify object
#' x <- CNV.segment(CNV.detail(CNV.bin(CNV.fit(query = d['GroupB_1'],
#'     ref = d[c('GroupA_1', 'GroupA_2', 'GroupA_3')], anno))))
#'
#' # output plots
#' CNV.genomeplot(x)
#' CNV.genomeplot(x, chr = 'chr6')
#' CNV.detailplot(x, name = 'PTEN', output = "pdf", directory = dir)
#' CNV.detailplot_wrap(x)
#' CNV.summaryplot(x)
#'
#' # output text files
#' CNV.write(x, what = 'segments')
#' CNV.write(x, what = 'detail')
#' CNV.write(x, what = 'bins')
#' CNV.write(x, what = 'probes')
#' @author Volker Hovestadt, Bjarne Daenekas \email{conumee@@hovestadt.bio}
#' @export
setGeneric("CNV.detailplot", function(object, ...) {
  standardGeneric("CNV.detailplot")
})

#' @rdname CNV.detailplot
setMethod("CNV.detailplot", signature(object = "CNV.analysis"),
          function(object, name, yaxt = "l", ylim = c(-1.25, 1.25), set_par = TRUE, output = "local", columns = NULL, main = NULL,
                   directory = getwd(), width = 12, height = 8, res = 720, cols = c("darkblue","darkblue", "lightgrey", "#F16729", "#F16729")) {

            if (!is.element(name, values(object@anno@detail)$name))
              stop("detail_name not in list of detail regions.")

            if (length(object@fit) == 0)
              stop("fit unavailable, run CNV.fit")
            if (length(object@bin) == 0)
              stop("bin unavailable, run CNV.bin")
            if (length(object@detail) == 0)
              stop("bin unavailable, run CNV.detail")
            # if(length(object@seg) == 0) stop('bin unavailable, run CNV.seg')
            if (is.null(columns)){
              columns = seq(ncol(object@fit$ratio))
            }
            if (is.null(main)) {
              main = paste(colnames(object@fit$ratio),"-",name)
            }
            if (set_par) {
              mfrow_original <- par()$mfrow
              mar_original <- par()$mar
              oma_original <- par()$oma
              par(mfrow = c(1, 1), mar = c(8, 4, 4, 4), oma = c(0, 0, 0, 0))
            }

            for (i in columns) {
              message(colnames(object@fit$ratio)[i])

              if(output == "pdf"){
                p_names <- paste(directory,"/",colnames(object@fit$ratio),"_",name,".pdf",sep="")
                pdf(p_names[i], width = width, height = height)
                par(mfrow = c(1, 1), mar = c(4, 4, 4, 4), oma = c(0, 0, 0, 0))
              }

              if(output == "png"){
                p_names <- paste(directory,"/",colnames(object@fit$ratio),"_",name,".png",sep="")
                png(p_names[i], units = "in", width = width, height = height, res = res)
                par(mfrow = c(1, 1), mar = c(4, 4, 4, 4), oma = c(0, 0, 0, 0))
              }

              detail.gene <- object@anno@detail[match(name, values(object@anno@detail)$name)]
              detail.region <- detail.gene
              ranges(detail.region) <- values(detail.gene)$thick

              plot(NA, xlim = c(start(detail.region), end(detail.region)), ylim = ylim,
                   xaxt = "n", yaxt = "n", xlab = NA, ylab = NA, main = main[i])
              axis(1, at = mean(c(start(detail.region), end(detail.region))), labels = as.vector(seqnames(detail.region)),
                   tick = 0, las = 1)
              axis(1, at = start(detail.region), labels = format(start(detail.region),
                                                                 big.mark = ",", scientific = FALSE), las = 2, padj = 1)
              axis(1, at = end(detail.region), labels = format(end(detail.region),
                                                               big.mark = ",", scientific = FALSE), las = 2, padj = 0)
              if (yaxt != "n")
                if (all(ylim == c(-1.25, 1.25))) {
                  axis(ifelse(yaxt == "r", 4, 2), at = round(seq(-1.2, 1.2, 0.4),
                                                             1), las = 2)
                } else {
                  axis(ifelse(yaxt == "r", 4, 2), las = 2)
                }

              axis(3, at = c(start(detail.gene), end(detail.gene)), labels = NA)

              detail.bins <- names(object@bin$ratio[[i]])[as.matrix(findOverlaps(detail.region,
                                                                                 object@anno@bins, maxgap = width(detail.region)))[, 2]]
              detail.probes <- names(object@anno@probes)[as.matrix(findOverlaps(detail.region,
                                                                                object@anno@probes, maxgap = width(detail.region)))[, 2]]

              detail.ratio <- object@fit$ratio[detail.probes,i] - object@bin$shift[i]
              detail.ratio[detail.ratio > ylim[2]] <- ylim[2]
              detail.ratio[detail.ratio < ylim[1]] <- ylim[1]
              detail.ratio.cols <- apply(colorRamp(cols)((detail.ratio + max(abs(ylim)))/(2 *
                                                                                            max(abs(ylim)))), 1, function(x) rgb(x[1], x[2], x[3], maxColorValue = 255))
              names(detail.ratio.cols) <- names(detail.ratio)
              lines(start(object@anno@probes[detail.probes]),detail.ratio,
                    type = "p", pch = 4, cex = 0.75, col = detail.ratio.cols)

              anno.bins.detail <- object@anno@bins[detail.bins]
              anno.bins.ratio <- object@bin$ratio[[i]][detail.bins] - object@bin$shift[i]
              anno.bins.ratio[anno.bins.ratio > ylim[2]] <- ylim[2]
              anno.bins.ratio[anno.bins.ratio < ylim[1]] <- ylim[1]
              lines(as.vector(rbind(rep(start(anno.bins.detail), each = 2), rep(end(anno.bins.detail),
                                                                                each = 2))), as.vector(rbind(NA, anno.bins.ratio, anno.bins.ratio, NA)),
                    col = "black", lwd = 1)

              segs <- GRanges(seqnames = object@seg$summary[[i]]$chrom, IRanges(start = object@seg$summary[[i]]$loc.start, end = object@seg$summary[[i]]$loc.end),
                              seqinfo = Seqinfo(genome = "hg19"))
              segs$seg.median <- object@seg$summary[[i]]$seg.median - object@bin$shift[i]
              segs <- segs[subjectHits(findOverlaps(query = detail.gene, subject = segs, type = "any"))]
              lines(as.vector(rbind(rep(start(segs), each = 2), rep(end(segs), each = 2))),
                    as.vector(rbind(NA, segs$seg.median, segs$seg.median, NA)),col = "darkblue", lwd = 2)

              if(is.element(output, c("pdf", "png"))){
                dev.off()
              }

            }

            if(is.element(output, c("pdf", "png"))){
              message(paste(ncol(object@fit$ratio)," files were created.", sep = ""))
            }

            if (set_par)
              par(mfrow = mfrow_original, mar = mar_original, oma = oma_original)

          })


#' CNV.detailplot_wrap
#' @description Create CNV plots for all detail regions. If the \code{CNV.analysis} object holds the information for multiple samples, the plots get either loaded individually in the graphical output or directly saved as .pdf or .png files.
#' @param object \code{CNV.analysis} object.
#' @param set_par logical. Use recommended graphical parameters for \code{oma} and \code{mar}? Defaults to \code{TRUE}. Original parameters are restored afterwards.
#' @param directory character. Export directory for saving the files
#' @param output character. Choose between \code{pdf} and \code{png}. Defaults to \code{local}
#' @param width numeric. Width in inches of the saved files. Defaults to \code{12}.
#' @param height numeric. Height in inches of the saved files. Defaults to \code{8}
#' @param res numeric. Resolution of the saved .png files. Defaults to \code{720}
#' @param main character. Used for \code{CNV.detailplot}. Do not manipulate
#' @param header character vector. Title of the plot(s). Defaults to sample names. Please provide a vector of the same length than the number of samples.
#' @param ... Additional paramters supplied to \code{CNV.detailplot}.
#' @return \code{NULL}.
#' @details This method is a wrapper of the \code{CNV.detailplot} method to plot all detail regions.
#' @examples
#' # prepare
#' library(minfiData)
#' data(MsetEx)
#' d <- CNV.load(MsetEx)
#' data(detail_regions)
#' anno <- CNV.create_anno(detail_regions = detail_regions)
#'
#' # create/modify object
#' x <- CNV.segment(CNV.detail(CNV.bin(CNV.fit(query = d['GroupB_1'],
#'     ref = d[c('GroupA_1', 'GroupA_2', 'GroupA_3')], anno))))
#'
#' # output plots
#' CNV.genomeplot(x)
#' CNV.genomeplot(x, chr = 'chr6')
#' CNV.detailplot(x, name = 'PTEN')
#' CNV.detailplot_wrap(x)
#' CNV.summaryplot(x)
#'
#' # output text files
#' CNV.write(x, what = 'segments')
#' CNV.write(x, what = 'detail')
#' CNV.write(x, what = 'bins')
#' CNV.write(x, what = 'probes')
#' @author Volker Hovestadt, Bjarne Daenekas \email{conumee@@hovestadt.bio}
#' @export
setGeneric("CNV.detailplot_wrap", function(object, ...) {
  standardGeneric("CNV.detailplot_wrap")
})

#' @rdname CNV.detailplot_wrap
setMethod("CNV.detailplot_wrap", signature(object = "CNV.analysis"), function(object,
                                                                              set_par = TRUE, main = NULL, header = NULL, output = "local", directory = getwd(), width = 12, height = 8, res = 720,...) {
  if (length(object@fit) == 0)
    stop("fit unavailable, run CNV.fit")
  if (length(object@bin) == 0)
    stop("bin unavailable, run CNV.bin")
  if (length(object@detail) == 0)
    stop("bin unavailable, run CNV.detail")
  # if(length(object@seg) == 0) stop('bin unavailable, run CNV.seg')

  if (set_par) {
    mfrow_original <- par()$mfrow
    mar_original <- par()$mar
    oma_original <- par()$oma
    par(mfrow = c(1, length(object@anno@detail) + 2), mar = c(8, 0,
                                                              4, 0), oma = c(0, 0, 4, 0))
  }

  if (is.null(header)) {
    header <- colnames(object@fit$ratio)
  } else if (length(header) != ncol(object@fit$ratio)) {
    stop("please provide names for all samples")
  }

  for (i in 1:ncol(object@fit$ratio)) {

    if (output == "pdf"){
      p_names <- paste(directory,"/",header,"_","detailplot_wrap",".pdf",sep="")
      pdf(p_names[i], width = width, height = height)
      par(mfrow = c(1, length(object@anno@detail) + 2), mar = c(8, 0, 4, 0), oma = c(0, 0, 4, 0))
    }

    if (output == "png"){
      p_names <- paste(directory,"/",header,"_","detailplot_wrap",".png",sep="")
      png(p_names[i], width = width, height = height, units = "in", res = res)
      par(mfrow = c(1, length(object@anno@detail) + 2), mar = c(8, 0, 4, 0), oma = c(0, 0, 4, 0))
    }

    frame()
    for (l in seq(length(object@anno@detail))) {
      if (l == 1) {
        CNV.detailplot(object, columns = i, name = object@anno@detail$name[l],
                       main = rep(object@anno@detail$name[l], ncol(object@fit$ratio)), yaxt = "l", set_par = FALSE, ...)
      } else if (l == length(object@anno@detail)) {
        CNV.detailplot(object, columns = i, name = object@anno@detail$name[l],
                       main = rep(object@anno@detail$name[l], ncol(object@fit$ratio)), yaxt = "r", set_par = FALSE, ...)
      } else {
        CNV.detailplot(object, columns = i, name = object@anno@detail$name[l],
                       main = rep(object@anno@detail$name[l], ncol(object@fit$ratio)),
                       yaxt = "n", set_par = FALSE, ...)
      }
    }

    frame()
    title(header[i], outer = TRUE)

    if(is.element(output, c("pdf", "png"))){
      dev.off()
    }
  }

  if(is.element(output, c("pdf", "png"))){
    message(paste(ncol(object@fit$ratio)," files were created.", sep = ""))
  }

  if (set_par) {
    par(mfrow = mfrow_original, mar = mar_original, oma = oma_original)
  }
})

#' CNV.summaryplot
#' @description Create a summaryplot that shows the CNVs in a all the query samples in the CNV.analysis object.
#' @param object \code{CNV.analysis} object.
#' @param set_par logical. Use recommended graphical parameters for \code{oma} and \code{mar}? Defaults to \code{TRUE}. Original parameters are restored afterwards.
#' @param main character. Specify the title of the plot. Defaults to \code{NULL}.
#' @param threshold numeric. Threshold for determining the copy number state. Defaults to \code{0.1}. See Details for details.
#' @param output character. Choose between \code{output} (to display it in the graphical output) or \code{pdf} and \code{png} to save it. Defaults to \code{local}.
#' @param directory character. Export directory for saving the files
#' @param output character. Choose between \code{output} (to display it in the graphical output) or \code{pdf} and \code{png} to save it. Defaults to \code{local}.
#' @param width numeric. Width in inches of the saved files. Defaults to \code{12}.
#' @param height numeric. Height in inches of the saved files. Defaults to \code{8}
#' @param res numeric. Resolution of the saved .png files. Defaults to \code{720}
#' @param ... Additional parameters (\code{CNV.write} generic, currently not used).
#' @details This function creates a plot that illustrates the changes in copy number states within the set of query samples that are stored in the \code{CNV.analysis} object. The y axis is showing the percentage of samples that are exhibiting a CNV at the genomic location shown on the x axis. The threshold for the log2-ratio to identify gains or losses is \code{0.1} by default.
#' @return \code{NULL}
#' @examples
#' # prepare
#' library(minfiData)
#' data(MsetEx)
#' d <- CNV.load(MsetEx)
#' data(detail_regions)
#' anno <- CNV.create_anno(detail_regions = detail_regions)
#'
#' # create/modify object
#' x <- CNV.segment(CNV.detail(CNV.bin(CNV.fit(query = d['GroupB_1'],
#'     ref = d[c('GroupA_1', 'GroupA_2', 'GroupA_3')], anno))))
#'
#' # output plots
#' CNV.genomeplot(x)
#' CNV.genomeplot(x, chr = 'chr6')
#' CNV.detailplot(x, name = 'PTEN')
#' CNV.detailplot_wrap(x)
#' CNV.summaryplot(x, threshold = 0.15)
#' CNV.heatmap(x)
#'
#' # output text files
#' CNV.write(x, what = 'segments')
#' CNV.write(x, what = 'detail')
#' CNV.write(x, what = 'bins')
#' CNV.write(x, what = 'probes')
#' @author Bjarne Daenekas \email{conumee@@hovestadt.bio}
#' @export
setGeneric("CNV.summaryplot", function(object, ...) {
  standardGeneric("CNV.summaryplot")
})

#' @rdname CNV.summaryplot
setMethod("CNV.summaryplot", signature(object = "CNV.analysis"), function(object,
                                                                          set_par = TRUE, main = NULL, output = "local", directory = getwd(), width = 12, height = 6, res = 720, threshold = 0.1,...) {

  if (set_par) {
    mfrow_original <- par()$mfrow
    mar_original <- par()$mar
    oma_original <- par()$oma
  }

  if(ncol(object@fit$ratio) <= 1) {
    stop("Please use multiple query samples to ceate a summaryplot")
  }

  if(output == "pdf"){
    p_names <- paste(directory,"/","genome_summaryplot",".pdf",sep="")
    pdf(p_names, width = width, height = height)
    par(mfrow = c(1, 1), mar = c(4, 4, 4, 4), oma = c(0, 0, 0, 0))
  }

  if(output == "png"){
    p_names <- paste(directory,"/", "genome_summaryplot",".png",sep="")
    png(p_names, units = "in", width = width, height = height, res = res)
    par(mfrow = c(1, 1), mar = c(4, 4, 4, 4), oma = c(0, 0, 0, 0))
  }

  y <- CNV.write(object, what = "threshold", threshold = threshold)

  message("creating summaryplot")

  segments.i <- GRanges(seqnames=y$Chromosome,ranges=IRanges(y$Start_Position, y$End_Position))
  segments_chromosomes <- GRanges(seqnames = object@anno@genome$chr, ranges = IRanges(start = 1, end = object@anno@genome$size))
  segments <- c(segments.i, segments_chromosomes)
  d_segments <-as.data.frame(GenomicRanges::disjoin(segments))

  overview <- as.data.frame(matrix(nrow = 0, ncol = 4))
  for (i in 1:nrow(d_segments)) {
    x <- d_segments[i,]
    involved_segements <- y[y$Chromosome == x$seqnames & y$Start_Position <= x$start & y$End_Position >= x$end,]
    balanced <- sum(involved_segements$Alteration == "balanced")
    gain <- sum(involved_segements$Alteration == "gain")
    loss <- sum(involved_segements$Alteration == "loss")
    c(as.character(x$seqnames), balanced, gain, loss)
    overview <- rbind(overview, c(as.character(x$seqnames), balanced, gain, loss))
  }

  colnames(overview) <- c("disjoined_segment", "count_balanced", "count_gains", "count_losses")
  overview$count_balanced <- as.numeric(overview$count_balanced)
  overview$count_gains <- as.numeric(overview$count_gains)
  overview$count_losses <- as.numeric(overview$count_losses)

  d_segments$gains <- overview$count_gains/length(unique(y$Sample))*100
  d_segments$losses <- overview$count_losses/length(unique(y$Sample))*100
  d_segments$balanced <- overview$count_balanced/length(unique(y$Sample))*100

  segments_pl <- d_segments[rep(1:nrow(d_segments),each=2),]
  odd_indexes<-seq(1,nrow(segments_pl),2)
  even_indexes<-seq(2,nrow(segments_pl),2)
  segments_pl$xpos<-NA
  segments_pl$xpos[odd_indexes]<-segments_pl$start[odd_indexes]
  segments_pl$xpos[even_indexes]<-segments_pl$end[even_indexes]

  segments_pl$seqnames <- factor(segments_pl$seqnames, levels = object@anno@genome$chr)
  segments_pl <- segments_pl[order(segments_pl$seqnames, segments_pl$start),]

  par(mfrow = c(1, 1), mgp=c(4,1,0), mar = c(4, 8, 4, 4), oma = c(0, 0, 0, 0))
  plot(NA, xlim = c(0, sum(as.numeric(object@anno@genome$size))), ylim = c(-100, 100), xaxs = "i", xaxt = "n", yaxt = "n",
       xlab = NA, ylab = "percentage of samples exhibiting the CNV [%]", main = main, cex=1.5, cex.lab=1, cex.axis=1, cex.main=1.5)

  abline(v = .cumsum0(object@anno@genome$size, right = TRUE),
         col = "grey")
  abline(v = .cumsum0(object@anno@genome$size) + object@anno@genome$pq,
         col = "grey", lty = 2)
  axis(1, at = .cumsum0(object@anno@genome$size) + object@anno@genome$size/2,
       labels = object@anno@genome$chr, las = 2)

  axis(2, las = 2, lty=1, at = seq(0, 100, 20), labels=abs(seq(0, 100, 20)),cex.axis=1)
  axis(2, las = 2, lty=1, at = seq(-100, 0, 20), labels=abs(seq(-100, 0, 20)),cex.axis=1)

  chr = object@anno@genome$chr
  chr.cumsum0 <- .cumsum0(object@anno@genome[chr, "size"], n = chr)

  polygon(as.numeric(chr.cumsum0[match(segments_pl$seqnames, names(chr.cumsum0))]) + segments_pl$xpos, segments_pl$gains, col="#F16729",lwd=1.5)
  polygon(as.numeric(chr.cumsum0[match(segments_pl$seqnames, names(chr.cumsum0))]) + segments_pl$xpos, -segments_pl$loss, col="darkblue",lwd=1.5)

  if(is.element(output, c("pdf", "png"))){
    dev.off()
    message("file was created")
  }

  if (set_par)
    par(mfrow = mfrow_original, mar = mar_original, oma = oma_original)
})

#' CNV.heatmap
#' @description Create a heatmap to illustrate CNVs in a set of query samples. Colors correspond to the other plots.
#' @param object \code{CNV.analysis} object.
#' @param set_par logical. Use recommended graphical parameters for \code{oma} and \code{mar}? Defaults to \code{TRUE}. Original parameters are restored afterwards.
#' @param main character. Specify the title of the plot. Defaults to \code{NULL}.
#' @param hclust logical. Should hierarchical clustering be performed? Default to \code{TRUE}.
#' @param hclust_method character. The agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#' @param dist_method character. The distance measure to be used. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski". Any unambiguous substring can be given.
#' @param cexRow numeric. Adjust the font size for row labeling. Default to 0.5.
#' @param zlim numeric. The minimum and maximum z values for which colors should be plotted. Default to \code{c(-0.5,0.5)}.
#' @param useRaster logical. Should a bitmap raster be used to create the plot instead of polygons? Default to \code{TRUE}.
#' @param output character. Choose between \code{output} (to display it in the graphical output) or \code{pdf} and \code{png} to save it. Defaults to \code{local}.
#' @param directory character. Export directory for saving the files
#' @param output character. Choose between \code{output} (to display it in the graphical output) or \code{pdf} and \code{png} to save it. Defaults to \code{local}.
#' @param width numeric. Width in inches of the saved files. Defaults to \code{12}.
#' @param height numeric. Height in inches of the saved files. Defaults to \code{8}
#' @param res numeric. Resolution of the saved .png files. Defaults to \code{720}
#' @param ... Additional parameters
#' @examples
#' #' # prepare
#' library(minfiData)
#' data(MsetEx)
#' d <- CNV.load(MsetEx)
#' data(detail_regions)
#' anno <- CNV.create_anno(detail_regions = detail_regions)
#'
#' # create/modify object
#' x <- CNV.segment(CNV.detail(CNV.bin(CNV.fit(query = d['GroupB_1'],
#'     ref = d[c('GroupA_1', 'GroupA_2', 'GroupA_3')], anno))))
#'
#' # output plots
#' CNV.genomeplot(x)
#' CNV.genomeplot(x, chr = 'chr6')
#' CNV.detailplot(x, name = 'PTEN')
#' CNV.detailplot_wrap(x)
#' CNV.summaryplot(x)
#' CNV.heatmap(x)
#'
#' # output text files
#' CNV.write(x, what = 'segments')
#' CNV.write(x, what = 'detail')
#' CNV.write(x, what = 'bins')
#' CNV.write(x, what = 'probes')
#' @author Bjarne Daenekas \email{conumee@@hovestadt.bio}
#' @export
setGeneric("CNV.heatmap", function(object, ...) {
  standardGeneric("CNV.heatmap")
})

#' @rdname CNV.heatmap
setMethod("CNV.heatmap", signature(object = "CNV.analysis"), function(object,
                                                                      set_par = TRUE, main = NULL, hclust = TRUE, hclust_method = "average", dist_method = "euclidian", cexRow = 1/2, zlim = c(-0.5,0.5), useRaster = TRUE, output = "local",
                                                                      directory = getwd(), width = 8, height = 8, res = 720,...) {

  if (set_par) {
    mfrow_original <- par()$mfrow
    mar_original <- par()$mar
    oma_original <- par()$oma
  }

  options(max.print = 1000)
  options(stringsAsFactors = FALSE)
  options(scipen = 999)

  if(output == "pdf"){
    p_names <- paste(directory,"/","genome_heatmap",".pdf",sep="")
    pdf(p_names, width = width, height = height)
    par(mfrow = c(1, 1), mar = c(4, 4, 4, 4), oma = c(0, 0, 0, 0))
  }

  if(output == "png"){
    p_names <- paste(directory,"/", "genome_heatmap",".png",sep="")
    png(p_names, units = "in", width = width, height = height, res = res)
    par(mfrow = c(1, 1), mar = c(4, 4, 4, 4), oma = c(0, 0, 0, 0))
  }

  bins <- CNV.write(object, what = "bins")
  annotation <- bins[,c(1:4)]
  annotation$X <- 1:nrow(annotation)
  bins <- bins[,-c(1:4)]
  bins <- t(bins)

  b <- which(head(annotation$Chromosome, -1) != tail(annotation$Chromosome, -1))  # between chromosomes
  l <- round(sapply(split(annotation$X, annotation$Chromosome), mean))
  ll <- rep(NA, nrow(annotation))
  ll[l] <- names(l)

  my_palette <- colorRampPalette(c("darkblue","darkblue", "white", "#F16729", "#F16729"))(n = 1000)

  if(hclust){

    bins.dist <- dist(bins, method = dist_method)
    bins.hc <- hclust(bins.dist, method = hclust_method)


    heatmap(bins, Colv = NA, Rowv = as.dendrogram(bins.hc), scale="n", useRaster = useRaster, main = main,
            col = my_palette, cexRow = cexRow, zlim = zlim, add.expr = abline(v=b), labCol = ll, cexCol = 1)

  } else {

    heatmap(bins, Colv = NA, Rowv = NA, scale="n", useRaster = useRaster, main = main,
            col = my_palette, cexRow = cexRow, zlim = zlim, add.expr = abline(v=b), labCol = ll, cexCol = 1)

  }

  if(is.element(output, c("pdf", "png"))){
    dev.off()
    message("file was created")

    if (set_par)
      par(mfrow = mfrow_original, mar = mar_original, oma = oma_original)
  }

})

#' CNV.write
#' @description Output CNV analysis results as table.
#' @param object \code{CNV.analysis} object.
#' @param file Path where output file should be written to. Defaults to \code{NULL}: No file is written, table is returned as data.frame object.
#' @param what character. This should be (an unambiguous abbreviation of) one of \code{'probes'}, \code{'bins'}, \code{'detail'}, \code{'segments'}, \code{gistic}, \code{threshold} (for CNV.summaryplot) or \code{focal}. Defaults to \code{'segments'}.
#' @param threshold numeric. Threshold for determining the copy number state. Defaults to \code{0.1}. See Description for details.
#' @param ... Additional parameters (\code{CNV.write} generic, currently not used).
#' @details  Function shows the output of the CNV analysis with conumee 2. To use the results as input for GISTIC choose \code{what = 'gistic'}. To assign the resulting segments to their copy number state and their size (focal, arm-level or whole chromosome) choose \code{what = 'overview'}. The threshold for the log2-ratio to identify gains or losses is \code{0.1} by default. To access the results from the Segmented Block Bootstrapping, use \code{what = focal}.
#' @examples
#' # prepare
#' library(minfiData)
#' data(MsetEx)
#' d <- CNV.load(MsetEx)
#' data(detail_regions)
#' anno <- CNV.create_anno(detail_regions = detail_regions)
#'
#' # create/modify object
#' x <- CNV.segment(CNV.detail(CNV.bin(CNV.fit(query = d['GroupB_1'],
#'     ref = d[c('GroupA_1', 'GroupA_2', 'GroupA_3')], anno))))
#'
#' # output plots
#' CNV.genomeplot(x)
#' CNV.genomeplot(x, chr = 'chr6')
#' CNV.detailplot(x, name = 'PTEN')
#' CNV.detailplot_wrap(x)
#'
#' # output text files
#' CNV.write(x, what = 'segments')
#' CNV.write(x, what = 'detail')
#' CNV.write(x, what = 'bins')
#' CNV.write(x, what = 'probes')
#' CNV.write(x, what = 'gistic')
#' CNV.write(x, what = 'overview')
#' CNV.write(x, what = 'focal')
#' @return if parameter \code{file} is not supplied, the table is returned as a \code{data.frame} object.
#' @author Bjarne Daenekas, Volker Hovestadt \email{conumee@@hovestadt.bio}
#' @export
setGeneric("CNV.write", function(object, ...) {
  standardGeneric("CNV.write")
})

#' @rdname CNV.write
setMethod("CNV.write", signature(object = "CNV.analysis"), function(object, file = NULL, what = "segments", threshold = 0.1) {
  w <- pmatch(what, c("probes", "bins", "detail", "segments", "gistic", "threshold", "focal"))
  if (w == 1) {
    if (length(object@fit) == 0)
      stop("fit unavailable, run CNV.fit")
    if (!is.null(file))
      if (!grepl(".igv$", file))
        warning("filename does not end in .igv")

    x <- data.frame(Chromosome = as.vector(seqnames(object@anno@probes)),
                    Start = start(object@anno@probes) - 1, End = end(object@anno@probes),
                    Feature = rownames(object@fit$ratio), row.names = NULL)

    for (i in 1:ncol(object@fit$ratio)) {
      x <- cbind(x,round(object@fit$ratio[,i] - object@bin$shift[i], 3))
    }
    colnames(x) <- c("Chromosome","Start", "End","Feature", colnames(object@fit$ratio))
  } else if (w == 2) {
    if (length(object@bin) == 0)
      stop("bin unavailable, run CNV.bin")
    if (!is.null(file))
      if (!grepl(".igv$", file))
        warning("filename does not end in .igv")

    x <- data.frame(Chromosome = as.vector(seqnames(object@anno@bins)),
                    Start = start(object@anno@bins), End = end(object@anno@bins),
                    Feature = names(object@anno@bins), row.names = NULL)
    for (i in 1:ncol(object@fit$ratio)){
      x <- cbind(x, round(object@bin$ratio[[i]] - object@bin$shift[i], 3))
    }
    colnames(x) <- c("Chromosome","Start", "End","Feature", colnames(object@fit$ratio))
  } else if (w == 3) {
    if (length(object@detail) == 0)
      stop("detail unavailable, run CNV.bin")
    if (!is.null(file))
      if (!grepl(".txt$", file))
        warning("filename does not end in .txt")

    x <- vector(mode = "list", length = ncol(object@fit$ratio))

    for (i in 1:ncol(object@fit$ratio)) {
      x[[i]] <- data.frame(Chromosome= as.vector(seqnames(object@anno@detail)),
                           Start = start(object@anno@detail), End = end(object@anno@detail),
                           Name = names(object@detail$probes), Sample = colnames(object@fit$ratio)[i],
                           Value = round(object@detail$ratio[[i]] - object@bin$shift[i], 3), Probes = as.numeric(object@detail$probes), row.names = NULL)
    }
    names(x) <- colnames(object@fit$ratio)
  } else if (w == 4) {
    if (length(object@seg) == 0)
      stop("seg unavailable, run CNV.bin")
    if (!is.null(file))
      if (!grepl(".seg$", file))
        warning("filename does not end in .seg")

    x <- vector(mode = "list", length = ncol(object@fit$ratio))

    for (i in 1:ncol(object@fit$ratio)) {
      x[[i]] <- data.frame(ID = colnames(object@fit$ratio)[i], chrom = object@seg$summary[[i]]$chrom,
                           loc.start = object@seg$summary[[i]]$loc.start, loc.end = object@seg$summary[[i]]$loc.end,
                           num.mark = object@seg$summary[[i]]$num.mark, bstat = object@seg$p[[i]]$bstat,
                           pval = object@seg$p[[i]]$pval, seg.mean = round(object@seg$summary[[i]]$seg.mean -
                                                                             object@bin$shift[i], 3), seg.median = round(object@seg$summary[[i]]$seg.median -
                                                                                                                           object@bin$shift[i], 3), row.names = NULL)
    }
    names(x) <- colnames(object@fit$ratio)
  } else if (w == 5) {
    if (length(object@seg) == 0)
      stop("seg unavailable, run CNV.bin")
    if (!is.null(file))
      if (!grepl(".seg$", file))
        warning("filename does not end in .seg")
    # seg format, last numeric column is used in igv
    if(nrow(object@anno@genome) == 19) {
      warning("GISTIC is not compatible with Illumina Methylation arrays for mice.")
    }
    x <- data.frame(matrix(ncol = 0, nrow = 0))
    for (i in 1:ncol(object@fit$ratio)) {
      y <- object@seg$summary[[i]]
      y$seg.median <- round(object@seg$summary[[i]]$seg.median-object@bin$shift[i], 3)
      y <- y[-which(is.element(y$chrom, c("chrX", "chrY"))),]
      x <- rbind(x, y)
    }
    x <- x[,-c(6,7,9)]
    colnames(x) <- c("Sample", "Chromosome", "Start_Position", "End_Position", "Num_Markers", "Seg.CN")
    x$Chromosome <- as.numeric(gsub("chr", "", x$Chromosome))
  } else if (w == 6) {
    if (length(object@seg) == 0)
      stop("seg unavailable, run CNV.segment")
    if (!is.null(file))
      if (!grepl(".seg$", file))
        warning("filename does not end in .seg")
    # seg format, last numeric column is used in igv
    x <- data.frame(matrix(ncol = 0, nrow = 0))
    for (i in 1:ncol(object@fit$ratio)) {
      y <- object@seg$summary[[i]]
      y$seg.median <- round(object@seg$summary[[i]]$seg.median -object@bin$shift[i], 3)
      x <- rbind(x, y)
    }
    x <- x[,-c(6,7,9)]
    colnames(x) <- c("Sample", "Chromosome", "Start_Position", "End_Position", "Num_Markers", "Seg.CN")
    new_col <- replicate(nrow(x), "balanced")
    new_col[which(x$Seg.CN >= threshold)] <- "gain"
    new_col[which(x$Seg.CN <= -threshold)] <- "loss"
    x$Alteration <- new_col
  } else if (w == 7){
    stop("Please run CNV.focal")
    x <- vector(mode='list', length = 6)
    x[[1]] <- object@detail$amp.bins
    x[[2]] <- object@detail$del.bins
    x[[3]] <- object@detail$amp.detail.regions
    x[[4]] <- object@detail$del.detail.regions
    x[[5]] <- object@detail$amp.cancer.genes
    x[[6]] <- object@detail$del.cancer.genes
    names(x) <- c("bins within amplified regions", "bins within lost regions", "amplified detail regions",
                  "deleted detail regions", "amplified genes from the Cancer Gene Census", "deleted genes from the Cancer Gene Census")
  } else{
    stop("value for what is ambigious.")
  }
  if (is.null(file)) {
    return(x)
  } else {
    write.table(x, file = file, quote = FALSE, sep = "\t", row.names = FALSE)
  }
})


#' CNV.plotly
#'
#' \code{CNV.plotly} plots an interactive copy number profile.
#'
#' @param x A \code{CNV.analysis} object after \code{CNV.segment} and \code{CNV.detail} is performed.
#' @param sample character. Name of the single sample that should be plotted. Default to first sample in the set of query samples. Check sample names with \code{colnames(object@@fit$ratio)}
#' @export
#' @import ggplot2
#' @import plotly

CNV.plotly <- function(x, sample = colnames(x@fit$ratio)[1]){

  if (!any(colnames(x@fit$coef) == sample)){
    stop(message("Please provide the correct sample name."))
  }

  if(is.null(x@anno@detail)){
    stop("Please use CNV.detail prior to CNV.plotly.")
  }

  sample_n <- which(colnames(x@fit$coef) == sample)

  ylim = c(-1.25, 1.25)
  bin.ratio <- x@bin$ratio[[sample_n]] - x@bin$shift[sample_n]
  bin.ratio[bin.ratio < ylim[1]] <- ylim[1]
  bin.ratio[bin.ratio > ylim[2]] <- ylim[2]
  cols2 = c("darkblue","darkblue", "lightgrey", "#F16729", "#F16729")

  chr <- x@anno@genome$chr
  chr.cumsum0 <- .cumsum0(x@anno@genome[chr, "size"], n = chr)

  y <- chr.cumsum0[as.vector(seqnames(x@anno@bins))] + values(x@anno@bins)$midpoint

  chrs <- .cumsum0(x@anno@genome[chr, "size"], right = TRUE)
  chr.cumsum0 <- .cumsum0(x@anno@genome[chr, "size"], n = chr)

  if (ncol(x@anno@genome) == 3){
    chrspq <- .cumsum0(x@anno@genome[chr, "size"]) + x@anno@genome[chr,"pq"]
  }

  tickl <- .cumsum0(x@anno@genome[chr, "size"]) + x@anno@genome[chr,"size"]/2

  cols = c("darkblue","darkblue", "lightgrey", "#F16729", "#F16729")
  bin.ratio.cols <- apply(colorRamp(cols)((bin.ratio + max(abs(ylim)))/(2 *max(abs(ylim)))),
                          1, function(x) rgb(x[1], x[2], x[3], maxColorValue = 255))

  df <- data.frame(y,bin.ratio,bin.ratio.cols)

  xs<-  x@seg$summary[[sample_n]]$loc.start + chr.cumsum0[x@seg$summary[[sample_n]]$chrom]
  xe <- x@seg$summary[[sample_n]]$loc.end + chr.cumsum0[x@seg$summary[[sample_n]]$chrom]
  ys <- x@seg$summary[[sample_n]]$seg.median - x@bin$shift[sample_n]
  ye <- x@seg$summary[[sample_n]]$seg.median - x@bin$shift[sample_n]

  df2 <- data.frame(xs,xe,ys,ye)

  detail.ratio <- x@detail$ratio[[sample_n]] - x@bin$shift[sample_n]
  detail.ratio[detail.ratio < ylim[1]] <- ylim[1]
  detail.ratio[detail.ratio > ylim[2]] <- ylim[2]
  detail.ratio.above <- (detail.ratio > 0 & detail.ratio < 0.85) | detail.ratio < -0.85

  detail.x <- start(x@anno@detail) + (end(x@anno@detail) - start(x@anno@detail)) /2 + chr.cumsum0[as.vector(seqnames(x@anno@detail))]

  df3 <- data.frame(detail.ratio,detail.x,names=values(x@anno@detail)$name)

  if (ncol(x@anno@genome) == 3){
    p <- ggplot(df,aes(x=y, y=bin.ratio)) +
      geom_point(colour=bin.ratio.cols,size=.5) + geom_vline(xintercept = chrs,color="black",size=0.1) +
      theme_bw()+
      ggtitle(names(x@fit$coef[sample_n]))+
      theme(text = element_text(family = "Arial"))+
      geom_vline(xintercept = chrspq,color="black",size=.1,linetype="dotted")+ ylim(-1.25, 1.25)+
      geom_segment(aes(x = xs, y = ys, xend = xe, yend = ye),size=.5, data = df2,color="darkblue")+
      xlab("")+
      ylab("")+
      geom_point(aes(x=detail.x,y=detail.ratio),size=1.15,alpha=0.9,data=df3,color="red") +
      scale_x_continuous(breaks=tickl,labels = c(chr))+#,expand = c(0, 0),limits = c(0, max(x)))+
      theme(axis.text.x= element_text(size=10,angle = 90))

    ggp <- ggplotly(p)
    ggpb <- plotly_build(ggp)

    ggpb$x$data[[1]]$text <- paste0(seqnames(x@anno@bins),"<br>","start: ",
                                    start(x@anno@bins),"<br>","end: ",end(x@anno@bins),"<br>",
                                    "probes: ",values(x@anno@bins)$probes, "<br>", "genes: ", x@anno@bins$genes)
    ggpb$x$data[[2]]$text <- ""
    ggpb$x$data[[3]]$text <- ""
    ggpb$x$data[[4]]$text <- paste0(x@seg$summary[[sample_n]]$chrom,"<br>","start: ",
                                    x@seg$summary[[sample_n]]$loc.start,"<br>","end: ",x@seg$summary[[sample_n]]$loc.end,"<br>",
                                    "median: ",x@seg$summary[[sample_n]]$seg.median)
    ggpb$x$data[[5]]$text <- values(x@anno@detail)$name
    ggpb%>%suppressWarnings(toWebGL())
  }

  if (ncol(x@anno@genome) == 2){
    p <- ggplot(df,aes(x=y, y=bin.ratio)) +
      geom_point(colour=bin.ratio.cols,size=.5) + geom_vline(xintercept = chrs,color="black",size=0.1) +
      theme_bw()+
      ggtitle(names(x@fit$coef[sample_n]))+
      theme(text = element_text(family = "Arial"))+
      ylim(-1.25, 1.25)+
      geom_segment(aes(x = xs, y = ys, xend = xe, yend = ye),size=.5, data = df2,color="darkblue")+
      xlab("")+
      ylab("")+
      geom_point(aes(x=detail.x,y=detail.ratio),size=1.15,alpha=0.9,data=df3,color="red") +
      scale_x_continuous(breaks=tickl,labels = c(chr))+#,expand = c(0, 0),limits = c(0, max(x)))+
      theme(axis.text.x= element_text(size=10,angle = 90))


    ggp <- ggplotly(p)
    ggpb <- plotly_build(ggp)

    ggpb$x$data[[1]]$text <- paste0(seqnames(x@anno@bins),"<br>","start: ",
                                    start(x@anno@bins),"<br>","end: ",end(x@anno@bins),"<br>",
                                    "probes: ",values(x@anno@bins)$probes, "<br>", "genes: ", x@anno@bins$genes)
    ggpb$x$data[[2]]$text <- paste0(x@seg$summary[[sample_n]]$chrom,"<br>","start: ",
                                    x@seg$summary[[sample_n]]$loc.start,"<br>","end: ",x@seg$summary[[sample_n]]$loc.end,"<br>",
                                    "median: ",x@seg$summary[[sample_n]]$seg.median)
    ggpb$x$data[[4]]$text <- values(x@anno@detail)$name
    ggpb%>%suppressWarnings(toWebGL())
  }
  return(ggpb)
}





