##### PROCESSING methods #####



#' CNV.fit
#' @description Normalize query sample intensities by fitting intensities to reference set using a linear regression model.
#' @param query \code{CNV.data} object of query sample (multiple samples).
#' @param ref \code{CNV.data} object of reference set.
#' @param anno \code{CNV.anno} object. Use \code{CNV.create_anno} do create.
#' @param intercept logical. Should intercept be considered? Defaults to \code{TRUE}.
#' @param ... Additional parameters (\code{CNV.fit} generic, currently not used).
#' @return \code{CNV.analysis} object.
#' @details The log2 ratio of query intensities versus a linear combination of reference set intensities that best reflects query intensities is calculated (as determined by linear regression). Every query sample in the \code{CNV.analysis} object is compared to the linear combination fo control samples individually. The annotations provided to \code{CNV.fit} are saved within the returned \code{CNV.analysis} object and used for subsequent analysis steps.
#' @examples
#' # prepare
#' library(minfiData)
#' data(MsetEx)
#' d <- CNV.load(MsetEx)
#' data(detail_regions)
#' anno <- CNV.create_anno(detail_regions = detail_regions)
#'
#' # create object
#' x <- CNV.fit(query = d, ref = d[c('GroupA_1', 'GroupA_2', 'GroupA_3')], anno)
#'
#' # modify object
#' #x <- CNV.bin(x)
#' #x <- CNV.detail(x)
#' #x <- CNV.segment(x)
#'
#' # general information
#' x
#' show(x)
#'
#' # coefficients from linear regression
#' x@@fit$coef
#'
#' @author Volker Hovestadt, Bjarne Daenekas \email{conumee@@hovestadt.bio}
#' @export
setGeneric("CNV.fit", function(query, ref, anno, ...) {
    standardGeneric("CNV.fit")
})

#' @rdname CNV.fit
setMethod("CNV.fit", signature(query = "CNV.data", ref = "CNV.data", anno = "CNV.anno"),
          function(query, ref, anno, intercept = TRUE) {
            if (ncol(query@intensity) == 0)
              stop("query intensities unavailable, run CNV.load")
            if (ncol(ref@intensity) == 0)
              stop("reference set intensities unavailable, run CNV.load")

            if (ncol(query@intensity) != 1)
              message("using multiple query samples")
            if (ncol(ref@intensity) == 1)
              warning("reference set contains only a single sample. use more samples for better results.")

            p <- unique(names(anno@probes))  # ordered by location
            if (!all(is.element(p, rownames(query@intensity))))
              stop("query intensities not given for all probes.")
            if (!all(is.element(p, rownames(ref@intensity))))
              stop("reference set intensities not given for all probes.")

            object <- new("CNV.analysis")
            object@date <- date()
            object@fit$args <- list(intercept = intercept)

            object@anno <- anno

            object@fit$coef <- data.frame(matrix(ncol = 0, nrow = ncol(ref@intensity)))
            object@fit$ratio <- data.frame(matrix(ncol = 0, nrow = length(p)))
            for (i in 1:ncol(query@intensity)) {

              message(paste(colnames(query@intensity)[i]), " (",round(i/ncol(query@intensity)*100, digits = 3), "%", ")", sep = "")
              r <- cor(query@intensity[p, ], ref@intensity[p, ])[i, ] < 0.99
              if (any(!r)) message("query sample seems to also be in the reference set. not used for fit.")
              if (intercept) {
                ref.fit <- lm(y ~ ., data = data.frame(y = log2(query@intensity[p,i]), X = log2(ref@intensity[p, r])))
              } else {
                ref.fit <- lm(y ~ . - 1, data = data.frame(y = log2(query@intensity[p,i]), X = log2(ref@intensity[p, r])))
              }
              object@fit$coef <- cbind(object@fit$coef,as.numeric(ref.fit$coefficients[-1]))

              ref.predict <- predict(ref.fit)
              ref.predict[ref.predict < 0] <- 0

              object@fit$ratio <- cbind(object@fit$ratio, log2(query@intensity[p,i]) - ref.predict[p])
            }


            colnames(object@fit$coef) <- colnames(query@intensity)
            rownames(object@fit$coef) <- colnames(ref@intensity)
            colnames(object@fit$ratio) <- colnames(query@intensity)
            rownames(object@fit$ratio) <- p

            object@fit$noise <- as.numeric()
            for (i in 1:ncol(query@intensity)) {
              object@fit$noise <- c(object@fit$noise, sqrt(mean((object@fit$ratio[-1,i] - object@fit$ratio[-nrow(object@fit$ratio),i])^2,na.rm = TRUE)))
            }

            names(object@fit$noise) <- colnames(query@intensity)
            return(object)
          })


#' CNV.bin
#' @description Combine single probe intensitiy values into predefined bins.
#' @param object \code{CNV.analysis} object.
#' @param ... Additional parameters (\code{CNV.bin} generic, currently not used).
#' @return \code{CNV.analysis} object.
#' @details The median intensity per bin and its variance are calculated. Bins are defined using \code{CNV.create_anno}. A value by which all probe and bin intensity values are shifted in subsequent analysis steps is calculated by minimizing the median absolute deviation from all bins to zero (ideally shifting the copy-neutral state to 0).
#' @examples
#' # prepare
#' library(minfiData)
#' data(MsetEx)
#' d <- CNV.load(MsetEx)
#' data(detail_regions)
#' anno <- CNV.create_anno(detail_regions = detail_regions)
#'
#' # create object
#' x <- CNV.fit(query = d['GroupB_1'], ref = d[c('GroupA_1', 'GroupA_2', 'GroupA_3')], anno)
#'
#' # modify object
#' x <- CNV.bin(x)
#' #x <- CNV.detail(x)
#' #x <- CNV.segment(x)
#'
#' # general information
#' x
#' show(x)
#'
#' # coefficients of linear regression
#' coef(x)
#'
#' # show or replace sample name
#' names(x)
#' names(x) <- 'Sample 1'
#' @author Volker Hovestadt, Bjarne Daenekas \email{conumee@@hovestadt.bio}
#' @export
setGeneric("CNV.bin", function(object, ...) {
    standardGeneric("CNV.bin")
})

#' @rdname CNV.bin
setMethod("CNV.bin", signature(object = "CNV.analysis"), function(object) {
  if (length(object@fit) == 0)
    stop("fit unavailable, run CNV.fit")

  o1 <- as.matrix(findOverlaps(query = object@anno@bins, subject = object@anno@probes))
  o2 <- data.frame(bin = names(object@anno@bins)[o1[, "queryHits"]],
                   probe = names(object@anno@probes)[o1[, "subjectHits"]], stringsAsFactors = FALSE)

  object@bin$ratio <- vector(mode = "list", length = ncol(object@fit$ratio))
  object@bin$variance <- vector(mode = "list", length = ncol(object@fit$ratio))
  object@bin$shift <- as.numeric()
  for (i in 1:ncol(object@fit$ratio)) {

    message(paste(colnames(object@fit$ratio)[i]), " (",round(i/ncol(object@fit$ratio)*100, digits = 3), "%", ")", sep = "")

    object@bin$ratio[[i]] <- sapply(split(object@fit$ratio[o2[, "probe"],i], o2[,"bin"]),
                                    median, na.rm = TRUE)[names(object@anno@bins)]

    object@bin$variance[[i]] <- sapply(split(object@fit$ratio[o2[, "probe"],i], o2[,"bin"]),
                                       var, na.rm = TRUE)[names(object@anno@bins)]

    object@bin$shift <- c(object@bin$shift, optim(0, function(s) median(abs(object@bin$ratio[[i]] -s),na.rm = TRUE),
                                                  method = "Brent", lower = -100, upper = 100)$par)
  }

  names(object@bin$shift) <- colnames(object@fit$ratio)
  names(object@bin$ratio) <- colnames(object@fit$ratio)
  names(object@bin$variance) <- colnames(object@fit$ratio)

  return(object)
})


#' CNV.detail
#' @description Combine single probe values within detail regions.
#' @param object \code{CNV.analysis} object.
#' @param ... Additional parameters (\code{CNV.detail} generic, currently not used).
#' @return \code{CNV.analysis} object.
#' @details The median intensity per detail region is calculated. Detail regions are defined using \code{CNV.create_anno(detail_bed=)}
#' @examples
#' # prepare
#' library(minfiData)
#' data(MsetEx)
#' d <- CNV.load(MsetEx)
#' data(detail_regions)
#' anno <- CNV.create_anno(detail_regions = detail_regions)
#'
#' # create object
#' x <- CNV.fit(query = d['GroupB_1'], ref = d[c('GroupA_1', 'GroupA_2', 'GroupA_3')], anno)
#'
#' # modify object
#' x <- CNV.bin(x)
#' x <- CNV.detail(x)
#' #x <- CNV.segment(x)
#'
#' # general information
#' x
#' show(x)
#'
#' # coefficients of linear regression
#' coef(x)
#'
#' # show or replace sample name
#' names(x)
#' names(x) <- 'Sample 1'
#' @author Volker Hovestadt, Bjarne Daenekas \email{conumee@@hovestadt.bio}
#' @export
setGeneric("CNV.detail", function(object, ...) {
    standardGeneric("CNV.detail")
})

#' @rdname CNV.detail
setMethod("CNV.detail", signature(object = "CNV.analysis"), function(object) {
  if (length(object@fit) == 0)
    stop("fit unavailable, run CNV.fit")
  # if(length(object@bin) == 0) stop('bin unavailable, run CNV.bin')

  if (length(object@anno@detail) == 0) {
    message("no detail regions provided, define using CNV.create_anno")
  } else {
    d1 <- as.matrix(findOverlaps(query = object@anno@detail, subject = object@anno@probes))
    d2 <- data.frame(detail = values(object@anno@detail)$name[d1[,"queryHits"]], probe = names(object@anno@probes[d1[, "subjectHits"]]),stringsAsFactors = FALSE)


    object@detail$ratio <- vector(mode = "list", length = ncol(object@fit$ratio))
    for (i in 1:ncol(object@fit$ratio)) {
      object@detail$ratio [[i]]<- sapply(split(object@fit$ratio[d2[, "probe"],i],
                                               d2[, "detail"]), median, na.rm = TRUE)[values(object@anno@detail)$name]
    }
    names(object@detail$ratio) <- colnames(object@fit$ratio)
    object@detail$probes <- table(d2[, 1])[values(object@anno@detail)$name]

  }
  return(object)
})



#' @import nullranges
NULL

#' CNV.focal
#' @description This function provides a statistical assessment for focal CNVs.
#' @param object \code{CNV.analysis} object.
#' @param conf numeric. Confidence level to calculate to determine the log2-threshold for high-level alterations. Choose between \code{0.95} and \code{0.99}. Default to \code{0.95}.
#' @param R numeric. Parameter for the \code{bootRanges} function. The number of bootstrap samples to generate. Default to \code{100}.
#' @param blockLength numeric. Parameter for the \code{bootRanges} function. The length (in basepairs) of the blocks for segmented block bootstrapping. Default to \code{500000}.
#' @param proportionLength logical. From the \code{nullranges} package: for the segmented block bootstrap, whether to use scaled block lengths, (scaling by the proportion of the segmentation state out of the total genome length)
#' @param minoverlap integer. The function determines the bins that overlap with the genes of interest. Which minimum number of basepairs should be considered for an overlap? Defaul to \code{1000L}.
#' @param ... Additional parameters (\code{CNV.detailplot} generic, currently not used).
#' @return A \code{CNV.analysis} object with significantly altered bins and genes from the Cancer Gene Census (curated by the Sanger Institute).
#' @details This function should facilitate the detection of diagnostically relevant CNVs that affect single genes. Segmented Block Bootstrapping is performed to define sample-specific Log2ratio thresholds for high-level alterations.
#' @examples
#'
#' x <- CNV.focal(x, conf = 0.95, R = 100, blockLength = 500000, minoverlap = 10000L)
#' x@@detail$del.bins  #bins that are part of high-level deletions.
#' x@@detail$amp.bins  #bins that are part of high-level amplifications.
#' x@@detail$sig.genes #significantly altered genes from the Cancer Gene Census and the predefined detail_regions.
#'
#' @author Bjarne Daenekas \email{conumee@@hovestadt.bio}
#' @export
setGeneric("CNV.focal", function(object, ...) {
  standardGeneric("CNV.focal")
})

#' @rdname CNV.focal
setMethod("CNV.focal", signature(object = "CNV.analysis"), function(object, conf = 0.95, R = 100, blockLength = 500000, proportionLength = TRUE, ...){
  if(ncol(object@anno@genome) == 2) {
    stop("CNV.focal is not compatible with mouse arrays.")
  }

  data("consensus_cancer_genes_hg19")
  data("genes")

  amp_bins <- vector(mode='list', length=ncol(object@fit$ratio))
  del_bins <- vector(mode='list', length=ncol(object@fit$ratio))
  amp_genes <- vector(mode='list', length=ncol(object@fit$ratio))
  del_genes <- vector(mode='list', length=ncol(object@fit$ratio))
  for(i in 1:ncol(object@fit$ratio)){

    message(paste(colnames(object@fit$ratio)[i]), " (",round(i/ncol(object@fit$ratio)*100, digits = 3), "%", ")", sep = "")

    segs <- object@seg$summary[[i]]
    segs <- GRanges(seqnames = segs$chrom, IRanges(start = segs$loc.start, end = segs$loc.end), seqinfo = Seqinfo(genome = "hg19"))
    seqlevels(segs) <- object@anno@genome$chr
    segs$seg.median <- object@seg$summary[[i]]$seg.median - object@bin$shift[i]
    segs$num.markers <- object@seg$summary[[i]]$num.mark
    segs <- sort(segs)

    seq <- rep(segs$seg.median, segs$num.markers)
    km <- kmeans(seq, centers=3, nstart=25)

    bins.log2 <- object@bin$ratio[[i]] - object@bin$shift[i]
    bins <- object@anno@bins[names(bins.log2)]
    bins$weight <- 1/object@bin$variance[[i]][names(bins.log2)]
    bins$log2 <- as.numeric(bins.log2)
    bins$state <- km$cluster
    seqinfo(bins) <- Seqinfo(genome = "hg19")


    seg <- lapply(seq_len(3), function(s){ # Combine nearby regions within same states
      x <- reduce(bins[bins$state == s])
      mcols(x)$state <- s
      x
    })

    seg <- c(seg[[1]], seg[[2]], seg[[3]])

    # if(min(km$size)>=20){
    #   boots <- bootRanges(bins, blockLength = blockLength, R = R, seg = seg, exclude = object@anno@exclude, proportionLength = FALSE)
    # } else{ #bootRanges can't deal with clusters that are too small (<20), proportionLength = TRUE then
    #   message("One kmeans cluster is very small (<20), setting proportionLength=TRUE.")
    #   boots <- bootRanges(bins, blockLength = blockLength, R = R, seg = seg, exclude = object@anno@exclude, proportionLength = TRUE)
    # }

    boots <- bootRanges(bins, blockLength = blockLength, R = R, seg = seg, exclude = object@anno@exclude, proportionLength = proportionLength)

    t <- round(length(boots) * (1-conf), digits = 0)/2
    del_t <- round(sort(boots$log2)[t], digits = 3)
    amp_t <- round(sort(boots$log2, decreasing =  TRUE)[t], digits = 3)

    bins.total <- object@bin$ratio[[i]] - object@bin$shift[i]
    bins.amp <- sort(bins.total[bins.total >= amp_t], decreasing = TRUE)
    bins.del <- sort(bins.total[bins.total <= del_t])
    sig.bins <- sort(c(abs(bins.amp), abs(bins.del)), decreasing = TRUE)

    amp_bins[[i]] <- bins.amp
    del_bins[[i]] <- bins.del

    genes.n <- setdiff(cancer_genes$SYMBOL, genes$SYMBOL)
    genes.c <- sort(c(genes, cancer_genes[genes.n]))

    if(length(object@anno@detail) >= 1){
      dif <- setdiff(genes.c$SYMBOL, object@anno@detail$name)
      details <- object@anno@detail
      mcols(details) <- data.frame(SYMBOL = object@anno@detail$name)
      cgenes <- c(details, genes.c[dif])
    }

    if(length(object@anno@detail) == 0){
      cgenes <- genes.c
    }

    d1 <- as.matrix(findOverlaps(query = cgenes, subject = anno@probes))
    d2 <- data.frame(gene = values(cgenes)$SYMBOL[d1[,"queryHits"]], probe = names(anno@probes[d1[, "subjectHits"]]),stringsAsFactors = FALSE)
    cgenes.ratio <- sapply(split(object@fit$ratio[d2[, "probe"],i], d2[, "gene"]), median, na.rm = TRUE)[values(cgenes)$SYMBOL]
    cgenes.ratio <- cgenes.ratio - object@bin$shift[i]

    amp.genes <- names(which(cgenes.ratio >= amp_t))
    del.genes <- names(which(cgenes.ratio <= del_t))


    amp_genes[[i]] <- amp.genes
    del_genes[[i]] <- del.genes
    }

  names(del_bins) <- colnames(object@fit$ratio)
  names(amp_bins) <- colnames(object@fit$ratio)
  names(amp_genes) <- colnames(object@fit$ratio)
  names(del_genes) <- colnames(object@fit$ratio)


  object@detail$del.bins <- del_bins
  object@detail$amp.bins <- amp_bins
  object@detail$amp.genes <- amp_genes
  object@detail$del.genes <- del_genes


  return(object)

})


#' @import DNAcopy
NULL

#' CNV.segment
#' @description Segment bin values (wrapper of \code{DNAcopy} package). Each bin is assigned to a weight that is inversely proprotional to its variance.
#' @param object \code{CNV.analysis} object.
#' @param alpha See details. Defaults to 0.001.
#' @param nperm See details. Defaults to 50000.
#' @param min.width See details. Defaults to 5.
#' @param undo.splits See details. Defaults to 'sdundo'.
#' @param undo.SD See details. Defaults to 2.2.
#' @param verbose See details. Defaults to 0.
#' @param ... Additional parameters supplied to the \code{segment} method of the \code{DNAcopy} package.
#' @return \code{CNV.analysis} object.
#' @details This method is a wrapper of the CNA, segment, segments.summary and segments.p methods of the DNAcopy package. Please refer to the respective man pages for more detailed information. The default parameters of \code{CNV.segment} override some of the default parameters of segment and are optimized for 450k data CNV analysis.
#' @examples
#' # prepare
#' library(minfiData)
#' data(MsetEx)
#' d <- CNV.load(MsetEx)
#' data(detail_regions)
#' anno <- CNV.create_anno(detail_regions = detail_regions)
#'
#' # create object
#' x <- CNV.fit(query = d['GroupB_1'], ref = d[c('GroupA_1', 'GroupA_2', 'GroupA_3')], anno)
#'
#' # modify object
#' x <- CNV.bin(x)
#' x <- CNV.detail(x)
#' x <- CNV.segment(x)
#'
#' # general information
#' x
#' show(x)
#'
#' # coefficients of linear regression
#' coef(x)
#'
#' # show or replace sample name
#' names(x)
#' names(x) <- 'Sample 1'
#' @author Volker Hovestadt, Bjarne Daenekas \email{conumee@@hovestadt.bio}
#' @export
setGeneric("CNV.segment", function(object, ...) {
    standardGeneric("CNV.segment")
})

#' @rdname CNV.segment
setMethod("CNV.segment", signature(object = "CNV.analysis"), function(object,
                                                                      alpha = 0.001, nperm = 50000, min.width = 5, undo.splits = "sdundo",
                                                                      undo.SD = 2.2, verbose = 0, ...) {
  if(length(object@fit) == 0){
    stop('fit unavailable, run CNV.fit')
  }
  if (length(object@bin) == 0){
    stop("bin unavailable, run CNV.bin")
  }

  a1 <- formals()
  a2 <- as.list(match.call())[-1]
  object@seg$args <- as.list(sapply(setdiff(unique(names(c(a1, a2))),
                                            c("object", "verbose")), function(an) if (is.element(an, names(a2)))
                                              a2[[an]] else a1[[an]], simplify = FALSE))

  object@seg$summary <- vector(mode = "list", length = ncol(object@fit$ratio))
  object@seg$p <- vector(mode = "list", length = ncol(object@fit$ratio))

  for (i in 1:ncol(object@fit$ratio)) {

    message(paste(colnames(object@fit$ratio)[i]), " (",round(i/ncol(object@fit$ratio)*100, digits = 3), "%", ")", sep = "")

    x1 <- DNAcopy::CNA(genomdat = object@bin$ratio[[i]][names(object@anno@bins)],
                       chrom = as.vector(seqnames(object@anno@bins)), maploc = values(object@anno@bins)$midpoint,
                       data.type = "logratio", sampleid = "sampleid")

    x2 <- DNAcopy::segment(x = x1, weights = 1/object@bin$variance[[i]][names(object@anno@bins)], verbose = verbose, min.width = min.width,
                           nperm = nperm, alpha = alpha, undo.splits = undo.splits, undo.SD = undo.SD,
                           ...)

    object@seg$summary[[i]] <- DNAcopy::segments.summary(x2)
    object@seg$summary[[i]]$chrom <- as.vector(object@seg$summary[[i]]$chrom)
    object@seg$summary[[i]]$ID <- colnames(object@fit$ratio)[i]
    object@seg$p[[i]] <- DNAcopy::segments.p(x2)
    object@seg$p[[i]]$chrom <- as.vector(object@seg$p[[i]]$chrom)
  }
  names(object@seg$summary) <- colnames(object@fit$ratio)
  names(object@seg$p) <- colnames(object@fit$ratio)
  return(object)
})

