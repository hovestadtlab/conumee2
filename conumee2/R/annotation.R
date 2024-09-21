##### ANNOTATION methods #####

#' @import minfi
#' @import GenomicRanges
#' @import IRanges
#' @import GenomeInfoDb
#' @importFrom rtracklayer import
NULL

#' CNV.create_anno
#' @description Create annotations for CNV analysis.
#' @param bin_minprobes numeric. Minimum number of probes per bin. Bins are iteratively merged with neighboring bin until minimum number is reached.
#' @param bin_minsize numeric. Minimum size of a bin.
#' @param bin_maxsize numeric. Maximum size of a bin. Merged bins that are larger are filtered out.
#' @param array_type character. One of \code{450k}, \code{EPIC}, \code{EPICv2}, \code{mouse}. When analyzing data from multiple array types, choose between \code{overlap.1} for 450k and EPIC, \code{overlap.2} for EPIC and EPICv2 or \code{overlap.3} for 450k, EPIC, EPICv2.
#' @param genome character. This parameter can be set to \code{hg38} when working with EPICv2 data. It will be ignored if \code{array_type = mouse}.
#' @param features vector. Per default, all unique CpGs on the array are used for the analysis. Please provide a vector with Probe IDs to create a customized set of probes.
#' @param exclude_regions GRanges object or path to bed file containing genomic regions to be excluded.
#' @param detail_regions GRanges object or path to bed file containing genomic regions to be examined in detail.
#' @param chrXY logical. Should chromosome X and Y be included in the analysis?
#' @return \code{CNV.anno} object.
#' @details This function collects all annotations required for CNV analysis using Illumina 450k, EPIC or Mouse arrays. The output \code{CNV.anno} object is not editable. Rerun \code{CNV.create_anno} to change parameters.
#' @examples
#' # create annotation object
#' anno <- CNV.create_anno(array_type = "450k", detail_regions = detail_regions, exclude_regions = exclude_regions)
#' anno
#' @author Volker Hovestadt \email{conumee@@hovestadt.bio}, Bjarne Daenekas
#' @export
CNV.create_anno <- function(bin_minprobes = 15, bin_minsize = 50000, bin_maxsize = 5e+06,
                            array_type = "450k", genome = "hg19", features = "all", exclude_regions = NULL, detail_regions = NULL, chrXY = FALSE) {
  object <- new("CNV.anno")
  object@date <- date()

  a1 <- formals()
  a2 <- as.list(match.call())[-1]
  object@args <- as.list(sapply(unique(names(c(a1, a2))), function(an) if (is.element(an, names(a2)))
    a2[[an]] else a1[[an]], simplify = FALSE))

  if (is.null(array_type)) {
    array_type <- "450k"
  }
  if (!all(is.element(array_type, c("450k", "EPIC", "EPICv2", "mouse")))) {
    stop("array_type must be one/multiple of 450k, EPIC, EPICv2, or mouse")
  }
  if (is.element("mouse", array_type) & length(array_type) > 1) {
    stop("cannot mix mouse and human arrays")
  }

  if(all(is.element(array_type, c("450k", "EPIC", "EPICv2")))) {   # human annotations
    if (chrXY) {
      object@genome <- data.frame(chr = paste("chr", c(1:22, "X", "Y"),
                                              sep = ""), stringsAsFactors = FALSE)
    } else {
      object@genome <- data.frame(chr = paste("chr", 1:22, sep = ""),
                                  stringsAsFactors = FALSE)
    }
    rownames(object@genome) <- object@genome$chr

    if (genome == "hg19") {
      message("using hg19 genome annotations from UCSC")
      tbl <- tbl_ucsc
    }

    if (genome == "hg38") {
      if (!is.element("EPICv2", array_type)) {
        stop("only the EPICv2 array supports the hg38 annotation")
      }
      message("using hg38 genome annotations from UCSC")
      tbl <- tbl_ucsc_hg38
    }

    tbl.chromInfo <- tbl$chromInfo[match(object@genome$chr, tbl$chromInfo$chrom), "size"]
    object@genome$size <- tbl.chromInfo
    tbl.gap <- tbl$gap[is.element(tbl$gap$chrom, object@genome$chr),]
    object@gap <- sort(GRanges(as.vector(tbl.gap$chrom), IRanges(tbl.gap$chromStart + 1, tbl.gap$chromEnd),
                               seqinfo = Seqinfo(object@genome$chr, object@genome$size)))

    if(genome == "hg19"){
      tbl.cytoBand <- tbl$cytoBand[is.element(tbl$cytoBand$chrom, object@genome$chr), ]

      # find the gap that is overlapping with the end of the last p-band, use
      # center of that gap for indicating centromeres in the genome plots

      pq <- sapply(split(tbl.cytoBand$chromEnd[grepl("p", tbl.cytoBand$name)],
                         as.vector(tbl.cytoBand$chrom[grepl("p", tbl.cytoBand$name)])), max)
      object@genome$pq <- start(resize(subsetByOverlaps(object@gap, GRanges(names(pq),
                                                                            IRanges(pq, pq))), 1, fix = "center"))
    } else if(genome == "hg38"){
      cent.min <- sapply(split(tbl$centromeres$chromStart[which(is.element(tbl$centromeres$chrom, object@genome$chr))],
                               tbl$centromeres$chrom[which(is.element(tbl$centromeres$chrom, object@genome$chr))]), min)
      cent.max <- sapply(split(tbl$centromeres$chromStart[which(is.element(tbl$centromeres$chrom, object@genome$chr))],
                               tbl$centromeres$chrom[which(is.element(tbl$centromeres$chrom, object@genome$chr))]), max)
      object@genome$pq <- start(sort(resize(GRanges(as.vector(names(cent.min)), IRanges(start = cent.min, stop = cent.max),
                                                    seqinfo = Seqinfo(object@genome$chr, object@genome$size)), 1, fix = "center")))
    }

    if (is.element("450k", array_type)) {
      message("loading 450k annotations")
    }

    if (is.element("EPIC", array_type)) {
      message("loading EPIC annotations")
    }

    if (is.element("EPICv2", array_type)) {
      if(genome == "hg19"){
        message("loading EPICv2 annotations")
      } else if(genome == "hg38"){
        message("loading EPICv2 annotations (hg38)")
        probesEPICv2 <- probesEPICv2.hg38
      }
    }

    if (all(array_type=="450k")) {
      probes <- probes450k
    } else if (all(array_type %in% "EPIC")) {
      probes <- probesEPIC
    } else if (all(array_type %in% "EPICv2")) {
      probes <- probesEPICv2
      probes$EPICv1_Loci <- probes$Methyl450_Loci <- NULL
    } else if (all(array_type %in% c("450k", "EPIC"))) {
      probes <- probes450k[intersect(names(probes450k), names(probesEPIC))]
    } else if (all(array_type %in% c("450k", "EPICv2"))) {
      probes <- probesEPICv2[probesEPICv2$Methyl450_Loci != ""]
      probes$EPICv1_Loci <- probes$Methyl450_Loci <- NULL
    } else if (all(array_type %in% c("EPIC", "EPICv2"))) {
      probes <- probesEPICv2[probesEPICv2$EPICv1_Loci != ""]
      probes$EPICv1_Loci <- probes$Methyl450_Loci <- NULL
    } else if (all(array_type %in% c("450k", "EPIC", "EPICv2"))) {
      probes <- probesEPICv2[probesEPICv2$Methyl450_Loci != "" & probesEPICv2$EPICv1_Loci != ""]
      probes$EPICv1_Loci <- probes$Methyl450_Loci <- NULL
    }

    if(features == "all"){
      object@probes <- sort(probes)
      message(" - ", length(object@probes), " probes used")
    }else{
      object@probes <- sort(probes[features])
      message(" - ", length(object@probes), " probes used")
    }

    object@probes <- sort(probes)
    message(" - ", length(object@probes), " probes used")

    if (!is.null(exclude_regions)) {
      message("importing regions to exclude from analysis")
      if (class(exclude_regions) == "GRanges") {
        object@exclude <- GRanges(as.vector(seqnames(exclude_regions)),
                                  ranges(exclude_regions), seqinfo = Seqinfo(object@genome$chr,
                                                                             object@genome$size))
        values(object@exclude) <- values(exclude_regions)
        object@exclude <- sort(object@exclude)
      } else {
        object@exclude <- sort(rtracklayer::import(exclude_regions,
                                                   seqinfo = Seqinfo(object@genome$chr, object@genome$size)))
      }
    } else {
      object@exclude <- GRanges(seqinfo = Seqinfo(object@genome$chr,
                                                  object@genome$size))
    }

    if (!is.null(detail_regions)) {
      message("importing regions for detailed analysis")
      if (class(detail_regions) == "GRanges") {
        object@detail <- GRanges(as.vector(seqnames(detail_regions)),
                                 ranges(detail_regions), seqinfo = Seqinfo(object@genome$chr,
                                                                           object@genome$size))
        if (any(grepl("name", names(values(detail_regions))))) {
          values(object@detail)$name <- values(detail_regions)[[grep("name",
                                                                     names(values(detail_regions)))[1]]]
        }
        if (any(grepl("IRanges", sapply(values(detail_regions), class)))) {
          values(object@detail)$thick <- values(detail_regions)[[grep("IRanges",
                                                                      sapply(values(detail_regions), class))[1]]]
        }
        object@detail <- sort(object@detail)
      } else {
        object@detail <- sort(rtracklayer::import(detail_regions, seqinfo = Seqinfo(object@genome$chr,
                                                                                    object@genome$size)))
      }
      if (!is.element("name", names(values(object@detail)))) {
        stop("detailed region bed file must contain name column.")
      }
      if (!all(table(values(object@detail)$name) == 1)) {
        stop("detailed region names must be unique.")
      }
    } else {
      object@detail <- GRanges(seqinfo = Seqinfo(object@genome$chr, object@genome$size))
    }
    if (!is.element("thick", names(values(object@detail)))) {
      values(object@detail)$thick <- resize(ranges(object@detail), fix = "center",
                                            1e+06)
    }

    message("creating bins")
    anno.tile <- CNV.create_bins(genome.anno = object@genome, bin_minsize = bin_minsize,
                                 genome.gap = object@gap, genome.exclude = object@exclude)
    message(" - ", length(anno.tile), " bins created")

    message("merging bins")
    object@bins <- CNV.merge_bins(genome.anno = object@genome, genome.tile = anno.tile,
                                  bin_minprobes = bin_minprobes, genome.probes = object@probes, bin_maxsize = bin_maxsize)
    message(" - ", length(object@bins), " bins remaining")

    message("getting the gene annotations for each bin")

    o <- findOverlaps(object@probes, object@bins)
    bin_genes <- sapply(lapply(sapply(split(object@probes$genes[queryHits(o)], names(object@bins)[subjectHits(o)]),
                                      function(x) na.omit(unlist(strsplit(x,split = ";")))), unique), paste, collapse = ";")

    object@bins$genes <- bin_genes

    return(object)

  } else if (array_type == "mouse") { # mouse annotations

    if (chrXY) {
      object@genome <- data.frame(chr = paste("chr", c(1:19, "X", "Y"), sep = ""), stringsAsFactors = FALSE)
    } else {
      object@genome <- data.frame(chr = paste("chr", 1:19, sep = ""), stringsAsFactors = FALSE)
    }

    rownames(object@genome) <- object@genome$chr
    message("using mm10 genome annotations from UCSC")

    object@genome$size <- mouse_annotation[[1]]$size[1:nrow(object@genome)]
    tbl.gap <- mouse_annotation[[2]][is.element(mouse_annotation[[2]]$chrom, object@genome$chr),]
    object@gap <- sort(GRanges(as.vector(tbl.gap$chrom), IRanges(tbl.gap$chromStart + 1,
                                                                 tbl.gap$chromEnd), seqinfo = Seqinfo(object@genome$chr, object@genome$size)))

    ind.chr <- is.element(mouse_annotation[[3]]$CHR, unlist(lapply(strsplit(object@genome$chr, "r"), function(x)x[2])))

    mouse_probes <- GRanges(as.vector(paste("chr", mouse_annotation[[3]]$CHR[ind.chr], sep = "")),
                            IRanges(start = mouse_annotation[[3]]$MAPINFO[ind.chr],
                                    end = mouse_annotation[[3]]$MAPINFO[ind.chr]), seqinfo = Seqinfo(object@genome$chr, object@genome$size, genome = "mm10"))

    names(mouse_probes) <- mouse_annotation[[3]]$Name[ind.chr]
    mouse_probes$IlmnID <- mouse_annotation[[3]]$IlmnID[ind.chr]
    mouse_probes$genes <- mouse_annotation[[3]]$genes[ind.chr]

    # CpG probes only
    mouse_probes <- sort(mouse_probes[substr(names(mouse_probes),1, 2) == "cg"])
    object@probes <- mouse_probes

    message(" - ", length(object@probes), " probes used")

    if (!is.null(exclude_regions)) {
      message("importing regions to exclude from analysis")
      if (class(exclude_regions) == "GRanges") {
        object@exclude <- GRanges(as.vector(seqnames(exclude_regions)),
                                  ranges(exclude_regions), seqinfo = Seqinfo(object@genome$chr,
                                                                             object@genome$size))
        values(object@exclude) <- values(exclude_regions)
        object@exclude <- sort(object@exclude)
      } else {
        object@exclude <- sort(rtracklayer::import(exclude_regions,
                                                   seqinfo = Seqinfo(object@genome$chr, object@genome$size)))
      }
    } else {
      object@exclude <- GRanges(seqinfo = Seqinfo(object@genome$chr,
                                                  object@genome$size))
    }

    if (!is.null(detail_regions)) {
      message("importing regions for detailed analysis")
      if (class(detail_regions) == "GRanges") {
        object@detail <- GRanges(as.vector(seqnames(detail_regions)),
                                 ranges(detail_regions), seqinfo = Seqinfo(object@genome$chr,
                                                                           object@genome$size))
        if (any(grepl("name", names(values(detail_regions))))) {
          values(object@detail)$name <- values(detail_regions)[[grep("name",
                                                                     names(values(detail_regions)))[1]]]
        }
        if (any(grepl("IRanges", sapply(values(detail_regions), class)))) {
          values(object@detail)$thick <- values(detail_regions)[[grep("IRanges",
                                                                      sapply(values(detail_regions), class))[1]]]
        }
        object@detail <- sort(object@detail)
      } else {
        object@detail <- sort(rtracklayer::import(detail_regions, seqinfo = Seqinfo(object@genome$chr,
                                                                                    object@genome$size)))
      }
      if (!is.element("name", names(values(object@detail)))) {
        stop("detailed region bed file must contain name column.")
      }
      if (!all(table(values(object@detail)$name) == 1)) {
        stop("detailed region names must be unique.")
      }
    } else {
      object@detail <- GRanges(seqinfo = Seqinfo(object@genome$chr, object@genome$size))
    }
    if (!is.element("thick", names(values(object@detail)))) {
      values(object@detail)$thick <- resize(ranges(object@detail), fix = "center",
                                            1e+06)
    }

    message("creating bins")
    anno.tile <- CNV.create_bins(genome.anno = object@genome, bin_minsize = bin_minsize,
                                 genome.gap = object@gap, genome.exclude = object@exclude)
    message(" - ", length(anno.tile), " bins created")

    message("merging bins")
    object@bins <- CNV.merge_bins(genome.anno = object@genome, genome.tile = anno.tile,
                                  bin_minprobes = bin_minprobes, genome.probes = object@probes, bin_maxsize = bin_maxsize)
    message(" - ", length(object@bins), " bins remaining")

    message("getting the gene annotations for each bin")

    o <- findOverlaps(object@probes, object@bins)
    bin_genes <- sapply(lapply(sapply(split(object@probes$genes[queryHits(o)], names(object@bins)[subjectHits(o)]),
                                      function(x) na.omit(unlist(strsplit(x,split = ";")))), unique), paste, collapse = ";")

    object@bins$genes <- bin_genes

    return(object)
  }
}


#' CNV.create_bins
#' @description Split genome into bins of defined size.
#' @param genome.anno foo
#' @param bin_minsize foo
#' @param genome.gap foo
#' @param genome.exclude foo
#' @return \code{GRanges} object.
CNV.create_bins <- function(genome.anno, bin_minsize = 50000, genome.gap, genome.exclude) {
  genome.tile <- sort(tileGenome(Seqinfo(genome.anno$chr, genome.anno$size),
                                 tilewidth = bin_minsize, cut.last.tile.in.chrom = TRUE))
  # setdiff for gaps (on every second window to avoid merging)
  genome.tile <- sort(c(setdiff(genome.tile[seq(1, length(genome.tile), 2)],
                                genome.gap), setdiff(genome.tile[seq(2, length(genome.tile), 2)], genome.gap)))
  # setdiff for exluded regions
  genome.tile <- sort(c(setdiff(genome.tile[seq(1, length(genome.tile), 2)],
                                genome.exclude), setdiff(genome.tile[seq(2, length(genome.tile), 2)],
                                                         genome.exclude)))
  return(genome.tile)
}

#' CNV.merge_bins
#' @description Merge bins containing less than the defined number probes with neighboring bin containing fewer probes.
#' @param genome.anno foo
#' @param genome.tile foo
#' @param bin_minprobes foo
#' @param genome.probes foo
#' @param bin_maxsize foo
#' @param verbose foo
#' @return \code{GRanges} object.
CNV.merge_bins <- function(genome.anno, genome.tile, bin_minprobes = 15, genome.probes,
                           bin_maxsize = 5e+06, verbose = FALSE) {
  values(genome.tile)$probes <- countOverlaps(genome.tile, genome.probes)
  genome.tile.df <- as.data.frame(genome.tile)[, c("seqnames", "start", "end",
                                                   "probes")]
  genome.tile.df$seqnames <- as.vector(genome.tile.df$seqnames)  # not factor

  genome.tile.df.bin <- do.call(rbind, lapply(split(genome.tile.df, genome.tile.df$seqnames),
                                              function(mdf) {
                                                while (min(mdf$probes) < bin_minprobes) {
                                                  mw <- which(mdf$probes == min(mdf$probes))[1]
                                                  mwn <- NA
                                                  mwns <- Inf
                                                  # left
                                                  if (is.element(mdf[mw, "start"] - 1, mdf[, "end"])) {
                                                    mwn <- mw - 1
                                                    mwns <- mdf[mw - 1, "probes"]
                                                    # }
                                                  }
                                                  # right
                                                  if (is.element(mdf[mw, "end"] + 1, mdf[, "start"])) {
                                                    if (mdf[mw + 1, "probes"] < mwns) {
                                                      mwn <- mw + 1
                                                      mwns <- mdf[mw + 1, "probes"]
                                                    }
                                                  }
                                                  if (is.na(mwn)) {
                                                    if (verbose)
                                                      message(paste(mdf[mw, 1:3], collapse = "-"), " has only ",
                                                              mdf[mw, "probes"], " probes and cannot be merged - remove")
                                                    mdf <- mdf[-mw, ]
                                                  } else {
                                                    # merge
                                                    mdf[mwn, "start"] <- min(mdf[c(mwn, mw), "start"])
                                                    mdf[mwn, "end"] <- max(mdf[c(mwn, mw), "end"])
                                                    mdf[mwn, "probes"] <- sum(mdf[c(mwn, mw), "probes"])
                                                    mdf <- mdf[-mw, ]
                                                  }
                                                }
                                                return(mdf)
                                              }))

  genome.tile.bin <- sort(GRanges(genome.tile.df.bin$seqnames, IRanges(genome.tile.df.bin$start,
                                                                       genome.tile.df.bin$end), seqinfo = seqinfo(genome.tile)))
  genome.tile.bin <- genome.tile.bin[width(genome.tile.bin) <= bin_maxsize]

  values(genome.tile.bin)$probes <- countOverlaps(genome.tile.bin, genome.probes)
  values(genome.tile.bin)$midpoint <- as.integer(start(genome.tile.bin) +
                                                   (end(genome.tile.bin) - start(genome.tile.bin))/2)

  names(genome.tile.bin) <- paste(as.vector(seqnames(genome.tile.bin)), formatC(unlist(lapply(table(seqnames(genome.tile.bin)),
                                                                                              seq)), width = nchar(max(table(seqnames(genome.tile.bin)))), format = "d",
                                                                                flag = "0"), sep = "-")
  return(genome.tile.bin)
}
