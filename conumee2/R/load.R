##### LOADING methods #####
#' @import illuminaio
#' @import RnBeads
NULL

#' CNV.load
#' @description Prepare combined intensities from various input objects.
#' @param input Object of MethylSet class (minfi package), data.frame class, matrix class or numeric class.
#' @param names Vector specifying sample names. If not supplied, colnames are used. For MethylSet input, the first column of pData(input) matching 'name' (grep) is used.
#' @param ... Additional parameters (\code{CNV.load} generic, currently not used).
#' @return \code{CNV.data} object.
#' @details This method gathers combined intensities of the Methylated and Unmethylated signals for all supplied probes. Probe IDs must be supplied as row names or in a seperate column named `ID_REF` or `TargetID`.
#' If column names match 'intensity', only those columns are used. Else, if column names match 'signal' or 'methylated', only those columns are used. Otherwise, all columns are used.
#' @examples
#' library(minfiData)
#' d <- CNV.load(MsetEx)
#' d
#' @author Volker Hovestadt \email{conumee@@hovestadt.bio}
#' @export
setGeneric("CNV.load", function(input, ...) {
    standardGeneric("CNV.load")
})

#' @rdname CNV.load
setMethod("CNV.load", signature(input = "MethylSet"), function(input, names = NULL) {
    object <- new("CNV.data")

    object@intensity <- as.data.frame(minfi::getMeth(input) + minfi::getUnmeth(input))

    input.names <- grep("Name", setdiff(colnames(minfi::pData(input)),
        c("Basename", "filenames")), ignore.case = TRUE)
    if (length(input.names) > 0)
        names(object) <- minfi::pData(input)[, grep("name", setdiff(colnames(minfi::pData(input)),
            c("Basename", "filenames")), ignore.case = TRUE)[1]]
    if (!is.null(names))
        names(object) <- names

    object <- CNV.check(object)

    return(object)
})

#' @rdname CNV.load
setMethod("CNV.load", signature(input = "data.frame"), function(input,
    names = NULL) {
    object <- new("CNV.data")
    object@date <- date()

    if (any(grepl("TargetID", colnames(input))))
        rownames(input) <- input[, "TargetID"]
    if (any(grepl("ID_REF", colnames(input))))
        rownames(input) <- input[, "ID_REF"]
    if (is.null(rownames(input)))
        stop("intensities not given for all probes.")

    if (any(grepl("intensity", colnames(input), ignore.case = TRUE))) {
        input.i <- grep("intensity", colnames(input), ignore.case = TRUE)
        input.n <- sapply(strsplit(colnames(input), "\\.", colnames(input))[input.i],
            function(x) paste(x[!grepl("intensity", x, ignore.case = TRUE)],
                collapse = "."))
        object@intensity <- as.data.frame(input[, input.i])
        colnames(object@intensity) <- make.names(input.n, unique = TRUE)
    } else if (any(grepl("signal", colnames(input), ignore.case = TRUE))) {
        input.i <- grep("signal", colnames(input), ignore.case = TRUE)
        input.n <- sapply(strsplit(colnames(input), "\\.", colnames(input))[input.i],
            function(x) paste(x[!grepl("signal|methylated", x, ignore.case = TRUE)],
                collapse = "."))
        if (!all(input.n[seq(1, length(input.i), 2)] == input.n[seq(2,
            length(input.i), 2)]))
            stop("names of both signal columns do not match.")
        object@intensity <- as.data.frame(input[, input.i[seq(1, length(input.i),
            2)]] + input[, input.i[seq(2, length(input.i), 2)]])
        colnames(object@intensity) <- make.names(input.n[seq(1, length(input.i),
            2)], unique = TRUE)
    } else {
        object@intensity <- as.data.frame(input)
    }
    if (!is.null(names))
        names(object) <- names

    object <- CNV.check(object)

    return(object)
})

#' @rdname CNV.load
setMethod("CNV.load", signature(input = "matrix"), function(input, names = NULL,...) {
    CNV.load(as.data.frame(input), names)
})

#' @rdname CNV.load
setMethod("CNV.load", signature(input = "numeric"), function(input, names = NULL) {
    object <- new("CNV.data")

    if (is.null(names(input)))
        stop("intensities not given for all probes.")
    object@intensity <- data.frame(sampleid = input)
    if (!is.null(names))
        names(object) <- names

    object <- CNV.check(object)

    return(object)
})

#' @rdname CNV.load
setMethod("CNV.load", signature(input = "RnBeadRawSet"), function(input, names = NULL) {
  object <- new("CNV.data")
  object@date <- date()

  M <- input@M
  U <- input@U

  combined_probes <- round(M + U, digits = 0)
  rownames(combined_probes) <- substr(rownames(input@sites),1, 10)
  colnames(combined_probes) <- input@pheno$Sample_Name

  object@intensity <- as.data.frame(combined_probes)

  object <- CNV.check(object)

  return(object)
})


#' CNV.check
#' @description Check intensity values.
#' @param object \code{CNV.data} object.
#' @return \code{CNV.data} object.
#' @details This method checks if intensities are positive and not NA. If not, they are set to 1. Warnings are given if intensities are abnormally high or low (> 50000 or < 5000, respectively).
#' @author Volker Hovestadt \email{conumee@@hovestadt.bio}
setGeneric("CNV.check", function(object) {
    standardGeneric("CNV.check")
})

#' @rdname CNV.check
setMethod("CNV.check", signature(object = "CNV.data"), function(object) {
    if (any(is.na(object@intensity))) {
        warning("some intensities are NA, now set to 1.")
        object@intensity[is.na(object@intensity)] <- 1
    }
    if (any(object@intensity < 0))
        warning("some intensities are smaller than 0, now set to 1.")
    object@intensity[object@intensity < 1] <- 1
    if (any(colMeans(object@intensity) < 5000))
        warning("intensities are abnormally low (< 5000).")
    if (any(colMeans(object@intensity) > 50000))
        warning("intensities are abnormally high (> 50000).")

    return(object)
})


#' read.450k.url
#' @description Read IDAT files from the web.
#' @param url URL of the directory in which the IDAT files are located.
#' @param idat Vector of IDAT names. \code{url} and \code{idat} default to the TCGA example described in the vignette.
#' @return \code{RGChannelSet} object.
#' @details This method downloads the provided list of IDAT files to a temporary folder (using the \code{RCurl} package). It then uses the `read.450k.exp` method of the `minfi` package.
#' @examples
#' RGsetTCGA <- read.450k.url()
#' @author Volker Hovestadt \email{conumee@@hovestadt.bio}
#' @export
read.450k.url <- function(url = NULL, idat = NULL) {
    if (is.null(url))
        url <- "https://github.com/hovestadt/conumeeData/raw/master/"
    if (is.null(idat))
        idat <- c("6042324037_R05C02", "6042324037_R06C01")
    tmp <- paste0(tempdir(), .Platform$file.sep)
    if (!grepl("/$", url))
        url <- paste0(url, "/")
    if (!any(grepl("_Grn.idat", idat)))
        idat <- unlist(lapply(idat, paste0, c("_Grn.idat", "_Red.idat")))
    for (i in idat) .curl(url = paste0(url, i), file = paste0(tmp, i))
    idatRG <- read.metharray.exp(base = tmp)
    for (i in idat) file.remove(paste0(tmp, i))
    return(idatRG)
}

.curl <- function(url, file, verbose = TRUE) {
    if (.Platform$OS.type == "unix") {
        if (!RCurl::url.exists(url))
            stop("url does not exist.")
    } else {
        if (!RCurl::url.exists(url, .opts = list(ssl.verifypeer = FALSE)))
            stop("url does not exist.")
    }
    if (verbose)
        message("downloading ", tail(strsplit(url, "/")[[1]], 1), appendLF = FALSE)
    f <- RCurl::CFILE(file, mode = "wb")
    if (.Platform$OS.type == "unix") {
        r <- RCurl::curlPerform(url = url, writedata = f@ref, noprogress = TRUE,
                                .opts = list(followlocation = TRUE))
    } else {
        r <- RCurl::curlPerform(url = url, writedata = f@ref, noprogress = TRUE,
                                .opts = list(followlocation = TRUE, ssl.verifypeer = FALSE))
    }
    RCurl::close(f)
    if (verbose)
        message(" - done.")
}


#' CNV.import
#' @description Load combined signal intensities from .idat-Files generated with the Illumina Mouse array or the EPICv2 array. Next, use the resulting \code{data.frame} for \code{CNV.load}.
#' @param array_type character. Choose either \code{"EPICv2"} or \code{"mouse"} as array_type. Default to \code{"EPICv2"}.
#' @param directory Specify the folder that stores the .idat-Files.
#' @param sample_sheet dataframe. Provide a sample sheet with at least three columns: \code{Sample_Name}, \code{Sentrix_ID} and \code{Sentrix_Position}. The spelling of the colnames must be exactly as shown.
#' @param ... Additional parameters (\code{CNV.load} generic, currently not used).
#' @return \code{dataframe} object.
#' @details This method loads the unmethylated and methylated signal intensities for each probe by using the \code{illuminaio} package and calculates the combined signal intensities. It is designed to be used for the Illumina EPICv2 arrays and the Illumina Mouse arrays. Subsequently, the resulting \code{data.frame} should be used for \code{CNV.load}
#' @author Bjarne Daenekas \email{conumee@@hovestadt.bio}
#' @export
setGeneric("CNV.import", function(array_type, directory, sample_sheet, ...) {
  standardGeneric("CNV.import")
})

#' @rdname CNV.import
setMethod("CNV.import", signature(array_type = "character", directory = "character", sample_sheet = "data.frame"),
          function(array_type = "EPICv2", directory = getwd(), sample_sheet = NULL) {
           if(is.null(sample_sheet)){
             stop("please provide a sample sheet")
           }

            object <- new("data.frame")

            if(array_type == "mouse"){
              data("CNV.import_mouse_data")
              anno <- CNV.import_mouse_data
            }

            if(array_type == "EPICv2"){
              data("CNV.import_EPICv2")
              anno <- CNV.import_EPICv2
            }

            if(!array_type %in% c("mouse", "EPICv2")){
              stop("Please choose EPICv2 or mouse as array_type.")
            }

            lf <- list.files(directory, pattern=".idat$", full.names = TRUE)

            object <- suppressWarnings(do.call(cbind, lapply(sample_sheet$Sample_Name, function(i) {
              message(i)
              f <- lf[grepl(paste0(sample_sheet[match(i, sample_sheet$Sample_Name), "Sentrix_ID"], "_", sample_sheet[match(i, sample_sheet$Sample_Name), "Sentrix_Position"], "_Grn"), lf)]
              g <- readIDAT(f)[["Quants"]]
              r <- readIDAT(sub("_Grn", "_Red", f))[["Quants"]]
              t1g <- g[anno[["type1g"]]$AddressA_ID, "Mean", drop = FALSE] + g[anno[["type1g"]]$AddressB_ID, "Mean"]
              t1r <- r[anno[["type1r"]]$AddressA_ID, "Mean", drop = FALSE] + r[anno[["type1r"]]$AddressB_ID, "Mean"]
              t2 <- g[anno[["type2"]]$AddressA_ID, "Mean", drop = FALSE] + r[anno[["type2"]]$AddressA_ID, "Mean"]
              names(t1g) <- anno[["type1g"]]$Name
              names(t1r) <- anno[["type1r"]]$Name
              names(t2) <- anno[["type2"]]$Name
              c(t1g, t1r, t2)
            })))
            colnames(object) <- sample_sheet$Sample_Name
            object <- as.data.frame(object)
            return(object)
            })

#' CNV.define_detail
#' @description Create a \code{GRanges} object for detail regions. Downstream plotting functions will create individual plots for each detail region.
#' @param array_type character. One of \code{450k}, \code{EPIC}, \code{EPICv2} or \code{mouse}. Defaults to \code{450k}.
#' @param gene character vector. Specify the genes by using their \code{SYMBOL} (case sensitive). Default to \code{predefined} which loads a predefined set of onco- and tumor suppressor genes.
#' @param ... Additional parameters (\code{CNV.load} generic, currently not used).
#' @return \code{GRanges} object
#' @examples
#' #create a GRanges object of detail regions
#' detail_regions <- CNV.define_detail(array_type = "450k", symbol = c("CDK6", "PTEN", "MYCN"))
#' detail_regions
#' @author Bjarne Daenekas \email{conumee@@hovestadt.bio}
#' @export
        CNV.define_detail <- function(array_type = "450k", gene = "predefined") {
            if(array_type %in% c("450k", "EPIC", "EPICv2")) {
            if(gene[1] == "predefined"){
              message("using set of predefined regions (hg19 genome)")
              object <- new("GRanges")
              data("detail_regions")
              object <- detail_regions
              return(object)
            }
            data("genes")

            if(any(gene %in% genes$SYMBOL == FALSE)){
              ind <- which(symbol %in% genes$SYMBOL == FALSE)
              message(paste(symbol[ind]," is not part of the gene annotation. ", sep = ""))
            }

            subset <- genes[which(genes$SYMBOL %in% gene)]
            subset$thick <- IRanges(start = start(subset) - 100000, end = end(subset) + 100000)
            colnames(mcols(subset)) <- c("name", "thick")
            names(subset) <- 1:length(subset)
            return(subset)
            }
          if(array_type == "mouse") {
            if(gene[1] == "predefined"){
              message("using set of predefined regions (mm10 genome)")
              object <- new("GRanges")
              data("detail_regions_mouse")
              object <- detail_regions_mouse
              return(object)
            }
            data("genes_mm10")

            if(any(gene %in% genes_mm10$SYMBOL == FALSE)){
              ind <- which(gene %in% genes_mm10$SYMBOL == FALSE)
              message(paste(gene[ind]," is not part of the gene annotation. ", sep = ""))
            }

            subset <- genes_mm10[which(genes_mm10$SYMBOL %in% gene)]
            subset$thick <- subset$thick <- IRanges(start = start(subset) - 100000, end = end(subset) + 100000)
            colnames(mcols(subset)) <- c("name", "thick")
            names(subset) <- 1:length(subset)
            return(subset)
          } else {
            message("please choose 450k, EPIC, EPICv2 or mouse as array_type")
          }
          }


