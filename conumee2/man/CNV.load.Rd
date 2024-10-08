% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load.R
\name{CNV.load}
\alias{CNV.load}
\alias{CNV.load,MethylSet-method}
\alias{CNV.load,data.frame-method}
\alias{CNV.load,matrix-method}
\alias{CNV.load,numeric-method}
\alias{CNV.load,RnBeadRawSet-method}
\title{CNV.load}
\usage{
CNV.load(input, ...)

\S4method{CNV.load}{MethylSet}(input, names = NULL)

\S4method{CNV.load}{data.frame}(input, names = NULL)

\S4method{CNV.load}{matrix}(input, names = NULL, ...)

\S4method{CNV.load}{numeric}(input, names = NULL)

\S4method{CNV.load}{RnBeadRawSet}(input, names = NULL)
}
\arguments{
\item{input}{Object of MethylSet class (minfi package), data.frame class, matrix class or numeric class.}

\item{...}{Additional parameters (\code{CNV.load} generic, currently not used).}

\item{names}{Vector specifying sample names. If not supplied, colnames are used. For MethylSet input, the first column of pData(input) matching 'name' (grep) is used.}
}
\value{
\code{CNV.data} object.
}
\description{
Prepare combined intensities from various input objects.
}
\details{
This method gathers combined intensities of the Methylated and Unmethylated signals for all supplied probes. Probe IDs must be supplied as row names or in a seperate column named `ID_REF` or `TargetID`.
If column names match 'intensity', only those columns are used. Else, if column names match 'signal' or 'methylated', only those columns are used. Otherwise, all columns are used.
}
\examples{
library(minfiData)
d <- CNV.load(MsetEx)
d

library(sesame)
sdfs = openSesame(idat_dir, prep = "QCDPB", func = NULL)
combined.intensities <- lapply(sdfs, totalIntensities)
d <- CNV.load(combined.intensities)
d
}
\author{
Volker Hovestadt \email{conumee@hovestadt.bio}
}
