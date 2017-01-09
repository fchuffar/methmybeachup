% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analyse_meth.R
\name{analyse_meth}
\alias{analyse_meth}
\title{A Function Analyses Methylome data}
\usage{
analyse_meth(gene, sunexp_data, sunexp_platform, cols, PLOT = TRUE,
  up_str = 5000, dwn_str = 5000, win_size = 1500, wig_size = 50)
}
\arguments{
\item{gene}{A vector describing the gene (line of a bed file).}

\item{sunexp_data}{A matrix of beta values.}

\item{sunexp_platform}{A data frame describing CpG positions.}

\item{cols}{A color vectors indexed by by samples names.}

\item{PLOT}{A boolean defining if graphical output must be dispayed on the graphical output.}

\item{up_str}{An integer specifying up stream size (in bp).}

\item{dwn_str}{An integer specifying down stream size (in bp).}

\item{win_size}{An integer specifying slidding window size (in bp).}

\item{wig_size}{An integer specifying wiggle size (in bp).}
}
\value{
A matrix of convolved probes signal around the TSS of the selected gene.
}
\description{
This function analyses methylome data.
}
\examples{
cols = as.numeric(sunexp_design$sex) * 2
names(cols) = rownames(sunexp_design)
gene = genes[1,]
res1 = analyse_meth(gene, sunexp_data, sunexp_platform, cols, PLOT=TRUE)
legend("topright", col=as.numeric(unique(sunexp_design$sex)) * 2, legend=unique(sunexp_design$sex), lty=1)
}
