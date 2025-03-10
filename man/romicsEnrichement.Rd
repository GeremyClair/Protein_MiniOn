% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2_romicsEnrichment.R
\name{romicsEnrichement}
\alias{romicsEnrichement}
\title{romicsEnrichment()
allows to perform proteinMinion enrichments directly on a romics_object. The Query list can be generated using the filters.}
\usage{
romicsEnrichement(
  romics_object,
  organismID = "9606",
  type = c("GO", "KEGG", "REACTOME"),
  ANOVA_filter = c("none", "p", "padj"),
  p = 0.05,
  cluster_filter = "none",
  statCol_filter = "none",
  statCol_text = "<=0.05",
  statCol_filter2 = "none",
  statCol_text2 = "<=0.05",
  statCol_filter3 = "none",
  statCol_text3 = "<=0.05",
  enrichment_function = "EASE",
  ...
)
}
\arguments{
\item{romics_object}{A romics_object created using romicsCreateObject() and transformed using the function log2transform() or log10transform()}
}
\value{
This function will print the pmartR filtering details and will return 2 plots the first one is a scatter plot of the pvalue by log2(Robust Mahalanobis Distance) the second one is a scatter plot of the log2(Robust Mahalanobis Distance) per sample.
}
\description{
romicsEnrichment()
allows to perform proteinMinion enrichments directly on a romics_object. The Query list can be generated using the filters.
}
\details{
This function requires the package 'pmartR' to be installed and loaded to be excecuted. It will calculate and plot the samples to be filtered out using the function Romics_outlier_eval().
}
\references{
Matzke, M., Waters, K., Metz, T., Jacobs, J., Sims, A., Baric, R., Pounds, J., and Webb-Robertson, B.J. (2011), Improved quality control processing of peptide-centric LC-MS proteomics data. Bioinformatics. 27(20): 2866-2872.
}
\author{
Geremy Clair
}
