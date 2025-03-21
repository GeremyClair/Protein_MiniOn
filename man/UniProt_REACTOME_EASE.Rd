% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/1_Proteinminion-functions.R
\name{UniProt_REACTOME_EASE}
\alias{UniProt_REACTOME_EASE}
\title{UniProt_REACTOME_EASE()}
\usage{
UniProt_REACTOME_EASE(
  query,
  universe,
  organismID = "9606",
  reviewed = TRUE,
  preloaded_UniProtTable = TRUE
)
}
\arguments{
\item{query}{should be a character vector corresponding to a query list (UniProt_Accession or UniProt_ID)}

\item{universe}{should be a character vector corresponding to a universe list (Uniprot Entry), if missing the all table from the UniProtTable will be used as universe}

\item{organismID}{Should be a UniProt organism ID, by default Homo sapiens will be used ("9606")}

\item{reviewed}{this indicate if the proteins extracted have to be reviewed or not (typically the ones present on swissProt), by default TRUE (the non reviewed proteins are discarded)}

\item{preloaded_UniProtTable}{this indicate wether or not the preloaded UniProtTable and UniProtTable_GO should be used, if FALSE the UniProtTable and UniProtTable_GO will be generated again}
}
\value{
a data.frame containing enrichment results.
}
\description{
Allows to use UniprotDB REACTOME information to perform a modified Fisher Exact test corresponding to the one performed by the DAVID Bioinformatics tool for enrichment analysis
}
\author{
Geremy Clair
}
