% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/1_Proteinminion-functions.R
\name{Enrich_EASE}
\alias{Enrich_EASE}
\title{UniProt_Fisher()}
\usage{
Enrich_EASE(
  query,
  universe,
  type = c("GO", "KEGG", "REACTOME"),
  organismID = "9606",
  reviewed = TRUE,
  preloaded_UniProtTable = TRUE
)
}
\arguments{
\item{query}{should be a character vector corresponding to a query list (UniProt_Accession or UniProt_ID)}

\item{universe}{should be a character vector corresponding to a universe list (Uniprot Entry), if missing the all table from the UniProtTable will be used as universe}

\item{type}{should be acharacter vector comprising KEGG,REACTOME, and GO : by default c('GO','KEGG','REACTOME')}

\item{organismID}{Should be a UniProt organism ID, by default Homo sapiens will be used ("9606")}

\item{reviewed}{this indicate if the proteins extracted have to be reviewed or not (typically the ones present on swissProt), by default TRUE (the non reviewed proteins are discarded)}

\item{preloaded_UniProtTable}{this indicate wether or not the preloaded UniProtTable and UniProtTable_REACTOME should be used, if FALSE the UniProtTable and UniProtTable_REACTOME will be generated again}
}
\value{
a data.frame containing enrichment results
}
\description{
Allows to use UniprotDB information to perform a FisherExact test for multiple database enrichment analysis
}
\author{
Geremy Clair
}
