% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/1_Proteinminion-functions.R
\name{UniProt_REACTOME_KS}
\alias{UniProt_REACTOME_KS}
\title{UniProt_REACTOME_KS()}
\usage{
UniProt_REACTOME_KS(
  rankingTable,
  order = "ascending",
  organismID = "9606",
  reviewed = TRUE,
  preloaded_UniProtTable = TRUE
)

UniProt_REACTOME_KS(
  rankingTable,
  order = "ascending",
  organismID = "9606",
  reviewed = TRUE,
  preloaded_UniProtTable = TRUE
)
}
\arguments{
\item{rankingTable}{should be a two column data.frame containing the Uniprot identifiers list as first column (UniProt_Accession or UniProt_ID) and the ranking values as second column}

\item{order}{should either be "ascending" or "descending" to indicate the order of the ranking values to use.}

\item{organismID}{Should be a UniProt organism ID, by default Homo sapiens will be used ("9606")}

\item{reviewed}{this indicate if the proteins extracted have to be reviewed or not (typically the ones present on swissProt), by default TRUE (the non reviewed proteins are discarded)}

\item{preloaded_UniProtTable}{this indicate wether or not the preloaded UniProtTable and UniProtTable_GO should be used, if FALSE the UniProtTable and UniProtTable_GO will be generated again}
}
\value{
a data.frame containing enrichment results from the GO retrieve from UniProt.

a data.frame containing enrichment results from the GO retrieve from UniProt.
}
\description{
Allows to use UniprotDB REACTOME information to perform a modified FisherExact test corresponding to the one performed by the DAVID Bioinformatics tool for enrichment analysis

Allows to use UniprotDB REACTOME information to perform a modified FisherExact test corresponding to the one performed by the DAVID Bioinformatics tool for enrichment analysis
}
