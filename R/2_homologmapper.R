#Organism.Unihomolog<- function(organism.id){
#  baseURL<-"http://www.uniprot.org/uniprot/?query=organism:"
#  URL<-paste0(baseURL,organism.id,"&columns=id,entry name,genes,reviewed,go,organism,database(refseq),database(KEGG),database(GeneID),database(HOVERGEN),database(eggNOG),database(InParanoid)&format=tab")
#  infos<-read.csv(URL,sep="\t",header=T)
#  return<-infos
#}

#example
#start_time <- Sys.time()
#organism.id<-"43179"
#squirrelhomologs<-Organism.Unihomolog(organism)
#end_time <- Sys.time()
#end_time-start_time
