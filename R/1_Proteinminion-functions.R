#' GO_details()
#' @description Retrieve GO_details from the GO.db from Bioconductor
#' @return a df containing the following details : GO.ID, Description, and GO.ontology.type
#' @examples GO_details<-GO_details()
#' @examples GO_details
GO_details<-function(){
  if(!"GO.db" %in% rownames(installed.packages())){stop("The package 'GO.db' is required to execute this function")}
  if(!"GO.db" %in% (.packages())){library("GO.db")}

GO_details = select(GO.db, keys(GO.db, "GOID"), c("TERM", "ONTOLOGY"),)
colnames(GO_details)<-c("GO.ID","Description","GO.ontology.type")
GO_details
}

#' download_UniProtTable()
#' @description Retrieve UniProt info for the proteins of a given organism.
#' @param organismID Should be a UniProt organism ID, by default Homo sapiens will be used ("9606").
#' @param reviewed indicate if the proteins extracted have to be reviewed or not (typically the ones present on swissProt), by default TRUE (the non reviewed proteins are discarded), if FALSE both reviewed and not reviewed will be conserved.
#' @param export indicates if the table should be exported for future use (using a previously downloaded table dramatically reduces the analysis time).
#' @param file should be a string indicating the location and name of the destination file, by default the file will be named with the format "UniProtTable_organismID_date.txt".
#'
#' @return This function returns a data.frame named UniProtTable containing the detail on the proteins for the queried organism, if the export is set on TRUE, the function will also save the table in the .rda format at the specified location
#'
#' @examples download_UniProtTable(organismID="9606", reviewed=TRUE, export=FALSE, file="UniProtTable_organismID_date.txt")
#'
download_UniProtTable<- function(organismID="9606",proteomeID="proteomeID",reviewed=FALSE, export=FALSE, file="UniProtTable_organismID_date.rda"){
  #general checkings
  if(missing(organismID)&missing(proteomeID)){warning("The 'organismID' and 'proteomeID' were left empty, 'homo sapiens' (9006) will be used by default")
    organismID<-"9606"}
  if(missing(reviewed)){reviewed<-FALSE}
  if(missing(reviewed)){reviewed<-FALSE}
  if(!is.logical(reviewed)){stop("<reviewed> should be either TRUE or FALSE")}
  if(missing(export)){export<-FALSE}
  if(!is.logical(reviewed)){stop("<reviewed> should be either TRUE or FALSE")}
  if(missing(file)){file=paste0("UniProtTable_OrganismID_",organismID,"_",format(Sys.time(),"%b_%d_%Y"),".txt")}


  #Base URL
  URL<-"http://www.uniprot.org/uniprot/?query="
  #reviewed?
  if(reviewed==TRUE){URL<-paste0(URL,"reviewed:yes+AND+")}
  #organismID
  if(!missing(organismID)){
    URL<-paste0(URL,"organism:",organismID)
    if(!missing(proteomeID)){URL<-paste0(URL,"+AND+")}
  }
  #proteomeID
  if(!missing(proteomeID)){
    URL<-paste0(URL,"proteome:",proteomeID)}
    #query type
    URL<-paste0(URL,"&columns=id,entry%20name,genes,reviewed,go,organism,database(refseq),database(KEGG),database(GeneID)&format=tab")

  #read the table from the internet
  UniProtTable<-read.csv(URL,sep="\t",header=T)

  #Adjust column names
  colnames(UniProtTable)<-c("Uniprot_Accession", "Uniprot_ID", "Gene_Name", "Status", "GO_List", "Organism", "Crossref_RefSeq","Crossref_KEGG", "Crossref_geneID")

  #export the global environment
  UniProtTable<<-UniProtTable

  if(export==TRUE){
    save(UniProtTable,file = file)
    #write.table(UniProtTable,file = file, row.names = FALSE, col.names = TRUE, quote=FALSE, sep="\t")
  }

}

#' download_UniProtFasta()
#' @description Retrieve UniProt fasta file for a given organism
#' @param organismID Should be a UniProt organism ID
#' @param proteome Should be a UniProt proteome ID
#' @param reviewed this indicate if the proteins extracted have to be reviewed or not (typically the ones present on swissProt), by default FALSE.
#' @param export indicates if the table should be exported for future use (using a previously downloaded table dramatically reduces the analysis time).
#' @param file should be a string indicating the location and name of the destination file, by default the file will be named with the format "UniProtFasta_organismID_date.fasta".
#'
#' @return This function returns a data.frame named UniProtFasta containing the detail on the proteins for the queried organism, if the export is set on TRUE, the function will also return the table in the fasta file
#'
#' @examples download_UniProtFasta(organismID="9606", reviewed=TRUE, export=TRUE,file="UniProtFasta_organismID_date.fasta")
#'
download_UniProtFasta<-function(organismID="9606", proteomeID="proteomeID", reviewed=FALSE, export=FALSE, file="UniProtFasta_organismID_date.fasta"){
  #general checkings
  if(missing(organismID)&missing(proteomeID)){warning("The 'organismID' and 'proteomeID' were left empty, 'homo sapiens' (9006) will be used by default")
    organismID<-"9606"}
  if(missing(reviewed)){reviewed<-FALSE}
  if(!is.logical(reviewed)){warning("<reviewed> should be either TRUE or FALSE, by default reviewed=FALSE was used")
  reviewed<- FALSE}
  if(missing(export)){export<-FALSE}
  if(!is.logical(reviewed)){stop("<reviewed> should be either TRUE or FALSE")}
  if(missing(file)){file=paste0("UniProtFasta_OrganismID_",organismID,"_",format(Sys.time(),"%b_%d_%Y"),".fasta")}

  #Base URL
  URL<-"http://www.uniprot.org/uniprot/?query="
  #reviewed?
  if(reviewed==TRUE){URL<-paste0(URL,"reviewed:yes+AND+")}
  #organism ID ?
  if(!missing(organismID)){
    URL<-paste0(URL,"organism:",organismID)
    if(!missing(proteomeID)){URL<-paste0(URL,"+AND+")}
  }
  #proteome
  if(!missing(proteomeID)){
    URL<-paste0(URL,"proteome:",proteomeID)}

     #query type
  URL<-paste0(URL,"&format=fasta")

  #read the table from the internet
  UniProtFasta<-read.csv(URL,sep="\t",header=F)

  #print your fasta file number of sequences
  print(paste("The UniProtFasta generated contains",sum(grepl(">",as.character(t(UniProtFasta)))),"sequences"))

  #export the global environment
  UniProtFasta<<-UniProtFasta

  #write fasta file if export=TRUE
  if(export==TRUE){
  write.table(UniProtFasta,file = file, row.names = FALSE, col.names = FALSE, quote=FALSE)
  }
}

#' UniprotFastaParser()
#' @description allow to extract the information from a fasta file in a table format
#' @param fasta.file should be a character vector corresponding to a query list (UniProt_Accession or UniProt_ID)
#'
#' @return a data.frame containing enrichment results
#'
#' @examples UniProtGOfisher(QueryExample, UniverseExample, organism.id=organismExample)
#'
UniprotFastaParser<-function(file="file.fasta"){
  name_file<-print(file)
  if(!grepl(".fasta|.faa",name_file)){stop("Your file shoud be a <.fasta> or a <.faa> file.")}
  fasta<-data.frame(read.delim(file = name_file))
  fasta<-as.character(fasta[grepl(">",fasta[,1]),])

  #the following is only valid for UniProtKB generated fasta files
  #create the db object
  db<-sub("\\|.*", "", fasta)
  db<-sub(">","",db)
  db<-sub("sp","UniProtKB/Swiss-Prot",db)
  db<-sub("tr","UniProtKB/TrEMBL",db)

  #get the Uniprot_accession
  Uniprot_accession<-sub(".*\\|(.*)\\|.*", "\\1", fasta, perl=TRUE)

  #get the Uniprot_name
  Uniprot_name<-sub("^(.*? | .*?) | .*", "\\1", fasta)
  Uniprot_name<-sub(".*\\|","",Uniprot_name)

  #get the uniprot description
  all_description<-sub("^\\S+\\s+", "", fasta)

  #get only the general description
  description<-sub(" OS\\=.*$", "\\1", all_description)

  #get the OrganismName
  all_description<-sub(".*OS\\=", "\\1", all_description)
  OrganismName<-sub(" OX\\=.*$", "\\1", all_description)

  #get the SequenceVersion
  SequenceVersion<-sub(".* SV=", "\\1", all_description)

  #get the ProteinExistence
  all_description<-sub(" SV=.*$", "\\1", all_description)
  ProteinExistence<-sub(".* PE=", "\\1", all_description)

  #get the GeneName
  all_description<-sub(" PE=.*$", "\\1", all_description)
  GeneName<-sub(".* GN=", "\\1", all_description)

  #get the OrganismIdentifier
  all_description<-sub(" GN=.*$", "\\1", all_description)
  OrganismIdentifier<-sub(".* OX=", "\\1", all_description)

  output<-data.frame(Uniprot_accession=Uniprot_accession,Uniprot_name=Uniprot_name,description=description,OrganismName=OrganismName,OrganismIdentifier=OrganismIdentifier,GeneName=GeneName,ProteinExistence=ProteinExistence,SequenceVersion=SequenceVersion)
}

#' generate_UniProtTable_GO()
#' @description Retrieve UniProt GOterms for the proteins of a given organism
#' @param organismID Should be a UniProt organism ID, by default Homo sapiens will be used ("9606")
#' @param reviewed this indicate if the proteins extracted have to be reviewed or not (typically the ones present on swissProt), by default TRUE (the non reviewed proteins are discarded)
#' @param preloaded_UniProtTable this indicate wether or not the preloaded UniProtTable and UniProtTable_GO should be used, if FALSE the UniProtTable and UniProtTable_GO will be generated again
#'
#' @return This file return a data.frame named UniProtTable_GO containing the list of GOs associated with the different proteins from a UniProtTable. If the export is set on, the UniProtTable_Go will also be saved in the .rda format at the specified location.
#'
#' @examples generate_UniProtTable_GO(organismID="9606",reviewed=TRUE,preloaded_UniProtTable=TRUE,export=TRUE,file="UniProtTable_GO_organismID_date.rda")
#'
generate_UniProtTable_GO<- function(organismID="9606",reviewed=TRUE,preloaded_UniProtTable=TRUE,export=FALSE,file="UniProtTable_GO_organismID_date.rda"){
  #general checkings
  if(missing(organismID)){organismID<-"9606"}
  if(missing(reviewed)){reviewed<-TRUE}
  if(!is.logical(reviewed)){stop("<reviewed> should be either TRUE or FALSE")}
  if(missing(preloaded_UniProtTable)){preloaded_UniProtTable<-TRUE}
  if(!is.logical(preloaded_UniProtTable)){stop("<preloaded_UniProtTable> should be either TRUE or FALSE")}
  if(missing(export)){export<-FALSE}
  if(!is.logical(reviewed)){stop("<reviewed> should be either TRUE or FALSE")}
  if(missing(file)){file=paste0("UniProtTable_GO_OrganismID_",organismID,"_",format(Sys.time(),"%b_%d_%Y"),".txt")}


  #if the object UniProtTable already exists, use it only if preloaded_UniProtTable is TRUE
  if(preloaded_UniProtTable==TRUE&&exists("UniProtTable")){
    UPT<- UniProtTable
  }else{
    download_UniProtTable(organismID=organismID, reviewed=reviewed)
    UPT<-UniProtTable
  }

  #extract the required columns from the UniProtTable (UPT) and ensure the content is considered as character by R
  UPT<-UPT[,c(1,2,5)]
  for(i in 1:3){
    UPT[,i]<-as.character(t(UPT[,i]))
  }

  #remove the end bracket
  UPT[,3]<- gsub("\\]","",UPT[,3])

  #create a list to store the different GO terms associated with each protein
  uniList<-list()
  #create objects to store the Uniprot_accession and the Uniprot_ID
  Uniprot_Accession<-character()
  Uniprot_ID<-character()

  #fill the three objects created above
  for (i in 1:nrow(UPT)){
    uniList[[i]]<-trimws(unlist(strsplit(UPT[i,3],";")))
    Uniprot_Accession<-c(Uniprot_Accession,rep(UPT[i,1],length(uniList[[i]])))
    Uniprot_ID<-c(Uniprot_ID,rep(UPT[i,2],length(uniList[[i]])))
    }

  #create the Data frame containing these details
  UniProtTable_GO<-data.frame(Uniprot_Accession=Uniprot_Accession, Uniprot_ID=Uniprot_ID, GO=unlist(uniList))

  #split the last column using the "[GO"
  UniProtTable_GO <- UniProtTable_GO %>% separate(GO, into= c('GO_description','GO_accession'),"\\[GO\\:" )

  #Place the reordered table in the Global Environment
  UniProtTable_GO <<- UniProtTable_GO[,c(1,2,4,3)]

  #export file if export==TRUE
  if(export==TRUE){
    save(UniProtTable_GO,file = file)
  }

  }

#' generate_UniProtTable_KEGG()
#' @description Retrieve UniProt GOterms for the proteins of a given organism
#' @param organismID Should be a UniProt organism ID, by default Homo sapiens will be used ("9606")
#' @param reviewed this indicate if the proteins extracted have to be reviewed or not (typically the ones present on swissProt), by default TRUE (the non reviewed proteins are discarded)
#' @param preloaded_UniProtTable this indicate wether or not the preloaded UniProtTable and UniProtTable_GO should be used, if FALSE the UniProtTable and UniProtTable_GO will be generated again
#'
#' @return This function generate a KEGG enrichment table named UniProtTable_KEGG usable with the enrichment functions Uniprot_KEGG_*. If the export is set on, the UniProtTable_KEGG will also be saved in the .rda format at the specified location.
#'
#' @examples generate_UniProtTable_KEGG(organismID="9606",reviewed=TRUE,preloaded_UniProtTable=TRUE,export=FALSE,file="UniProtTable_KEGG_organismID_date.rda")
#'
generate_UniProtTable_KEGG<- function(organismID="9606",reviewed=TRUE,preloaded_UniProtTable=TRUE,export=FALSE,file="UniProtTable_KEGG_organismID_date.rda"){
  #general checkings
  if(missing(organismID)){organismID<-"9606"}
  if(missing(reviewed)){reviewed<-TRUE}
  if(!is.logical(reviewed)){stop("<reviewed> should be either TRUE or FALSE")}
  if(missing(preloaded_UniProtTable)){preloaded_UniProtTable<-TRUE}
  if(!is.logical(preloaded_UniProtTable)){stop("<preloaded_UniProtTable> should be either TRUE or FALSE")}
  if(missing(export)){export<-FALSE}
  if(!is.logical(reviewed)){stop("<reviewed> should be either TRUE or FALSE")}
  if(missing(file)){file=paste0("UniProtTable_GO_OrganismID_",organismID,"_",format(Sys.time(),"%b_%d_%Y"),".txt")}


  #if the object UniProtTable already exists, use it only if preloaded_UniProtTable is TRUE
  if(preloaded_UniProtTable==TRUE&&exists("UniProtTable")){
    UPT<- UniProtTable
  }else{
    download_UniProtTable(organismID=organismID, reviewed=reviewed)
    UPT<-UniProtTable
  }

  #extract the required columns from the UniProtTable (UPT) and ensure the content is considered as character by R
  UPT<-UPT[,c(1,2,8)]
  for(i in 1:3){
    UPT[,i]<-as.character(t(UPT[,i]))
  }

  #create a list to store the different GO terms associated with each protein
  uniList<-list()
  #create objects to store the Uniprot_accession and the Uniprot_ID
  Uniprot_Accession<-character()
  Uniprot_ID<-character()

  #fill the three objects created above
  for (i in 1:nrow(UPT)){
    uniList[[i]]<-trimws(unlist(strsplit(UPT[i,3],";")))
    Uniprot_Accession<-c(Uniprot_Accession,rep(UPT[i,1],length(uniList[[i]])))
    Uniprot_ID<-c(Uniprot_ID,rep(UPT[i,2],length(uniList[[i]])))
  }

  #create the Data frame containing these details
  UniProtTable_KEGG<-data.frame(Uniprot_Accession=Uniprot_Accession, Uniprot_ID=Uniprot_ID, KEGG=unlist(uniList))

  #split the last column using the ":"
  UniProtTable_KEGG <- UniProtTable_KEGG %>% separate(KEGG, into= c('KEGG_organism','KEGG_ID'),"\\:" )
  UniProtTable_KEGG$KEGG_ID <- paste0(UniProtTable_KEGG$KEGG_organism,":", UniProtTable_KEGG$KEGG_ID)

  #Get the lists of correspondance for the different organisms
  #first extract the unique organisms IDs
  KEGG_orglist <- unique(UniProtTable_KEGG$KEGG_organism)

  #Now retrieve the organisms pathway lists from KEGG using their API
  KEGG_URL<-"http://rest.kegg.jp/"

  KEGG_path_ID<- data.frame()
  KEGG_path_names<- data.frame()
  for(i in 1:length(KEGG_orglist)){
    KEGG_path_ID<-rbind(KEGG_path_ID,read.delim(file=paste0(KEGG_URL,"link/pathway/",KEGG_orglist[i]),header=FALSE,col.names = c("KEGG_ID","Pathway_ID")))
    KEGG_path_names<- rbind(KEGG_path_names,read.delim(file=paste0(KEGG_URL,"list/pathway/",KEGG_orglist[i]),header=FALSE,col.names = c("Pathway_ID","Pathway_name")))
    }

  KEGG_path_ID<-merge(KEGG_path_ID, KEGG_path_names, by = "Pathway_ID")

  #now merge the table with the other info from the UniProtTable_KEGG table
  UniProtTable_KEGG<-merge(KEGG_path_ID,UniProtTable_KEGG, by="KEGG_ID")

  #Place the reordered table in the Global Environment
  UniProtTable_KEGG <<- UniProtTable_KEGG[,c(4,5,1,2,3)]

  #count the unique proteins mapped to KEGG pathways and the number of unique pathways
  print(paste0(length(unique(UniProtTable_KEGG$Uniprot_Accession))," proteins were mapped onto ",length(unique(UniProtTable_KEGG$Pathway_ID))," KEGG pathways  (listed in the object KEGG_path_names)"))

  #Also place the pathway list in the global Environment
  KEGG_path_names <<- KEGG_path_names

  #export file if export==TRUE
  if(export==TRUE){
    save(UniProtTable_KEGG,file = file)

  }
}

#' UniProt_GO_Fisher()
#' @description Allows to use UniprotDB Gene Ontology information to perform a FisherExact test for enrichment analysis
#' @param query should be a character vector corresponding to a query list (UniProt_Accession or UniProt_ID)
#' @param universe should be a character vector corresponding to a universe list (Uniprot Entry), if missing the all table from the UniProtTable will be used as universe
#' @param organismID Should be a UniProt organism ID, by default Homo sapiens will be used ("9606")
#' @param reviewed this indicate if the proteins extracted have to be reviewed or not (typically the ones present on swissProt), by default TRUE (the non reviewed proteins are discarded)
#' @param preloaded_UniProtTable this indicate wether or not the preloaded UniProtTable and UniProtTable_GO should be used, if FALSE the UniProtTable and UniProtTable_GO will be generated again
#'
#' @return a data.frame containing enrichment results
#'
#' @examples UniProt_GO_Fisher(QueryExample, UniverseExample, organism.id=organismExample)
#'

UniProt_GO_Fisher<-function(query, universe, organismID="9606", reviewed=TRUE, preloaded_UniProtTable=TRUE){
  #general checkings
  if(missing(organismID)){organismID<-"9606"}
  if(missing(reviewed)){reviewed<-TRUE}
  if(!is.logical(reviewed)){stop("<reviewed> should be either TRUE or FALSE")}
  if(missing(preloaded_UniProtTable)){preloaded_UniProtTable<-TRUE}
  if(!is.logical(preloaded_UniProtTable)){stop("<preloaded_UniProtTable> should be either TRUE or FALSE")}

  #load the UniProtTable and UniProtTable_GO (not if it is preloaded)
  if(preloaded_UniProtTable==TRUE&&exists("UniProtTable")){
  UPT<- UniProtTable
  }else {
    download_UniProtTable(organismID=organismID,reviewed=reviewed)
    UPT<- UniProtTable
    }

  if(preloaded_UniProtTable==TRUE&&exists("UniProtTable_GO")){
    UPT_GO<- UniProtTable_GO
  }else {
    generate_UniProtTable_GO(organismID=organismID,reviewed=reviewed,preloaded_UniProtTable==TRUE)
    UPT_GO<- UniProtTable_GO
    }

  #check query
  if(missing(query)){stop("Your <query> is missing")}

  #check how many of the query are in the lists UniProtTable
  query_ID <- query[query %in% UniProtTable$Uniprot_ID]
  query_Accession <- query[query %in% UniProtTable$Uniprot_Accession]

  print(paste0("Your <query> contained ", length(query_ID), " UniProt_IDs and ", length(query_Accession), " UniProt_Accession (some might be redundant)." ))

  #if no values in the UniProtTable stop
  if (length(query_ID)==0&&length(query_Accession)==0){stop("Your <query> should be Uniprot_IDs or a Uniprot_Accessions, ensure that the organism of your query is present in the UniProtTable.")}

  #convert uniprot_Accession in Uniprot
  if(length(query_Accession)>0){
    print("The uniprot_Accession of your query were converted in Uniprot_IDs.")
    query_Accession<- as.character(t(UniProtTable[match(query_Accession,UniProtTable$Uniprot_Accession),2]))
    query_ID<-c(query_ID,query_Accession)
  }

  #remove duplicates
  duplicated_ID<-duplicated(query_ID)
  if(sum(duplicated_ID)>0){
    print(paste0("Your list of Uniprot_IDs was containing ", sum(duplicated_ID)," duplicated value(s), duplicates were removed, the query now comprises ", length(query_ID)-sum(duplicated_ID), " unique IDs."))
    query_ID<-query_ID[!duplicated_ID]
    }

  #check universe
  if(missing(universe)){
  print("Your <universe> was missing the global list of ID present in UniProtTable was used as background.")
  universe<-as.character(t(UniProtTable$Uniprot_ID))
  }

  #check how many of the universe are in the lists UniProtTable
  universe_ID <- universe[universe %in% UniProtTable$Uniprot_ID]
  universe_Accession <- universe[universe %in% UniProtTable$Uniprot_Accession]

  print(paste0("Your universe contained ", length(universe_ID), " UniProt_IDs and ", length(universe_Accession), " UniProt_Accession (some might be redundant)." ))

  #if no values in the UniProtTable stop
  if (length(universe_ID)==0&&length(universe_Accession)==0){stop("Your <universe> should be Uniprot_IDs or a Uniprot_Accessions, ensure that the organism of your universe is present in the UniProtTable.")}

  #convert uniprot_Accession in Uniprot
  if(length(universe_Accession)>0){
    print("The uniprot_Accession of your universe were converted in Uniprot_IDs.")
    universe_Accession<- as.character(t(UniProtTable[match(universe_Accession,UniProtTable$Uniprot_Accession),2]))
    universe_ID<-c(universe_ID,universe_Accession)
  }

  #remove duplicates
  duplicated_ID<-duplicated(universe_ID)
  if(sum(duplicated_ID)>0){
    print(paste0("Your list of Uniprot_IDs was containing ", sum(duplicated_ID)," duplicated value(s), duplicates were removed, the universe now comprises ", length(universe_ID)-sum(duplicated_ID), " unique IDs."))
    universe_ID<-universe_ID[!duplicated_ID]
  }

  #check if the universe contains all the proteins from the query
  if(sum(is.na(match(query_ID,universe_ID)))>0){stop("Your <universe> does not comprises all the IDs present in your <query>")}

  #generate the UniProtTable_GO for the query and the universe
  UPT_GO_query <- UPT_GO[UPT_GO$Uniprot_ID %in% query_ID,]
  UPT_GO_universe <-UPT_GO[UPT_GO$Uniprot_ID %in% universe_ID,]

  #create the unique_GO list based on the UPT_GO_query
  unique_GO<- unique(UPT_GO_query[,3:4])

  #create the Contingencies matrix
  contingency_list<-list()
  for (i in 1:nrow(unique_GO)){
    contingency_list[[i]]<-matrix(ncol=2,nrow=2)
    contingency_list[[i]][1,1]<- sum(UPT_GO_query[,3]==unique_GO[i,1])
    contingency_list[[i]][1,2]<- sum(UPT_GO_universe[,3]==unique_GO[i,1])
    contingency_list[[i]][2,1]<- length(query_ID)-contingency_list[[i]][1,1]
    contingency_list[[i]][2,2]<- length(universe)-contingency_list[[i]][1,2]
    }
  names(contingency_list)<-unique_GO[,1]

  #create Final table
  fisher_results<-data.frame(matrix(NA, nrow=nrow(unique_GO), ncol=11))
  colnames(fisher_results)<- c("Ontology_type","Ontology_accession","Ongology_description","Count_query","Count_universe","%_query","%_universe","pval","adjpval","fold_change","Proteins_in_query")
  fisher_results[,2:3]<-unique_GO[,1:2]

  #run the test and populate the table
  for (i in 1:nrow(unique_GO)){
    test<-fisher.test(contingency_list[[i]])
    fisher_results$`Ontology_type`[i]<-"GO"
    fisher_results$pval[i]<-test$p.value
    fisher_results$adjpval[i]<-p.adjust(fisher_results$pval[i] , method = "BH", nrow(unique_GO))
    fisher_results$Count_query[i]<- paste0(contingency_list[[i]][1,1],"/",length(query_ID))
    fisher_results$Count_universe[i]<- paste0(contingency_list[[i]][1,2],"/",length(universe_ID))
    fisher_results$`%_query` [i]<- contingency_list[[i]][1,1]/length(query_ID)*100
    fisher_results$`%_universe`[i]<- contingency_list[[i]][1,2]/length(universe_ID)*100
    fisher_results$fold_change[i]<- fisher_results$`%_query`[i]/fisher_results$`%_universe`[i]
    fisher_results$Proteins_in_query[i]<-paste(UPT_GO_query[UPT_GO_query$GO_accession==fisher_results$Ontology_accession[i],]$Uniprot_ID,collapse =";")
  }

  fisher_results[order(fisher_results$pval),]
}

#' UniProt_GO_EASE()
#' @description Allows to use UniprotDB Gene Ontology information to perform a modified FisherExact test corresponding to the one performed by the DAVID Bioinformatics tool for enrichment analysis
#' @param query should be a character vector corresponding to a query list (UniProt_Accession or UniProt_ID)
#' @param universe should be a character vector corresponding to a universe list (Uniprot Entry), if missing the all table from the UniProtTable will be used as universe
#' @param organismID Should be a UniProt organism ID, by default Homo sapiens will be used ("9606")
#' @param reviewed this indicate if the proteins extracted have to be reviewed or not (typically the ones present on swissProt), by default TRUE (the non reviewed proteins are discarded)
#' @param preloaded_UniProtTable this indicate wether or not the preloaded UniProtTable and UniProtTable_GO should be used, if FALSE the UniProtTable and UniProtTable_GO will be generated again
#'
#' @return a data.frame containing enrichment results from the GO retrieve from UniProt.
#'
#' @examples UniProt_GO_EASE(QueryExample, UniverseExample, organism.id=organismExample)
#'

UniProt_GO_EASE<-function(query, universe, organismID="9606", reviewed=TRUE, preloaded_UniProtTable=TRUE){
  #general checkings
  if(missing(organismID)){organismID<-"9606"}
  if(missing(reviewed)){reviewed<-TRUE}
  if(!is.logical(reviewed)){stop("<reviewed> should be either TRUE or FALSE")}
  if(missing(preloaded_UniProtTable)){preloaded_UniProtTable<-TRUE}
  if(!is.logical(preloaded_UniProtTable)){stop("<preloaded_UniProtTable> should be either TRUE or FALSE")}

  #load the UniProtTable and UniProtTable_GO (not if it is preloaded)
  if(preloaded_UniProtTable==TRUE&&exists("UniProtTable")){
    UPT<- UniProtTable
  }else {
    download_UniProtTable(organismID=organismID,reviewed=reviewed)
    }

  if(preloaded_UniProtTable==TRUE&&exists("UniProtTable_GO")){
    UPT_GO<- UniProtTable_GO
  }else {generate_UniProtTable_GO(organismID=organismID,reviewed=reviewed,preloaded_UniProtTable==TRUE)}

  #check query
  if(missing(query)){stop("Your <query> is missing")}

  #check how many of the query are in the lists UniProtTable
  query_ID <- query[query %in% UniProtTable$Uniprot_ID]
  query_Accession <- query[query %in% UniProtTable$Uniprot_Accession]

  print(paste0("Your <query> contained ", length(query_ID), " UniProt_IDs and ", length(query_Accession), " UniProt_Accession (some might be redundant)." ))

  #if no values in the UniProtTable stop
  if (length(query_ID)==0&&length(query_Accession)==0){stop("Your <query> should be Uniprot_IDs or a Uniprot_Accessions, ensure that the organism of your query is present in the UniProtTable.")}

  #convert uniprot_Accession in Uniprot
  if(length(query_Accession)>0){
    print("The uniprot_Accession of your query were converted in Uniprot_IDs.")
    query_Accession<- as.character(t(UniProtTable[match(query_Accession,UniProtTable$Uniprot_Accession),2]))
    query_ID<-c(query_ID,query_Accession)
  }

  #remove duplicates
  duplicated_ID<-duplicated(query_ID)
  if(sum(duplicated_ID)>0){
    print(paste0("Your list of Uniprot_IDs was containing ", sum(duplicated_ID)," duplicated value(s), duplicates were removed, the query now comprises ", length(query_ID)-sum(duplicated_ID), " unique IDs."))
    query_ID<-query_ID[!duplicated_ID]
  }

  #check universe
  if(missing(universe)){
    print("Your <universe> was missing the global list of ID present in UniProtTable was used as background.")
    universe<-as.character(t(UniProtTable$Uniprot_ID))
  }

  #check how many of the universe are in the lists UniProtTable
  universe_ID <- universe[universe %in% UniProtTable$Uniprot_ID]
  universe_Accession <- universe[universe %in% UniProtTable$Uniprot_Accession]

  print(paste0("Your universe contained ", length(universe_ID), " UniProt_IDs and ", length(universe_Accession), " UniProt_Accession (some might be redundant)." ))

  #if no values in the UniProtTable stop
  if (length(universe_ID)==0&&length(universe_Accession)==0){stop("Your <universe> should be Uniprot_IDs or a Uniprot_Accessions, ensure that the organism of your universe is present in the UniProtTable.")}

  #convert uniprot_Accession in Uniprot
  if(length(universe_Accession)>0){
    print("The uniprot_Accession of your universe were converted in Uniprot_IDs.")
    universe_Accession<- as.character(t(UniProtTable[match(universe_Accession,UniProtTable$Uniprot_Accession),2]))
    universe_ID<-c(universe_ID,universe_Accession)
  }

  #remove duplicates
  duplicated_ID<-duplicated(universe_ID)
  if(sum(duplicated_ID)>0){
    print(paste0("Your list of Uniprot_IDs was containing ", sum(duplicated_ID)," duplicated value(s), duplicates were removed, the universe now comprises ", length(universe_ID)-sum(duplicated_ID), " unique IDs."))
    universe_ID<-universe_ID[!duplicated_ID]
  }

  #check if the universe contains all the proteins from the query
  if(sum(is.na(match(query_ID,universe_ID)))>0){stop("Your <universe> does not comprises all the IDs present in your <query>")}

  #generate the UniProtTable_GO for the query and the universe
  UPT_GO_query <- UPT_GO[UPT_GO$Uniprot_ID %in% query_ID,]
  UPT_GO_universe <-UPT_GO[UPT_GO$Uniprot_ID %in% universe_ID,]

  #create the unique_GO list based on the UPT_GO_query
  unique_GO<- unique(UPT_GO_query[,3:4])

  #create the Contingencies matrix
  contingency_list<-list()
  for (i in 1:nrow(unique_GO)){
    contingency_list[[i]]<-matrix(ncol=2,nrow=2)
    contingency_list[[i]][1,1]<- sum(UPT_GO_query[,3]==unique_GO[i,1])-1
    contingency_list[[i]][1,2]<- sum(UPT_GO_universe[,3]==unique_GO[i,1])
    contingency_list[[i]][2,1]<- length(query_ID)-contingency_list[[i]][1,1]
    contingency_list[[i]][2,2]<- length(universe)-contingency_list[[i]][1,2]
  }
  names(contingency_list)<-unique_GO[,1]

  #create Final table
  EASE_results<-data.frame(matrix(NA, nrow=nrow(unique_GO), ncol=10))
  colnames(EASE_results)<- c("Ontology_type","Ontology_accession","Ongology_description","Count_universe","%_query","%_universe","pval","adjpval","fold_change","Proteins_in_query")
  EASE_results[,2:3]<-unique_GO[,1:2]

  #run the test and populate the table
  for (i in 1:nrow(unique_GO)){
    test<-fisher.test(contingency_list[[i]])
    EASE_results$`Ontology_type`[i]<-"GO"
    EASE_results$pval[i]<-test$p.value
    EASE_results$adjpval[i]<-p.adjust(EASE_results$pval[i] , method = "BH", nrow(unique_GO))
    EASE_results$Count_query[i]<- paste0((contingency_list[[i]][1,1]+1),"/",length(query_ID))
    EASE_results$Count_universe[i]<- paste0(contingency_list[[i]][1,2],"/",length(universe_ID))
    EASE_results$`%_query` [i]<- contingency_list[[i]][1,1]/length(query_ID)*100
    EASE_results$`%_universe`[i]<- contingency_list[[i]][1,2]/length(universe_ID)*100
    EASE_results$fold_change[i]<- EASE_results$`%_query`[i]/EASE_results$`%_universe`[i]
    EASE_results$Proteins_in_query[i]<-paste(UPT_GO_query[UPT_GO_query$GO_accession==EASE_results$Ontology_accession[i],]$Uniprot_ID,collapse =";")
  }

  EASE_results[order(EASE_results$pval),]


}

#' UniProt_GO_KS()
#' @description Allows to use UniprotDB Gene Ontology information to perform a modified FisherExact test corresponding to the one performed by the DAVID Bioinformatics tool for enrichment analysis
#' @param rankingTable should be a two column data.frame containing the Uniprot identifiers list as first column (UniProt_Accession or UniProt_ID) and the ranking values as second column
#' @param order should either be "ascending" or "descending" to indicate the order of the ranking values to use.
#' @param organismID Should be a UniProt organism ID, by default Homo sapiens will be used ("9606")
#' @param reviewed this indicate if the proteins extracted have to be reviewed or not (typically the ones present on swissProt), by default TRUE (the non reviewed proteins are discarded)
#' @param preloaded_UniProtTable this indicate wether or not the preloaded UniProtTable and UniProtTable_GO should be used, if FALSE the UniProtTable and UniProtTable_GO will be generated again
#'
#' @return a data.frame containing enrichment results from the GO retrieve from UniProt.
#'
#' @examples UniProt_GO_EASE(QueryExample, UniverseExample, organism.id=organismExample)
#'

UniProt_GO_KS<-function(rankingTable, order= "ascending", organismID="9606", reviewed=TRUE, preloaded_UniProtTable=TRUE){
  #general checkings
  rankingTable<-data.frame(rankingTable)
  if(ncol(rankingTable)!=2){stop("df should contain two columns")}
  rankingTable[rankingTable==""]<-NA
  if(sum(is.na(rankingTable))>0){stop("Your ranking table should not contain missing values")}
  if(missing(order)){order="ascending"}
  if(!order %in% c("ascending", "descending")){stop("order should either be 'ascending' or 'descending'")}
  if(missing(organismID)){organismID<-"9606"}
  if(missing(reviewed)){reviewed<-TRUE}
  if(!is.logical(reviewed)){stop("<reviewed> should be either TRUE or FALSE")}
  if(missing(preloaded_UniProtTable)){preloaded_UniProtTable<-TRUE}
  if(!is.logical(preloaded_UniProtTable)){stop("<preloaded_UniProtTable> should be either TRUE or FALSE")}

  #load the UniProtTable and UniProtTable_GO (not if it is preloaded)
  if(preloaded_UniProtTable==TRUE&&exists("UniProtTable")){
    UPT<- UniProtTable
  }else {
    download_UniProtTable(organismID=organismID,reviewed=reviewed)
    UPT<- UniProtTable
  }

  if(preloaded_UniProtTable==TRUE&&exists("UniProtTable_GO")){
    UPT_GO<- UniProtTable_GO
  }else{
    generate_UniProtTable_GO(organismID=organismID,reviewed=reviewed,preloaded_UniProtTable==TRUE)
    UPT_GO<- UniProtTable_GO
    }

  #order based on the ascending/descending order parameter
  if(order=="ascending"){
    rankingTable<-rankingTable[order(rankingTable[,2]),]
  }else{
      rankingTable<-rankingTable[order(rankingTable[,2],decreasing=TRUE),]
      }

  #add a KS column in this table
  rankingTable$KSRank<-1:nrow(rankingTable)

  #check how many of the rankingTable are in the lists UniProtTable
  rankingTable_ID <- rankingTable[rankingTable[,1] %in% UniProtTable$Uniprot_ID,1]
  rankingTable_Accession <- rankingTable[rankingTable[,1] %in% UniProtTable$Uniprot_Accession,1]

  print(paste0("Your <rankingTable> contained ", length(rankingTable_ID), " UniProt_IDs and ", length(rankingTable_Accession), " UniProt_Accession (some might be redundant)." ))

  #if no values in the UniProtTable stop
  if (length(rankingTable_ID)==0&&length(rankingTable_Accession)==0){stop("Your <query> should be Uniprot_IDs or a Uniprot_Accessions, ensure that the organism of your query is present in the UniProtTable.")}

  #convert uniprot_Accession in Uniprot_ID
  if(length(rankingTable_Accession)>0){
    print("The uniprot_Accession of your rankingTable were converted in Uniprot_IDs.")
    rankingTable_Accession<- as.character(t(UniProtTable[match(rankingTable_Accession,UniProtTable$Uniprot_Accession),2]))
    rankingTable_ID<-c(rankingTable_ID,rankingTable_Accession)
  }

  #remove duplicates
  duplicated_ID<-duplicated(rankingTable[,1])
  if(sum(duplicated_ID)>0){
    print(paste0("Your list of Uniprot_IDs was containing ", sum(duplicated_ID)," duplicated value(s), duplicates were removed, the rankingTable now comprises ", length(rankingTable_ID)-sum(duplicated_ID), " unique IDs."))
    rankingTable<-rankingTable[!duplicated_ID,]
    rankingTable$KSRank<-1:nrow(rankingTable)
    }

  #generate the UniProtTable_GO for the rankingTable
  UPT_GO_rankingTable <- UPT_GO[UPT_GO$Uniprot_ID %in% rankingTable_ID,]

  #create the unique_GO list based on the UPT_GO_query
  unique_GO<- unique(UPT_GO_rankingTable[,3:4])

  #create a list of sub RankingTables
  RT_by_GO<-list()
  IDs_by_GO<-list()

  for (i in 1:nrow(unique_GO)){
    IDs<-UPT_GO_rankingTable[unique_GO$GO_accession[i]==UPT_GO_rankingTable$GO_accession,]
    IDs_by_GO[[i]]<-paste(IDs$Uniprot_ID,collapse=";")
    IDs<-IDs$Uniprot_ID
    Rank<-rankingTable[rankingTable[,1] %in% IDs,]
    RT_by_GO[[i]]<-Rank
    }

  names(RT_by_GO)<-paste0(unique_GO$GO_accession,"@",unique_GO$GO_description)

  #realise the KS.test
  p<-numeric()
  for (i in 1:nrow(unique_GO)){
    if(nrow(rankingTable)>length(RT_by_GO[[i]])){
      p[i]<-ks.test(RT_by_GO[[i]]$KSRank,seq_len(nrow(rankingTable))[-RT_by_GO[[i]]$KSRank],alternative="greater")$p.value}else{p[i]<-1}
  }

  #calculate the adjusted pvalues
  adjp<-p.adjust(p)

  #create the final table
  final_table<-data.frame(matrix(ncol=4,nrow= nrow(unique_GO)))
  colnames(final_table) <- c("Count.in.list", "p.value", "adj.p","IDs.in.list")
  final_table<-cbind(unique_GO,final_table)
  for(i in 1:nrow(unique_GO)){
    final_table$Count.in.list[i]<-nrow(RT_by_GO[[i]])
    final_table$IDs.in.list[i]<-IDs_by_GO[[i]]
    }
  final_table$`p.value`<-p
  final_table$adj.p<-adjp

  return(final_table)

}

#' UniProt_KEGG_Fisher()
#' @description Allows to use UniprotDB and KEGG db information to perform a FisherExact test for enrichment analysis
#' @param query should be a character vector corresponding to a query list (UniProt_Accession or UniProt_ID)
#' @param universe should be a character vector corresponding to a universe list (Uniprot Entry), if missing the all table from the UniProtTable will be used as universe
#' @param organismID Should be a UniProt organism ID, by default Homo sapiens will be used ("9606")
#' @param reviewed this indicate if the proteins extracted have to be reviewed or not (typically the ones present on swissProt), by default TRUE (the non reviewed proteins are discarded)
#' @param preloaded_UniProtTable this indicate wether or not the preloaded UniProtTable and UniProtTable_GO should be used, if FALSE the UniProtTable and UniProtTable_GO will be generated again
#'
#' @return a data.frame containing enrichment results
#'
#' @examples UniProtGOfisher(QueryExample, UniverseExample, organism.id=organismExample)
#'
UniProt_KEGG_Fisher<-function(query, universe, organismID="9606", reviewed=TRUE, preloaded_UniProtTable=TRUE){
  #general checkings
  if(missing(organismID)){organismID<-"9606"}
  if(missing(reviewed)){reviewed<-TRUE}
  if(!is.logical(reviewed)){stop("<reviewed> should be either TRUE or FALSE")}
  if(missing(preloaded_UniProtTable)){preloaded_UniProtTable<-TRUE}
  if(!is.logical(preloaded_UniProtTable)){stop("<preloaded_UniProtTable> should be either TRUE or FALSE")}

  #load the UniProtTable and UniProtTable_GO (not if it is preloaded)
  if(preloaded_UniProtTable==TRUE&&exists("UniProtTable")){
    UPT<- UniProtTable
  }else {
    download_UniProtTable(organismID=organismID,reviewed=reviewed)
    UPT<- UniProtTable
    }

  if(preloaded_UniProtTable==TRUE&&exists("UniProtTable_KEGG")){
    UPT_KEGG<- UniProtTable_KEGG
  }else {
    generate_UniProtTable_KEGG(organismID=organismID,reviewed=reviewed,preloaded_UniProtTable==TRUE)
    UPT_KEGG<- UniProtTable_KEGG
    }

  #check query
  if(missing(query)){stop("Your <query> is missing")}

  #check how many of the query are in the lists UniProtTable
  query_ID <- query[query %in% UniProtTable$Uniprot_ID]
  query_Accession <- query[query %in% UniProtTable$Uniprot_Accession]

  print(paste0("Your <query> contained ", length(query_ID), " UniProt_IDs and ", length(query_Accession), " UniProt_Accession (some might be redundant)." ))

  #if no values in the UniProtTable stop
  if (length(query_ID)==0&&length(query_Accession)==0){stop("Your <query> should be Uniprot_IDs or a Uniprot_Accessions, ensure that the organism of your query is present in the UniProtTable.")}

  #convert uniprot_Accession in Uniprot
  if(length(query_Accession)>0){
    print("The uniprot_Accession of your query were converted in Uniprot_IDs.")
    query_Accession<- as.character(t(UniProtTable[match(query_Accession,UniProtTable$Uniprot_Accession),2]))
    query_ID<-c(query_ID,query_Accession)
  }

  #remove duplicates
  duplicated_ID<-duplicated(query_ID)
  if(sum(duplicated_ID)>0){
    print(paste0("Your list of Uniprot_IDs was containing ", sum(duplicated_ID)," duplicated value(s), duplicates were removed, the query now comprises ", length(query_ID)-sum(duplicated_ID), " unique IDs."))
    query_ID<-query_ID[!duplicated_ID]
  }

  #check universe
  if(missing(universe)){
    print("Your <universe> was missing the global list of ID present in UniProtTable was used as background.")
    universe<-as.character(t(UniProtTable$Uniprot_ID))
  }

  #check how many of the universe are in the lists UniProtTable
  universe_ID <- universe[universe %in% UniProtTable$Uniprot_ID]
  universe_Accession <- universe[universe %in% UniProtTable$Uniprot_Accession]

  print(paste0("Your universe contained ", length(universe_ID), " UniProt_IDs and ", length(universe_Accession), " UniProt_Accession (some might be redundant)." ))

  #if no values in the UniProtTable stop
  if (length(universe_ID)==0&&length(universe_Accession)==0){stop("Your <universe> should be Uniprot_IDs or a Uniprot_Accessions, ensure that the organism of your universe is present in the UniProtTable.")}

  #convert uniprot_Accession in Uniprot
  if(length(universe_Accession)>0){
    print("The uniprot_Accession of your universe were converted in Uniprot_IDs.")
    universe_Accession<- as.character(t(UniProtTable[match(universe_Accession,UniProtTable$Uniprot_Accession),2]))
    universe_ID<-c(universe_ID,universe_Accession)
  }

  #remove duplicates
  duplicated_ID<-duplicated(universe_ID)
  if(sum(duplicated_ID)>0){
    print(paste0("Your list of Uniprot_IDs was containing ", sum(duplicated_ID)," duplicated value(s), duplicates were removed, the universe now comprises ", length(universe_ID)-sum(duplicated_ID), " unique IDs."))
    universe_ID<-universe_ID[!duplicated_ID]
  }

  #check if the universe contains all the proteins from the query
  if(sum(is.na(match(query_ID,universe_ID)))>0){stop("Your <universe> does not comprises all the IDs present in your <query>")}

  #generate the UniProtTable_KEGG for the query and the universe
  UPT_KEGG_query <- UPT_KEGG[UPT_KEGG$Uniprot_ID %in% query_ID,]
  UPT_KEGG_universe <-UPT_KEGG[UPT_KEGG$Uniprot_ID %in% universe_ID,]

  #create the unique_GO list based on the UPT_GO_query
  unique_PATH<- unique(UPT_KEGG_query[,4:5])

  #create the Contingencies matrix
  contingency_list<-list()
  for (i in 1:nrow(unique_PATH)){
    contingency_list[[i]]<-matrix(ncol=2,nrow=2)
    contingency_list[[i]][1,1]<- sum(UPT_KEGG_query[,4]==unique_PATH[i,1])
    contingency_list[[i]][1,2]<- sum(UPT_KEGG_universe[,4]==unique_PATH[i,1])
    contingency_list[[i]][2,1]<- length(query_ID)-contingency_list[[i]][1,1]
    contingency_list[[i]][2,2]<- length(universe)-contingency_list[[i]][1,2]
  }
  names(contingency_list)<-unique_PATH[,1]

  #create Final table
  fisher_results<-data.frame(matrix(NA, nrow=nrow(unique_PATH), ncol=11))
  colnames(fisher_results)<- c("Ontology_type","Ontology_accession","Ongology_description","Count_query","Count_universe","%_query","%_universe","pval","adjpval","fold_change","Proteins_in_query")
  fisher_results[,2:3]<-unique_PATH[,1:2]

  #run the test and populate the table
  for (i in 1:nrow(unique_PATH)){
    fisher_results$`Ontology_type`[i]<-"KEGG"
    test<-fisher.test(contingency_list[[i]])
    fisher_results$pval[i]<-test$p.value
    fisher_results$adjpval[i]<-p.adjust(fisher_results$pval[i] , method = "BH", nrow(unique_PATH))
    fisher_results$Count_query[i]<- paste0(contingency_list[[i]][1,1],"/",length(query_ID))
    fisher_results$Count_universe[i]<- paste0(contingency_list[[i]][1,2],"/",length(universe_ID))
    fisher_results$`%_query` [i]<- contingency_list[[i]][1,1]/length(query_ID)*100
    fisher_results$`%_universe`[i]<- contingency_list[[i]][1,2]/length(universe_ID)*100
    fisher_results$fold_change[i]<- fisher_results$`%_query`[i]/fisher_results$`%_universe`[i]
    fisher_results$Proteins_in_query[i]<-paste(UPT_KEGG_query[UPT_KEGG_query$Pathway_ID==fisher_results$Ontology_accession[i],]$Uniprot_ID,collapse =";")
  }

  fisher_results[order(fisher_results$pval),]
}

#' UniProt_KEGG_Fisher()
#' @description Allows to use UniprotDB and KEGG db information to perform a modified FisherExact test corresponding to the one performed by the DAVID Bioinformatics tool for enrichment analysis.
#' @param query should be a character vector corresponding to a query list (UniProt_Accession or UniProt_ID)
#' @param universe should be a character vector corresponding to a universe list (Uniprot Entry), if missing the all table from the UniProtTable will be used as universe
#' @param organismID Should be a UniProt organism ID, by default Homo sapiens will be used ("9606")
#' @param reviewed this indicate if the proteins extracted have to be reviewed or not (typically the ones present on swissProt), by default TRUE (the non reviewed proteins are discarded)
#' @param preloaded_UniProtTable this indicate wether or not the preloaded UniProtTable and UniProtTable_GO should be used, if FALSE the UniProtTable and UniProtTable_GO will be generated again
#'
#' @return a data.frame containing enrichment results
#'
#' @examples UniProtGOfisher(QueryExample, UniverseExample, organism.id=organismExample)
#'

UniProt_KEGG_EASE<-function(query, universe, organismID="9606", reviewed=TRUE, preloaded_UniProtTable=TRUE){
  #general checkings
  if(missing(organismID)){organismID<-"9606"}
  if(missing(reviewed)){reviewed<-TRUE}
  if(!is.logical(reviewed)){stop("<reviewed> should be either TRUE or FALSE")}
  if(missing(preloaded_UniProtTable)){preloaded_UniProtTable<-TRUE}
  if(!is.logical(preloaded_UniProtTable)){stop("<preloaded_UniProtTable> should be either TRUE or FALSE")}

  #load the UniProtTable and UniProtTable_GO (not if it is preloaded)
  if(preloaded_UniProtTable==TRUE&&exists("UniProtTable")){
    UPT<- UniProtTable
  }else {download_UniProtTable(organismID=organismID,reviewed=reviewed)}

  if(preloaded_UniProtTable==TRUE&&exists("UniProtTable_KEGG")){
    UPT_KEGG<- UniProtTable_KEGG
  }else {generate_UniProtTable_KEGG(organismID=organismID,reviewed=reviewed,preloaded_UniProtTable==TRUE)}

  #check query
  if(missing(query)){stop("Your <query> is missing")}

  #check how many of the query are in the lists UniProtTable
  query_ID <- query[query %in% UniProtTable$Uniprot_ID]
  query_Accession <- query[query %in% UniProtTable$Uniprot_Accession]

  print(paste0("Your <query> contained ", length(query_ID), " UniProt_IDs and ", length(query_Accession), " UniProt_Accession (some might be redundant)." ))

  #if no values in the UniProtTable stop
  if (length(query_ID)==0&&length(query_Accession)==0){stop("Your <query> should be Uniprot_IDs or a Uniprot_Accessions, ensure that the organism of your query is present in the UniProtTable.")}

  #convert uniprot_Accession in Uniprot
  if(length(query_Accession)>0){
    print("The uniprot_Accession of your query were converted in Uniprot_IDs.")
    query_Accession<- as.character(t(UniProtTable[match(query_Accession,UniProtTable$Uniprot_Accession),2]))
    query_ID<-c(query_ID,query_Accession)
  }

  #remove duplicates
  duplicated_ID<-duplicated(query_ID)
  if(sum(duplicated_ID)>0){
    print(paste0("Your list of Uniprot_IDs was containing ", sum(duplicated_ID)," duplicated value(s), duplicates were removed, the query now comprises ", length(query_ID)-sum(duplicated_ID), " unique IDs."))
    query_ID<-query_ID[!duplicated_ID]
  }

  #check universe
  if(missing(universe)){
    print("Your <universe> was missing the global list of ID present in UniProtTable was used as background.")
    universe<-as.character(t(UniProtTable$Uniprot_ID))
  }

  #check how many of the universe are in the lists UniProtTable
  universe_ID <- universe[universe %in% UniProtTable$Uniprot_ID]
  universe_Accession <- universe[universe %in% UniProtTable$Uniprot_Accession]

  print(paste0("Your universe contained ", length(universe_ID), " UniProt_IDs and ", length(universe_Accession), " UniProt_Accession (some might be redundant)." ))

  #if no values in the UniProtTable stop
  if (length(universe_ID)==0&&length(universe_Accession)==0){stop("Your <universe> should be Uniprot_IDs or a Uniprot_Accessions, ensure that the organism of your universe is present in the UniProtTable.")}

  #convert uniprot_Accession in Uniprot
  if(length(universe_Accession)>0){
    print("The uniprot_Accession of your universe were converted in Uniprot_IDs.")
    universe_Accession<- as.character(t(UniProtTable[match(universe_Accession,UniProtTable$Uniprot_Accession),2]))
    universe_ID<-c(universe_ID,universe_Accession)
  }

  #remove duplicates
  duplicated_ID<-duplicated(universe_ID)
  if(sum(duplicated_ID)>0){
    print(paste0("Your list of Uniprot_IDs was containing ", sum(duplicated_ID)," duplicated value(s), duplicates were removed, the universe now comprises ", length(universe_ID)-sum(duplicated_ID), " unique IDs."))
    universe_ID<-universe_ID[!duplicated_ID]
  }

  #check if the universe contains all the proteins from the query
  if(sum(is.na(match(query_ID,universe_ID)))>0){stop("Your <universe> does not comprises all the IDs present in your <query>")}

  #generate the UniProtTable_KEGG for the query and the universe
  UPT_KEGG_query <- UPT_KEGG[UPT_KEGG$Uniprot_ID %in% query_ID,]
  UPT_KEGG_universe <-UPT_KEGG[UPT_KEGG$Uniprot_ID %in% universe_ID,]

  #create the unique_GO list based on the UPT_GO_query
  unique_PATH<- unique(UPT_KEGG_query[,4:5])

  #create the Contingencies matrix
  contingency_list<-list()
  for (i in 1:nrow(unique_PATH)){
    contingency_list[[i]]<-matrix(ncol=2,nrow=2)
    contingency_list[[i]][1,1]<- sum(UPT_KEGG_query[,4]==unique_PATH[i,1])-1
    contingency_list[[i]][1,2]<- sum(UPT_KEGG_universe[,4]==unique_PATH[i,1])
    contingency_list[[i]][2,1]<- length(query_ID)-contingency_list[[i]][1,1]
    contingency_list[[i]][2,2]<- length(universe)-contingency_list[[i]][1,2]
  }
  names(contingency_list)<-unique_PATH[,1]

  #create Final table
  fisher_results<-data.frame(matrix(NA, nrow=nrow(unique_PATH), ncol=11))
  colnames(fisher_results)<- c("Ontology_type","Ontology_accession","Ongology_description","Count_query","Count_universe","%_query","%_universe","pval","adjpval","fold_change","Proteins_in_query")
  fisher_results[,2:3]<-unique_PATH[,1:2]

  #run the test and populate the table
  for (i in 1:nrow(unique_PATH)){
    fisher_results$`Ontology_type`[i]<-"KEGG"
    test<-fisher.test(contingency_list[[i]])
    fisher_results$pval[i]<-test$p.value
    fisher_results$adjpval[i]<-p.adjust(fisher_results$pval[i] , method = "BH", nrow(unique_PATH))
    fisher_results$Count_query[i]<- paste0(contingency_list[[i]][1,1],"/",length(query_ID))
    fisher_results$Count_universe[i]<- paste0(contingency_list[[i]][1,2],"/",length(universe_ID))
    fisher_results$`%_query` [i]<- contingency_list[[i]][1,1]/length(query_ID)*100
    fisher_results$`%_universe`[i]<- contingency_list[[i]][1,2]/length(universe_ID)*100
    fisher_results$fold_change[i]<- fisher_results$`%_query`[i]/fisher_results$`%_universe`[i]
    fisher_results$Proteins_in_query[i]<-paste(UPT_KEGG_query[UPT_KEGG_query$Pathway_ID==fisher_results$Ontology_accession[i],]$Uniprot_ID,collapse =";")
  }

  fisher_results[order(fisher_results$pval),]
}



#' CustomDB_Fisher()
#' @description Allows to perform enrichment analysis from a custom database using a query list, a universe list and an annotation table (containing the identifiers and associated ontology terms)
#' @param query should be a character vector corresponding to a query list (UniProt_Accession or UniProt_ID)
#' @param universe should be a character vector corresponding to a universe list (Uniprot Entry), if missing the all table from the UniProtTable will be used as universe
#' @param custom_db should be a data.frame containing at least 2 columns one with the identifiers matching to the query/universe identifiers and one containing the associated ontology terms/groups. For each ID multiple ontology terms can be associated, multiple rows containing the same ID and different ontology terms should be present.
#' @param id_colname character vector corresponding to the colname of the column to use for the ids
#' @param ontology_colname
#'
#' @return a data.frame containing enrichment results
#'
#' @examples UniProtGOfisher(QueryExample, UniverseExample, organism.id=organismExample)
#'
CustomDB_Fisher<-function(query, universe, custom_db, id_colname, ontology_colname){
  #general checkings
  if(missing(query)){stop("Your <query> is missing.")}
  if(missing(universe)){stop("Your <universe> is missing.")}
  if(missing(custom_db)){stop("Your <custom_db> is missing.")}
  if(!is.character(query)){stop("Your <query> should be a character vector")}
  if(!is.character(universe)){stop("Your <universe> should be a character vector")}
  if(sum(query %in% universe)!= length(query)){stop("Not all the IDs present in the <query> are in the <universe> list")}
  if(missing(id_colname)){
    warning("Your <id_colname> is missing the first column will be used as id column")
    id_colname<-colnames(custom_db)[1]
    }
  if(missing(ontology_colname)){
    warning("Your <ontology_colname> is missing the second column will be used as id column")
    ontology_colname<-colnames(custom_db)[2]
  }

  #remove duplicates
  duplicated_query<-duplicated(query)
  if(sum(duplicated_query)>0){
    print(paste0("Your <query> list  was containing ", sum(duplicated_query)," duplicated IDs, duplicates were removed, the query now comprises ", length(query_ID)-sum(duplicated_ID), " unique IDs."))
    query<-query[!duplicated_query]
  }

  duplicated_universe<-duplicated(universe)
  if(sum(duplicated_universe)>0){
    print(paste0("Your <universe> list  was containing ", sum(duplicated_ID)," duplicated IDs, duplicates were removed, the query now comprises ", length(query_ID)-sum(duplicated_ID), " unique IDs."))
    universe<-universe[!duplicated_universe]
  }

  duplicated_ontologies<-duplicated(paste0(custom_db[,colnames(custom_db)==id_colname],"@",custom_db[,colnames(custom_db)==ontology_colname]))
  if(sum(duplicated_ontologies)>0){
    print(paste0("Your <custom_db> list  was containing ", sum(duplicated_ontologies)," duplicated IDs/ontology pairs. Duplicates were removed, it now comprises ", nrow(custom_db)-sum(duplicated_ontologies), " unique IDs."))
    custom_db<-custom_db[!duplicated_ontologies,]
  }


  #check how many of the universe IDs are in the custom_db
  universe_in_custom_db <- universe[universe %in% custom_db[,colnames(custom_db)==id_colname]]


  #generate the UniProtTable_GO for the query and the universe
  ontology_query <- custom_db[custom_db[,colnames(custom_db)==id_colname] %in% query,]
  ontology_universe <- custom_db[custom_db[,colnames(custom_db)==id_colname] %in% universe,]

  #create the unique_GO list based on the UPT_GO_query
  unique_ontologies<- unique(as.character(t(ontology_query[,colnames(ontology_query)==ontology_colname])))

  #create the Contingencies matrix
  contingency_list<-list()
  for (i in 1:length(unique_ontologies)){
    contingency_list[[i]]<-matrix(ncol=2,nrow=2)
    contingency_list[[i]][1,1]<- sum(as.character(t(ontology_query[,colnames(ontology_query)==ontology_colname]))==unique_ontologies[i])
    contingency_list[[i]][1,2]<- sum(as.character(t(ontology_universe[,colnames(ontology_universe)==ontology_colname]))==unique_ontologies[i])
    contingency_list[[i]][2,1]<- length(query)-contingency_list[[i]][1,1]
    contingency_list[[i]][2,2]<- length(universe)-contingency_list[[i]][1,2]
  }
  names(contingency_list)<-unique_ontologies

  #create Final table
  fisher_results<-data.frame(matrix(NA, nrow=length(unique_ontologies), ncol=9))
  colnames(fisher_results)<- c("Ontology_accession","Count_query","Count_universe","%_query","%_universe","pval","adjpval","fold_change","IDs_in_query")
  fisher_results[,1]<-unique_ontologies

  #run the test and populate the table
  for (i in 1:nrow(fisher_results)){
    test<-fisher.test(contingency_list[[i]])
    fisher_results$pval[i]<-test$p.value
    fisher_results$adjpval[i]<-p.adjust(fisher_results$pval[i] , method = "BH", length(unique_ontologies))
    fisher_results$Count_query[i]<- paste0(contingency_list[[i]][1,1],"/",length(query))
    fisher_results$Count_universe[i]<- paste0(contingency_list[[i]][1,2],"/",length(universe))
    fisher_results$`%_query` [i]<- contingency_list[[i]][1,1]/length(query)*100
    fisher_results$`%_universe`[i]<- contingency_list[[i]][1,2]/length(universe)*100
    fisher_results$fold_change[i]<- fisher_results$`%_query`[i]/fisher_results$`%_universe`[i]
    fisher_results$IDs_in_query[i]<-paste(ontology_query[ontology_query[,colnames(ontology_query)==ontology_colname]==unique_ontologies[i],colnames(ontology_query)==id_colname],collapse =";")
  }

  fisher_results<-fisher_results[order(fisher_results$pval),]
}

