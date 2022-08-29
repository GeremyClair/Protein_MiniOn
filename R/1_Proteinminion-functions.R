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
  if(!is.logical(reviewed)){warning("<reviewed> should be either TRUE or FALSE, by default reviewed=FALSE was used")
    reviewed<- FALSE}
  if(missing(export)){export<-FALSE}
  if(!is.logical(reviewed)){stop("<reviewed> should be either TRUE or FALSE")}
  if(missing(file)){file=paste0("UniProtFasta_OrganismID_",organismID,"_",format(Sys.time(),"%b_%d_%Y"),".fasta")}

  URL<-"https://rest.uniprot.org/uniprotkb/stream?compressed=false"

  #query type
  URL<-paste0(URL,"&fields=accession%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Clength%2Cgo_p%2Cgo_c%2Cgo_f%2Cgo%2Cgo_id%2Cxref_kegg%2Cxref_reactome")
  URL<-paste0(URL,"&format=tsv&query=")

  # organism ID or proteome ID missing
  if(missing(organismID)&missing(proteomeID)){
    warning("Both the OrganismID and the ProteomeID were missing, the human proteome was used")
    URL<-paste0(URL,"%28","organism_id%3A","9606","%29")
  }else{
    if(missing(organismID)==FALSE & missing(proteomeID)==FALSE){
      warning("Both the OrganismID and the ProteomeID were filled, ONLY the OrganismID will be used")
      URL<-paste0(URL,"%28","organism_id%3A",organismID,"%29")
    }else{
      if(!missing(organismID)){URL<-paste0(URL,"%28","organism_id%3A",organismID,"%29")}
      if(!missing(proteomeID)){URL<-paste0(URL,"%28","proteome%3A",proteomeID,"%29")}
    }
  }

  #reviewed?
  if(reviewed==TRUE){URL<-paste0(URL,"%20","AND","%20","%28","reviewed","%3A","true","%29")}

  #read the table from the internet
  UniProtTable<-read.csv(URL,sep="\t",header=T)

  #Adjust column names
  colnames(UniProtTable)<-c("Uniprot_Accession", "Uniprot_ID","Description", "Gene_Name", "Organism", "Length", "GO_BP", "GO_CC","GO_MF","GO_List","GO_IDs","Crossref_KEGG","Reactome_List")

  #export the global environment
  UniProtTable<<-UniProtTable

  if(export==TRUE){
    save(UniProtTable,file = file)
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

  URL<-"https://rest.uniprot.org/uniprotkb/stream?compressed=false"

  #query type
  URL<-paste0(URL,"&format=fasta&query=")

  # organism ID or proteome ID missing
  if(missing(organismID)&missing(proteomeID)){
    warning("Both the OrganismID and the ProteomeID were missing, the human proteome was used")
    URL<-paste0(URL,"%28","organism_id%3A","9606","%29")
  }else{
    if(missing(organismID)==FALSE & missing(proteomeID)==FALSE){
    warning("Both the OrganismID and the ProteomeID were filled, ONLY the OrganismID will be used")
      URL<-paste0(URL,"%28","organism_id%3A",organismID,"%29")
    }else{
      if(!missing(organismID)){URL<-paste0(URL,"%28","organism_id%3A",organismID,"%29")}
      if(!missing(proteomeID)){URL<-paste0(URL,"%28","proteome%3A",proteomeID,"%29")}
    }
  }

  #reviewed?
  if(reviewed==TRUE){URL<-paste0(URL,"%20","AND","%20","%28","reviewed","%3A","true","%29")}

  #read the table from the internet
  f<-read.csv(URL,sep="\t",header=F)

  #print your fasta file number of sequences
  print(paste("The UniProtFasta generated contains",sum(grepl(">",as.character(t(f)))),"sequences"))

  #export the global environment
  UniProtFasta<<-f

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
  if(missing(preloaded_UniProtTable)){preloaded_UniProtTable<-TRUE}
  if(!is.logical(preloaded_UniProtTable)){stop("<preloaded_UniProtTable> should be either TRUE or FALSE")}
  if(preloaded_UniProtTable==TRUE&&exists("UniProtTable")){
    UPT<- UniProtTable
  }else{
    #general checking
    if(missing(organismID)&missing(proteomeID)){warning("The 'organismID' and 'proteomeID' were left empty, 'homo sapiens' (9006) will be used by default")
      organismID<-"9606"}
    if(missing(reviewed)){reviewed<-FALSE}
    if(!is.logical(reviewed)){warning("<reviewed> should be either TRUE or FALSE, by default reviewed=FALSE was used")
      reviewed<- FALSE}
    if(!is.logical(reviewed)){stop("<reviewed> should be either TRUE or FALSE")}
    #download the table
    download_UniProtTable(organismID=organismID, reviewed=reviewed)
    UPT<-UniProtTable
  }

  #extract the required columns from the UniProtTable (UPT) and ensure the content is considered as character by R
  UPT<-UPT[,colnames(UPT) %in% c("Uniprot_Accession","Uniprot_ID","GO_BP","GO_CC","GO_MF")]

  for(i in 1:5){
    UPT[,i]<-as.character(t(UPT[,i]))
  }

  #remove the end bracket
  UPT[,3]<- gsub("\\]","",UPT[,3])
  UPT[,4]<- gsub("\\]","",UPT[,4])
  UPT[,5]<- gsub("\\]","",UPT[,5])

  #create a list to store the different GO terms associated with each protein
  uniList_BP<-list()
  uniList_CC<-list()
  uniList_MF<-list()

  #create objects to store the Uniprot_accession and the Uniprot_ID
  Uniprot_Accession_BP<-character()
  Uniprot_Accession_CC<-character()
  Uniprot_Accession_MF<-character()

  Uniprot_ID_BP<-character()
  Uniprot_ID_CC<-character()
  Uniprot_ID_MF<-character()

  #fill the three objects created above
  for (i in 1:nrow(UPT)){
    uniList_BP[[i]]<-trimws(unlist(strsplit(UPT[i,3],";")))
    Uniprot_Accession_BP<-c(Uniprot_Accession_BP,rep(UPT[i,1],length(uniList_BP[[i]])))
    Uniprot_ID_BP<-c(Uniprot_ID_BP,rep(UPT[i,2],length(uniList_BP[[i]])))
  }

  for (i in 1:nrow(UPT)){
    uniList_CC[[i]]<-trimws(unlist(strsplit(UPT[i,4],";")))
    Uniprot_Accession_CC<-c(Uniprot_Accession_CC,rep(UPT[i,1],length(uniList_CC[[i]])))
    Uniprot_ID_CC<-c(Uniprot_ID_CC,rep(UPT[i,2],length(uniList_CC[[i]])))
  }

  for (i in 1:nrow(UPT)){
    uniList_MF[[i]]<-trimws(unlist(strsplit(UPT[i,5],";")))
    Uniprot_Accession_MF<-c(Uniprot_Accession_MF,rep(UPT[i,1],length(uniList_MF[[i]])))
    Uniprot_ID_MF<-c(Uniprot_ID_MF,rep(UPT[i,2],length(uniList_MF[[i]])))
  }

  #create the Data frames containing these details
  UniProtTable_GO_BP<-data.frame(Uniprot_Accession=Uniprot_Accession_BP, Uniprot_ID=Uniprot_ID_BP, GO=unlist(uniList_BP))
  UniProtTable_GO_BP<-cbind(UniProtTable_GO_BP,GO_type="GO_BP")
  UniProtTable_GO_BP <- UniProtTable_GO_BP %>% separate(GO, into= c('GO_description','GO_accession'),"\\[GO\\:" )

  UniProtTable_GO_CC<-data.frame(Uniprot_Accession=Uniprot_Accession_CC, Uniprot_ID=Uniprot_ID_CC, GO=unlist(uniList_CC))
  UniProtTable_GO_CC<-cbind(UniProtTable_GO_CC,GO_type="GO_CC")
  UniProtTable_GO_CC <- UniProtTable_GO_CC %>% separate(GO, into= c('GO_description','GO_accession'),"\\[GO\\:" )

  UniProtTable_GO_MF<-data.frame(Uniprot_Accession=Uniprot_Accession_MF, Uniprot_ID=Uniprot_ID_MF, GO=unlist(uniList_MF))
  UniProtTable_GO_MF<-cbind(UniProtTable_GO_MF,GO_type="GO_MF")
  UniProtTable_GO_MF <- UniProtTable_GO_MF %>% separate(GO, into= c('GO_description','GO_accession'),"\\[GO\\:" )

  UniProtTable_GO<-rbind(UniProtTable_GO_BP,UniProtTable_GO_CC,UniProtTable_GO_MF)


  #Place the reordered table in the Global Environment
  UniProtTable_GO <<- UniProtTable_GO[,c(1,2,4,3,5)]

  #export file if export==TRUE
  if(export==TRUE){
    save(UniProtTable_GO,file = file)
     }
  }

#' generate_UniProtTable_KEGG()
#' @description Retrieve UniProt KEGG pathways for the proteins of a given organism
#' @param organismID Should be a UniProt organism ID, by default Homo sapiens will be used ("9606")
#' @param reviewed this indicate if the proteins extracted have to be reviewed or not (typically the ones present on swissProt), by default TRUE (the non reviewed proteins are discarded)
#' @param preloaded_UniProtTable this indicate whether or not the preloaded UniProtTable and UniProtTable_GO should be used, if FALSE the UniProtTable and UniProtTable_GO will be generated again
#'
#' @return This function generate a KEGG enrichment table named UniProtTable_KEGG usable with the enrichment functions Uniprot_KEGG_*. If the export is set on, the UniProtTable_KEGG will also be saved in the .rda format at the specified location.
#'
#' @examples generate_UniProtTable_KEGG(organismID="9606",reviewed=TRUE,preloaded_UniProtTable=TRUE,export=FALSE,file="UniProtTable_KEGG_organismID_date.rda")
#'
generate_UniProtTable_KEGG<- function(organismID="9606",reviewed=TRUE,preloaded_UniProtTable=TRUE,export=FALSE,file="UniProtTable_KEGG_organismID_date.rda"){
  if(missing(preloaded_UniProtTable)){preloaded_UniProtTable<-TRUE}
  if(!is.logical(preloaded_UniProtTable)){stop("<preloaded_UniProtTable> should be either TRUE or FALSE")}

  if(preloaded_UniProtTable==TRUE&&exists("UniProtTable")){
    UPT<- UniProtTable
  }else{
    #general checking
    if(missing(organismID)&missing(proteomeID)){warning("The 'organismID' and 'proteomeID' were left empty, 'homo sapiens' (9006) will be used by default")
      organismID<-"9606"}
    if(missing(reviewed)){reviewed<-FALSE}
    if(!is.logical(reviewed)){warning("<reviewed> should be either TRUE or FALSE, by default reviewed=FALSE was used")
      reviewed<- FALSE}
    if(!is.logical(reviewed)){stop("<reviewed> should be either TRUE or FALSE")}
    #download the table
    download_UniProtTable(organismID=organismID, reviewed=reviewed)
    UPT<-UniProtTable
  }


  #extract the required columns from the UniProtTable (UPT) and ensure the content is considered as character by R
  UPT<-UPT[,colnames(UPT) %in% c("Uniprot_Accession","Uniprot_ID","Crossref_KEGG" )]

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

#' generate_UniProtTable_REACTOME()
#' @description Retrieve UniProt Reactome pathways for the proteins of a given organism
#' @param organismID Should be a UniProt organism ID, by default Homo sapiens will be used ("9606")
#' @param reviewed this indicate if the proteins extracted have to be reviewed or not (typically the ones present on swissProt), by default TRUE (the non reviewed proteins are discarded)
#' @param preloaded_UniProtTable this indicate wether or not the preloaded UniProtTable and UniProtTable_GO should be used, if FALSE the UniProtTable and UniProtTable_GO will be generated again
#'
#' @return This function generate a reactome table
#'
#' @examples generate_UniProtTable_REACTOME(organismID="9606",reviewed=TRUE,preloaded_UniProtTable=TRUE,export=FALSE,file="UniProtTable_KEGG_organismID_date.rda")
#'
generate_UniProtTable_REACTOME<-function(organismID="9606",reviewed=TRUE,preloaded_UniProtTable=TRUE,export=FALSE,file="UniProtTable_KEGG_organismID_date.rda"){
if(missing(preloaded_UniProtTable)){preloaded_UniProtTable<-TRUE}
  if(!is.logical(preloaded_UniProtTable)){stop("<preloaded_UniProtTable> should be either TRUE or FALSE")}

  if(preloaded_UniProtTable==TRUE&&exists("UniProtTable")){
    UPT<- UniProtTable
  }else{
    #general checking
    if(missing(organismID)&missing(proteomeID)){warning("The 'organismID' and 'proteomeID' were left empty, 'homo sapiens' (9006) will be used by default")
      organismID<-"9606"}
    if(missing(reviewed)){reviewed<-FALSE}
    if(!is.logical(reviewed)){warning("<reviewed> should be either TRUE or FALSE, by default reviewed=FALSE was used")
      reviewed<- FALSE}
    if(!is.logical(reviewed)){stop("<reviewed> should be either TRUE or FALSE")}
    #download the table
    download_UniProtTable(organismID=organismID, reviewed=reviewed)
    UPT<-UniProtTable
  }

  #extract the required columns from the UniProtTable (UPT) and ensure the content is considered as character by R
  UPT<-UPT[,colnames(UPT) %in% c("Uniprot_Accession","Uniprot_ID","Reactome_List")]

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
  UniProtTable_REACTOME<-data.frame(Uniprot_Accession=Uniprot_Accession, Uniprot_ID=Uniprot_ID, REACTOME=unlist(uniList))

  #download the REACTOME pathway names
  #to be added

  #Place the reordered table in the Global Environment
  UniProtTable_REACTOME <<- UniProtTable_REACTOME

  #export file if export==TRUE
  if(export==TRUE){
    save(UniProtTable_GO,file = file)
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
  #general checking
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

  #check universe
  if(missing(universe)){
    print("Your <universe> was missing the global list of ID present in UniProtTable was used as background.")
    universe<-as.character(t(UniProtTable$Uniprot_ID))
  }

  #if query contained multiple IDs separated by semicolons only the first one was conserved
  if(sum(grepl(";",query))>0){
    print(paste0(sum(grepl(";",query))," proteins of your <query> list had more than one identifier, only the first one was conserved"))
    query<-gsub("\\;.*","",query)
  }

  #if universe contained multiple IDs separated by semicolons only the first one was conserved
  if(sum(grepl(";",universe))>0){
    print(paste0(sum(grepl(";",universe))," proteins of your <universe> list had more than one identifier, only the first one was conserved"))
    universe<-gsub("\\;.*","",universe)
  }

  #check how many of the query are in the lists UniProtTable
  query_ID <- query[query %in% UniProtTable$Uniprot_ID]
  query_Accession <- query[query %in% UniProtTable$Uniprot_Accession]
  print(paste0("Your <query> contained ", length(query_ID), " UniProt_IDs and ", length(query_Accession), " UniProt_Accession (some might be redundant)." ))

  #if no values in the UniProtTable stop
  if (length(query_ID)==0&&length(query_Accession)==0){stop("Your <query> should be Uniprot_IDs or a Uniprot_Accessions, ensure that the organism of your query is present in the UniProtTable.")}

  #check how many of the query are in the lists UniProtTable
  universe_ID <- universe[universe %in% UniProtTable$Uniprot_ID]
  universe_Accession <- universe[universe %in% UniProtTable$Uniprot_Accession]
  print(paste0("Your <universe> contained ", length(universe_ID), " UniProt_IDs and ", length(universe_Accession), " UniProt_Accession (some might be redundant)." ))

  #if no values in the UniProtTable stop
  if (length(query_ID)==0&&length(query_Accession)==0){stop("Your <query> should be Uniprot_IDs or a Uniprot_Accessions, ensure that the organism of your query is present in the UniProtTable.")}
  if (length(universe_ID)==0&&length(universe_Accession)==0){stop("Your <universe> should be Uniprot_IDs or a Uniprot_Accessions, ensure that the organism of your universe is present in the UniProtTable.")}

  #convert uniprot_Accession in Uniprot
  if(length(query_Accession)>0){
    print("The uniprot_Accession of your query were converted in Uniprot_IDs.")
    query_Accession<- as.character(t(UniProtTable[match(query_Accession,UniProtTable$Uniprot_Accession),2]))
    query_ID<-c(query_ID,query_Accession)
  }

  if(length(universe_Accession)>0){
    print("The uniprot_Accession of your universe were converted in Uniprot_IDs.")
    universe_Accession<- as.character(t(UniProtTable[match(universe_Accession,UniProtTable$Uniprot_Accession),2]))
    universe_ID<-c(universe_ID,universe_Accession)
  }

  #remove duplicates
  duplicated_ID<-duplicated(query_ID)
  if(sum(duplicated_ID)>0){
    print(paste0("Your list of Uniprot_IDs was containing ", sum(duplicated_ID)," duplicated value(s), duplicates were removed, the query now comprises ", length(query_ID)-sum(duplicated_ID), " unique IDs."))
    query_ID<-query_ID[!duplicated_ID]
    }

  duplicated_ID<-duplicated(universe_ID)
  if(sum(duplicated_ID)>0){
    print(paste0("Your list of Uniprot_IDs was containing ", sum(duplicated_ID)," duplicated value(s), duplicates were removed, the universe now comprises ", length(universe_ID)-sum(duplicated_ID), " unique IDs."))
    universe_ID<-universe_ID[!duplicated_ID]
  }

  #check if the universe contains all the proteins from the query
  if(sum(is.na(match(query_ID,universe_ID)))>0){warning("Your <universe> does not comprises all the IDs present in your <query>")}

  #generate the UniProtTable_GO for the query and the universe
  UPT_GO_query <- UPT_GO[UPT_GO$Uniprot_ID %in% query_ID,]
  UPT_GO_universe <-UPT_GO[UPT_GO$Uniprot_ID %in% universe_ID,]

  #create the unique_GO list based on the UPT_GO_query
  unique_GO<- unique(UPT_GO_query[,3:5])

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
  fisher_results<-data.frame(matrix(NA, nrow=nrow(unique_GO), ncol=13))
  colnames(fisher_results)<- c("Category","Term_accession","Term_description","Count_query","Pop_query","Count_universe","Pop_universe","%_query","%_universe","pval","adjpval","fold_change","Proteins_in_query")
  fisher_results$Category<-unique_GO$GO_type
  fisher_results$Term_accession<-unique_GO$GO_accession
  fisher_results$Term_description<-unique_GO$GO_description
  fisher_results$Pop_query<-rep(length(query_ID),nrow(fisher_results))
  fisher_results$Pop_universe<-rep(length(universe_ID),nrow(fisher_results))

  #run the test and populate the table
  for (i in 1:nrow(unique_GO)){
    test<-fisher.test(contingency_list[[i]])
    fisher_results$pval[i]<-test$p.value
    fisher_results$Count_query[i]<-contingency_list[[i]][1,1]
    fisher_results$Count_universe[i]<-contingency_list[[i]][1,2]
    fisher_results$`%_query` [i]<- contingency_list[[i]][1,1]/length(query_ID)*100
    fisher_results$`%_universe`[i]<- contingency_list[[i]][1,2]/length(universe_ID)*100
    fisher_results$fold_change[i]<- fisher_results$`%_query`[i]/fisher_results$`%_universe`[i]
    fisher_results$Proteins_in_query[i]<-paste(UPT_GO_query[UPT_GO_query$GO_accession==fisher_results$Term_accession[i],]$Uniprot_ID,collapse =";")
  }

  fisher_results_BP<-fisher_results[fisher_results$Category=="GO_BP",]
  fisher_results_BP$adjpval<-p.adjust(p=fisher_results_BP$pval,method = "BH",length(fisher_results_BP$pval))

  fisher_results_CC<-fisher_results[fisher_results$Category=="GO_CC",]
  fisher_results_CC$adjpval<-p.adjust(p=fisher_results_CC$pval,method = "BH",length(fisher_results_CC$pval))

  fisher_results_MF<-fisher_results[fisher_results$Category=="GO_MF",]
  fisher_results_MF$adjpval<-p.adjust(p=fisher_results_MF$pval,method = "BH",length(fisher_results_MF$pval))

  fisher_results<-rbind(fisher_results_BP,fisher_results_CC,fisher_results_MF)
  fisher_results<-fisher_results[order(fisher_results$pval),]
  return(fisher_results)
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
  #general checking
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

  #check universe
  if(missing(universe)){
    print("Your <universe> was missing the global list of ID present in UniProtTable was used as background.")
    universe<-as.character(t(UniProtTable$Uniprot_ID))
  }

  #if query contained multiple IDs separated by semicolons only the first one was conserved
  if(sum(grepl(";",query))>0){
    print(paste0(sum(grepl(";",query))," proteins of your <query> list had more than one identifier, only the first one was conserved"))
    query<-gsub("\\;.*","",query)
  }

  #if universe contained multiple IDs separated by semicolons only the first one was conserved
  if(sum(grepl(";",universe))>0){
    print(paste0(sum(grepl(";",universe))," proteins of your <universe> list had more than one identifier, only the first one was conserved"))
    universe<-gsub("\\;.*","",universe)
  }

  #check how many of the query are in the lists UniProtTable
  query_ID <- query[query %in% UniProtTable$Uniprot_ID]
  query_Accession <- query[query %in% UniProtTable$Uniprot_Accession]
  print(paste0("Your <query> contained ", length(query_ID), " UniProt_IDs and ", length(query_Accession), " UniProt_Accession (some might be redundant)." ))

  #if no values in the UniProtTable stop
  if (length(query_ID)==0&&length(query_Accession)==0){stop("Your <query> should be Uniprot_IDs or a Uniprot_Accessions, ensure that the organism of your query is present in the UniProtTable.")}

  #check how many of the query are in the lists UniProtTable
  universe_ID <- universe[universe %in% UniProtTable$Uniprot_ID]
  universe_Accession <- universe[universe %in% UniProtTable$Uniprot_Accession]
  print(paste0("Your <universe> contained ", length(universe_ID), " UniProt_IDs and ", length(universe_Accession), " UniProt_Accession (some might be redundant)." ))

  #if no values in the UniProtTable stop
  if (length(query_ID)==0&&length(query_Accession)==0){stop("Your <query> should be Uniprot_IDs or a Uniprot_Accessions, ensure that the organism of your query is present in the UniProtTable.")}
  if (length(universe_ID)==0&&length(universe_Accession)==0){stop("Your <universe> should be Uniprot_IDs or a Uniprot_Accessions, ensure that the organism of your universe is present in the UniProtTable.")}

  #convert uniprot_Accession in Uniprot
  if(length(query_Accession)>0){
    print("The uniprot_Accession of your query were converted in Uniprot_IDs.")
    query_Accession<- as.character(t(UniProtTable[match(query_Accession,UniProtTable$Uniprot_Accession),2]))
    query_ID<-c(query_ID,query_Accession)
  }

  if(length(universe_Accession)>0){
    print("The uniprot_Accession of your universe were converted in Uniprot_IDs.")
    universe_Accession<- as.character(t(UniProtTable[match(universe_Accession,UniProtTable$Uniprot_Accession),2]))
    universe_ID<-c(universe_ID,universe_Accession)
  }

  #remove duplicates
  duplicated_ID<-duplicated(query_ID)
  if(sum(duplicated_ID)>0){
    print(paste0("Your list of Uniprot_IDs was containing ", sum(duplicated_ID)," duplicated value(s), duplicates were removed, the query now comprises ", length(query_ID)-sum(duplicated_ID), " unique IDs."))
    query_ID<-query_ID[!duplicated_ID]
  }

  duplicated_ID<-duplicated(universe_ID)
  if(sum(duplicated_ID)>0){
    print(paste0("Your list of Uniprot_IDs was containing ", sum(duplicated_ID)," duplicated value(s), duplicates were removed, the universe now comprises ", length(universe_ID)-sum(duplicated_ID), " unique IDs."))
    universe_ID<-universe_ID[!duplicated_ID]
  }

  #check if the universe contains all the proteins from the query
  if(sum(is.na(match(query_ID,universe_ID)))>0){warning("Your <universe> does not comprises all the IDs present in your <query>")}

  #generate the UniProtTable_GO for the query and the universe
  UPT_GO_query <- UPT_GO[UPT_GO$Uniprot_ID %in% query_ID,]
  UPT_GO_universe <-UPT_GO[UPT_GO$Uniprot_ID %in% universe_ID,]

  #create the unique_GO list based on the UPT_GO_query
  unique_GO<- unique(UPT_GO_query[,3:5])

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
  EASE_results<-data.frame(matrix(NA, nrow=nrow(unique_GO), ncol=13))
  colnames(EASE_results)<- c("Category","Term_accession","Term_description","Count_query","Pop_query","Count_universe","Pop_universe","%_query","%_universe","pval","adjpval","fold_change","Proteins_in_query")
  EASE_results$Category<-unique_GO$GO_type
  EASE_results$Term_accession<-unique_GO$GO_accession
  EASE_results$Term_description<-unique_GO$GO_description
  EASE_results$Pop_query<-rep(length(query_ID),nrow(EASE_results))
  EASE_results$Pop_universe<-rep(length(universe_ID),nrow(EASE_results))

  #run the test and populate the table
  for (i in 1:nrow(unique_GO)){
    test<-fisher.test(contingency_list[[i]])
    EASE_results$pval[i]<-test$p.value
    EASE_results$Count_query[i]<-contingency_list[[i]][1,1]+1
    EASE_results$Count_universe[i]<-contingency_list[[i]][1,2]
    EASE_results$`%_query` [i]<- contingency_list[[i]][1,1]/length(query_ID)*100
    EASE_results$`%_universe`[i]<- contingency_list[[i]][1,2]/length(universe_ID)*100
    EASE_results$fold_change[i]<- EASE_results$`%_query`[i]/EASE_results$`%_universe`[i]
    EASE_results$Proteins_in_query[i]<-paste(UPT_GO_query[UPT_GO_query$GO_accession==EASE_results$Term_accession[i],]$Uniprot_ID,collapse =";")
  }

  EASE_results_BP<-EASE_results[EASE_results$Category=="GO_BP",]
  EASE_results_BP$adjpval<-p.adjust(p=EASE_results_BP$pval,method = "BH",length(EASE_results_BP$pval))

  EASE_results_CC<-EASE_results[EASE_results$Category=="GO_CC",]
  EASE_results_CC$adjpval<-p.adjust(p=EASE_results_CC$pval,method = "BH",length(EASE_results_CC$pval))

  EASE_results_MF<-EASE_results[EASE_results$Category=="GO_MF",]
  EASE_results_MF$adjpval<-p.adjust(p=EASE_results_MF$pval,method = "BH",length(EASE_results_MF$pval))

  EASE_results<-rbind(EASE_results_BP,EASE_results_CC,EASE_results_MF)
  EASE_results<-EASE_results[order(EASE_results$pval),]
  return(EASE_results)
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
  #general checking
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
  }else {
    generate_UniProtTable_GO(organismID=organismID,reviewed=reviewed,preloaded_UniProtTable==TRUE)
    UPT_GO<- UniProtTable_GO
  }

  #if rankingTable_IDs contained multiple IDs separated by semicolons only the first one was conserved
  if(sum(grepl(";",rankingTable[,1]))>0){
    print(paste0(sum(grepl(";",rankingTable[,1]))," proteins of your <rankingTable_IDs> list had more than one identifier, only the first one was conserved"))
    rankingTable[,1]<-gsub("\\;.*","",rankingTable[,1])
  }

  #check how many of the rankingTable_IDs are in the lists UniProtTable
  rankingTable_ID <- rankingTable[rankingTable[,1] %in% UniProtTable$Uniprot_ID,1]
  rankingTable_Accession <- rankingTable[rankingTable[,1] %in% UniProtTable$Uniprot_Accession,1]

   #if no values in the UniProtTable stop
  if (length(rankingTable_ID)==0&&length(rankingTable_Accession)==0){stop("Your <rankingTable> should be Uniprot_IDs or a Uniprot_Accessions, ensure that the organism of your rankingTable is present in the UniProtTable.")}

  #convert uniprot_Accession in Uniprot
  if(length(rankingTable_Accession)>0){
    print("The uniprot_Accession of your rankingTable were converted in Uniprot_IDs.")
    rankingTable[rankingTable[,1] %in% UniProtTable$Uniprot_Accession,1]<-as.character(t(UniProtTable[match(rankingTable[rankingTable[,1] %in% UniProtTable$Uniprot_Accession,1],UniProtTable$Uniprot_Accession),2]))
  }

  #remove duplicates
  duplicated_ID<-duplicated(rankingTable[,1])
  if(sum(duplicated_ID)>0){
    print(paste0("Your list of Uniprot_IDs was containing ", sum(duplicated_ID)," duplicated value(s), duplicates were removed, the query now comprises ", length(query_ID)-sum(duplicated_ID), " unique IDs."))
    rankingTable<-rankingTable[!duplicated_ID,]
  }

    #order based on the ascending/descending order parameter
  if(order=="ascending"){
    rankingTable<-rankingTable[order(as.numeric(rankingTable[,2])),]
  }else{
      rankingTable<-rankingTable[order(as.numeric(rankingTable[,2]),decreasing=TRUE),]
      }

  #generate the UniProtTable_GO for the rankingTable
  UPT_GO_rankingTable <- UPT_GO[UPT_GO$Uniprot_ID %in% rankingTable[,1],]

  #create the unique_GO list based on the UPT_GO_query
  unique_GO<- unique(UPT_GO_rankingTable[,3:5])

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

  names(RT_by_GO)<-paste0(unique_GO$GO_accession,"@",unique_GO$GO_description,"@",unique_GO$GO_type)

  #realise the KS.test
  p<-numeric()
  for (i in 1:nrow(unique_GO)){
    if(nrow(rankingTable)>length(RT_by_GO[[i]])){
      p[i]<-ks.test(RT_by_GO[[i]]$KSRank,seq_len(nrow(rankingTable))[-RT_by_GO[[i]]$KSRank],alternative="greater")$p.value}else{p[i]<-1}
  }

  #create the final table
  final_table<-data.frame(matrix(ncol=4,nrow= nrow(unique_GO)))
  colnames(final_table) <- c("Count.in.list", "pval", "adjpval","IDs.in.list")
  final_table<-cbind(unique_GO,final_table)
  for(i in 1:nrow(unique_GO)){
    final_table$Count.in.list[i]<-nrow(RT_by_GO[[i]])
    final_table$IDs.in.list[i]<-IDs_by_GO[[i]]
    }
  final_table$pval<-p

  final_table_BP<-final_table[final_table$GO_type=="GO_BP",]
  final_table_BP$adjpval<-p.adjust(p=final_table_BP$pval,method = "BH",length(final_table_BP$pval))

  final_table_CC<-final_table[final_table$GO_type=="GO_CC",]
  final_table_CC$adjpval<-p.adjust(p=final_table_CC$pval,method = "BH",length(final_table_CC$pval))

  final_table_MF<-final_table[final_table$GO_type=="GO_MF",]
  final_table_MF$adjpval<-p.adjust(p=final_table_MF$pval,method = "BH",length(final_table_MF$pval))

  final_table<-rbind(final_table_BP,final_table_CC,final_table_MF)

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
    #general checking
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
      UPT_KEGG<- UniProtTable_KEGG
    }else {
      generate_UniProtTable_GO(organismID=organismID,reviewed=reviewed,preloaded_UniProtTable==TRUE)
      UPT_KEGG<- UniProtTable_KEGG
    }

    #check query
    if(missing(query)){stop("Your <query> is missing")}

    #check universe
    if(missing(universe)){
      print("Your <universe> was missing the global list of ID present in UniProtTable was used as background.")
      universe<-as.character(t(UniProtTable$Uniprot_ID))
    }

    #if query contained multiple IDs separated by semicolons only the first one was conserved
    if(sum(grepl(";",query))>0){
      print(paste0(sum(grepl(";",query))," proteins of your <query> list had more than one identifier, only the first one was conserved"))
      query<-gsub("\\;.*","",query)
    }

    #if universe contained multiple IDs separated by semicolons only the first one was conserved
    if(sum(grepl(";",universe))>0){
      print(paste0(sum(grepl(";",universe))," proteins of your <universe> list had more than one identifier, only the first one was conserved"))
      universe<-gsub("\\;.*","",universe)
    }

    #check how many of the query are in the lists UniProtTable
    query_ID <- query[query %in% UniProtTable$Uniprot_ID]
    query_Accession <- query[query %in% UniProtTable$Uniprot_Accession]
    print(paste0("Your <query> contained ", length(query_ID), " UniProt_IDs and ", length(query_Accession), " UniProt_Accession (some might be redundant)." ))

    #if no values in the UniProtTable stop
    if (length(query_ID)==0&&length(query_Accession)==0){stop("Your <query> should be Uniprot_IDs or a Uniprot_Accessions, ensure that the organism of your query is present in the UniProtTable.")}

    #check how many of the query are in the lists UniProtTable
    universe_ID <- universe[universe %in% UniProtTable$Uniprot_ID]
    universe_Accession <- universe[universe %in% UniProtTable$Uniprot_Accession]
    print(paste0("Your <universe> contained ", length(universe_ID), " UniProt_IDs and ", length(universe_Accession), " UniProt_Accession (some might be redundant)." ))

    #if no values in the UniProtTable stop
    if (length(query_ID)==0&&length(query_Accession)==0){stop("Your <query> should be Uniprot_IDs or a Uniprot_Accessions, ensure that the organism of your query is present in the UniProtTable.")}
    if (length(universe_ID)==0&&length(universe_Accession)==0){stop("Your <universe> should be Uniprot_IDs or a Uniprot_Accessions, ensure that the organism of your universe is present in the UniProtTable.")}

    #convert uniprot_Accession in Uniprot
    if(length(query_Accession)>0){
      print("The uniprot_Accession of your query were converted in Uniprot_IDs.")
      query_Accession<- as.character(t(UniProtTable[match(query_Accession,UniProtTable$Uniprot_Accession),2]))
      query_ID<-c(query_ID,query_Accession)
    }

    if(length(universe_Accession)>0){
      print("The uniprot_Accession of your universe were converted in Uniprot_IDs.")
      universe_Accession<- as.character(t(UniProtTable[match(universe_Accession,UniProtTable$Uniprot_Accession),2]))
      universe_ID<-c(universe_ID,universe_Accession)
    }

    #remove duplicates
    duplicated_ID<-duplicated(query_ID)
    if(sum(duplicated_ID)>0){
      print(paste0("Your list of Uniprot_IDs was containing ", sum(duplicated_ID)," duplicated value(s), duplicates were removed, the query now comprises ", length(query_ID)-sum(duplicated_ID), " unique IDs."))
      query_ID<-query_ID[!duplicated_ID]
    }

    duplicated_ID<-duplicated(universe_ID)
    if(sum(duplicated_ID)>0){
      print(paste0("Your list of Uniprot_IDs was containing ", sum(duplicated_ID)," duplicated value(s), duplicates were removed, the universe now comprises ", length(universe_ID)-sum(duplicated_ID), " unique IDs."))
      universe_ID<-universe_ID[!duplicated_ID]
    }

    #check if the universe contains all the proteins from the query
    if(sum(is.na(match(query_ID,universe_ID)))>0){warning("Your <universe> does not comprises all the IDs present in your <query>")}

    #generate the UniProtTable_GO for the query and the universe
    UPT_KEGG_query <- UPT_KEGG[UPT_KEGG$Uniprot_ID %in% query_ID,]
    UPT_KEGG_universe <-UPT_KEGG[UPT_KEGG$Uniprot_ID %in% universe_ID,]

    #create the unique_KEGGO list based on the UPT_KEGG_query
    unique_KEGG<- unique(UPT_KEGG_query[,4:5])
    unique_KEGG<- unique_KEGG[!is.na(unique_KEGG$Pathway_ID),]

    #create the Contingencies matrix
    contingency_list<-list()
    for (i in 1:nrow(unique_KEGG)){
      contingency_list[[i]]<-matrix(ncol=2,nrow=2)
      contingency_list[[i]][1,1]<- sum(UPT_KEGG_query[,4]==unique_KEGG[i,1])
      contingency_list[[i]][1,2]<- sum(UPT_KEGG_universe[,4]==unique_KEGG[i,1])
      contingency_list[[i]][2,1]<- length(query_ID)-contingency_list[[i]][1,1]
      contingency_list[[i]][2,2]<- length(universe)-contingency_list[[i]][1,2]
    }
    names(contingency_list)<-unique_KEGG[,1]

    #create Final table
    fisher_results<-data.frame(matrix(NA, nrow=nrow(unique_KEGG), ncol=13))
    colnames(fisher_results)<- c("Category","Term_accession","Term_description","Count_query","Pop_query","Count_universe","Pop_universe","%_query","%_universe","pval","adjpval","fold_change","Proteins_in_query")
    fisher_results$Category<-"KEGG_pathway"
    fisher_results$Term_accession<-unique_KEGG$Pathway_ID
    fisher_results$Term_description<-unique_KEGG$Pathway_name
    fisher_results$Pop_query<-rep(length(query_ID),nrow(fisher_results))
    fisher_results$Pop_universe<-rep(length(universe_ID),nrow(fisher_results))

    #run the test and populate the table
    for (i in 1:nrow(unique_KEGG)){
      test<-fisher.test(contingency_list[[i]])
      fisher_results$pval[i]<-test$p.value
      fisher_results$Count_query[i]<-contingency_list[[i]][1,1]
      fisher_results$Count_universe[i]<-contingency_list[[i]][1,2]
      fisher_results$`%_query` [i]<- contingency_list[[i]][1,1]/length(query_ID)*100
      fisher_results$`%_universe`[i]<- contingency_list[[i]][1,2]/length(universe_ID)*100
      fisher_results$fold_change[i]<- fisher_results$`%_query`[i]/fisher_results$`%_universe`[i]
      fisher_results$Proteins_in_query[i]<-paste(UPT_KEGG_query[UPT_KEGG_query$KEGG_ID==fisher_results$Term_accession[i],]$Uniprot_ID,collapse =";")
    }

    fisher_results<-fisher_results[order(fisher_results$pval),]
    fisher_results$adjpval<-p.adjust(fisher_results$pval)
    return(fisher_results)
  }

#' UniProt_KEGG_EASE()
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
  #general checking
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
    UPT_KEGG<- UniProtTable_KEGG
  }else {
    generate_UniProtTable_GO(organismID=organismID,reviewed=reviewed,preloaded_UniProtTable==TRUE)
    UPT_KEGG<- UniProtTable_KEGG
  }

  #check query
  if(missing(query)){stop("Your <query> is missing")}

  #check universe
  if(missing(universe)){
    print("Your <universe> was missing the global list of ID present in UniProtTable was used as background.")
    universe<-as.character(t(UniProtTable$Uniprot_ID))
  }

  #if query contained multiple IDs separated by semicolons only the first one was conserved
  if(sum(grepl(";",query))>0){
    print(paste0(sum(grepl(";",query))," proteins of your <query> list had more than one identifier, only the first one was conserved"))
    query<-gsub("\\;.*","",query)
  }

  #if universe contained multiple IDs separated by semicolons only the first one was conserved
  if(sum(grepl(";",universe))>0){
    print(paste0(sum(grepl(";",universe))," proteins of your <universe> list had more than one identifier, only the first one was conserved"))
    universe<-gsub("\\;.*","",universe)
  }

  #check how many of the query are in the lists UniProtTable
  query_ID <- query[query %in% UniProtTable$Uniprot_ID]
  query_Accession <- query[query %in% UniProtTable$Uniprot_Accession]
  print(paste0("Your <query> contained ", length(query_ID), " UniProt_IDs and ", length(query_Accession), " UniProt_Accession (some might be redundant)." ))

  #if no values in the UniProtTable stop
  if (length(query_ID)==0&&length(query_Accession)==0){stop("Your <query> should be Uniprot_IDs or a Uniprot_Accessions, ensure that the organism of your query is present in the UniProtTable.")}

  #check how many of the query are in the lists UniProtTable
  universe_ID <- universe[universe %in% UniProtTable$Uniprot_ID]
  universe_Accession <- universe[universe %in% UniProtTable$Uniprot_Accession]
  print(paste0("Your <universe> contained ", length(universe_ID), " UniProt_IDs and ", length(universe_Accession), " UniProt_Accession (some might be redundant)." ))

  #if no values in the UniProtTable stop
  if (length(query_ID)==0&&length(query_Accession)==0){stop("Your <query> should be Uniprot_IDs or a Uniprot_Accessions, ensure that the organism of your query is present in the UniProtTable.")}
  if (length(universe_ID)==0&&length(universe_Accession)==0){stop("Your <universe> should be Uniprot_IDs or a Uniprot_Accessions, ensure that the organism of your universe is present in the UniProtTable.")}

  #convert uniprot_Accession in Uniprot
  if(length(query_Accession)>0){
    print("The uniprot_Accession of your query were converted in Uniprot_IDs.")
    query_Accession<- as.character(t(UniProtTable[match(query_Accession,UniProtTable$Uniprot_Accession),2]))
    query_ID<-c(query_ID,query_Accession)
  }

  if(length(universe_Accession)>0){
    print("The uniprot_Accession of your universe were converted in Uniprot_IDs.")
    universe_Accession<- as.character(t(UniProtTable[match(universe_Accession,UniProtTable$Uniprot_Accession),2]))
    universe_ID<-c(universe_ID,universe_Accession)
  }

  #remove duplicates
  duplicated_ID<-duplicated(query_ID)
  if(sum(duplicated_ID)>0){
    print(paste0("Your list of Uniprot_IDs was containing ", sum(duplicated_ID)," duplicated value(s), duplicates were removed, the query now comprises ", length(query_ID)-sum(duplicated_ID), " unique IDs."))
    query_ID<-query_ID[!duplicated_ID]
  }

  duplicated_ID<-duplicated(universe_ID)
  if(sum(duplicated_ID)>0){
    print(paste0("Your list of Uniprot_IDs was containing ", sum(duplicated_ID)," duplicated value(s), duplicates were removed, the universe now comprises ", length(universe_ID)-sum(duplicated_ID), " unique IDs."))
    universe_ID<-universe_ID[!duplicated_ID]
  }

  #check if the universe contains all the proteins from the query
  if(sum(is.na(match(query_ID,universe_ID)))>0){warning("Your <universe> does not comprises all the IDs present in your <query>")}

  #generate the UniProtTable_GO for the query and the universe
  UPT_KEGG_query <- UPT_KEGG[UPT_KEGG$Uniprot_ID %in% query_ID,]
  UPT_KEGG_universe <-UPT_KEGG[UPT_KEGG$Uniprot_ID %in% universe_ID,]

  #create the unique_KEGG list based on the UPT_KEGG_query
  unique_KEGG<- unique(UPT_KEGG_query[,4:5])
  unique_KEGG<- unique_KEGG[!is.na(unique_KEGG$Pathway_ID),]

  #create the Contingencies matrix
  contingency_list<-list()
  for (i in 1:nrow(unique_KEGG)){
    contingency_list[[i]]<-matrix(ncol=2,nrow=2)
    contingency_list[[i]][1,1]<- sum(UPT_KEGG_query[,4]==unique_KEGG[i,1])-1
    contingency_list[[i]][1,2]<- sum(UPT_KEGG_universe[,4]==unique_KEGG[i,1])
    contingency_list[[i]][2,1]<- length(query_ID)-contingency_list[[i]][1,1]
    contingency_list[[i]][2,2]<- length(universe)-contingency_list[[i]][1,2]
  }
  names(contingency_list)<-unique_KEGG[,1]

  #create Final table
  fisher_results<-data.frame(matrix(NA, nrow=nrow(unique_KEGG), ncol=13))
  colnames(fisher_results)<- c("Category","Term_accession","Term_description","Count_query","Pop_query","Count_universe","Pop_universe","%_query","%_universe","pval","adjpval","fold_change","Proteins_in_query")
  fisher_results$Category<-"KEGG_pathway"
  fisher_results$Term_accession<-unique_KEGG$Pathway_ID
  fisher_results$Term_description<-unique_KEGG$Pathway_name
  fisher_results$Pop_query<-rep(length(query_ID),nrow(fisher_results))
  fisher_results$Pop_universe<-rep(length(universe_ID),nrow(fisher_results))

  #run the test and populate the table
  for (i in 1:nrow(unique_KEGG)){
    test<-fisher.test(contingency_list[[i]])
    fisher_results$pval[i]<-test$p.value
    fisher_results$Count_query[i]<-contingency_list[[i]][1,1]+1
    fisher_results$Count_universe[i]<-contingency_list[[i]][1,2]
    fisher_results$`%_query` [i]<- contingency_list[[i]][1,1]/length(query_ID)*100
    fisher_results$`%_universe`[i]<- contingency_list[[i]][1,2]/length(universe_ID)*100
    fisher_results$fold_change[i]<- fisher_results$`%_query`[i]/fisher_results$`%_universe`[i]
    fisher_results$Proteins_in_query[i]<-paste(UPT_KEGG_query[UPT_KEGG_query$KEGG_ID==fisher_results$Term_accession[i],]$Uniprot_ID,collapse =";")
  }

  fisher_results<-fisher_results[order(fisher_results$pval),]
  fisher_results$adjpval<-p.adjust(fisher_results$pval)
  return(fisher_results)
}

