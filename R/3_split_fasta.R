#' splitFasta()
#' @description split a fasta in a sequence of fasta files containing a set number of
#' @param organismID Should be a UniProt organism ID, by default Homo sapiens will be used ("9606").
#' @param reviewed indicate if the proteins extracted have to be reviewed or not (typically the ones present on swissProt), by default TRUE (the non reviewed proteins are discarded), if FALSE both reviewed and not reviewed will be conserved.
#' @param export indicates if the table should be exported for future use (using a previously downloaded table dramatically reduces the analysis time).
#' @param file should be a string indicating the location and name of the destination file, by default the file will be named with the format "UniProtTable_organismID_date.txt".
#'
#' @return This function returns a data.frame named UniProtTable containing the detail on the proteins for the queried organism, if the export is set on TRUE, the function will also save the table in the .rda format at the specified location
#'
#' @examples download_UniProtTable(organismID="9606", reviewed=TRUE, export=FALSE, file="UniProtTable_organismID_date.txt")
#'
splitFasta<- function(file="file.fasta",seq_number=10000){
  #
  if(missing(file)){stop("The location of the fasta file is missing. Please indicate the fasta file to split")}
  if(missing(seq_number)){
    seq_number=10000
    warning("The number of sequences to use for the split fasta file was not set, your file will be splitted in files smaller than 10000 sequences.")
  }
  if(!is.numeric(seq_number)){stop("seq_number has to be numerical")}

  #attempt to read the file
  fasta<-read.delim(file=file,header= F)
  header <- grepl(">",fasta[,1])

  #create a vector containing the sequence number position
  sequence_number<-numeric()
  if(header[1]!=TRUE){stop("The first sequence of the fasta files does not start with the fasta header delimiter '>'")}else{sequence_number[1]<-1}
  for(i in 2:length(header)){
  if(header[i]==T){sequence_number[i] <- sequence_number[i-1]+1}else{sequence_number[i]<-sequence_number[i-1]}
  }

  #create a list of cutting points
  cutpoint<-seq(0,max(sequence_number),by=seq_number)
  cutpoint<-c(cutpoint, max(sequence_number))

  #create the sub fasta files and place them in the same folder as the original fasta
    #find the folder location
    folder<- sub("\\/[^\\/]+$","/",file)
    filename<- sub(".*\\/", "", file)
    filename<- sub("\\..*","",filename)
    extension<- sub(".*\\.",".",file)

  for(i in 2:length(cutpoint)){
    seq_start <- cutpoint[i-1]+1
    seq_end <- cutpoint[i]
    selection_fasta <- fasta[sequence_number>=seq_start&sequence_number<=seq_end,]
    write.table(selection_fasta,file = paste0(folder,filename,"_seq_",seq_start,"_to_",seq_end,extension),row.names = F,col.names = F,quote=F)
  }
  }

