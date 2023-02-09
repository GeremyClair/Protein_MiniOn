#' romicsEnrichment()
#' allows to perform proteinMinion enrichments directly on a romics_object. The Query list can be generated using the filters.
#'
#' @param romics_object A romics_object created using romicsCreateObject() and transformed using the function log2transform() or log10transform()
#'
#' @details This function requires the package 'pmartR' to be installed and loaded to be excecuted. It will calculate and plot the samples to be filtered out using the function Romics_outlier_eval().
#' @return This function will print the pmartR filtering details and will return 2 plots the first one is a scatter plot of the pvalue by log2(Robust Mahalanobis Distance) the second one is a scatter plot of the log2(Robust Mahalanobis Distance) per sample.
#'
#' @references Matzke, M., Waters, K., Metz, T., Jacobs, J., Sims, A., Baric, R., Pounds, J., and Webb-Robertson, B.J. (2011), Improved quality control processing of peptide-centric LC-MS proteomics data. Bioinformatics. 27(20): 2866-2872.
#'
#' @author Geremy Clair
#' @export
#'
romicsEnrichement<-function(romics_object, organismID="9606", type= c("GO", "KEGG","REACTOME"), ANOVA_filter= c("none", "p", "padj"), p = 0.05,  cluster_filter="none",statCol_filter="none", statCol_text="<=0.05", statCol_filter2="none", statCol_text2="<=0.05", statCol_filter3="none", statCol_text3="<=0.05", enrichment_function="EASE", ...){
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(missing(ANOVA_filter)){ANOVA_filter="none"}
  if(!ANOVA_filter %in% c("none","p","padj")){stop("Your ANOVA filter should be either none, p or padj")}
  if(missing(p)){p<-0.05}
  if(missing(cluster_filter)){cluster_filter="none"}
  if(!cluster_filter %in% c("none",romicsCalculatedStats(romics_object))){stop("<cluster_filter> is not a calculated stat column of the romics_object. Use the function romicsCalculatedStats() to identify the usable columns.")}
  if(cluster_filter!="none" & !grepl("_cluster",cluster_filter)){stop("<cluster_filter> has to be a '_cluster' statistic")}
  if(missing(statCol_filter)){statCol_filter="none"}
  if(!statCol_filter %in% c("none",romicsCalculatedStats(romics_object))){stop("<statCol_filter> is not a calculated stat column of the romics_object. Use the function romicsCalculatedStats() to identify the usable columns.")}
  if(missing(statCol_filter2)){statCol_filter2="none"}
  if(!statCol_filter2 %in% c("none",romicsCalculatedStats(romics_object))){stop("<statCol_filter2> is not a calculated stat column of the romics_object. Use the function romicsCalculatedStats() to identify the usable columns.")}
  if(missing(statCol_filter3)){statCol_filter3="none"}
  if(!statCol_filter3 %in% c("none",romicsCalculatedStats(romics_object))){stop("<statCol_filter3> is not a calculated stat column of the romics_object. Use the function romicsCalculatedStats() to identify the usable columns.")}
  if(is.null(romics_object$statistics)){stop("your romics_object does not contain a statistical layer")}
  if(missing(type)){type= c("KEGG", "GO","REACTOME")}
  if(sum(type %in% c("KEGG", "GO","REACTOME"))!=length(type)){stop("<type> has to be c('GO','KEGG','REACTOME'), or to contain 'KEGG', 'REACTOME', or 'GO'")}
  if(sum(c(ANOVA_filter, statCol_filter, statCol_filter, statCol_filter)!="none")>1 & cluster_filter!="none"){stop("this function cannot do simultaneously a cluster filter AND an other type of statistic filter please either set cluster_filter to 'none' or all the other filters to 'none'.")}
  if(missing(enrichment_function)){enrichment_function<-"EASE"}
  if(!enrichment_function %in% c("EASE","Fisher")){stop("<enrichment_function> has to be either 'EASE' or 'Fisher'")}

  data<-romics_object$data
  statistics<-romics_object$statistics

  names<- rownames(romics_object$statistics)
  if(sum(grepl(";",names))>0){
    print("Your ID list was containing more than one identifier separated with semicolons, only the first IDs were conserved")
    names<-gsub(";.*", "", names)
  }
  rownames(data)<-names

  #define the universe
  universe <- names

  #create a groups list.
  groups<-list()

  #Filter based on ANOVA
  if(ANOVA_filter=="p"){
    if(is.null(romics_object$statistics$ANOVA_p)){
      warning("The ANOVA_p has not been calculated, no filtering was applied")}else{
        data <- data[statistics$ANOVA_p<p,]
        statistics<-statistics[statistics$ANOVA_p<p,]
        groups[[1]]<-rownames(data)
        names(groups)[1]<-paste0("ANOVA_p<",p)
    }}

  if(ANOVA_filter=="padj"){
    if(is.null(romics_object$statistics$ANOVA_p)){
      warning("The ANOVA_padj has not been calculated, no filtering was applied")}else{
        data <- data[statistics$ANOVA_p<p,]
        statistics<-statistics[statistics$ANOVA_padj<p,]
        groups[[1]]<-rownames(data)
        names(groups)[1]<-paste0("ANOVA_padj<",p)
    }}

  #Filter based on statCol_filter
  if(statCol_filter!="none"){
    if(missing(statCol_text)){
      warning("the stat <statCol_text> was missing the stat column was not filtered")
    }else{
      text<-paste0("statistics$`",statCol_filter,"`",statCol_text)
      message("the following filter was applied:")
      message(text)
      test<-eval(parse(text=text))
      test[is.na(test)]<-FALSE
      if(!is.logical(test) | length(test)!= nrow(statistics)){
        message("The result of following test is not a logical vector of the same length as the dataset :")
        message(test)
      }else{
      data <- data[test,]
      statistics <-statistics[test,]
      groups[[1]]<-rownames(data)
      names(groups)[1]<-text
      }}}

  if(statCol_filter2!="none"){
    if(missing(statCol_text2)){
      warning("the stat <statCol_text2> was missing the stat column was not filtered")
    }else{
      text<-paste0("statistics$`",statCol_filter2,"`",statCol_text2)
      message("the following filter was applied:")
      message(text)
      test<-eval(parse(text=text))
      test[is.na(test)]<-FALSE
      if(!is.logical(test) | length(test)!= nrow(statistics)){
        message("The result of following test is not a logical vector of the same length as the dataset :")
        message(test)
      }else{
        data <- data[test,]
        statistics <-statistics[test,]
        groups[[1]]<-rownames(data)
        names(groups)[1]<-paste0(names(groups)[1]," & ",text)
              }}}

  if(statCol_filter3!="none"){
    if(missing(statCol_text3)){
      warning("the stat <statCol_text3> was missing the stat column was not filtered")
    }else{
      text<-paste0("statistics$`",statCol_filter3,"`",statCol_text3)
      message("the following filter was applied:")
      message(text)
      test<-eval(parse(text=text))
      test[is.na(test)]<-FALSE
      if(!is.logical(test) | length(test)!= nrow(statistics)){
        message("The result of following test is not a logical vector of the same length as the dataset :")
        message(test)
      }else{
        data <- data[test,]
        statistics <-statistics[test,]
        groups[[1]]<-rownames(data)
        names(groups)[1]<-paste0(names(groups)[1]," & ",text)
        }}}

  #if cluster column
  if (cluster_filter!="none"){
    clusters<-as.numeric(t(statistics[,colnames(statistics)==cluster_filter]))
    data<-data[!is.na(clusters),]
    clusters<- clusters[!is.na(clusters)]
  for(i in 1:max(clusters)){
      groups[[i]]<- rownames(data)[clusters==i]
      names(groups)[i]<-paste0("Cluster_",i)
    }
  }

  if(enrichment_function=="Fisher"){
  #create an enrichment result table
  enrichment_results<-data.frame(matrix(ncol=12,nrow=0))

  #Enrichments using proteinMinion
  for (i in 1:length(groups)){
    query<-as.character(t(groups[[i]]))

    if("GO" %in% type){
      enrichment_table<-UniProt_GO_Fisher(query,universe,organismID = organismID , ...)
      enrichment_table<-cbind(enriched_in=names(groups)[i],enrichment_table)
      enrichment_results<-rbind(enrichment_results,enrichment_table)
      }

    if("KEGG" %in% type){
      enrichment_table<-UniProt_KEGG_Fisher(query,universe,organismID = organismID , ...)
      enrichment_table<-cbind(enriched_in=names(groups)[i],enrichment_table)
      enrichment_results<-rbind(enrichment_results,enrichment_table)
      }

    if("REACTOME" %in% type){
      enrichment_table<-UniProt_REACTOME_Fisher(query,universe,organismID = organismID , ...)
      enrichment_table<-cbind(enriched_in=names(groups)[i],enrichment_table)
      enrichment_results<-rbind(enrichment_results,enrichment_table)
    }
  }
  }

  if(enrichment_function=="EASE"){
    #create an enrichment result table
    enrichment_results<-data.frame(matrix(ncol=12,nrow=0))

    #Enrichments using proteinMinion
    for (i in 1:length(groups)){
      query<-as.character(t(groups[[i]]))

      if("GO" %in% type){
        enrichment_table<-UniProt_GO_EASE(query,universe,organismID = organismID , ...)
        enrichment_table<-cbind(enriched_in=names(groups)[i],enrichment_table)
        enrichment_results<-rbind(enrichment_results,enrichment_table)
      }

      if("KEGG" %in% type){
        enrichment_table<-UniProt_KEGG_EASE(query,universe,organismID = organismID , ...)
        enrichment_table<-cbind(enriched_in=names(groups)[i],enrichment_table)
        enrichment_results<-rbind(enrichment_results,enrichment_table)
      }

      if("REACTOME" %in% type){
        enrichment_table<-UniProt_REACTOME_EASE(query,universe,organismID = organismID , ...)
        enrichment_table<-cbind(enriched_in=names(groups)[i],enrichment_table)
        enrichment_results<-rbind(enrichment_results,enrichment_table)
      }
    }
  }


  return(enrichment_results)
  }
