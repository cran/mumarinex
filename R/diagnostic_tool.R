#' diagnostic_tool
#'
#' @title Diagnostic tool to identify the key components that best explain the MUMARINEX sub-indices
#'
#' @seealso [decomplot()]
#'
#' @description
#' Identifies, for each sub-index, the species or taxa that contribute most to its variation.
#'
#' @note
#' To reduce the risk of misleading interpretations, a one-sided significance test (t-test) is applied to sub-indices against the reference condition. Nevertheless, taxa identified as contributing to sub-index and MUMARINEX variations may not always be ecologically relevant, and the results are provided without any guarantee. This tool is not a substitute for a thorough ecological knowledge of the studied site and careful examination of the data, although it may help guide users toward potential hypotheses. The significance tests can be disabled, but the resulting outputs should then be interpreted with extreme caution.
#'
#' @param x A data frame or a matrix organized with samples in rows and species in columns.
#' @param g A vector of length \code{nrow(x)} indicating how the samples should be grouped (e.g., stations, treatments).
#' @param ref A logical or numeric vector identifying the reference row positions.
#' @param signif_test Logical; if TRUE, only sub-indices significantly lower than the reference conditions (t-test, p < 0.05) are returned. Conditions that are not statistically significant are indicated by 'ns'.
#' @param mute A logical indicating whether the results are displayed in the console.
#'
#' @returns A data frame summarizing the key information explaining variations in CSR, CBCD, and CPI.
#' @export
#'
#' @importFrom stats sd
#' @importFrom stats t.test
#' @importFrom utils head
#'
#' @include mumarinex.R
#'
#' @examples
#' data("Simulated_data")
#' ref_idx<-41:50
#' stations<-matrix(unlist(strsplit(rownames(Simulated_data),".",fixed=TRUE)),ncol=2,byrow=TRUE)[,1]
#' diagnostic_tool(x=Simulated_data,g=stations,ref=ref_idx)
diagnostic_tool<-function(x,g,ref,signif_test=TRUE,mute=FALSE){

  if(!(is.data.frame(x)||is.matrix(x))){
    stop("'x' must be a data frame or a matrix.")
  }

  if(!all(apply(x,2,is.numeric))){
    stop("All columns in 'x' must be numeric.")
  }

  if(!(is.factor(g)||is.vector(g))||length(g)!=nrow(x)){
    stop("'g' must be a vector of length equal to nrow(x).")
  }

  if(!(is.logical(ref)||is.numeric(ref))){
    stop("'ref' must be either a logical or numeric vector.")
  }

  if(is.logical(ref)){
    if(length(ref) != nrow(x)){
      stop("If 'ref' is logical, it must have length equal to nrow(x).")
    }
  }

  if(is.numeric(ref)){
    if(any(ref<1|ref>nrow(x))){
      stop("Numeric 'ref' indices must be between 1 and nrow(x).")
    }
  }

  if(!is.logical(signif_test)||length(signif_test)!= 1){
    stop("'signif_test' must be a single logical value (TRUE or FALSE).")
  }

  if(!is.logical(mute)||length(mute)!=1){
    stop("'mute' must be a single logical value (TRUE or FALSE).")
  }

  subind<-mumarinex(x,ref,T)$subindices

  ref_fact<-rep(NA,nrow(subind));ref_fact[ref]<-"R";ref_fact[-ref]<-as.character(g[-ref]);ref_fact<-factor(ref_fact)
  ref_fact<-ref_fact[-ref]

  test_idx<-seq(nrow(subind))[-ref]
  ref_idx<-seq(nrow(subind))[ref]

  if(signif_test==T){

    pCSR<-data.frame(NULL)
    for(i in seq(unique(ref_fact))){
      station<-as.character(unique(ref_fact)[i])
      if(sd(subind$CSR[ref_fact==unique(ref_fact)[i]])!=0)
      {p<-t.test(subind$CSR[ref_fact==unique(ref_fact)[i]],subind$CSR[ref],alternative = "less")$p.value}else{p<-1}
      lin<-c(station,p)
      pCSR<-rbind(pCSR,lin)
      colnames(pCSR)<-c("station","p")
      pCSR$p<-as.numeric(pCSR$p)
    }

    pCBCD<-data.frame(NULL)
    for(i in seq(unique(ref_fact))){
      station<-as.character(unique(ref_fact)[i])
      if(sd(subind$CBCD[ref_fact==unique(ref_fact)[i]])!=0)
      {p<-t.test(subind$CBCD[ref_fact==unique(ref_fact)[i]],subind$CBCD[ref],alternative = "less")$p.value}else{p<-1}
      lin<-c(station,p)
      pCBCD<-rbind(pCBCD,lin)
      colnames(pCBCD)<-c("station","p")
      pCBCD$p<-as.numeric(pCBCD$p)
    }

    pCPI<-data.frame(NULL)
    for(i in seq(unique(ref_fact))){
      station<-as.character(unique(ref_fact)[i])
      if(sd(subind$CPI[ref_fact==unique(ref_fact)[i]])!=0)
      {p<-t.test(subind$CPI[ref_fact==unique(ref_fact)[i]],subind$CPI[ref],alternative = "less")$p.value}else{p<-1}
      lin<-c(station,p)
      pCPI<-rbind(pCPI,lin)
      colnames(pCPI)<-c("station","p")
      pCPI$p<-as.numeric(pCPI$p)
    }

  }else{
    pCSR<-data.frame(station=unique(ref_fact),p=0)
    pCBCD<-data.frame(station=unique(ref_fact),p=0)
    pCPI<-data.frame(station=unique(ref_fact),p=0)
  }

  CSR_df<-NULL

  if(min(pCSR$p)<0.05){

    xref<-x[ref_idx,which(colSums(x[ref_idx,])>0)]
    ref_specif<-apply(xref,2,function(x){sum(x>0)/length(x)});ref_fidel<-colSums(xref)/sum(xref);indval_ref<-ref_specif*ref_fidel
    xref<-xref[,names(sort(indval_ref,decreasing=T))]
    ref_sp<-colnames(xref)

    for(i in seq(unique(g[test_idx]))){

      station<-as.character(unique(g[test_idx])[i])

      if(pCSR[which(pCSR$station==station),"p"]<0.05|signif_test==F){

        df<-x[which(g==station),]

        station_sp<-colnames(df)[which(colSums(df)>0)]
        sp_comm<-station_sp[which(station_sp%in%ref_sp)]
        sp_new<-setdiff(station_sp,ref_sp)

        if(length(sp_new)>1){
          new_specif<-apply(df[,sp_new],2,function(x){sum(x>0)/length(x)})
          new_fidel<-colSums(df[,sp_new])/sum(df[,sp_new])
          indval_new<-new_specif*new_fidel
          sp_new<-names(sort(indval_new,decreasing=T))
        }

        sp_missing_raw<-setdiff(ref_sp,station_sp)

        missing<-NULL
        for(j in 1:length(unique(g[ref_idx]))){
          df_ref<-x[which(g==unique(g[ref])[j]),]
          df_ref<-df_ref[,which(colSums(df_ref)>0)]
          ref_temp_sp<-names(df_ref)
          sp_missing_temp<-setdiff(ref_temp_sp,station_sp)
          missing<-c(missing,length(sp_missing_temp))
        }

        new_species<-head(sp_new,5);new_species[5]<-NA;new_species[6]<-"----------------------";new_species[is.na(new_species)]<-"/"
        Station<-rep("",5);Station[3]<-station;Station[6]<-"------"
        Mean_missing<-rep("",5);Mean_missing[3]<-round(mean(missing),1);Mean_missing[6]<-"------";Mean_missing[is.na(Mean_missing)]<-"/"
        N_missing<-rep("",5);N_missing[3]<-length(sp_missing_raw);N_missing[6]<-"------";N_missing[is.na(N_missing)]<-"/"
        Missing_species<-head(sp_missing_raw,5);Missing_species[6]<-"----------------------";Missing_species[is.na(Missing_species)]<-"/"

        lin<-data.frame(Sample=Station,Raw=N_missing,Mean=Mean_missing,Missing_species=Missing_species,New_species=new_species)
        CSR_df<-rbind(CSR_df,lin);rownames(CSR_df)<-NULL

      }else{

        new_species<-"ns";new_species[2]<-"----------------------"
        Station<-station;Station[2]<-"------"
        Mean_missing<-"ns";Mean_missing[2]<-"------"
        N_missing<-"ns";N_missing[2]<-"------"
        Missing_species<-"ns";Missing_species[2]<-"----------------------"

        lin<-data.frame(Sample=Station,Raw=N_missing,Mean=Mean_missing,Missing_species=Missing_species,New_species=new_species)
        CSR_df<-rbind(CSR_df,lin)
      }
    }

  }else{

    CSR_df<-data.frame(Sample="ns",Raw="ns",Mean="ns",Missing_species="ns",New_species="ns")
  }

  if(mute==F){
    cat("\n","|-----------------------------------------------------------------------------------|","\n",sep="")
    cat("|--------------------------------- CSR  diagnostic ---------------------------------|","\n")
    cat("|-----------------------------------------------------------------------------------|","\n")
    cat("> Raw: Raw taxa difference between sample and reference pool","\n")
    cat("> Mean: Mean taxa difference between sample and reference pool","\n")
    cat("> Missing_species: Top 5 missing species (sorted by IndVal of the reference)","\n")
    cat("> New_species: Top 5 new species (sorted by IndVal of the sample)","\n")
    print(kable(CSR_df,caption="CSR diagnostic",format="simple"))
  }

  if(min(pCBCD$p)<0.05){

    xref<-x[ref_idx,]

    CBCD_df<-NULL

    for(i in seq(unique(g[test_idx]))){

      station<-as.character(unique(g[test_idx])[i])

      if(pCBCD[which(pCBCD$station==station),"p"]<0.05|signif_test==F){

        df<-x[which(g==station),]

        moy_diff_ab<-NULL
        for(c in 1:ncol(df)){
          moy_diff_c<-mean(outer(df[,c],xref[,c],function(x,y)x-y))
          moy_diff_ab<-c(moy_diff_ab,moy_diff_c)
        }
        names(moy_diff_ab)<-colnames(df)
        moy_diff_ab_increase<-sort(moy_diff_ab[which(moy_diff_ab>0)],decreasing = T)
        moy_diff_ab_decrease<-sort(moy_diff_ab[which(moy_diff_ab<0)])

        Station<-rep("",5);Station[3]<-station;Station[6]<-"------"
        lower_ab<-head(names(moy_diff_ab_decrease),5);lower_ab[6]<-"----------------------";lower_ab[is.na(lower_ab)]<-"/"
        Mean_abund_decrease<-round(head(moy_diff_ab_decrease,5),1);Mean_abund_decrease[6]<-"------";Mean_abund_decrease[is.na(Mean_abund_decrease)]<-"/"
        higher_ab<-head(names(moy_diff_ab_increase),5);higher_ab[6]<-"----------------------";higher_ab[is.na(higher_ab)]<-"/"
        Mean_abund_increase<-round(head(moy_diff_ab_increase,5),1);Mean_abund_increase[6]<-"------";Mean_abund_increase[is.na(Mean_abund_increase)]<-"/"

        lin<-data.frame(Sample=Station,Lower_abundance=lower_ab,Decrease=Mean_abund_decrease,Higher_abundance=higher_ab,Increase=Mean_abund_increase)
        CBCD_df<-rbind(CBCD_df,lin);rownames(CBCD_df)<-NULL

      }else{

        Station<-station<-station;Station[2]<-"------"
        lower_ab<-"ns";lower_ab[2]<-"----------------------"
        Mean_abund_decrease<-"ns";Mean_abund_decrease[2]<-"------"
        higher_ab<-"ns";higher_ab[2]<-"----------------------"
        Mean_abund_increase<-"ns";Mean_abund_increase[2]<-"------"

        lin<-data.frame(Sample=Station,Lower_abundance=lower_ab,Decrease=Mean_abund_decrease,Higher_abundance=higher_ab,Increase=Mean_abund_increase)
        CBCD_df<-rbind(CBCD_df,lin);rownames(CBCD_df)<-NULL
      }
    }

  }else{
    CBCD_df<-data.frame(Sample="ns",Lower_abundance="ns",Decrease="ns",Higher_abundance="ns",Increase="ns")
  }

  if(mute==F){
    cat("\n","|-----------------------------------------------------------------------------------|","\n",sep="")
    cat("|--------------------------------- CBCD diagnostic ---------------------------------|","\n")
    cat("|-----------------------------------------------------------------------------------|","\n")
    cat("> Lower_abundance: Important reference taxa which present lower abundances","\n")
    cat("> Decrease: Mean decrease (vs reference) of the corresponding taxa","\n")
    cat("> Higher_abundance: Important reference taxa which present higher abundances","\n")
    cat("> Increase: Mean increase (vs reference) of the corresponding taxa","\n")

    print(kable(CBCD_df,format="simple"))
  }

  if(min(pCPI$p)<0.05|signif_test==F){

    xref<-x[ref_idx,]
    xref_mean<-colMeans(xref)

    CPI_df<-NULL

    for(i in seq(unique(g[test_idx]))){

      station<-as.character(unique(g[test_idx])[i])

      if(pCPI[which(pCPI$station==station),"p"]<0.05|signif_test==F){
        df<-x[which(g==station),]
        df_corrected<-df-rep(xref_mean,each=nrow(df))
        df_corrected[df_corrected<0]<-0
        df_sorted<-df_corrected[,order(colSums(df_corrected),decreasing=T)]


        Station<-rep("",5);Station[3]<-station;Station[6]<-"------"
        domin_sp<-colnames(df_sorted)[1:5];domin_sp[6]<-"----------------------";domin_sp[is.na(domin_sp)]<-"/"
        contrib<-head((colSums(df_sorted)/sum(df_sorted))*100,5);contrib[6]<-"------";contrib[is.na(contrib)]<-"/"

        lin<-data.frame(Sample=Station,Dominant_species=domin_sp,Contribution=contrib)
        CPI_df<-rbind(CPI_df,lin);rownames(CPI_df)<-NULL

      }else{

        Station<-station<-station;Station[2]<-"------"
        domin_sp<-"ns";domin_sp[2]<-"----------------------"
        contrib<-"ns";contrib[2]<-"------"

        lin<-data.frame(Sample=Station,Dominant_species=domin_sp,Contribution=contrib)
        CPI_df<-rbind(CPI_df,lin)
      }
    }

  }else{

    CPI_df<-data.frame(Sample="ns",Dominant_species="ns",Contribution="ns")

  }
  if(mute==F){
    cat("\n","|-----------------------------------------------------------------------------------|","\n",sep="")
    cat("|--------------------------------- CPI  diagnostic ---------------------------------|","\n")
    cat("|-----------------------------------------------------------------------------------|","\n")

    cat("> Dominant_species: most abundant species (corrected by reference)","\n")
    cat("> Contribution: Taxa contribution (%) to total abundance (corrected by reference)","\n")

    print(kable(CPI_df,caption="CPI diagnostic",format="simple"))
  }

  all_diag_results<-list(CSR=CSR_df,CBCD=CBCD_df,CPI=CPI_df)
  invisible(all_diag_results)
}

