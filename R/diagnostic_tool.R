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
#' @param log A logical indicating whether the data must be log transformed.
#' @param signif_test Logical; if TRUE, only sub-indices significantly lower than the reference conditions (t-test, p < 0.05) are returned. Conditions that are not statistically significant are indicated by 'ns'.
#' @param mute A logical indicating whether the results are displayed in the console.
#'
#' @returns A data frame summarizing the key information explaining variations in SCSR, CBCS, and SPI.
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
#' ref_idx<-1:10
#' stations<-matrix(unlist(strsplit(rownames(Simulated_data),".",fixed=TRUE)),ncol=2,byrow=TRUE)[,1]
#' diagnostic_tool(x=Simulated_data,g=stations,ref=ref_idx,log=FALSE)
diagnostic_tool<-function(x,g,ref,log=TRUE,signif_test=TRUE,mute=FALSE){

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

  subind<-mumarinex(x,ref,T,log=log)$subindices

  ref_fact<-rep(NA,nrow(subind));ref_fact[ref]<-"R";ref_fact[-ref]<-as.character(g[-ref]);ref_fact<-factor(ref_fact)
  ref_fact<-ref_fact[-ref]

  test_idx<-seq(nrow(subind))[-ref]
  ref_idx<-seq(nrow(subind))[ref]

  if(signif_test==T){

    pSCSR<-data.frame(NULL)
    for(i in seq(unique(ref_fact))){
      station<-as.character(unique(ref_fact)[i])
      if(sd(subind$SCSR[ref_fact==unique(ref_fact)[i]])!=0)
      {p<-t.test(subind$SCSR[ref_fact==unique(ref_fact)[i]],subind$SCSR[ref],alternative = "less")$p.value}else{
        ifelse(unique(subind$SCSR[ref_fact==unique(ref_fact)[i]])==unique(subind$SCSR[ref]),p<-1,p<-0)}
      lin<-c(station,p)
      pSCSR<-rbind(pSCSR,lin)
      colnames(pSCSR)<-c("station","p")
      pSCSR$p<-as.numeric(pSCSR$p)
    }

    pCBCS<-data.frame(NULL)
    for(i in seq(unique(ref_fact))){
      station<-as.character(unique(ref_fact)[i])
      if(sd(subind$CBCS[ref_fact==unique(ref_fact)[i]])!=0)
      {p<-t.test(subind$CBCS[ref_fact==unique(ref_fact)[i]],subind$CBCS[ref],alternative = "less")$p.value}else{p<-1}
      lin<-c(station,p)
      pCBCS<-rbind(pCBCS,lin)
      colnames(pCBCS)<-c("station","p")
      pCBCS$p<-as.numeric(pCBCS$p)
    }

    pSPI<-data.frame(NULL)
    for(i in seq(unique(ref_fact))){
      station<-as.character(unique(ref_fact)[i])
      if(sd(subind$SPI[ref_fact==unique(ref_fact)[i]])!=0)
      {p<-t.test(subind$SPI[ref_fact==unique(ref_fact)[i]],subind$SPI[ref],alternative = "less")$p.value}else{p<-1}
      lin<-c(station,p)
      pSPI<-rbind(pSPI,lin)
      colnames(pSPI)<-c("station","p")
      pSPI$p<-as.numeric(pSPI$p)
    }

  }else{
    pSCSR<-data.frame(station=unique(ref_fact),p=0)
    pCBCS<-data.frame(station=unique(ref_fact),p=0)
    pSPI<-data.frame(station=unique(ref_fact),p=0)
  }

  SCSR_df<-NULL

  if(min(pSCSR$p)<0.05){

    xref<-x[ref_idx,which(colSums(x[ref_idx,])>0)]
    ref_specif<-apply(xref,2,function(x){sum(x>0)/length(x)});ref_fidel<-colSums(xref)/sum(xref);indval_ref<-ref_specif*ref_fidel
    xref<-xref[,names(sort(indval_ref,decreasing=T))]
    ref_sp<-colnames(xref)

    for(i in seq(unique(g[test_idx]))){

      station<-as.character(unique(g[test_idx])[i])

      if(pSCSR[which(pSCSR$station==station),"p"]<0.05|signif_test==F){

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
        sp_new_raw<-setdiff(station_sp,ref_sp)
        sp_diff_taxa_raw<-c(sp_missing_raw,sp_new_raw)

        diff_taxa<-NULL
        for(j in 1:length(unique(g[ref_idx]))){
          df_ref<-x[which(g==unique(g[ref])[j]),]
          df_ref<-df_ref[,which(colSums(df_ref)>0)]
          ref_temp_sp<-names(df_ref)
          sp_missing_temp<-setdiff(ref_temp_sp,station_sp)
          sp_new_temp<-setdiff(station_sp,ref_temp_sp)
          sp_diff<-c(sp_missing_temp,sp_new_temp)
          diff_taxa<-c(diff_taxa,length(sp_diff))
        }

        new_species<-head(sp_new,5);new_species[6]<-NA;new_species[6]<-"----------------------";new_species[is.na(new_species)]<-"/"
        Station<-rep("",5);Station[3]<-station;Station[6]<-"------"
        Mean_diff_taxa<-rep("",5);Mean_diff_taxa[3]<-round(mean(diff_taxa),1);Mean_diff_taxa[6]<-"------";Mean_diff_taxa[is.na(Mean_diff_taxa)]<-"/"
        N_diff_taxa<-rep("",5);N_diff_taxa[3]<-length(sp_diff_taxa_raw);N_diff_taxa[6]<-"------";N_diff_taxa[is.na(N_diff_taxa)]<-"/"
        Missing_species<-head(sp_missing_raw,5);Missing_species[6]<-"----------------------";Missing_species[is.na(Missing_species)]<-"/"

        lin<-data.frame(Sample=Station,Raw=N_diff_taxa,Mean=Mean_diff_taxa,Missing_species=Missing_species,New_species=new_species)
        SCSR_df<-rbind(SCSR_df,lin);rownames(SCSR_df)<-NULL

      }else{

        new_species<-"ns";new_species[2]<-"----------------------"
        Station<-station;Station[2]<-"------"
        Mean_diff_taxa<-"ns";Mean_diff_taxa[2]<-"------"
        N_diff_taxa<-"ns";N_diff_taxa[2]<-"------"
        Missing_species<-"ns";Missing_species[2]<-"----------------------"

        lin<-data.frame(Sample=Station,Raw=N_diff_taxa,Mean=Mean_diff_taxa,Missing_species=Missing_species,New_species=new_species)
        SCSR_df<-rbind(SCSR_df,lin)
      }
    }

  }else{

    SCSR_df<-data.frame(Sample="ns",Raw="ns",Mean="ns",Missing_species="ns",New_species="ns")
  }

  if(mute==F){
    cat("\n","|-----------------------------------------------------------------------------------|","\n",sep="")
    cat("|--------------------------------- SCSR  diagnostic --------------------------------|","\n")
    cat("|-----------------------------------------------------------------------------------|","\n")
    cat("> Raw: Raw taxa difference between sample and reference pool","\n")
    cat("> Mean: Mean taxa difference between sample and reference pool","\n")
    cat("> Missing_species: Top 5 diff_taxa species (sorted by IndVal of the reference)","\n")
    cat("> New_species: Top 5 new species (sorted by IndVal of the sample)","\n")
    print(kable(SCSR_df,caption="SCSR diagnostic",format="simple"))
  }

  if(min(pCBCS$p)<0.05){

    xref<-x[ref_idx,]

    CBCS_df<-NULL

    for(i in seq(unique(g[test_idx]))){

      station<-as.character(unique(g[test_idx])[i])

      if(pCBCS[which(pCBCS$station==station),"p"]<0.05|signif_test==F){

        df<-x[which(g==station),]

        moy_diff_ab<-NULL
        moy_diff_rel<-NULL

        for(c in 1:ncol(df)){
          moy_diff_c<-mean(outer(df[,c],xref[,c],function(x,y)x-y))
          moy_diff_rel_c<-(mean(df[,c])/mean(xref[,c])*100)-100

          moy_diff_ab<-c(moy_diff_ab,moy_diff_c)
          moy_diff_rel<-c(moy_diff_rel,moy_diff_rel_c)
        }
        names(moy_diff_ab)<-colnames(df)
        names(moy_diff_rel)<-colnames(df)

        moy_diff_rel_increase<-sort(moy_diff_rel[which(moy_diff_rel>0&moy_diff_rel<Inf)],decreasing = T)
        moy_diff_rel_decrease<-sort(moy_diff_rel[which(moy_diff_rel<0&moy_diff_rel!=(-100))])
        moy_diff_ab_increase<-moy_diff_ab[names(moy_diff_rel_increase)]
        moy_diff_ab_decrease<-sort(moy_diff_ab[names(moy_diff_rel_decrease)])

        Station<-rep("",5);Station[3]<-station;Station[6]<-"------"
        lower_ab<-head(names(moy_diff_ab_decrease),5);lower_ab[6]<-"----------------------";lower_ab[is.na(lower_ab)]<-"/"
        Mean_abund_decrease<-round(head(moy_diff_ab_decrease,5),1);Mean_abund_decrease[6]<-"------";Mean_abund_decrease[is.na(Mean_abund_decrease)]<-"/"
        Rel_abund_decrease<-round(head(moy_diff_rel_decrease,5),1);Rel_abund_decrease[6]<-"------";Rel_abund_decrease[is.na(Rel_abund_decrease)]<-"/"
        higher_ab<-head(names(moy_diff_ab_increase),5);higher_ab[6]<-"----------------------";higher_ab[is.na(higher_ab)]<-"/"
        Mean_abund_increase<-round(head(moy_diff_ab_increase,5),1);Mean_abund_increase[6]<-"------";Mean_abund_increase[is.na(Mean_abund_increase)]<-"/"
        Rel_abund_increase<-round(head(moy_diff_rel_increase,5),1);Rel_abund_increase[6]<-"------";Rel_abund_increase[is.na(Rel_abund_increase)]<-"/"

        lin<-data.frame(Sample=Station,Lower_abundance=lower_ab,Decrease=Mean_abund_decrease,Relative_D=Rel_abund_decrease,Higher_abundance=higher_ab,Increase=Mean_abund_increase,Relative_I=Rel_abund_increase)
        CBCS_df<-rbind(CBCS_df,lin);rownames(CBCS_df)<-NULL

      }else{

        Station<-station<-station;Station[2]<-"------"
        lower_ab<-"ns";lower_ab[2]<-"----------------------"
        Mean_abund_decrease<-"ns";Mean_abund_decrease[2]<-"------"
        Rel_abund_decrease<-"ns";Rel_abund_decrease[2]<-"------"
        higher_ab<-"ns";higher_ab[2]<-"----------------------"
        Mean_abund_increase<-"ns";Mean_abund_increase[2]<-"------"
        Rel_abund_increase<-"ns";Rel_abund_increase[2]<-"------"

        lin<-data.frame(Sample=Station,Lower_abundance=lower_ab,Decrease=Mean_abund_decrease,Relative_D=Rel_abund_decrease,Higher_abundance=higher_ab,Increase=Mean_abund_increase,Relative_I=Rel_abund_increase)
        CBCS_df<-rbind(CBCS_df,lin);rownames(CBCS_df)<-NULL
      }
    }

  }else{
    CBCS_df<-data.frame(Sample="ns",Lower_abundance="ns",Decrease="ns",Relative_D="ns",Higher_abundance="ns",Increase="ns",Relative_I="ns")
  }

  if(mute==F){
    cat("\n","|-----------------------------------------------------------------------------------|","\n",sep="")
    cat("|--------------------------------- CBCS diagnostic ---------------------------------|","\n")
    cat("|-----------------------------------------------------------------------------------|","\n")
    cat("> Lower_abundance: Important reference taxa which present lower abundances","\n")
    cat("> Decrease: Mean decrease (vs reference) of the corresponding taxa","\n")
    cat("> Relative_D: Relative mean decrease (vs reference) of the corresponding taxa (%)","\n")
    cat("> Higher_abundance: Important reference taxa which present higher abundances","\n")
    cat("> Increase: Mean increase (vs reference) of the corresponding taxa","\n")
    cat("> Relative_I: Relative mean increase (vs reference) of the corresponding taxa (%)","\n")

    print(kable(CBCS_df,format="simple"))
  }

  if(min(pSPI$p)<0.05|signif_test==F){

    xref<-x[ref_idx,]
    xref_mean<-colMeans(xref)

    SPI_df<-NULL

    for(i in seq(unique(g[test_idx]))){

      station<-as.character(unique(g[test_idx])[i])

      if(pSPI[which(pSPI$station==station),"p"]<0.05|signif_test==F){
        df<-x[which(g==station),]
        df_corrected<-df-rep(xref_mean,each=nrow(df))
        df_corrected[df_corrected<0]<-0
        df_sorted<-df_corrected[,order(colSums(df_corrected),decreasing=T)]


        Station<-rep("",5);Station[3]<-station;Station[6]<-"------"
        domin_sp<-colnames(df_sorted)[1:5];domin_sp[6]<-"----------------------";domin_sp[is.na(domin_sp)]<-"/"
        contrib<-round(head((colSums(df_sorted)/sum(df_sorted))*100,5),1);contrib[6]<-"------";contrib[is.na(contrib)]<-"/"

        lin<-data.frame(Sample=Station,Dominant_species=domin_sp,Contribution=contrib)
        SPI_df<-rbind(SPI_df,lin);rownames(SPI_df)<-NULL

      }else{

        Station<-station<-station;Station[2]<-"------"
        domin_sp<-"ns";domin_sp[2]<-"----------------------"
        contrib<-"ns";contrib[2]<-"------"

        lin<-data.frame(Sample=Station,Dominant_species=domin_sp,Contribution=contrib)
        SPI_df<-rbind(SPI_df,lin)
      }
    }

  }else{

    SPI_df<-data.frame(Sample="ns",Dominant_species="ns",Contribution="ns")

  }
  if(mute==F){
    cat("\n","|-----------------------------------------------------------------------------------|","\n",sep="")
    cat("|--------------------------------- SPI  diagnostic ---------------------------------|","\n")
    cat("|-----------------------------------------------------------------------------------|","\n")

    cat("> Dominant_species: most abundant species (corrected by reference)","\n")
    cat("> Contribution: Taxa contribution (%) to total abundance (corrected by reference)","\n")

    print(kable(SPI_df,caption="SPI diagnostic",format="simple"))
  }

  all_diag_results<-list(SCSR=SCSR_df,CBCS=CBCS_df,SPI=SPI_df)
  invisible(all_diag_results)
}

