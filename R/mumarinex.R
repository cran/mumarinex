#' mumarinex
#'@description
#'Computes the MUMARINEX index and its sub-indices (SCSR, CBCS, and SPI), following the method described in Chauvel et al. (2026).
#'
#' @references Chauvel, N., Pezy, J.P., Grall, J., Thiébaut, E. (2026). A general-purpose Multivariate Marine Recovery Index for quantifying the influence of human activities on benthic habitat ecological status. Ecological Indicator.
#'
#' @title MUMARINEX and subindices computation
#'
#' @param x A data frame or a matrix organized with samples in rows and species in columns.
#' @param ref A logical or numeric vector identifying the reference row positions.
#' @param subindices A logical indicating whether the sub-indices (SCSR, CBCS, and SPI) should be returned.
#' @param log A logical indicating whether the data must be log transformed.
#'
#' @returns A data frame with the MUMARINEX values. When subindices = TRUE, an additional data frame is returned containing the results of the sub-indices (SCSR, CBCS, and SPI).
#' @importFrom vegan vegdist
#' @importFrom stats median
#'
#' @export
#'
#' @examples
#' data("Simulated_data")
#' ref_idx<-1:10
#' mumarinex(x=Simulated_data,ref=ref_idx,log=FALSE)
#' mumarinex(x=Simulated_data,ref=ref_idx,subindices=TRUE,log=FALSE)
mumarinex<-function(x,ref,subindices=FALSE,log=TRUE){

  if(!(is.data.frame(x)||is.matrix(x))){
    stop("'x' must be a data frame or a matrix.")
  }

  if(!all(apply(x, 2, is.numeric))){
    stop("All columns in 'x' must be numeric.")
  }

  if(!(is.logical(ref)||is.numeric(ref))){
    stop("'ref' must be either a logical or a numeric vector.")
  }

  if(is.logical(ref)){
    if(length(ref)!=nrow(x)){
      stop("If 'ref' is logical, it must have length equal to nrow(x).")
    }
  }

  if (is.numeric(ref)){
    if(any(ref<1|ref>nrow(x))){
      stop("Numeric 'ref' indices must be between 1 and nrow(x).")
    }
  }

  if(!is.logical(subindices)||length(subindices)!=1){
    stop("'subindices' must be a single logical value (TRUE or FALSE).")
  }

  if(log==T){
    xBC<-log1p(x)
  }else{
    xBC<-x
  }
  BC_dist_df<-data.frame(as.matrix(vegdist(xBC,method="bray",diag=T,upper=T)))
  colnames(BC_dist_df)<-rownames(x);rownames(BC_dist_df)<-rownames(x)
  diag(BC_dist_df)<-NA

  REF_dist_df<-BC_dist_df[ref,ref]
  mDR<-median(unlist(REF_dist_df),na.rm=T)

  REF_vs_TEST_df_corre<-BC_dist_df[ref,]-mDR
  REF_vs_TEST_df_corre[REF_vs_TEST_df_corre<0]<-0
  CBCS<-colMeans(REF_vs_TEST_df_corre,na.rm=T)

  df_sp_REF<-data.frame(x[,which(colSums(x[ref,])>0)])
  df_sp_NEW<-data.frame(x[,which(colSums(x[ref,])==0)])
  SR<-apply(df_sp_REF,1,function(x){sum(x>0)})
  NR<-apply(df_sp_NEW,1,function(x){sum(x>0)})
  SCSR<-SR-NR

  raw_J<-apply(x,1, function(x){
    pi<-x/sum(x)
    -sum(pi*log(pi,2),na.rm=TRUE)/log(sum(x>0),2)
  })

  data_bio<-cbind.data.frame(SCSR=SCSR,SPI=raw_J,CBCS=CBCS)
  data_bio[is.na(data_bio)]<-0

  mS<-median(data_bio$SCSR[ref])
  mJ<-median(data_bio$SPI[ref])

  data_bio_standard<-data.frame(SCSR=data_bio$SCSR/mS,SPI=data_bio$SPI/mJ,CBCS=(1-data_bio$CBCS))
  data_bio_standard$SCSR[data_bio_standard$SCSR>1]<-1;data_bio_standard$CBCS[data_bio_standard$CBCS>1]<-1;data_bio_standard$SPI[data_bio_standard$SPI>1]<-1
  data_bio_standard$SCSR[data_bio_standard$SCSR<0]<-0;data_bio_standard$CBCS[data_bio_standard$CBCS<0]<-0;data_bio_standard$SPI[data_bio_standard$SPI<0]<-0

  MUMARINEX<-data_bio_standard$SCSR*data_bio_standard$SPI*data_bio_standard$CBCS

  if(subindices==T){
    return(list(MUMARINEX=MUMARINEX,subindices=data_bio_standard))
  }else{
    return(MUMARINEX)
  }
}
