#' Decomposition of the MUMARINEX value into its three sub-indices
#' and representation as boxplots
#'
#' @title Decomposition of the MUMARINEX value into its three sub-indices
#' and representation as boxplots
#'
#' @description
#' Generates a graphical representation (boxplot) of MUMARINEX sub-indices to assess which component(s) contribute most strongly to the overall MUMARINEX score.
#'
#'
#' @seealso [diagnostic_tool()]
#'
#' @param x A data frame organized with samples in rows and species in columns.
#' @param g A vector of length \code{nrow(x)} indicating how the samples should be grouped (e.g., stations, treatments).
#' @param ref A logical or numeric vector identifying the reference row positions.
#' @param log A logical indicating whether the data must be log transformed.
#' @param fill Fill color of the boxplots.
#' @param border Border color of the boxplots.
#' @param main Main title of the plot.
#'
#' @returns A boxplot of length \code{g} representing the variations in the different MUMARINEX sub-indices (SCSR, CBCS, and SPI).
#' @export
#' @importFrom knitr kable
#' @importFrom graphics par
#' @importFrom graphics boxplot
#' @importFrom graphics legend
#' @importFrom graphics text
#' @importFrom graphics title
#' @importFrom graphics axis
#'
#' @examples
#' data("Simulated_data")
#' ref_idx<-1:10
#' stations<-matrix(unlist(strsplit(rownames(Simulated_data),".",fixed=TRUE)),ncol=2,byrow=TRUE)[,1]
#' decomplot(x=Simulated_data,g=stations,ref=ref_idx,log=FALSE,main="Simulated data")
decomplot<-function(x,g,ref,log=TRUE,fill=c("lightblue","#FFFFE0DF","#90F0909E"),border=c("#0080AB","#C7C700DF","#0091009E"),main=NULL){

  if(!(is.data.frame(x)||is.matrix(x))){
    stop("'x' must be a data frame or a matrix.")
  }
  if(!all(apply(x,2,is.numeric))){
    stop("All columns in 'x' must be numeric.")
  }

  if(!(is.factor(g)||is.vector(g))||length(g)!=nrow(x)){
    stop("'g' must be either a factor or a vector with length equal to nrow(x).")
  }

  if(!(is.logical(ref)||is.numeric(ref))){
    stop("'ref' must be either a logical or numeric vector.")
  }
  if(is.logical(ref)&&length(ref)!=nrow(x)) {
    stop("If 'ref' is logical, it must have length equal to nrow(x).")
  }
  if(is.numeric(ref)&&any(ref<1|ref>nrow(x))){
    stop("Numeric 'ref' indices must be between 1 and nrow(x).")
  }

  if(length(fill)!=3||length(border)!=3){
    stop("'fill' and 'border' must each be a character vector of length 3 (SCSR, CBCS, SPI).")
  }

  if(!is.null(main) && !is.character(main)){
    stop("'main' must be either NULL or a character string.")
  }

  subind<-mumarinex(x,ref,T,log=log)$subindices

  g<-factor(g,levels=unique(g))

  oldpar<-par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(fig = c(0, 1, 0, 1), mar=c(1,4.1,0.1,0))

  boxplot(subind$SCSR~g,
          boxfill=fill[1],border=border[1],
          whisklty = 3,pch=16,
          frame=F,xaxt="n",yaxt="n",ylab="",xlab="",ylim=c(0,1.1))

  boxplot(subind$CBCS~g,
          boxfill=fill[2],border=border[2],
          whisklty = 3,pch=16,
          frame=F,xaxt="n",yaxt="n",ylab="",xlab="",ylim=c(0,1.1),
          add=T)

  boxplot(subind$SPI~g,
          boxfill=fill[3],border=border[3],
          whisklty = 3,pch=16,
          frame=F,xaxt="n",yaxt="n",ylab="",xlab="",ylim=c(0,1.1),
          add=T)


  legend("topleft",c("SCSR","CBCS","SPI"),
         col=c(fill[1],fill[2],fill[3]),pch=15,text.col="white",
         horiz=T,bty="n",xpd=T,adj=c(0,0.3),text.font=2,cex=1.4)->leg_cord
  legend("topleft",c("SCSR","CBCS","SPI"),
         col=c(border[1],border[2],border[3]),pch=0,
         horiz=T,bty="n",xpd=T,adj=c(0,0.4),text.font=2,cex=1.4)

  text(par("usr")[2]*0.95,leg_cord$text$y[1],labels=main,font=2,cex=1.4,xpd=T,adj=c(1,0.4))

  title(ylab=expression(bold("MUMARINEX")),line=2.8,cex.lab=1.8)
  axis(1,font=2,at=seq(1:length(unique(g))),labels = unique(g),lty=0,line=-1.2,cex.axis=1.4)
  axis(2,font=2,cex.axis=1.4,las=2)

}
