#' Format scientific exponent
#'
#'  rstatix::p_format does good job but still leaves scientific
#'  notation with an "e". This function converts the e to standard
#'  x10 and superscript exponent. html verion not tested.
#'  @param pvals vector of p values
#'  @param threshold threshold below which pvals should be shown in scientific
#'  notation (e.g. 1e-4 rather than 0.0001). Default is 0.001.
#'  @param html boolean to say whether to convert to html markdown (T) or plotmath notation (F)
#'  @return vector of formatted pvals
#'  @export
prettyExponents<-function(pvals,threshold=0.001,html=F){
  idx<-pvals<threshold
  pvals.f<-pvals
  pvals.f[idx]<-scales::scientific(pvals[pvals<threshold], digits=2)
  pvals.f[!idx]<-rstatix::p_format(pvals[pvals>=threshold], digits=2, accuracy=1e-32)
  if(html==T){
    pvals.f<-gsub("e","x10<sup>",pvals)
    pvals.f<-gsub("x10<sup>(.*)$","x10<sup>\\1</sup>", pvals)
  } else {
    pvals.f<-gsub("e","%*%10",pvals.f)
    pvals.f<-gsub("10(.*)$","10^\\1", pvals.f)
    pvals.f<-gsub("<","p<",pvals.f,fixed=T)
  }
  return(pvals.f)
}



