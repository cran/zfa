#' Zoom-Focus Algorithm
#'
#' @description This function performs the Zoom-Focus Algorithm (ZFA) to locate optimal testing regions for rare variant association tests and performs the test based on the optimized regions. The package is suitable to be applied on sequencing data set that is composed of variants with minor allele frequency less than 0.01 (rare variants). The package calls existing rare variant test functions to conduct rare variant test.
#'
#' ZFA consists of two steps: Zooming and Focusing.  In the first step Zooming, a given genomic region is partitioned by an order of two, and the best partition is located using multiple testing corrected p-values returned by desired rare variant test. In the second step Focusing, the boundaries of the zoomed region are refined by allowing them to expand or shrink at micro-level. The computation complexity is linear to the number of variants for Zooming (when the option fast.path=FALSE); a fast-Zoom version can further reduce the complexity to the logarithm of data size (fast.path=TRUE, default).
#'
#' @param data a data frame or numeric matrix. Genotypes should be coded as 0,1 or 2.
#' @param y a numeric vector with two levels. Phenotype values are coded as 0 or 1.
#' @param bin a numeric integer taking value of power of two, namely, 2, 4, 8, 16, 32, 64, 128, 256, 512 etc. The bin size specifies the initial window size P to perform the Zoom-Focus Algorithm. Default bin =256.
#' @param fast.path a logical value indicating whether or not to use the fast-Zoom approach. The fast-Zoom performs a binary search instead of exhaustive search, such that at each partition order, the region is divided into two parts, only the part with smaller p-value is continued for the next level search. Default = TRUE.
#' @param filter.pval a p-value threshold to select zoomed regions for conducting focusing step. When specified, only zoomed regions with p-value smaller than the threshold will be passed to focusing step. Default=0.01. Set filter.pval=NULL for conducting focusing step with all the zoomed regions.
#' @param output.pval a p-value threshold for filtering the output. If set NULL, all the results will be listed; otherwise, the function will only output the regions with p-values smaller than output.pval. Default=0.05.
#' @param CommonRare_Cutoff MAF cutoff to define common and rare variants. Default=0.01.
#' @param test a character to choose the rare variant method that combines with the Zoom-Focus Algorithm.If test = "SKAT", the SKAT of variance component test is applied. If test = "SKATO", the SKAT-O of combination method test is applied. If test = "burden", the weighted burden test is applied. If test = "wtest", the W-test of burden test category is applied.
#' @return The \code{"zfa"} function returns a list with the following components:
#'
#' \item{n.regions}{Total number of regions to which the input genotype data is divided by initial bin size P.}
#'
#' \item{n.rare}{Total number of rare variants used for the analysis.}
#'
#' \item{n.common}{Total number of common variants excluded in the analysis.}
#'
#' \item{results}{The testing results consist of several elements: 1) \code{lower.bound and upper.bound} represent variant information which indicates the lower and upper bound of optimized testing region; 2) \code{opt.region.size} denotes the size of optimal testing region after performing ZFA; 3) \code{corrected.pvalue} displays the multiple testing (Bonferroni) corrected p-value of the optimal testing region.}
#'
#' \item{Bon.sig.level}{Suggested Bonferroni corrected significance level for the input data at threshold alpha=0.05, which equals to 0.05 / # of regions.}
#'
#' \item{variants}{The variants contained in each output optimal region.}
#'
#' Note that the \code{variants} in the optimized region will not be printed in default. User can can extract  the information by calling details of results. See an example in \strong{Examples} section.
#'
#' The zfa optimizes the testing region according to input variant sequence and assumes that they are arranged by chromosome positions. If variants from non-adjacent genomic regions are input as one data, the zfa will still treat them as adjacent. In this case, the user should be careful in interpreting the results: when an optimized region consists of distant variants, the region may not be biologically meaningful; when an optimized region consists of variants from two neighboring genes, the results may be meaningful.
#'
#' @details The algorithm divides sequencing data into multiple fixed genomic regions with a certain initial bin size. ZFA is conducted in each bin. Current version includes 4 existing rare variant tests. The SKAT, SKAT-O and weighted burden test are called from \code{'SKAT'} package, and the W-test Collapsing method is self-contained.
#'
#' @author Haoyi Weng & Maggie Wang
#'
#' @references M. H. Wang., H. Weng., et al. (2017) A Zoom-Focus algorithm (ZFA) to locate the optimal testing region for rare variant association tests. Bioinformatics.
#'
#' M. H. Wang., R. Sun., et al. (2016) A fast and powerful W-test for pairwise epistasis testing. Nucleic Acids Research.doi:10.1093/nar/gkw347.
#'
#' R. Sun., H. Weng., et al. (2016) A W-test collapsing method for rare-variant association testing in exome sequencing data. Genetic Epidemiology, 40(7): 591-596.
#'
#' Wu, M. C., Lee, S., et al. (2011) Rare Variant Association Testing for Sequencing Data Using the Sequence Kernel Association Test (SKAT). American Journal of Human Genetics, 89, 82-93.
#'
#' Lee, S., Emond, M.J., et al. (2012) Optimal unified approach for rare variant association testing with application to small sample case-control whole-exome sequencing studies. American Journal of Human Genetics, 91, 224-237.
#'
#' Lee, S., Fuchsberger, C., et al. (2015) An efficient resampling method for calibrating single and gene-based rare variant association analysis in case-control studies. Biostatistics, kxv033.
#'
#' @examples
#' data(zfa.example)
#' attach(zfa.example)
#'
#' # fast-zoom with wtest, all zoomed regions passed to focusing step, and output all results
#' zfa.result1<-zfa(X,y,bin = 32,fast.path = TRUE,filter.pval=NULL,output.pval=NULL,test = "wtest")
#'
#' # zooming with wtest, select zoomed regions for focusing and output regions with both p-value<0.01
#' zfa.result2<-zfa(X,y,bin = 32,fast.path = FALSE,filter.pval=0.01,output.pval=0.01,test = "wtest")
#'
#' ## an example to view the detail of variants in each output optimal region
#' result1.detail<-zfa.result1$variants
#'
#' @export
#' @importFrom stats pchisq

zfa<-function(data,y,bin=256,fast.path=TRUE,filter.pval=0.01,output.pval=0.05,CommonRare_Cutoff=0.01,test=c("SKAT","SKATO","burden","wtest")) {
  if (is.data.frame(y))
    y<-as.matrix(y)
  if (any(is.na(data)))
    stop ("NA occurs in data")
  if (!all(data %in% c(0,1,2)))
    stop ("all the genotypes in 'data' must be 0, 1 or 2")
  n.col<-ncol(data)
  if (n.col<2)
    stop ("'data' must contain at least two variants")
  if (any(is.na(y)))
    stop ("NA occurs in y")
  if (!all(y %in% c(0,1)))
    stop ("all the phenotypes in 'y' must be 0 or 1")
  if (length(y)!=nrow(data))
    stop ("'data' and 'y' must have the same length")
  if (!bin %in% 2^c(1:10))
    stop ("initial fixed window size must be an integer which belongs to the power of two, such as 4, 18 , 16, 32, 64, 128, 256, ...")
  if (CommonRare_Cutoff>0.5)
    stop ("The MAF cutoff should be a numeric value between 0 and 0.5")
  maf.data<-maf(data)
  if (any(maf.data>CommonRare_Cutoff)) {
    data<-data[,-c(which(maf.data>CommonRare_Cutoff))]
    n.common<-length(which(maf.data>CommonRare_Cutoff))
    n.rare<-n.col-n.common
    warning ("Common variants found in 'data', which were excluded from the analysis")
  } else {
    n.common<-0
    n.rare<-n.col
  }
  if (ncol(data)<bin) {
    bin<-2^c(1:10)
    bin<-bin[which(ncol(data)-bin>0)]
    bin<-bin[which.min(ncol(data)-bin)]
  }
  num.fix.window<-floor(ncol(data)/bin)
  result.final<-NULL
  if (test!="wtest") {obj<-SKAT::SKAT_Null_Model(y ~ 1, out_type="D")}
  if (fast.path) {
    for (i in 1:num.fix.window) {
      data1<-data[,c((bin*(i-1)+1):(bin*i))]
      result<-zooming.fast(data1,y,obj,test)
      result.final<-data.frame(rbind(result.final,focusing(result,data1,y,fast.path,filter.pval,bin,i,obj,test)))
      rm(data1,result)
    }
  } else {
    for (i in 1:num.fix.window) {
      data1<-data[,c((bin*(i-1)+1):(bin*i))]
      result<-zooming.full(data1,y,obj,test)
      result.final<-data.frame(rbind(result.final,focusing(result,data1,y,fast.path,filter.pval,bin,i,obj,test)))
      rm(data1,result)
    }
  }
  if (bin*num.fix.window!=ncol(data)) {
    total.region<-num.fix.window+1
    tail<-c(bin*num.fix.window+1,ncol(data))
    if (test=="wtest") {
      result.tail<-cbind(tail[1],tail[2],tail[2]-tail[1]+1,wtest(data[,tail[1]:tail[2]],y))
    } else if (test=="SKAT") {
      result.tail<-cbind(tail[1],tail[2],tail[2]-tail[1]+1,SKAT::SKATBinary(as.matrix(data[,tail[1]:tail[2]]),obj)$p.value)
    } else if (test=="SKATO") {
      result.tail<-cbind(tail[1],tail[2],tail[2]-tail[1]+1,SKAT::SKATBinary(as.matrix(data[,tail[1]:tail[2]]),obj,method="SKATO")$p.value)
    } else if (test=="burden") {
      result.tail<-cbind(tail[1],tail[2],tail[2]-tail[1]+1,SKAT::SKATBinary(as.matrix(data[,tail[1]:tail[2]]),obj,method="Burden")$p.value)
    }
    colnames(result.tail)<-c("lower.bound","upper.bound","opt.region.size","corrected.pvalue")
    result.final<-data.frame(rbind(result.final,result.tail))
    rm(tail,result.tail)
  } else {
    total.region<-num.fix.window
  }
  if(!is.null(output.pval)) {
    l.output.pval<-which(result.final[,"corrected.pvalue"]<=output.pval)
    result.final<-result.final[l.output.pval,]
  }
  if (nrow(result.final)==0) {
    out<-list(n.regions=total.region,n.rare=n.rare,n.common=n.common,results=NULL,Bon.sig.level=0.05/total.region,variants=NULL)
  } else {
    marker.info<-colnames(data.frame(data))
    region.num<-paste("region",rownames(result.final),sep = " ")
    variant.detail<-apply(result.final,1,function(x) marker.info[c(x[1]:x[2])])
    if (length(region.num)==1) {
      variant.detail<-list(variant.detail)
    }
    names(variant.detail)<-c(region.num)
    result.final[,1]<-marker.info[result.final[,1]]
    result.final[,2]<-marker.info[result.final[,2]]
    out<-list(n.regions=total.region,n.rare=n.rare,n.common=n.common,results=result.final,Bon.sig.level=0.05/total.region,variants=variant.detail)
  }
  class(out)<-"zfa"
  out
}


