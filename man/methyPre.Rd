\name{methyPre}
\alias{methyPre}
\title{
Wrapper function for assembling comprehensive upstream pipelines for methylation array data.
}
\description{
This function is for conveniently specifying a pipeline to preprocess methylation array data using common parameters from the minfi package and from common manipulations of experiment summary objects. It is recommended this be used to supplement, but not replace, effective analysis workflows, ie. to compare the efficacies of different normalization methods. Methylation experiment objects that have already been preprocessed can still be conveniently run through this workflow, ie. to map probes to the genome, convert to a new experiment object class, or conveniently filter on recognized lists of cross-reactive probes.

}
\usage{
methyPre(meExpObj, workflow = c("norm", "pfilt", "map", "minfifilt", "ufilt", "crxcgfilt", "batch"),
  normfun = c("noob", "SWAN", "illumina", "funnorm"), pfilt = TRUE, detPcutoff = 0.05,
  minfiFilt = TRUE, uFilt = NULL, crxcgFilt = c(NULL, "hm450", "epic"), batchCorrect = FALSE,
  bcCovariate = NULL, bcBatch = NULL, returnObj = c("preObj1", "preObj2", "preObj3", "gexp1", "gexp2", "gexp3"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{meExpObj}{
 Experiment summary object containing methylation, annotation, and sample data.
}
  \item{workflow}{
 User-specified options (as character strings) outlining workflow steps to be taken. All or part of the possible options may be used, with the implicit order of operations being:

 1. norm (normalize, see normfun arg);

 2. pfilt (intensity p-val filter, see pfilt and detPcutoff args);

 3. map (map data to genome);

 4. minfifilt (use minfi filtering functions and annotation, see minfiFilt arg);

 5. ufilt (user-specified CpG IDs to filter, see uFilt);

 6. crxcgfilt(filter low-quality probes, see crxcgFilt arg);

 7. batch (perform ComBat batch correction, see batchCorrect, bcCovariate, and bcBatch variables).
}
  \item{normfun}{
 Normalization function(s) to be used. User can specify any combination of options: Illumina; SWAN; noob; funnorm. The implicit order of operations is background (Illumina or noob) followed by quantile (SWAN or funnorm) normalizations.
}
  \item{pfilt}{
 Whether to obtain intensity p-values for filtering. Requires that meExpObj be an RGChannelSet.
}
  \item{detPcutoff}{
 Mean p-value cutoff for intensity p-value. ie. CpGs removed if mean of p-value across samples is above the cutoff.
}
  \item{minfiFilt}{
 Whether to perform a series of filters using minfi functions and annotation: CpGs are filtered if they map to or are associated with SNPs, are CH/non-CpG probes, or map to X or Y chromosomes.
}
  \item{uFilt}{
 User-specified list of CpGs to filter out. Should be a vector or list of characters strings.
}
  \item{crxcgFilt}{
  Filter cross-reactive or poor-quality probes. Options are "hm450" or "epic", which filter either the HM450 CpGs specified in Chen et al 2013 (see ?chen_crxcg for info) or EPIC CpGs identified by Illumina (?illumina_crxcg) and Pidsley et al 2016 (?pidsley_crxcg) as being cross reactive, sensitive to genetic variation, or of otherwise poor quality.
}
  \item{batchCorrect}{
 Should batch correction with ComBat be performed?
}
  \item{bcCovariate}{
 Covariate variable (as vector or list) to be included in batch correction (ie. tissue type).
}
  \item{bcBatch}{
 Variable (vector or list) used to specify batches for correction (ie. array plate ID).
}
  \item{returnObj}{
 Objects to return, where gcombat is automatically returned if batch correction is specified. Options are:

 1. preObj1: Illumina, noob, or unnormalized experiment object.

 2. preObj2: preObj1 after SWAN, funnorm, or no additional normalization.

 3. preObj3: preObj2 after detection-p value or no filter applied.

 4. gexp1: preObj3 after mapping or no mapping performed.

 5. gexp2: gexp1 after minfi CpG filters applied or not applied.

 6. gexp3: gexp2 after user-defined and/or crxcg filters are applied or aren't applied.


}
}
\details{
This is a wrapper function for common workflows to preprocess Illumina methylation arrays using pipelines in the minfi package.

Note, this function is intended enable not just apply a pipeline efficiently, but a convenient way to compare different workflows.

}
\value{
List of experiment summary objects containing methylation, annotation, and sample data.

  \item{(preObj1|preObj2|preObj3|gexp1|gexp2|gexp3|gcombat)}{Methylation experiment summary object.}
 ...
}
\references{

}
\author{
Sean Maden, smaden@fredhutch.org
}
\note{
}

\seealso{
}
\examples{
## The function is currently defined as
function(meExpObj, workflow=c("norm","pfilt","map","minfifilt","ufilt","crxcgfilt","batch"),
                     normfun=c("noob","SWAN","illumina","funnorm"),
                     pfilt=TRUE, detPcutoff=0.05,minfiFilt=TRUE,uFilt=NULL,
                     crxcgFilt=c(NULL,"hm450","epic"), batchCorrect=FALSE,
                     bcCovariate=NULL,bcBatch=NULL,
                     returnObj=c("preObj1","preObj2","preObj3","gexp1","gexp2","gexp3")){
  require(minfi)
  return.list <- list()

  if(length(workflow[grepl("norm",workflow)])>0){
    if(length(normfun[grepl("illumina",normfun)])>0){
      message("starting Illumina normalization..")
      preObj1 <- preprocessIllumina(meExpObj)
      message("finished Illumina normalization!")
    } else if(length(normfun[grepl("noob",normfun)])>0){
      message("starting Noob normalization..")
      preObj1 <- preprocessNoob(meExpObj)
      message("finished Noob normalization!")
    } else{
      message("no background normalizations performed. Continuing..")
      preObj1 <- meExpObj
    }
    if(length(returnObj[grepl("preObj1",returnObj)])>0){
      return.list <- append(return.list,list(preObj1))
    }

    if(length(normfun[grepl("swan",normfun)])>0){
      message("beginning SWAN normalization..")
      preObj2 <- preprocessSWAN(meExpObj,mSet=preObj1)
      message("finished SWAN normalization!")
    } else if(length(normfun[grepl("funnorm",normfun)])>0){
      message("beginning Functional Normalization..")
      preObj2 <- preprocessFunnorm(meExpObj,bgCorr=FALSE)
      message("finished Functional Normalization!")
    } else{
      preObj2 <- preObj1
    }
  } else{
    preObj2 <- meExpObj
  }
  if(length(returnObj[grepl("preObj2",returnObj)])>0){
    return.list <- append(return.list,list(preObj2))
  }

  if(length(workflow[grepl("pfilt",workflow)])>0 & class(meExpObj)=="RGChannelSet"){
    message("getting intensity p-values from RG set...")
    detP <- detectionP(meExpObj)
    message(paste0("filtering on mean detection cutoff =",detPcutoff))
    failed <-rowMeans(detP) > detPcutoff
    preObj3 <- preObj2[!failed,]
  } else{
    preObj3 <- preObj2
  }
  if(length(returnObj[grepl("preObj3",returnObj)])>0){
    return.list <- append(return.list,list(preObj3))
  }

  if(length(workflow[grepl("map",workflow)])>0){
    gexp1 <- mapToGenome(preObj3)
  } else{
    gexp1 <- preObj3
  }
  if(length(returnObj[grepl("gexp1",returnObj)])>0){
    return.list <- append(return.list,list(gexp1))
  }

  if(length(workflow[grepl("minfifilt",workflow)])>0){
    message("Applying minfi filters..")
    message("Removing cg probes containing SNPs..")
    gexp2 <- dropLociWithSnps(gexp1, snps=c("SBE", "CpG"))
    message("Removing CH and SNP-assoc. probes...")
    gexp2 <- dropMethylationLoci(gexp2,dropRS=TRUE,dropCH=TRUE)
    message("Removing chrX and chrY-assoc. probes..")
    gexp2 <- gexp2[!getAnnotation(gexp2)$chr %in% c("chrY","chrX"),]
    message("After applying minfi filters, ",nrow(gexp2)," CpGs remain.")
  } else{
    gexp2 <- gexp1
  }
  if(length(returnObj[grepl("gexp2",returnObj)])>0){
    return.list <- append(return.list,list(gexp2))
  }

  if(length(workflow[grepl("ufilt",workflow)])>0|
     length(workflow[grepl("crxcgfilt",workflow)])>0){
    if(crxcgFilt=="hm450"){
      message("Filtering cross-reactive HM450 probes...")
      data(chen_crxcg)
      gexp3 <- gexp2[!rownames(gexp2) %in% chen.crxcg,]
    }
    if(crxcgFilt=="epic"){
      message("Filtering cross-reactive EPIC probes")
      data(illumina_crxcg);data(pidsley_crxcg)
      gexp3 <- gexp2[!rownames(gexp2) %in% c(illumina.crxcg,pidsley.crxcg),]
    }
    if(!is.null(uFilt)){
      message("Filtering user-specified cg list...")
      gexp3 <- gexp2[!rownames(gexp2) %in% uFilt,]
    }
    } else{
    gexp3 <- gexp2
    }
  if(length(returnObj[grepl("gexp3",returnObj)])>0){
    return.list <- append(return.list,list(gexp3))
  }

  if(batchCorrect){
    bc.check <- readline(message("You selected to batch correct. Do you want to save now or continue? Enter 'save' or 'continue'"))
    if(length(bc.check[grepl("save",bc.check)])>0){
      return(gexp3)
    }
    if(length(bc.check[grepl("continue",bc.check)])>0){
      require("sva")
      if(!is.null(bcCovariate)){
        message("Applying ComBat correction to M-values..")
        combat_M <- ComBat(getM(gexp3),mod=model.matrix(~bcCovariate),batch=bcBatch)
        message("Assembling corrected experiment object..")
        gcombat <- GenomicRatioSet(gr=granges(gexp3),
                                   Beta=ilogit2(combat_M),
                                   M=combat_M,
                                   annotation=annotation(gexp3),
                                   preprocessMethod=preprocessMethod(gexp3),
                                   pData=pData(gexp3))
        return.list <- append(return.list,gcombat)
        message("Completed all specified steps! Returning..")
        return(return.list)

      } else{
        message("Applying ComBat correction to M-values..")
        combat_M <- ComBat(getM(gexp3),batch=bcBatch)
        message("Assembling corrected experiment object...")
        gcombat <- GenomicRatioSet(gr=granges(gexp3),
                                   Beta=ilogit2(combat_M),
                                   M=combat_M,
                                   annotation=annotation(gexp3),
                                   preprocessMethod=preprocessMethod(gexp3),
                                   pData=pData(gexp3))
        return.list <- append(return.list,gcombat)
        message("Completed all specified steps! Returning..")
        return(return.list)
      }
    }

  } else{
    message("Completed all specified steps! Returning...")
    return(return.list)
  }
}

