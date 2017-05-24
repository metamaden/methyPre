methyPre <- function(meExpObj, workflow=c("norm","pfilt","map","minfifilt","ufilt","crxcgfilt","batch"),
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
