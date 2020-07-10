load("rawData.Rdata")
source("rhythmic.functions.R")
library(rain)
library(reshape2)
library(data.table)
library(ggplot2)
library(scales)
library(matrixTests)
library(impute)

###############################################################
#impute metabolomics. No features have more than 20% NAs 

for(i in c("serum", "muscle")){
  x <- impute.knn(as.matrix(dataRaw[[i]][,-ncol(dataRaw[[i]]), with = F]))$data
  dataRaw[[i]] <- data.table(x, dataRaw[[i]][,ncol(dataRaw[[i]]), with = F])
}

##################################################################################
#Differential detection using rain

res <- NULL
for(i in c("muscle", "genes", "serum")){
  res[[paste0(i, ".U")]] <- rainAnalyze(data = dataRaw[[i]][,grepl(".U.", colnames(dataRaw[[i]]), fixed = T), with = F],
                                        meta = dataGroups[[i]][grepl("U", Food, fixed = T),],
                                        deltaPer = 12, period = 24)
  res[[paste0(i, ".R")]] <- rainAnalyze(data = dataRaw[[i]][,grepl(".R.", colnames(dataRaw[[i]]), fixed = T), with = F],
                                        meta = dataGroups[[i]][grepl("R", Food, fixed = T),],
                                        deltaPer = 12, period = 24)
  
  dataRaw[[paste0(i, ".R")]] <- dataRaw[[i]][,grepl("\\.R\\.|feature", colnames(dataRaw[[i]])), with = F]
  dataRaw[[paste0(i, ".U")]] <- dataRaw[[i]][,grepl("\\.U\\.|feature", colnames(dataRaw[[i]])), with = F]
  dataGroups[[paste0(i, ".R")]] <- dataGroups[[i]][grepl("R", Food, fixed = T),]
  dataGroups[[paste0(i, ".U")]] <- dataGroups[[i]][grepl("U", Food, fixed = T),]
  
}

temp <- names(res)
res <- lapply(names(res), function(x){
  
  temp <- res[[x]]
  temp$feature <- dataAnnotation[[gsub("\\.R|\\.U","",x)]]$feature
  temp
  
})
names(res) <- temp
rm(temp)

####################################################################################################################################################
#subset the data to match the periods as circadian genes have known periods (and its needed to run comparative stats)
#20-25% are not 24h.

coreCirc <- c("Clock", "ARNTL", "Npas2", "Per1", "Per2", "Per3", "Cry1", "Cry2", "Nr1d1", "Nr1d2", "Dbp", "RORA")
coreCirc <- dataAnnotation$genes[dataAnnotation$genes$Symbol %in% toupper(coreCirc), "feature"]

res$genes.U[feature %in% coreCirc, "period"] <- 24
res$genes.R[feature %in% coreCirc, "period"] <- 24

resSig <- lapply(res, function(x){
  x <- x[pValAdj<0.05,]
  x$Tags <- as.numeric(as.character(x$Tags))
  x[period %in% c(12,24),]#keep only 12 h and 24 hour as those are only intresting (other will phase out)
})

resSig$muscle <- data.table(period = 12, feature = intersect(resSig$muscle.U[period == 12,feature], resSig$muscle.R[period == 12,feature]))
resSig$muscle <- rbind(resSig$muscle, 
                       data.table(period = 24, feature = intersect(resSig$muscle.U[period == 24,feature], resSig$muscle.R[period == 24,feature])))

resSig$serum <- data.table(period = 12, feature = intersect(resSig$serum.U[period == 12,feature], resSig$serum.R[period == 12,feature]))
resSig$serum <- rbind(resSig$serum, 
                      data.table(period = 24, feature = intersect(resSig$serum.U[period == 24,feature], resSig$serum.R[period == 24,feature])))

resSig$genes <- data.table(period = 12, feature = intersect(resSig$genes.U[period == 12,feature], resSig$genes.R[period == 12,feature]))
resSig$genes <- rbind(resSig$genes, 
                      data.table(period = 24, feature = intersect(resSig$genes.U[period == 24,feature], resSig$genes.R[period == 24,feature])))

####################################################################################################################################################
#circa_compare for differences only in periodic
#sligh modification of circacompare function to output better column names
library(MASS)
library(tidyr)

resDiff <- NULL

for(tiss in c("muscle", "genes", "serum")){
  
  featAll <- resSig[[tiss]]
  
  data <- as.data.frame(dataRaw[[tiss]])
  data$feature <- dataRaw[[tiss]]$feature
  
  data <- melt(data, id.vars = "feature")
  temp <- strsplit2(data$variable, split = "\\.")
  data <- data.table(data, temp)
  data$V3 <- timeNum(data$V3)
  colnames(data)[4:6] <- c("participant", "Food", "Time")
  
  resDiff[[tiss]] <- NULL
  # dataCurves[[tiss]] <- NULL
  
  for(feat in featAll$feature){
    
    #transform data
    temp <- data[feature %in% feat]
    temp <- temp[temp$value != 0,]#necessary for boxcox
    
    #power transformation since not all features fullfill regression assumtions
    lambda <- geoR::boxcoxfit(temp$value)[[1]]
    lambda <- round(lambda, 2)
    
    if(lambda == 0|abs(lambda)>2){
      temp$value <- log10(temp$value)
    }else{
      if(lambda>0){
        temp$value <- temp$value^lambda
      }else{
        temp$value <- -temp$value^lambda
      }
    }
    
    per <- featAll[feature == feat, period]
    
    #circacompare calculation
    fit <- circacompare(x = temp, 
                        col_time = "Time", 
                        col_outcome = "value", 
                        col_group = "Food", 
                        period = per, 
                        alpha_threshold = 1)
    #threshold set to one in order to actually perform the analysis. It is later filtered.
    #rain has precalculated both the period but also the rythmicity.
    
    resDiff[[tiss]][[feat]] <- rbind(fit$results, data.frame(parameter = c("lamda", "period"), value  = c(lambda, per)))
    
    print(feat)
  }
  
  temp <- lapply(resDiff[[tiss]], function(x) x[,2])
  resDiff[[tiss]] <- temp[unlist(lapply(resDiff[[tiss]], nrow))>2]
  resDiff[[tiss]] <- t(matrix(unlist(resDiff[[tiss]]), nrow = 17, byrow=F))
  colnames(resDiff[[tiss]]) <- c("Both_rhythmic", 
                                 "Rhythmicity_pValR", 
                                 "Rhythmicity_pValU", 
                                 "R_MESOR", 
                                 "U_MESOR", 
                                 "MESOR_diff", 
                                 "MESOR_pVal", 
                                 "R_Amp", "U_Amp", 
                                 "Amp_diff", 
                                 "Amp_pVal", 
                                 "R_Acr", "U_Acr", 
                                 "Acr_diff", 
                                 "Acr_pVal", 
                                 "lamdba", 
                                 "period")
  resDiff[[tiss]] <- as.data.table(resDiff[[tiss]])
  resDiff[[tiss]]$features <- names(temp)[unlist(lapply(temp, length))>3]
  
  resDiff[[tiss]][,fdr_R:=p.adjust(Rhythmicity_pValR, "fdr")]
  resDiff[[tiss]][,fdr_U:=p.adjust(Rhythmicity_pValU, "fdr")]
  resDiff[[tiss]][,fdr_Acr:=p.adjust(Acr_pVal, "fdr")]
  resDiff[[tiss]][,fdr_MESOR:=p.adjust(MESOR_pVal, "fdr")]
  resDiff[[tiss]][,fdr_amp:=p.adjust(Amp_pVal, "fdr")]
  
}

####################################################################################################################################################
#calculate the single paramenters

resSingle <- NULL

for(tiss in names(dataRaw)[-1:-3]){
  
  featAll <- resSig[[tiss]]
  
  data <- as.data.frame(dataRaw[[tiss]])
  data$feature <- dataRaw[[tiss]]$feature
  
  data <- melt(data, id.vars = "feature")
  temp <- strsplit2(data$variable, split = "\\.")
  data <- data.table(data, temp)
  data$V3 <- timeNum(data$V3)
  colnames(data)[4:6] <- c("participant", "Food", "Time")
  
  resSingle[[tiss]] <- NULL
  
  for(feat in featAll$feature){
    
    #transform data
    temp <- data[feature %in% feat]
    temp <- temp[temp$value != 0,]#necessary for boxcox
    
    lambda <- geoR::boxcoxfit(temp$value)[[1]]
    lambda <- round(lambda, 2)
    
    if(lambda == 0|abs(lambda)>2){
      temp$value <- log10(temp$value)
    }else{
      if(lambda>0){
        temp$value <- temp$value^lambda
      }else{
        temp$value <- -temp$value^lambda
      }
    }
    
    per <- featAll[feature == feat, period]
    
    #circacompare calculation
    fit <- circa_single(x = temp, col_time = "Time", col_outcome = "value", period = per, alpha_threshold = 1)
    resSingle[[tiss]][[feat]] <- unlist(c(fit[[2]][1,], lambda = lambda, period = per))
    resSingle[[tiss]][[feat]][1:2] <- backTransformParameters(resSingle[[tiss]][[feat]]["mesor"], 
                                                              resSingle[[tiss]][[feat]]["amplitude"], 
                                                              lambda)
    print(feat)
  }
  
  
  temp <- as.data.table(matrix(unlist(resSingle[[tiss]]), 
                               ncol = length(resSingle[[tiss]][[1]]), 
                               byrow = T) )
  colnames(temp) <- names(resSingle[[tiss]][[1]])
  temp$feature <- names(resSingle[[tiss]])
  resSingle[[tiss]] <- temp
  
}

####################################################################################################################################################
#datacurves only for cortisol serum and core genes, and resDiff
dataCurvesFeatures <- c("Clock", 
                        "ARNTL", 
                        "Npas2", 
                        "Per1", 
                        "Per2", 
                        "Per3", 
                        "Cry1", 
                        "Cry2", 
                        "Nr1d1", 
                        "Nr1d2", 
                        "Dbp", 
                        "RORA", 
                        dataAnnotation$genes[dataAnnotation$genes$feature %in% resDiff$genes[fdr_Acr<0.05|fdr_amp<0.05|fdr_Acr<0.05]$features,"Symbol"])
dataCurvesFeatures <- toupper(dataCurvesFeatures)
#to get them in the same order
dataCurvesFeatures <- dataAnnotation$genes$Symbol[dataAnnotation$genes$Symbol %in% na.omit(dataCurvesFeatures)] 
names(dataCurvesFeatures) <- dataAnnotation$genes$genes[dataAnnotation$genes$Symbol %in% na.omit(dataCurvesFeatures)] 

dataCurves <- NULL
dataCurves$genes <- NULL
for(feat in names(dataCurvesFeatures)){
  
  data <- dataRaw$genes[feature == feat]
  
  data <- melt(data, id.vars = "feature")
  temp <- strsplit2(data$variable, split = "\\.")
  data <- data.table(data, temp)
  data$V3 <- timeNum(data$V3)
  colnames(data)[4:6] <- c("participant", "Food", "Time")
  
  #transform data
  temp <- data[feature %in% feat]
  temp <- temp[temp$value != 0,]#necessary for boxcox
  
  lambda <- geoR::boxcoxfit(temp$value)[[1]]
  lambda <- round(lambda, 2)
  
  if(lambda == 0|abs(lambda)>2){
    temp$value <- log10(temp$value)
  }else{
    if(lambda>0){
      temp$value <- temp$value^lambda
    }else{
      temp$value <- -temp$value^lambda
    }
  }
  fitU <- circa_single(x = temp[Food == "U"], col_time = "Time", col_outcome = "value", period = 24, alpha_threshold = 1)
  fitR <- circa_single(x = temp[Food == "R"], col_time = "Time", col_outcome = "value", period = 24, alpha_threshold = 1)
  
  time <- astroFns::hms2rad(h = seq(0,23.9,.1))
  fitU <- fitU[[3]]$m$getAllPars()["k"] + (fitU[[3]]$m$getAllPars()["alpha"] * cos(time - fitU[[3]]$m$getAllPars()["phi"]))
  fitR <- fitR[[3]]$m$getAllPars()["k"] + fitR[[3]]$m$getAllPars()["alpha"] * cos(time - fitR[[3]]$m$getAllPars()["phi"])
  
  dataCurves$genes[[feat]] <- data.table(value = c(fitU, fitR), 
                                         diet = c(rep("EXF", 240), rep("TRF", 240)),
                                         time = rep(seq(0,23.9,.1), 2))
  
  #backtransform the lamda value since reviewer 3 had concerns on the amplitude being non linear
  if(lambda == 0|abs(lambda)>2){
    dataCurves$genes[[feat]]$value <- 10^dataCurves$genes[[feat]]$value
  }else{
    dataCurves$genes[[feat]]$value <- sapply(dataCurves$genes[[feat]]$value, function(x) backTransformSingle(x, round(lambda, 2)))
  }
  
  print(feat)
  
}

dataCurvesFeatures <- c("cortisol", "cortisone", "corticosterone",
                        resDiff$serum[fdr_Acr<0.05|fdr_amp<0.05|fdr_Acr<0.05]$features)

dataCurves$serum <- NULL
for(feat in dataCurvesFeatures){
  
  data <- dataRaw$serum[feature == feat]
  
  data <- melt(data, id.vars = "feature")
  temp <- strsplit2(data$variable, split = "\\.")
  data <- data.table(data, temp)
  data$V3 <- timeNum(data$V3)
  colnames(data)[4:6] <- c("participant", "Food", "Time")
  
  #transform data
  temp <- data[feature %in% feat]
  temp <- temp[temp$value != 0,]#necessary for boxcox
  
  lambda <- geoR::boxcoxfit(temp$value)[[1]]
  lambda <- round(lambda, 2)
  
  if(lambda == 0|abs(lambda)>2){
    temp$value <- log10(temp$value)
  }else{
    if(lambda>0){
      temp$value <- temp$value^lambda
    }else{
      temp$value <- -temp$value^lambda
    }
  }
  fitU <- circa_single(x = temp[Food == "U"], col_time = "Time", col_outcome = "value", period = 24, alpha_threshold = 1)
  fitR <- circa_single(x = temp[Food == "R"], col_time = "Time", col_outcome = "value", period = 24, alpha_threshold = 1)
  
  time <- astroFns::hms2rad(h = seq(0,23.9,.1))
  fitU <- fitU[[3]]$m$getAllPars()["k"] + (fitU[[3]]$m$getAllPars()["alpha"] * cos(time - fitU[[3]]$m$getAllPars()["phi"]))
  fitR <- fitR[[3]]$m$getAllPars()["k"] + fitR[[3]]$m$getAllPars()["alpha"] * cos(time - fitR[[3]]$m$getAllPars()["phi"])
  
  dataCurves$serum[[feat]] <- data.table(value = c(fitU, fitR), 
                                         diet = c(rep("EXF", 240), rep("TRF", 240)),
                                         time = rep(seq(0,23.9,.1), 2))
  
  if(lambda == 0|abs(lambda)>2){
    dataCurves$serum[[feat]]$value <- 10^dataCurves$serum[[feat]]$value
  }else{
    dataCurves$serum[[feat]]$value <- sapply(dataCurves$serum[[feat]]$value, function(x) backTransformSingle(x, round(lambda, 2)))
  }
  
  print(feat)
  
}

dataCurvesFeatures <- resDiff$muscle[fdr_Acr<0.05|fdr_amp<0.05|fdr_Acr<0.05]$features

dataCurves$muscle <- NULL
for(feat in dataCurvesFeatures){
  
  data <- dataRaw$serum[feature == feat]
  
  data <- melt(data, id.vars = "feature")
  temp <- strsplit2(data$variable, split = "\\.")
  data <- data.table(data, temp)
  data$V3 <- timeNum(data$V3)
  colnames(data)[4:6] <- c("participant", "Food", "Time")
  
  #transform data
  temp <- data[feature %in% feat]
  temp <- temp[temp$value != 0,]#necessary for boxcox
  
  lambda <- geoR::boxcoxfit(temp$value)[[1]]
  lambda <- round(lambda, 2)
  
  if(lambda == 0|abs(lambda)>2){
    temp$value <- log10(temp$value)
  }else{
    if(lambda>0){
      temp$value <- temp$value^lambda
    }else{
      temp$value <- -temp$value^lambda
    }
  }
  fitU <- circa_single(x = temp[Food == "U"], col_time = "Time", col_outcome = "value", period = 24, alpha_threshold = 1)
  fitR <- circa_single(x = temp[Food == "R"], col_time = "Time", col_outcome = "value", period = 24, alpha_threshold = 1)
  
  time <- astroFns::hms2rad(h = seq(0,23.9,.1))
  fitU <- fitU[[3]]$m$getAllPars()["k"] + (fitU[[3]]$m$getAllPars()["alpha"] * cos(time - fitU[[3]]$m$getAllPars()["phi"]))
  fitR <- fitR[[3]]$m$getAllPars()["k"] + fitR[[3]]$m$getAllPars()["alpha"] * cos(time - fitR[[3]]$m$getAllPars()["phi"])
  
  dataCurves$muscle[[feat]] <- data.table(value = c(fitU, fitR), 
                                         diet = c(rep("EXF", 240), rep("TRF", 240)),
                                         time = rep(seq(0,23.9,.1), 2))
  
  if(lambda == 0|abs(lambda)>2){
    dataCurves$muscle[[feat]]$value <- 10^dataCurves$muscle[[feat]]$value
  }else{
    dataCurves$muscle[[feat]]$value <- sapply(dataCurves$muscle[[feat]]$value, function(x) backTransformSingle(x, round(lambda, 2)))
  }
  
  print(feat)
  
}

dataCurves <- lapply(dataCurves, function(dataCurvestemp){
  temp <- lapply(dataCurvestemp, function(x) x[,1])
  matrix(unlist(temp), nrow = 480, byrow=F)[,2] == temp[[2]]
  temp <- as.data.table(matrix(unlist(temp), nrow = 480, byrow=F))
  colnames(temp) <- names(dataCurvestemp)
  temp$group <- paste0(dataCurvestemp[[1]][,time],"_", dataCurvestemp[[1]][,diet])
  temp
})

#########################################################################################################################################################
#integrate the different results

resDiff
if(!all(all(resSingle$muscle.U$feature == resSig$muscle.U$feature),
        all(resSingle$muscle.R$feature == resSig$muscle.R$feature),
        all(resSingle$genes.U$feature == resSig$genes.U$feature),
        all(resSingle$genes.R$feature == resSig$genes.R$feature),
        all(resSingle$serum.U$feature == resSig$serum.U$feature),
        all(resSingle$serum.R$feature == resSig$serum.R$feature))){
  stop("stop!")
  }

resSig <- lapply(names(resSig)[-7:-9], function(x){
  x <- cbind(resSig[[x]], resSingle[[x]][,c("mesor", "amplitude", "peak_time_hours", "lambda")])
  colnames(x)[10] <- "Acrophase"
  x$Tags <- as.numeric(x$Tags)
  x
})
names(resSig) <- names(res)

resSigHalf <- lapply(resSig, function(x){
  x[period == 12]
})
resSig <- lapply(resSig, function(x){
  x[period == 24]
})

save(file = "rhythmic.metabolites.Rdata", 
     list = c("dataAnnotation",
              "dataCurves",
              "dataRaw",
              "dataGroups",
              "resSig",
              "resSigHalf",
              "res",
              "resDiff"))
