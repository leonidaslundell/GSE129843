library(reshape2)
library(ggplot2)
library(data.table)
library(clusterProfiler)
library(patchwork)
library(scales)
library(UpSetR)
library(VennDiagram)
library(impute)
library(Rtsne)
library(waffle)
library(ggupset)
library(openxlsx)

###############################################################
#functions to load

#very lightly modified to output more managable variable names...
circacompare <- function(x,
                         col_time,
                         col_group,
                         col_outcome,
                         period = 24,
                         alpha_threshold = 0.05,
                         timeout_n = 10000){
  
  if(!"ggplot2" %in% installed.packages()[, "Package"]){
    return(message("Please install 'ggplot2'"))
  }
  
  library(ggplot2)
  colnames(x)[grep(col_group, colnames(x))] <- "group"
  
  if(length(levels(as.factor(x$group))) != 2){
    return(message("Your grouping variable had more or less than 2 levels! \nThis function is used to compare two groups of data. \nTo avoid me having to guess, please send data with only two possible values in your grouping variable to this function."))
  }
  
  group_1_text <- levels(as.factor(x$group))[1]
  group_2_text <- levels(as.factor(x$group))[2]
  colnames(x)[grep(col_time, colnames(x))] <- "time"
  
  if(!class(x$time) %in% c("numeric", "integer")){
    return(message(paste("The time variable which you gave was a '",
                         class(x$time),
                         "' \nThis function expects time to be given as hours and be of class 'integer' or 'numeric'.",
                         "\nPlease convert the time variable in your dataframe to be of one of these classes",
                         sep = "")))
  }
  
  colnames(x)[grep(col_outcome, colnames(x))] <- "measure"
  
  if(!class(x$measure) %in% c("numeric", "integer")){
    return(message(paste("The measure variable which you gave was a '",
                         class(x$measure),
                         "' \nThis function expects measure to be number and be of class 'integer' or 'numeric'.",
                         "\nPlease convert the measure variable in your dataframe to be of one of these classes",
                         sep = "")))
  }
  
  x$time_r <- (x$time/24)*2*pi*(24/period)
  x$x_group <- ifelse(x$group == group_1_text, 0, 1)
  
  comparison_model_success <- 0
  comparison_model_timeout <- FALSE
  g1_success <- 0
  g2_success <- 0
  g1_alpha_p <- NA
  g2_alpha_p <- NA
  n <- 0
  dat_group_1 <- x[x$group == group_1_text,]
  dat_group_2 <- x[x$group == group_2_text,]
  
  while(g1_success !=1){
    g1_alpha_start <- (max(dat_group_1$measure, na.rm = TRUE) - min(dat_group_1$measure, na.rm = TRUE)) * runif(1)
    g1_phi_start <- runif(1)*6.15 - 3.15
    g1_k_start <- mean(dat_group_1$measure, na.rm = TRUE)*2*runif(1)
    
    fit.nls_group_1 <- try({nls(measure~k + alpha*cos(time_r-phi),
                                data = dat_group_1,
                                start = list(k=g1_k_start,alpha=g1_alpha_start,phi=g1_phi_start))},
                           silent = TRUE)
    if(class(fit.nls_group_1) == "try-error"){
      n <- n + 1
    }
    else{
      g1_k_out <- summary(fit.nls_group_1)$coef[1,1]
      g1_alpha_out <- summary(fit.nls_group_1)$coef[2,1]
      g1_alpha_p <- summary(fit.nls_group_1)$coef[2,4]
      g1_phi_out <- summary(fit.nls_group_1)$coef[3,1]
      g1_success <- ifelse(g1_alpha_out > 0,1,0)
      n <- n + 1
    }
    if(n >= timeout_n){
      return(message("Failed to converge group 1 data prior to timeout. \nYou may try to increase the allowed attempts before timeout by increasing the value of the 'timeout_n' argument or setting a new seed before this function.\nIf you have repeated difficulties, please contact me (via github) or Oliver Rawashdeh (contact details in manuscript)."))
    }
  }
  n <- 0
  while(g2_success !=1){
    g2_alpha_start <- (max(dat_group_2$measure, na.rm = TRUE) - min(dat_group_2$measure, na.rm = TRUE)) * runif(1)
    g2_phi_start <- runif(1)*6.15 - 3.15
    g2_k_start <- mean(dat_group_2$measure, na.rm = TRUE)*2*runif(1)
    
    fit.nls_group_2 <- try({nls(measure~k + alpha*cos(time_r-phi),
                                data = dat_group_2,
                                start = list(k=g2_k_start,alpha=g2_alpha_start,phi=g2_phi_start))},
                           silent = TRUE)
    if(class(fit.nls_group_2) == "try-error"){
      n <- n + 1
    }
    else{
      g2_k_out <- summary(fit.nls_group_2)$coef[1,1]
      g2_alpha_out <- summary(fit.nls_group_2)$coef[2,1]
      g2_alpha_p <- summary(fit.nls_group_2)$coef[2,4]
      g2_phi_out <- summary(fit.nls_group_2)$coef[3,1]
      g2_success <- ifelse(g2_alpha_out > 0,1,0)
      n <- n + 1
    }
    if(n >= timeout_n){
      return(message("Failed to converge group 2 data prior to timeout. \nYou may try to increase the allowed attempts before timeout by increasing the value of the 'timeout_n' argument or setting a new seed before this function.\nIf you have repeated difficulties, please contact me (via github) or Oliver Rawashdeh (contact details in manuscript)."))
    }
  }
  
  g1_rhythmic <- ifelse(g1_alpha_p < alpha_threshold, TRUE, FALSE)
  g2_rhythmic <- ifelse(g2_alpha_p < alpha_threshold, TRUE, FALSE)
  both_groups_rhythmic <- ifelse(g1_rhythmic ==TRUE & g2_rhythmic==TRUE, TRUE, FALSE)
  
  if(both_groups_rhythmic==FALSE){
    if(g1_rhythmic == FALSE & g2_rhythmic == FALSE){
      return(message("Both groups of data were arrhythmic (to the power specified by the argument 'alpha_threshold').\nThe data was, therefore, not used for a comparison between the two groups."))
    }
    if(g1_rhythmic == FALSE){
      return(message(group_1_text, " was arrhythmic (to the power specified by the argument 'alpha_threshold').\nThe data was, therefore, not used for a comparison between the two groups."))
    }else{
      return(message(group_2_text, " was arrhythmic (to the power specified by the argument 'alpha_threshold').\nThe data was, therefore, not used for a comparison between the two groups."))
    }
  }
  n <- 0
  if(both_groups_rhythmic == TRUE){
    while(comparison_model_success == 0 & comparison_model_timeout == FALSE){
      alpha_in <- g1_alpha_out*2*runif(1)
      alpha1_in <- (g2_alpha_out - g1_alpha_out)*2*runif(1)
      phi_in <- g1_phi_out*2*runif(1)
      phi1_in <- runif(1)*2*pi - pi
      k_in <- g1_k_out*2*runif(1)
      k1_in <- (g2_k_out - g1_k_out)*2*runif(1)
      
      fit.nls <- try({nls(measure~k+k1*x_group+(alpha+alpha1*x_group)*cos(time_r-(phi+phi1*x_group)),
                          data = x,
                          start = list(k=k_in, k1=k1_in, alpha=alpha_in, alpha1=alpha1_in, phi=phi_in, phi1=phi1_in),
                          nls.control(maxiter = 100, minFactor = 1/10000#, warnOnly = TRUE
                          ))},
                     silent = TRUE)
      
      if (class(fit.nls) == "try-error") {
        n <- n + 1
      }
      else{
        k_out <- coef(fit.nls)[1]
        k1_out <- coef(fit.nls)[2]
        k_out_p <- (summary(fit.nls)$coef)[1,4]
        k1_out_p <- (summary(fit.nls)$coef)[2,4]
        
        alpha_out <- coef(fit.nls)[3]
        alpha1_out <- coef(fit.nls)[4]
        alpha1_out_p <- (summary(fit.nls)$coef)[4,4]
        phi_out <- coef(fit.nls)[5]
        phi1_out <- coef(fit.nls)[6]
        phi1_out_p <- (summary(fit.nls)$coef)[6,4]
        
        comparison_model_success <- ifelse(alpha_out>0 & (alpha_out + alpha1_out) > 0 & phi1_out <pi & phi1_out >-pi, 1, 0)
        comparison_model_timeout <- ifelse(n>timeout_n, TRUE, FALSE)
        n <- n + 1
      }
    }
    
    if(comparison_model_timeout == TRUE){
      return(message("Both groups of data were rhythmic but the curve fitting procedure failed due to timing out. \nYou may try to increase the allowed attempts before timeout by increasing the value of the 'timeout_n' argument or setting a new seed before this function.\nIf you have repeated difficulties, please contact me (via github) or Oliver Rawashdeh (contact details in manuscript)."))
    }
    #loop curve fitting process (all data) until outputs are appropriate, or until looped more times than timeout_n
    if(comparison_model_timeout == FALSE){
      eq_1 <- function(time){k_out + alpha_out*cos((2*pi/period)*time - phi_out)}
      eq_2 <- function(time){k_out + k1_out + (alpha_out + alpha1_out)*cos((2*pi/period)*time - (phi_out + phi1_out))}
      
      fig_out <- ggplot2::ggplot(x, aes(time, measure)) +
        stat_function(fun = eq_1, colour = "blue", size=1) +
        stat_function(fun = eq_2, colour = "red", size=1) +
        geom_point(aes(colour = group)) +
        scale_colour_manual(breaks = c(group_1_text, group_2_text),
                            values = c("blue", "red")) +
        xlab("time (hours)") +
        xlim(min(floor(x$time/period) * period),
             max(ceiling(x$time/period) * period))
      
    }#if the nls was successful, create a graph to plot the data as well as curves of best fit, 'fig_out'
  }
  if(both_groups_rhythmic==TRUE & comparison_model_success==1){
    if(phi_out > pi){
      while(phi_out > pi){
        phi_out <- phi_out - 2*pi
      }
    }
    if(phi_out < -pi){
      while(phi_out < -pi){
        phi_out <- phi_out + 2*pi
      }
    }#adjust phi_out so that -pi < phi_out < pi
    baseline_diff_abs <- k1_out
    baseline_diff_pc <- ((k_out + k1_out)/k_out)*100 - 100
    amplitude_diff_abs <- alpha1_out
    amplitude_diff_pc <-  ((alpha_out+alpha1_out)/alpha_out)*100 - 100
    g1_peak_time <- phi_out*period/(2*pi)
    g2_peak_time <- (phi_out+phi1_out)*period/(2*pi)
    
    while(g1_peak_time >period | g1_peak_time < 0){
      if(g1_peak_time >period){
        g1_peak_time <- g1_peak_time - period
      }
      if(g1_peak_time<0){
        g1_peak_time <- g1_peak_time + period
      }
    }
    while(g2_peak_time >period| g2_peak_time <0){
      if(g2_peak_time>period){
        g2_peak_time <- g2_peak_time - period
      }
      if(g2_peak_time<0){
        g2_peak_time <- g2_peak_time + period
      }
    }
    peak_time_diff <- phi1_out*period/(2*pi)
  }
  
  output_parms <- data.frame(parameter = c("Both_rhythmic",
                                           paste("Rhythmicity_pVal", group_1_text, sep = ""),
                                           paste("Rhythmicity_pVal", group_2_text, sep = ""),
                                           paste(group_1_text, "_MESOR", sep = ""),
                                           paste(group_2_text, "_MESOR", sep = ""),
                                           "MESOR_diff",
                                           "MESOR_pVal",
                                           paste(group_1_text, "_Amp", sep = ""),
                                           paste(group_2_text, "_Amp", sep = ""),
                                           "Amp_diff",
                                           "Amp_pVal",
                                           paste(group_1_text, "_Acr", sep = ""),
                                           paste(group_2_text, "_Acr", sep = ""),
                                           "Acr_diff",
                                           "Acr_pVal"),
                             value = c(both_groups_rhythmic, 
                                       g1_alpha_p, 
                                       g2_alpha_p, 
                                       k_out, 
                                       (k_out + k1_out), 
                                       k1_out,
                                       k1_out_p, 
                                       alpha_out, 
                                       alpha_out + alpha1_out, 
                                       alpha1_out, 
                                       alpha1_out_p,
                                       g1_peak_time, 
                                       g2_peak_time, 
                                       peak_time_diff, 
                                       phi1_out_p))
  
  
  if(exists("fig_out")){
    return(list(fig = list(group1 = eq_1, group2 = eq_2), results = output_parms, fitBoth = fit.nls, fitGroup1 = fit.nls_group_1, fitGroup2 = fit.nls_group_2))
  }
}

circa_single <- function (x, col_time, col_outcome, period = 24, alpha_threshold = 0.05, 
                          timeout_n = 10000) 
{
  if (!"ggplot2" %in% installed.packages()[, "Package"]) {
    return(message("Please install 'ggplot2'"))
  }
  library(ggplot2)
  colnames(x)[grep(col_time, colnames(x))] <- "time"
  if (!class(x$time) %in% c("numeric", "integer")) {
    return(message(paste("The time variable which you gave was a '", 
                         class(x$time), "' \nThis function expects time to be given as hours and be of class 'integer' or 'numeric'.", 
                         "\nPlease convert the time variable in your dataframe to be of one of these classes", 
                         sep = "")))
  }
  colnames(x)[grep(col_outcome, colnames(x))] <- "measure"
  if (!class(x$measure) %in% c("numeric", "integer")) {
    return(message(paste("The measure variable which you gave was a '", 
                         class(x$measure), "' \nThis function expects measure to be number and be of class 'integer' or 'numeric'.", 
                         "\nPlease convert the measure variable in your dataframe to be of one of these classes", 
                         sep = "")))
  }
  x$time_r <- (x$time/24) * 2 * pi * (24/period)
  comparison_model_success <- 0
  comparison_model_timeout <- FALSE
  success <- 0
  n <- 0
  while (success != 1) {
    alpha_start <- (max(x$measure, na.rm = TRUE) - min(x$measure, 
                                                       na.rm = TRUE)) * runif(1)
    phi_start <- runif(1) * 6.15 - 3.15
    k_start <- mean(x$measure, na.rm = TRUE) * 2 * runif(1)
    fit.nls <- try({
      nls(measure ~ k + alpha * cos(time_r - phi), data = x, 
          start = list(k = k_start, alpha = alpha_start, 
                       phi = phi_start))
    }, silent = TRUE)
    if (class(fit.nls) == "try-error") {
      n <- n + 1
    }
    else {
      k_out <- summary(fit.nls)$coef[1, 1]
      alpha_out <- summary(fit.nls)$coef[2, 1]
      alpha_p <- summary(fit.nls)$coef[2, 4]
      phi_out <- summary(fit.nls)$coef[3, 1]
      success <- ifelse(alpha_out > 0 & phi_out >= 0 & 
                          phi_out <= 2 * pi, 1, 0)
      n <- n + 1
    }
    if (n >= timeout_n) {
      return(message("Failed to converge data prior to timeout. \nYou may try to increase the allowed attempts before timeout by increasing the value of the 'timeout_n' argument or setting a new seed before this function.\nIf you have repeated difficulties, please contact me (via github) or Oliver Rawashdeh (contact details in manuscript)."))
    }
  }
  data_rhythmic <- ifelse(alpha_p < alpha_threshold, TRUE, 
                          FALSE)
  eq <- function(time) {
    k_out + alpha_out * cos((2 * pi/period) * time - phi_out)
  }
  if (data_rhythmic == TRUE) {
    fig_out <- ggplot2::ggplot(x, aes(time, measure)) + stat_function(fun = eq, 
                                                                      size = 1) + geom_point() + xlab("time (hours)") + 
      xlim(min(floor(x$time/period) * period), max(ceiling(x$time/period) * 
                                                     period)) + labs(subtitle = "Data is rhythmic")
  }
  else {
    fig_out <- ggplot2::ggplot(x, aes(time, measure)) + geom_point() + 
      xlab("time (hours)") + xlim(min(floor(x$time/period) * 
                                        period), max(ceiling(x$time/period) * period)) + 
      labs(subtitle = "Data is arrhythmic")
  }
  k_out
  alpha_out
  alpha_p
  phi_out
  peak_time <- phi_out * period/(2 * pi)
  output_parms <- data.frame(mesor = k_out, amplitude = alpha_out, 
                             amplitude_p = alpha_p, phase_radians = phi_out, peak_time_hours = phi_out * 
                               period/(2 * pi))
  return(list(fig_out, output_parms, fit.nls))
}

###############
#backtransform power transformed data

backTransformSingle <- function(x, lambdaInt){
  
  if (lambdaInt > 0) {
    (x)^(1/lambdaInt)
  } else {
    (-x)^(1/lambdaInt)
  } 
  
}

backTransformParameters <- function(mesor, amplitude, lambda){
  
  backTransformInternal <- function(x, lambdaInt){
    if(lambdaInt == 0||abs(lambdaInt)>2){
      10^x
    }else{
      if (lambdaInt > 0) {
        (x)^(1/lambdaInt)
      } else {
        (-x)^(1/lambdaInt)
      } 
    }
  }
  
  #recalculate the amplitude in original value scale.
  pos <- backTransformInternal((mesor + amplitude), lambda) -
    backTransformInternal(mesor, lambda = lambda)
  neg <- backTransformInternal(mesor, lambda) -
    backTransformInternal((mesor - amplitude), lambda = lambda)
  
  amp <- (neg + pos)/2
  
  mesor <- backTransformInternal(mesor, lambda)
  
  return(c(mesor = mesor, amplitude = amp))
  
}

###############
#nearest time of acrophase
nearestTime <- function(times, times.to.match= c(0,4,8,12,16,20), hoursPerCycle = 24){
  #if over 24h substract 24
  if(times>24){times <- times - 24}
  
  sinMeasuredTimes <- sin(times.to.match/hoursPerCycle * 2 * pi)
  cosMeasuredTimes <- cos(times.to.match/hoursPerCycle * 2 * pi)
  
  sinNewTimepoint <- sin(times/hoursPerCycle * 2 * pi)
  cosNewTimepoint <- cos(times/hoursPerCycle * 2 * pi)
  
  distances <- sqrt(
    outer(sinMeasuredTimes, sinNewTimepoint, FUN = `-`)^2 +
      outer(cosMeasuredTimes, cosNewTimepoint, FUN = `-`)^2
  )
  minVals <- apply(distances, MARGIN = 2, min)
  closestTimepoint <- t(distances) - minVals < sqrt(.Machine$double.eps)
  
  out <- apply(closestTimepoint, 1, function(x){times.to.match[which(x)]})
  
  if(length(out)>1){stop("time ambiguity")}
  
  return(out)
}

##########################
#phase advance delay

phaseShiftAnalyzer <- function(resU, resR){
  inter <- intersect(resU$feature, resR$feature)
  resU <- resU[match(inter, feature)]
  resR <- resR[match(inter, feature)]
  
  if(all(resR$feature == resU$feature)){
    phase <- resR$Acrophase - resU$Acrophase
    phase <- sapply(phase, nearestTime)
    inter <- data.frame(feature = inter, phase = phase)
    # inter <- inter[!inter$phase == 0,]
    inter <- list(advance = inter[inter$phase<0,],
                  delay = inter[inter$phase>0,],
                  all = inter)
    #cant have more than 12h distance. If above 12, substract 24
    inter$all[abs(inter$all$phase)>12, "phase"] <- inter$all[abs(inter$all$phase)>12, "phase"] - 24
    return(inter)
    
  }else{warning("features not same order")}
}

plotPhaseShift <- function(data, ymax){
  cols <- rep("gray40", length(table(data$phase)))
  cols[names(table(data$phase)) == 0] <- "gray"
  
  ggplot(data, aes(x = phase)) + 
    geom_bar(stat = "count", width = 3, fill = cols) +
    theme(panel.background = element_blank(),
          axis.line = element_line(color = "black"),
          axis.ticks = element_line(size = 1),
          axis.text = element_text(size =12, color = "black")) +
    scale_y_continuous(expand = c(0,0), limits = c(0, max(table(data$phase)) +1)) +
    scale_x_continuous(labels = names(table(data$phase)), 
                       breaks = as.numeric(names(table(data$phase)))) +
    xlab("Phase shift (h compared to EXF)")
}

##########################
#chronogram plot

plotChrono <- function(data, limits, labelsExclude, hjust){
  
  plotChronoInternalBase <- function(data){
    ggplot(data, aes(x=time, y = value, fill=diet)) +
      geom_bar(stat="identity", position = "dodge", width = .7) +
      scale_fill_manual(values=c("black","red1")) +
      
      theme_minimal() +
      theme(axis.title.y =  element_text(color = "black", size = 16, hjust = 0.7),
            axis.title.x =  element_blank(),
            axis.text.x = element_text(size = 22, color = "black", vjust = 1),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            axis.line.y.left = element_blank(),
            legend.text = element_text(size = 20),
            legend.title =  element_blank(),
            legend.position = "none",
            legend.key.size = unit(3,"line"),
            plot.margin = unit(c(3,3,3,3), "cm")) +
      coord_polar(start = -.04, clip = "off")
  }
  
  plotChronoInternalAxis <- function(gg, limits, hjust = 1.15, labelsExclude){
    #manual axis
    labelsExclude <- !paste0(seq(10, 100, 10), "%") %in% labelsExclude
    gg +
      #axis line
      scale_y_continuous(name = NULL, limits = limits, breaks = seq(0,.8,.1), position = "right") +
      #axis labels
      #modify so that remove labelsExclude
      annotate(geom = "text", 
               x = rep(.5, 10)[labelsExclude], 
               y = seq(0.1, 1, .1)[labelsExclude], 
               label = paste0(seq(10, 100, 10), "%")[labelsExclude],
               hjust = hjust, size = 7) +
      #main axis
      geom_line(data = data.frame(x = c(0.5, 0.5), y = c(0, limits[2])), aes(x = x, y = y), inherit.aes = F) +
      # axis ticks
      geom_line(data = data.frame(x = c(0.45, 0.55), y = c(0.05, 0.05)), aes(x = x, y = y), inherit.aes = F, size =0.5) +
      geom_line(data = data.frame(x = c(0.45, 0.55), y = c(0.1, 0.1)), aes(x = x, y = y), inherit.aes = F, size =0.5) +
      geom_line(data = data.frame(x = c(0.45, 0.55), y = c(0.2, 0.2)), aes(x = x, y = y), inherit.aes = F, size =0.5) +
      geom_line(data = data.frame(x = c(0.45, 0.55), y = c(0.3, 0.3)), aes(x = x, y = y), inherit.aes = F) +
      geom_line(data = data.frame(x = c(0.45, 0.55), y = c(0.4, 0.4)), aes(x = x, y = y), inherit.aes = F) +
      geom_line(data = data.frame(x = c(0.45, 0.55), y = c(0.5, 0.5)), aes(x = x, y = y), inherit.aes = F) +
      geom_line(data = data.frame(x = c(0.45, 0.55), y = c(0.6, 0.6)), aes(x = x, y = y), inherit.aes = F) +
      geom_line(data = data.frame(x = c(0.45, 0.55), y = c(0.7, 0.7)), aes(x = x, y = y), inherit.aes = F) +
      geom_line(data = data.frame(x = c(0.45, 0.55), y = c(0.8, 0.8)), aes(x = x, y = y), inherit.aes = F) +
      geom_line(data = data.frame(x = c(0.45, 0.55), y = c(0.9, 0.9)), aes(x = x, y = y), inherit.aes = F) 
  }
  
  plotChronoInternalFood <- function(gg, limits){
    limits <- limits[2]*1.0
    gg + 
      annotate("text", x = 1, limits*.95, label = sprintf('\u25bc'), color="black", size=10 , angle=-30, fontface="bold") + #U 7:00
      annotate("text", x = 2.75, limits, label = sprintf('\u25bc'), color="black", size=10 , angle=-135, fontface="bold") + #U 14:00
      annotate("text", x = 4.5, limits, label = sprintf('\u25bc'), color="black", size=10 , angle=120, fontface="bold") +  # U21
      annotate("text", x = 1.75, limits, label = sprintf('\u25bc'), color="red", size=10 , angle=-75, fontface="bold") + #R 10
      annotate("text", x = 2.5, limits, label = sprintf('\u25bc'), color="red", size=10 , angle=-120, fontface="bold") + #R 13
      annotate("text", x = 3.5, limits, label = sprintf('\u25bc'), color="red", size=10 , angle=-180, fontface="bold")
  }
  
  gg <- plotChronoInternalBase(data)
  gg <- plotChronoInternalAxis(gg, limits, labelsExclude = labelsExclude)
  gg <- plotChronoInternalFood(gg, limits)
  return(gg)
}

##########################
#Upset plot

dataUpset <- function(data.upset){
  
  data.upset <- as.data.table(melt(data.upset))
  data.upset$value <- as.character(data.upset$value)
  
  test <- lapply(unique(data.upset$value), function(x){
    data.upset[value==x, L1]
  })
  
  # names(test) <- unique(data.upset$value)
  data.upset <- data.frame(features = unique(data.upset$value))
  data.upset$value <- test
  return(data.upset)
}

###############
#strsplit2
strsplit2 <-function(x, split, fixed = F)
{
  x <- strsplit(as.character(x), split, fixed = fixed)
  l <- length(x[[1]])
  
  matrix(unlist(x), ncol = l, byrow = T)
  
}

###############
#Plot raw values and mean.

plotCircadian <- function(resU, resR, raw, curves, features, log = F){
  
  #raw
  plot <- melt(raw[feature %in% features,], id.vars = "feature")
  plot <- cbind(plot, strsplit2(plot$variable, "\\."))[,-2]
  colnames(plot) <- c("feature", "value", "Part", "Diet", "Time")
  plot$Time <- timeNum(plot$Time)
  
  #cosinor
  plotCurves <- melt(curves, id.vars = "group")
  plotCurves <- cbind(plotCurves, strsplit2(plotCurves$group, "\\_"))[,-1]
  colnames(plotCurves) <- c("feature", "value", "Time","Diet")
  plotCurves[,id:=paste0(Diet,".",feature)]
  plotCurves$Time <- as.numeric(plotCurves$Time)
  
  #mesor amp acr
  plotMeta <- rbind(cbind(Diet = "EXF", resU[feature %in% features, c("feature", "pValAdj")]),
                    cbind(Diet = "TRF", resR[feature %in% features, c("feature", "pValAdj")]))
  plotMeta[,id:=paste0(Diet,".",feature)]
  
  plotCurves <- merge(plotMeta[,c("id", "pValAdj")], plotCurves, by = "id")
  plotCurves$Star <- symnum(plotCurves$pValAdj, cutpoints = c(0, 0.005, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns"))
  plotCurves$lty <- symnum(plotCurves$pValAdj, cutpoints = c(0, 0.05, 1), symbols = c("sig", "ns"))
  
  plotCurves$Star <- factor(plotCurves$Star, levels = c("***", "**", "*", "ns")) 
  plotCurves$Star <- relevel(plotCurves$Star, ref = "ns")
  plotCurves$lty <- factor(plotCurves$lty, levels = c("sig", "ns"))
  plotCurves$lty <- relevel(plotCurves$lty, ref = "ns")
  
  ###############
  #emtpy factor for getting all aesth.
  dataEmt <- plotCurves[1:4]
  dataEmt[,"pValAdj"] <- c(0.05,0.01,0.001,1)
  dataEmt[,"Star"] <- c("*","**","***","ns")
  dataEmt[,"lty"] <- c("sig","sig", "sig", "ns")
  # dataEmt[,"feature"] <- as.factor(letters[1:4])
  
  plot <- list(points = plot, curves = plotCurves[,-1], empty = dataEmt)
  plot <- lapply(plot, function(x){
    x$Diet <- as.factor(x$Diet)
    levels(x$Diet) <- c("TRF", "EXF")
    x$Diet <- relevel(x$Diet, ref = "EXF")
    x
  })
  
  out <- ggplot(plot$points, aes(x = Time, y = value, color = Diet)) + 
    geom_jitter(width = .5) + 
    geom_line(data = plot$curves, aes(size = Star, lty = lty), inherit.aes = T) +
    geom_blank(data = plot$empty, aes(size = Star, lty = lty)) +
    facet_wrap(~feature, scale = "free_y", drop = T) + 
    theme(legend.position = "bottom",
          plot.title = element_text(size = 18), 
          axis.text.y = element_text(size = 16, colour = "black"),
          axis.title.y = element_text(size = 16),
          strip.text = element_text(size=25),
          legend.text = element_text(size=15), 
          legend.title = element_text(size=20),
          axis.text.x = element_text(size = 14, color = "black", angle = -45, hjust = -.17)) +
    scale_x_continuous(labels = c("07:00", "11:00", "15:00", "19:00", "23:00", "03:00"), 
                       name = "Time of the day",
                       breaks =  sort(unique(plot$points$Time))) +
    scale_color_manual(values=c("black", "red1"), name = "Diet") +
    #add food triangles
    annotate("text", x = c(1, 8, 15), y = -Inf, label = sprintf('\u25bc'),
             color="black", size = 5, angle =180, alpha = .8) +
    annotate("text", x = c(4, 7, 11), y = -Inf, label = sprintf('\u25bc'),
             color="red",  size = 5, angle =180, alpha = .9) +
    #add meta for pvals and stuff
    scale_linetype_manual(labels = c("sig","ns"),
                          values = c("dashed","solid"),
                          breaks = c("","ns"),
                          name = "Significance") +
    scale_size_manual(labels = c("*", "**", "***", ""),
                      breaks = c("*", "**", "***", ""),
                      values = c(.5, 1, 2, 0.5),
                      name = NULL) +
    guides(size = guide_legend(order = 3),
           color = guide_legend(order = 1),
           linetype = guide_legend(order = 2)) 
  out
  if(log){out <- out + scale_y_continuous(labels = scientific)}
  
  return(out)
}

##########################
#plot in heatmap instead

dataHeatmap <- function(resU, resR, dataU, dataR, targets, geneAnnotation = NULL){
  dataHeatmapInternal <- function(data, sig, geneAnnotation.){
    
    temp <- melt(data, id.vars = "feature")
    
    #time
    temp[,time:=as.factor(as.numeric(gsub(".*T", "", temp$variable)))]
    levels(temp$time) <- c("07:00", "11:00", "15:00", "19:00", "23:00", "03:00")
    #data trans
    temp <- temp[,mean(value, na.rm = T), by = c("feature","time")]
    temp[,scale:=scale(V1), by = feature]
    #add significant
    if(is.null(geneAnnotation.)){
      temp <- merge(temp, sig[,c("feature", "pValAdj", "phase", "SUPERPATHWAY")], by = "feature")
    }else{
      temp <- merge(temp, sig[,c("feature", "pValAdj", "phase")], by = "feature")
    }
    temp$sig <- symnum(x = temp$pValAdj, cutpoints = c(0,0.05,1), symbols = c("*", "ns"))
    
    return(temp)
  }
  
  if(!is.null(geneAnnotation)){
    dataR$feature <- geneAnnotation$Symbol
    dataU$feature <- geneAnnotation$Symbol
    resR$feature <- geneAnnotation$Symbol
    resU$feature <- geneAnnotation$Symbol
  }
  
  temp <- cbind(dataHeatmapInternal(data = dataR[feature %in% targets], sig = resR, geneAnnotation. = geneAnnotation), diet = "TRF")
  temp2 <- cbind(dataHeatmapInternal(data = dataU[feature %in% targets], sig = resU, geneAnnotation. = geneAnnotation), diet = "EXF")
  temp <- rbind(temp, temp2)
  temp$sig <- factor(temp$sig, levels = c("*", "ns"))
  return(temp)
}


plotHeatmap <- function(data, order, cutline, leg.title = "Relative intensity (z-score)"){
  
  ggplot(data, 
         aes(x = feature, y = time, alpha = sig, fill = scale)) +
    geom_tile() +
    facet_wrap(~diet) +
    geom_vline(xintercept = cutline, color = "gray", size = 2) +
    coord_flip() +
    scale_x_discrete(limits = order) +
    # scale_y_discrete(expand = c(0.0,-1)) +
    scale_alpha_discrete(range=c(1, .2), guide=FALSE) +
    scale_fill_viridis_c(option = "E", name = leg.title) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 20),
      plot.title = element_text(hjust = 0.5, vjust = .1, size = 30),
      panel.background = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title = element_blank(),
      axis.text.y.right =  element_blank(), 
      legend.position = "bottom",
      legend.direction="horizontal", 
      legend.text = element_text(hjust = .7),
      legend.title = element_text(vjust = 0.85),
      axis.text.x =  element_text(size = 14, color = "black", angle = -45, hjust = 0.2, vjust = 0.3), 
      plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
}


###############
#Time transformer

timeNum <- function(x){
  x <- as.numeric(gsub("T", "", x))
  x[x==19] <- 21
  x
}

###############
#Rain analyzer

rainAnalyze <- function(data, meta, deltaPer, period){
  
  temp <- t(data)
  tempMeta <- meta
  
  timeSeries <- as.numeric(table(tempMeta$TimepointNum))
  
  temp <- rain(x = temp,
               period = period,
               measure.sequence = timeSeries, 
               period.delta = deltaPer,
               deltat = 4,
               verbose = T,
               na.rm = T,
               method = "independent")
  
  temp <- data.frame(Tags = rownames(temp), temp)
  setDT(temp)
  temp[,trough := fifelse((phase + peak.shape)>period, (phase + peak.shape)-period, (phase + peak.shape))]
  temp <- temp[, pValAdj := p.adjust(pVal, method = "BH")]
  temp[,c("peak.shape"):=NULL]
  return(temp)
}
