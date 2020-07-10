rm(list= ls())
source("scriptsPublication/circadian.functions.R")
load("dataPublication/circadian.metabolites.Rdata")

##########################################################################################################################################################################
#reformat res data abit
res$genes.U$Symbol <- dataAnnotation$genes$Symbol[match(res$genes.U$feature, dataAnnotation$genes$genes)]
res$genes.R$Symbol <- dataAnnotation$genes$Symbol[match(res$genes.R$feature, dataAnnotation$genes$genes)]

resSig$genes.U$Symbol <- dataAnnotation$genes$Symbol[match(resSig$genes.U$feature, dataAnnotation$genes$genes)]
resSig$genes.R$Symbol <- dataAnnotation$genes$Symbol[match(resSig$genes.R$feature, dataAnnotation$genes$genes)]

for(x in c("muscle", "serum")){
  resSig[[paste0(x, ".R")]]$SUPERPATHWAY <- dataAnnotation[[x]]$SUPERPATHWAY[match(resSig[[paste0(x, ".R")]]$feature, dataAnnotation[[x]]$feature)]
  res[[paste0(x, ".R")]]$SUPERPATHWAY <- dataAnnotation[[x]]$SUPERPATHWAY[match(res[[paste0(x, ".R")]]$feature, dataAnnotation[[x]]$feature)]
  resSig[[paste0(x, ".U")]]$SUPERPATHWAY <- dataAnnotation[[x]]$SUPERPATHWAY[match(resSig[[paste0(x, ".U")]]$feature, dataAnnotation[[x]]$feature)]
  res[[paste0(x, ".U")]]$SUPERPATHWAY <- dataAnnotation[[x]]$SUPERPATHWAY[match(res[[paste0(x, ".U")]]$feature, dataAnnotation[[x]]$feature)]
}

##########################################################################################################################################################################
#periods distribution

periods <- lapply(res, function(x){
  x[pValAdj<0.05,period]
})

periods <- melt(periods)
periods <- cbind(periods, strsplit2(periods$L1, "\\."))
colnames(periods) <- c("value", "group", "tissue", "diet")
periods$diet <- as.factor(periods$diet)
levels(periods$diet) <- c("TRF", "EXF")
periods$diet <- relevel(periods$diet, ref = "EXF")

levels(periods$tissue) <- c("Skeletal muscle\ntranscripts", "Skeletal muscle\n metabolites", "Serum\nmetabolites") 

periodsSum <- melt(table(periods$tissue, periods$diet, periods$value))
colnames(periodsSum) <- c("tissue", "diet", "period", "value")
periodsSum <- periodsSum[1:3,]
periodsSum$value <- c(2500, 50, 300)
periodsSum$period <- NULL

ggplot(periods, aes(x = value, fill = diet)) + geom_bar(stat = "count", position = "dodge", width = 2.5) + 
  geom_blank(data = periodsSum, aes(y = value, fill = diet), inherit.aes = F) +
  facet_wrap(~tissue, scales = "free") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(breaks = sort(unique(periods$value))) +
  scale_fill_manual(values=c("black","red1")) + 
  ylab("Counts") + 
  xlab("Period (h)") +
  theme_bw() + 
  theme(
    strip.text = element_text(size=14, color="black"),
    axis.title.x = element_text(size=18, color="black"),
    axis.title.y = element_text(size=18, color="black"),
    axis.text.x = element_text(size=16, color="black"),
    axis.text.y = element_text(size=16, color="black"),
    panel.border = element_blank(), 
    panel.grid = element_blank(),
    axis.line = element_line(),
    legend.position="right")
ggsave("figures/hist.period.distribution.pdf", height = 4, width = 7)

#############
#export to excel for publication
out <- as.data.table(periods)
out <- out[,table(value), by = c("tissue", "diet")]
out$time <- rep(c(12,16,20,24,28), 6)
out$tissue <- gsub("\n", " ", out$tissue)
colnames(out)[3] <- "counts"
openxlsx::write.xlsx(out, file = "figures/S2.xlsx")

##################################################################################
#Upset
###############
#Overall metabolomics
#ggupse instead of upsetr since it looks simply better

data <- lapply(X = resSig[c("muscle.U","muscle.R","serum.U", "serum.R")], function(x){
  x$feature
})

names(data) <- c("EXF muscle", "TRF muscle","EXF serum", "TRF serum")

order <- c("EXF muscle", "TRF muscle","EXF serum", "TRF serum")

data <- dataUpset(data)

temp <- data$value

temp <- lapply(temp, function(x){
  ((length(grep("TRF muscle", x)) + length(grep("EXF muscle", x)) == 2) & length(x) == 2)|
    ((length(grep("TRF serum", x)) + length(grep("EXF serum", x)) == 2) & length(x) == 2)|
    ((length(grep("TRF serum", x)) + length(grep("TRF muscle", x)) == 2) & length(x) == 2)|
    ((length(grep("EXF serum", x)) + length(grep("EXF muscle", x)) == 2) & length(x) == 2)|
    ((length(grep("EXF serum", x)) +
        length(grep("EXF muscle", x)) +
        length(grep("TRF muscle", x)) +
        length(grep("TRF serum", x)) == 4)
     & length(x) == 4) |
    length(x) == 1
  
})
sum(unlist(temp))
data <- data[unlist(temp),]

pdf("figures/upset.metabolomics.pdf", width = 5.4, 7)
ggplot(data, aes(x = value)) +
  geom_bar() + 
  geom_text(stat='count', aes(label=..count..), 
            vjust=-1, size = 7, 
            position=position_jitter(width=0,height=.7)) + #nudge_y = -1,
  scale_x_upset(order_by = "freq") +
  scale_y_continuous(expand = c(0,0), limits = c(0,150)) +
  axis_combmatrix(levels = order) +
  theme(panel.background = element_blank(),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        axis.text.y = element_text(size = 18, color = "black"),
        axis.title.y = element_text(size = 24, color = "black", vjust = -12),
        axis.title.x = element_blank(),
        plot.margin = unit(c(0.5,0,0,0), "cm")) +
  theme_combmatrix(combmatrix.panel.point.size = 4,
                   combmatrix.label.extra_spacing = 5, 
                   combmatrix.label.text = element_text(color = "black", size = 20)) +
  ylab("Periodic metabolites")
dev.off()

#############
#export to excel for publication
data <- lapply(X = resSig[c("muscle.U","muscle.R","serum.U", "serum.R")], function(x){
  x$feature
})
upset(fromList(data))
#same output (sans unintersing comparisons)
filler <- function(x, total = max(unlist(lapply(data, length)))){
  c(sort(x),rep("", total-length(x)))
}
dataAll <- data.frame(Muscle_EXF = filler(data$muscle.U), 
                      Muscle_TRF = filler(data$muscle.R), 
                      Serum_EXF = filler(data$serum.U),
                      Serum_TRF = filler(data$serum.R))
data <- fromList(data)
colnames(data) <- c("Muscle_EXF", "Muscle_TRF", "Serum_EXF", "Serum_TRF")
for(x in colnames(data)){
  data[data[,x] == 1, x] <- x
}
data <- apply(data, 1, function(x){
  paste0(x[x!=0], collapse = "/")
})
data <- table(data)
write.xlsx(list(Counts = data, Features = dataAll), file = "figures/F2d.xlsx")
#some extra manual polishing.

##################################################################################
#Venn diagram for muscle, serum and genes
##########################
#muscle
x <- venn.diagram(list(TRF = resSig$muscle.R$feature,
                       EXF = resSig$muscle.U$feature),
                  fill = c("red", "black"),
                  col = "transparent",
                  cex = 7, reverse = T,
                  fontfamily = "sans", 
                  category.names = c("",""),
                  rotation.degree = 0,
                  cat.fontfamily = "sans", 
                  filename = NULL, )

pdf("figures/venn.muscle.pdf")
grid.draw(x)
dev.off()

##########################
#serum

x <- venn.diagram(list(EXF = resSig$serum.U$feature,
                       TRF = resSig$serum.R$feature),
                  fill = c("black", "red"),
                  col = "transparent",
                  cex = c(5,7,6),
                  fontfamily = "sans", 
                  category.names = c("",""),
                  rotation.degree = 170,
                  cat.fontfamily = "sans", 
                  filename = NULL, )
dev.off()
pdf("figures/venn.serum.pdf")
grid.draw(x)
dev.off()

##########################
#genes

x <- venn.diagram(list(EXF = resSig$genes.U$feature,
                       TRF = resSig$genes.R$feature),
                  fill = c("black", "red"),
                  col = "transparent",
                  cex = c(6,7,7),
                  fontfamily = "sans", 
                  category.names = c("",""),
                  rotation.degree = 170,
                  cat.fontfamily = "sans", 
                  filename = NULL, )

pdf("figures/venn.genes.pdf")
grid.draw(x)
dev.off()

#############
#export to excel for publication
vennCounter <- function(x){
  un <- unique(intersect(x[[1]], x[[2]]))
  x1 <- x[[1]][!x[[1]] %in% un]
  x2 <- x[[2]][!x[[2]] %in% un]
  l <- max(c(length(un), length(x1), length(x2)))
  data.frame(union = filler(un, l),
             EXF_unique = filler(x1,l),
             TRF_unique = filler(x2,l))
}

out <- lapply(c("genes","muscle","serum"), function(x){
  vennCounter(lapply(resSig[grepl(x,names(resSig))], function(x) x$feature))
})
names(out) <- c("genes","muscle","serum")

write.xlsx(out, file = "figures/F2a-c.xlsx")
##################################################################################
#phase shift
dataPhasePlot <- lapply(c("muscle", "serum", "genes"), function(tiss){
  U <- paste0(tiss, ".U")
  R <- paste0(tiss, ".R")
  
  phaseShiftAnalyzer(resU = resSig[[U]], resR = resSig[[R]])$all
})

names(dataPhasePlot) <- c("muscle", "serum", "genes")

plotGenes <- plotPhaseShift(dataPhasePlot$genes) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,500)) + 
  # scale_x_continuous(breaks = seq(-12,12,4), labels = seq(-12,12,4), limits = c(-14,14)) +
  theme(axis.text.x = element_text(size = 13)) +
  ylab("Skeletal Muscle\nTranscripts (counts)")
plotGenes

plotMuscle <- plotPhaseShift(dataPhasePlot$muscle) + 
  theme(axis.text.x = element_text(size = 13)) +
  # scale_x_continuous(breaks = seq(-12,12,4), labels = seq(-12,12,4), limits = c(-14,14)) +
  ylab("Skeletal Muscle\nMetabolite (counts)") 
plotMuscle

plotSerum <- plotPhaseShift(dataPhasePlot$serum) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,120)) + 
  # scale_x_continuous(breaks = seq(-12,12,4), labels = seq(-12,12,4), limits = c(-14,14)) +
  theme(axis.text.x = element_text(size = 13))  +
  ylab("Serum\nMetabolite (counts)")
plotSerum

plotGenes + theme(axis.text.x = element_blank(),
                  axis.title.x = element_blank()) +
  plotMuscle + theme(axis.text.x = element_blank(),
                     axis.title.x = element_blank()) +
  plotSerum + plot_layout(nrow = 3) &
  scale_x_continuous(breaks = seq(-8,8,4), labels = seq(-8,8,4), limits = c(-10,10))

ggsave(filename = "figures/final.phase.shift.vs.EXF.pdf", width = 3, height = 5)

#############
#export figure for publication

write.xlsx(dataPhasePlot, file = "figures/F4b.xlsx")


##################################################################################
#circadian misalignment

##########################
#serum

cortisol <- dataRaw$serum[feature == "cortisol",-130]
cortisol <- melt(cortisol)
cortisol <- data.table(strsplit2(cortisol$variable, "\\."), value = cortisol)

cortisol <- cortisol[,.SD[value.value == max(value.value), c("value.variable", "value.value", "V3")], by = c("V2","V1")][,c("V2", "V1", "V3")]
cortisol$V3 <- timeNum(cortisol$V3)
colnames(cortisol) <- c("Diet", "ID", "TimepointNum")

features <- dataRaw$serum[!feature == "cortisol",]
features <- melt(features)
features <- data.table(strsplit2(features$variable, "\\."), value = features)

features <- features[,.SD[value.value == max(value.value), c("value.variable", "value.value", "V3")], by = c("V2","V1", "value.feature")][,c("V2","V1","value.feature", "V3")]
colnames(features) <- c("Diet", "ID", "feature", "TimepointNum")
features$TimepointNum <- timeNum(features$TimepointNum)

features[,sampleID:=paste0(Diet,"_", ID)]
cortisol[,sampleID:=paste0(Diet,"_", ID)]

features <- merge(cortisol, features, by = "sampleID")
features <- features[,c(1:4,7:8)]
colnames(features) <- c("sampleID","Diet","ID","TimepointCortisol", "feature", "TimepointFeature")
features[,phase:=TimepointFeature-TimepointCortisol]

features$Diet <- as.factor(features$Diet)
levels(features$Diet) <- c("TRF", "EXF")
features$Diet <- relevel(features$Die, ref = "EXF")

features[phase>12]$phase <- features[phase>12]$phase -24
features[abs(phase)>12]$phase <- features[abs(phase)>12]$phase + 24

ks.test(x = features[Diet == "EXF",phase],
        y = features[Diet == "TRF",phase], exact = F)

featuresSerum <- features

plotMissalignSerum <- ggplot(features, aes(x = phase, col = Diet)) + 
  geom_line(stat = "density", adjust = 1.5, size =2) +
  geom_vline(xintercept = 0, lty = 2) +
  theme(panel.background = element_blank(),
        plot.title = element_text(vjust = 4),
        axis.title.x = element_text(size=22, color="black"),
        axis.title.y = element_text(size=18, color="black"),
        axis.text.x = element_text(size=20, color="black"),
        axis.text.y = element_text(size=20, color="black"),
        panel.border = element_blank(), 
        panel.grid = element_blank(),
        axis.line = element_line(),
        plot.margin = unit(c(.5,.5,.5,.5), "cm")) +
  scale_x_continuous(labels = seq(-12, 12, 4), 
                     breaks = seq(-12, 12, 4),
                     expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 0.06)) +
  scale_color_manual(values = c("black", "red")) +
  xlab("Relative phase shift \n (h compared cortisol peak)") +
  ylab("Density (all participants)")
plotMissalignSerum
ggsave("figures/circadian.serum.misalignment.density.counts.pdf", width = 5, height = 5)

#<>==================================<>#
#<>===========Double check===========<>#
#<>==================================<>#
temp <- dataRaw$serum[dataAnnotation$serum$Label == "valine", grepl("RF029",colnames(dataRaw$serum)), with = F]
which.max(temp[,seq(1,12,2), with = F])#T9
which.max(temp[,seq(2,12,2), with = F])#T21
tempCort <- dataRaw$serum[dataAnnotation$serum$Label == "cortisol", grepl("RF029",colnames(dataRaw$serum)), with = F]
which.max(tempCort[,seq(1,12,2), with = F])#T5
which.max(tempCort[,seq(2,12,2), with = F])#T1
#TRF valine shift of 4
#EXF valine shift of -4
features[Diet == "TRF"&ID=="RF029"&feature=="valine"]
features[Diet == "EXF"&ID=="RF029"&feature=="valine"]
#its good

##########################
#muscle

features <- dataRaw$muscle[!feature == "cortisol",]
features <- melt(features)
features <- data.table(strsplit2(features$variable, "\\."), value = features)

features <- features[,.SD[value.value == max(value.value), c("value.variable", "value.value", "V3")],
                     by = c("V2","V1", "value.feature")][,c("V2","V1","value.feature", "V3")]
colnames(features) <- c("Diet", "ID", "feature", "TimepointNum")
features$TimepointNum <- timeNum(features$TimepointNum)

features[,sampleID:=paste0(Diet,"_", ID)]
cortisol[,sampleID:=paste0(Diet,"_", ID)]

features <- merge(cortisol, features, by = "sampleID")
features <- features[,c(1:4,7:8)]
colnames(features) <- c("sampleID","Diet","ID","TimepointCortisol", "feature", "TimepointFeature")
features[,phase:=TimepointFeature-TimepointCortisol]

features$Diet <- as.factor(features$Diet)
levels(features$Diet) <- c("TRF", "EXF")
features$Diet <- relevel(features$Die, ref = "EXF")

#cant have more than 12 h. remove 24 from +12
features[phase>12]$phase <- features[phase>12]$phase -24
features[abs(phase)>12]$phase <- features[abs(phase)>12]$phase + 24

ks.test(x = features[Diet == "EXF",phase],
        y = features[Diet == "TRF",phase], exact = F)

plotMissalignMuscle <- ggplot(features, aes(x = phase, col = Diet)) + 
  geom_line(stat = "density", adjust = 1.5, size =2) +
  geom_vline(xintercept = 0, lty = 2) +
  theme(panel.background = element_blank(),
        plot.title = element_text(vjust = 4),
        axis.title.x = element_text(size=22, color="black"),
        axis.title.y = element_text(size=18, color="black"),
        axis.text.x = element_text(size=20, color="black"),
        axis.text.y = element_text(size=20, color="black"),
        panel.border = element_blank(), 
        panel.grid = element_blank(),
        axis.line = element_line(),
        legend.position="none",
        plot.margin = unit(c(.5,.5,.5,.5), "cm")) +
  scale_x_continuous(labels = seq(-12, 12, 4), 
                     breaks = seq(-12, 12, 4),
                     expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 0.08)) +
  scale_color_manual(values = c("black", "red")) +
  xlab("Relative phase shift \n(h compared serum cortisol)") +
  ylab("Density (all participants)")
plotMissalignMuscle
ggsave("figures/circadian.muscle.misalignment.density.counts.pdf", width = 5, height = 5)

plotMissalignMuscle + theme(axis.title.x = element_blank(), legend.position = "none") + 
  plotMissalignSerum + plot_layout(ncol = 1) + theme(legend.position = "bottom")
ggsave("figures/final.circadian.misalignment.density.counts.pdf", width = 5, height = 7)

#############
#export for publication

write.xlsx(list(serum = featuresSerum,
                muscle = features),
           file = "figures/F4a.xlsx")

##################################################################################
#Dimensionality reduction
###############
#TSNE

#TSNE expects complete data. Impute only for TSNE
dataTSNE <- NULL

for(i in names(resSig)){
  dataTSNE[[i]] <- as.matrix(dataRaw[[i]][,-ncol(dataRaw[[i]]), with = F])
  rownames(dataTSNE[[i]]) <- dataAnnotation[[gsub("\\.R|\\.U","", i)]]$feature
  dataTSNE[[i]] <- dataTSNE[[i]][resSig[[i]]$feature,]#
  dataTSNE[[i]] <- t(scale(t(dataTSNE[[i]])))
}

clusterTSNE <- lapply(dataTSNE[c("genes.U", "genes.R", "serum.U", "serum.R")], function(x) Rtsne(x, perplexity = 30, verbose=T, normalize = F))
clusterTSNE$muscle.R <- Rtsne(dataTSNE$muscle.R, verbose=T, perplexity = 10,normalize = F)
clusterTSNE$muscle.U <- Rtsne(dataTSNE$muscle.U, verbose=T, perplexity = 10, normalize = F)

temp <- names(clusterTSNE)
clusterTSNE <- lapply(names(clusterTSNE), function(x){
  data.frame(clusterTSNE[[x]]$Y, Peak = resSig[[x]][, Acrophase])
})
names(clusterTSNE) <- temp
rm(temp)

cols <- c("7:00" = "forestgreen", "11:00" = "yellowgreen", "15:00" = "yellow2", "19:00" = "blueviolet", "23:00" = "navyblue", "3:00" = "darkcyan")

for(i in names(clusterTSNE)){
  x <- ggplot(data = clusterTSNE[[i]], aes(x=X1,y=X2, color=Peak)) + 
    geom_point(size = 2.5, alpha = .8) + 
    scale_colour_gradientn(colors = cols,
                           guide = "colourbar", labels = c("7:00", "11:00", "15:00", "19:00", "23:00", "3:00"), breaks = c(0, 4, 8, 13, 16, 20)) +
    theme(panel.border = element_rect(fill = NA, size = 1),
          panel.background = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          legend.text = element_text(size=15),
          legend.position = "right")
  print(x)
}

#make final patchwork figure
plot <- NULL
for(i in names(clusterTSNE)){
  plot[[i]] <- ggplot(data = clusterTSNE[[i]], aes(x=X1,y=X2, color=Peak)) + 
    scale_colour_gradientn(colors = cols,
                           guide = "colourbar", labels = c("7:00", "11:00", "15:00", "19:00", "23:00", "3:00"), breaks = c(0, 4, 8, 13, 16, 20)) +
    theme(panel.border = element_rect(fill = NA, size = 1),
          panel.background = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          legend.text = element_text(size=15))
}


plot$genes.U + ggtitle("Muscle EXF transcripts") + geom_point(size = 1, alpha = .8) + guides(colour = FALSE) +
  plot$genes.R + ggtitle("Muscle TRF transcripts") + geom_point(size = 1, alpha = .8) + guides(colour = FALSE) +
  plot$muscle.U + ggtitle("Muscle EXF metabolites") + geom_point(size = 2, alpha = .8) + guides(colour = FALSE) +
  plot$muscle.R + ggtitle("Muscle TRF metabolites") + geom_point(size = 2, alpha = .8) + guides(colour = FALSE) +
  plot$serum.U + ggtitle("Serum EXF metabolites") + geom_point(size = 2, alpha = .8) + guides(colour = FALSE) +
  plot$serum.R + ggtitle("Serum TRF metabolites") + geom_point(size = 2, alpha = .8) +
  guide_area() + 
  plot_layout(nrow = 2, ncol = 4, guides = "collect")

ggsave("figures/final.tsne.pdf", width = 12, height = 7)

###############
#PCA

clusterPCA <- NULL

for(i in c("genes", "muscle", "serum")){
  clusterPCA[[i]] <- as.matrix(dataRaw[[i]][,- ncol(dataRaw[[i]]), with = F])
  clusterPCA[[i]] <- prcomp(t(clusterPCA[[i]]), center = T, scale. = T)
}

pcsImportance <- NULL

for(i in names(clusterPCA)){
  temp <- data.frame(PC1 = melt(clusterPCA[[i]]$x[,1]),
                     PC2 = melt(clusterPCA[[i]]$x[,2]),
                     PC3 = melt(clusterPCA[[i]]$x[,3]))
  colnames(temp) <- c("PC1", "PC2", "PC3")
  
  pcsImportance[[i]] <- paste(names(summary(clusterPCA[[i]])$importance[2,1:3]),
                              " (",
                              (round(summary(clusterPCA[[i]])$importance[2,1:3],2)*100),"%)", sep = "")
  
  temp$Time <- factor(dataGroups[[i]]$Timepoint, levels = c("T1", "T5", "T9", "T13", "T17", "T19"))
  levels(temp$Time) <- c("7:00","11:00","15:00","19:00","23:00","3:00")
  
  temp$Diet <- as.factor(dataGroups[[i]]$Food)
  temp$Diet <- relevel(temp$Diet, ref = "U")
  levels(temp$Diet) <- c("EXF", "TRF")
  
  clusterPCA[[i]] <- temp
  
}

plotPCA <- NULL
for(x in names(clusterPCA)){
  plotPCA[[x]] <- ggplot(clusterPCA[[x]], aes(x = PC1, y = PC2, color = Time, shape = Diet)) + 
    geom_point(size = 2) +
    scale_fill_manual(values = cols, aesthetics = "color") +
    theme_bw() + 
    theme(axis.title.y = element_text(size=20, color="black"),
          axis.title.x = element_text(size=20, color="black"),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_blank(),
          legend.text = element_text(size=15),
          legend.title = element_text(size=20)) +
    xlab(pcsImportance[[x]][1]) + 
    ylab(pcsImportance[[x]][2])
  
  print(plotPCA[[x]])
}

plotPCA$genes + plotPCA$muscle + plotPCA$serum + plot_layout(ncol = 3, nrow = 1, guides = "collect")
ggsave("figures/final.PCA.pdf", height = 3.5, width = 11)

##################################################################################
#Enrichment

##################################################################################
#genes
dataEnrichGenes <- rbind(cbind(Diet = "EXF", resSig$genes.U),
                         cbind(Diet = "TRF", resSig$genes.R))
dataEnrichGenes <- as.data.frame(dataEnrichGenes)
dataEnrichGenes$Entrez <- dataAnnotation$genes$entrez[match(dataEnrichGenes$feature,
                                                            dataAnnotation$genes$feature)]
dataEnrichGenes$phase <- sapply(dataEnrichGenes$Acrophase, nearestTime)

enrichGenesOverallGO <- compareCluster(Entrez ~ Diet, fun = "enrichGO", data = dataEnrichGenes, 
                                       OrgDb = "org.Hs.eg.db", ont = "MF", readable = T,
                                       universe = dataAnnotation$genes$entrez)

dataEnrichGenesHalf <- rbind(cbind(Diet = "EXF", resSigHalf$genes.U),
                             cbind(Diet = "TRF", resSigHalf$genes.R))
dataEnrichGenesHalf$Entrez <- dataAnnotation$genes$entrez[match(dataEnrichGenesHalf$feature, 
                                                                dataAnnotation$genes$feature)]
dataEnrichGenesHalf$phase <- sapply(dataEnrichGenesHalf$Acrophase, nearestTime)

enrichGenesOverallGOHalf <- compareCluster(Entrez ~ Diet, fun = "enrichGO", data = dataEnrichGenesHalf, 
                                           OrgDb = "org.Hs.eg.db", ont = "MF", readable = T,
                                           universe = dataAnnotation$genes$entrez)

enrichGenesOverallGOHalf <- simplify(enrichGenesOverallGOHalf, cutoff = 0.5)

plotGenesOverall <- dotplot(enrichGenesOverallGO, x= ~ Diet, showCategory = 60) +
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.x = element_blank())
plotGenesOverall
ggsave(file = "figures/dotplot.genes.OVERALL.Whole.pdf", height = 6, width = 8)

plotGenesOverallHalf <- dotplot(enrichGenesOverallGOHalf, x= ~ Diet, showCategory = 60) +
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.x = element_text(size = 18))
plotGenesOverallHalf
ggsave(file = "figures/dotplot.genes.OVERALL.Half.pdf", height = 6, width = 8)

##################################################################################
#muscle

dataEnrichMuscle <- rbind(cbind(Diet = "EXF", resSig$muscle.U),
                          cbind(Diet = "TRF", resSig$muscle.R))
term2gene <- data.frame(term = dataAnnotation$muscle$SUBPATHWAY, 
                        feature = dataAnnotation$muscle$Label)
enrichMuscle <- compareCluster(feature ~ Diet, data = dataEnrichMuscle,
                               fun = "enricher", 
                               TERM2GENE = term2gene, minGSSize = 5,
                               universe = dataAnnotation$muscle$feature)

plotMuscle <- dotplot(enrichMuscle, x= ~ Diet, showCategory = 80) +
  scale_y_discrete(limits = levels(enrichMuscle@compareClusterResult$Description)) +
  theme(axis.text.x =  element_blank())
plotMuscle
ggsave("figures/dotplot.muslce.Whole.pdf")

dataEnrichMuscleHalf <- rbind(cbind(Diet = "EXF", resSigHalf$muscle.U),
                              cbind(Diet = "TRF", resSigHalf$muscle.R))
term2gene <- data.frame(term = dataAnnotation$muscle$SUBPATHWAY, 
                        feature = dataAnnotation$muscle$Label)
dataEnrichMuscleHalf$phase <- sapply(dataEnrichMuscleHalf$Acrophase, nearestTime)

enrichMuscleHalf <- compareCluster(feature ~ Diet, data = dataEnrichMuscleHalf,
                                   fun = "enricher",
                                   TERM2GENE = term2gene,
                                   minGSSize = 5,
                                   universe = dataAnnotation$muscle$feature)
plotMuscleHalf <- dotplot(enrichMuscleHalf)
plotMuscleHalf
ggsave("figures/supplementary.dotplot.muslce.Half.pdf")

##################################################################################
#serum

dataEnrichSerum <- rbind(cbind(Diet = "EXF", resSig$serum.U),
                         cbind(Diet = "TRF", resSig$serum.R))
term2gene <- data.frame(term = dataAnnotation$serum$SUBPATHWAY, 
                        feature = dataAnnotation$serum$Label)
dataEnrichSerum$phase <- sapply(dataEnrichSerum$Acrophase, nearestTime)

enrichSerum <- compareCluster(feature ~ Diet, data = dataEnrichSerum,
                              fun = "enricher",
                              TERM2GENE = term2gene,
                              universe = dataAnnotation$serum$feature)

plotSerum <- dotplot(enrichSerum, x= ~ Diet) +
  theme(axis.text.x = element_text(size = 16),
        plot.margin = unit(c(0,0,0,0), "cm"), 
        plot.background = element_blank())
plotSerum

dataEnrichSerumHalf <- rbind(cbind(Diet = "EXF", resSigHalf$serum.U),
                             cbind(Diet = "TRF", resSigHalf$serum.R))
term2gene <- data.frame(term = dataAnnotation$serum$SUBPATHWAY, 
                        feature = dataAnnotation$serum$Label)
dataEnrichSerumHalf$phase <- sapply(dataEnrichSerumHalf$Acrophase, nearestTime)

enrichSerumHalf <- compareCluster(feature ~ Diet, data = dataEnrichSerumHalf,
                                  fun = "enricher", 
                                  TERM2GENE = term2gene,
                                  universe = dataAnnotation$serum$feature)

plotSerumHalf <- dotplot(enrichSerumHalf)
plotSerumHalf

###############
#plot final figure
plotGenesOverall + theme(legend.position = "none") + 
  plotMuscle+theme(legend.position = "none") + 
  plotSerum+theme(legend.position = "bottom", axis.text.x = element_text(size = 12)) + 
  plot_layout(nrow = 3, guides = 'collect', height = c(.75, .1, .1)) &
  scale_size(range = c(1,6), limits = c(0.005, .4), breaks = seq(0.005, .4, length.out = 3), labels = scales::percent_format(.1)) &
  scale_color_gradient2(high = '#e41a1c', mid = "#3399FF", low = '#0033CC', midpoint = 0.025, limits = c(0, 0.05), breaks = c(0,0.025,0.05), name = "FDR") &
  theme(plot.tag = element_text(size = 12),
        axis.text.y = element_text(size = 8),
        strip.text = element_text(size = 9.5),
        plot.margin = unit(c(-.1,0,-.1,0), "cm")) 

ggsave("figures/final.enrichments.pdf", height = 6, width = 8)

plotGenesOverallHalf + theme(legend.position = "none", axis.text = element_blank()) + 
  plotMuscleHalf + theme(legend.position = "none", axis.text = element_blank()) + 
  plotSerumHalf + theme(legend.position = "bottom", axis.text.x = element_text(size = 15)) + 
  plot_layout(nrow = 3, guides = 'collect', height = c(.75, .1, .1)) &
  scale_size(range = c(1,6), limits = c(0.001, .4), breaks = seq(0.001, .4, length.out = 3), labels = scales::percent_format(.1)) &
  scale_color_gradient2(high = '#e41a1c', mid = "#3399FF", low = '#0033CC', midpoint = 0.025, limits = c(0, 0.05), breaks = c(0,0.025,0.05), name = "FDR") &
  theme(plot.tag = element_text(size = 12),
        axis.text.y = element_text(size = 8),
        strip.text = element_text(size = 9.5),
        plot.margin = unit(c(-.1,0,-.1,0), "cm")) 
ggsave("figures/final.enrichments.half.pdf", height = 6, width = 8)

#############
#export data for publication

out <- rbind(cbind(tissue = "Muscle Transcripts",
                   enrichGenesOverallGO@compareClusterResult),
             cbind(tissue = "Muscle Metabolites",
                   enrichMuscle@compareClusterResult),
             cbind(tissue = "Serum Metabolites",
                   enrichSerum@compareClusterResult))

write.xlsx(out, file = "figures/F6a.xlsx")

out <- rbind(cbind(tissue = "Muscle Transcripts",
                   enrichGenesOverallGOHalf@compareClusterResult),
             cbind(tissue = "Muscle Metabolites",
                   enrichMuscleHalf@compareClusterResult),
             cbind(tissue = "Serum Metabolites",
                   enrichSerumHalf@compareClusterResult))

write.xlsx(out, file = "figures/S3a.xlsx")


##################################################################################
#Periodicity of transporters

trans <- unique(unlist(strsplit(enrichGenesOverallGO@compareClusterResult[enrichGenesOverallGO@compareClusterResult$ID %in% 
                                                                            c("GO:0046943","GO:0005342", "GO:0008028"), "geneID"], split = "\\/")))
temp <- dataRaw$genes
temp$feature <- dataAnnotation$genes$Symbol
temp2 <- res
temp2$genes.U$feature <- dataAnnotation$genes$Symbol
temp2$genes.R$feature <- dataAnnotation$genes$Symbol

plotHeatmapData <- dataHeatmap(resU = res$genes.U, 
                               resR = res$genes.R, 
                               dataU = dataRaw$genes.U, 
                               dataR = dataRaw$genes.R,
                               targets = c(trans),#, "SLC27A4", "CD36"
                               geneAnnotation = dataAnnotation$genes)


order <- plotHeatmapData[diet == "TRF", which(max(scale)==scale), by = c("feature")][order(V1), feature]
order <- c(order[order %in% trans], order[!order %in% trans])

plotHeatmap(plotHeatmapData, order = rev(order), cutline = c(0)) +
  theme(axis.text.y = element_text(size=12),
        legend.title.align = .5,
        plot.margin = unit(c(0,.5,0,0), "cm"))
ggsave("figures/heatmap.transporters.pdf", width = 6, height = 6.5)

#############
#export data for publication

plotHeatmapData <- plotHeatmapData[,-3]
colnames(plotHeatmapData)[3] <- "scaledValue"
write.xlsx(plotHeatmapData, file = "figures/F7a.xlsx")

##################################################################################
#proportion of total metabolites

pathCounts <- lapply(resSig[c("serum.U", "serum.R", "muscle.U", "muscle.R")], function(x){
  table(x$SUPERPATHWAY)/sum(table(x$SUPERPATHWAY))*100
})

pathCounts <- melt(pathCounts, stringsAsfactors = F)
colnames(pathCounts) <- c("Pathway", "Counts", "Condition")

pathCounts$Pathway <- as.character(pathCounts$Pathway)
pathCounts$Pathway[pathCounts$Pathway %in% c("Partially Characterized Molecules", "Xenobiotics", "Cofactors and Vitamins")] <- "Other" 
pathCounts$Pathway <- as.factor(pathCounts$Pathway)

pathCounts <- as.data.table(pathCounts)
pathCounts[Pathway == "Other", Counts := sum(Counts), by = Condition]
pathCounts <- pathCounts[!duplicated(paste0(Pathway, Counts))]

pathCounts$Pathway <- gdata::reorder.factor(x = pathCounts$Pathway,
                                            new.order = levels(pathCounts$Pathway)[c(1:3,5:7,4)])

pathCounts$Condition

plotTotal <- ggplot(pathCounts, aes(y = Counts, x = Condition, fill = Pathway)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_x_discrete(expand = c(0,0.2), 
                   position = "bottom", 
                   labels = rep(c("TRF\n", "EXF\n"), 2), 
                   limits = c("serum.R", "serum.U", "muscle.R", "muscle.U")) +
  scale_fill_manual(values = c('#1b9e77','#d95f02','#e7298a','steelblue4',"gray",'#e6ab02' ,'#7570b3')) +
  scale_y_continuous(name = "Proportion (%)", position = "right") +
  theme(panel.background = element_blank(),
        axis.text.x = element_text(size = 14,color = "black"),
        axis.text.y = element_text(size = 8, colour = "black", vjust = .8),
        axis.ticks.x = element_line(size = 1),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") + labs(size = 50)# + xlab("")
plotTotal

plotTotal + theme(axis.text.y = element_blank(), legend.position = "none")
ggsave("figures/barplot.propotion.of.total.periodic.pdf", height = 1, width = 6)

#############
#export for publication
pathCounts <- cbind(pathCounts[,-3], strsplit2(pathCounts$Condition, split = "\\."))
colnames(pathCounts)[3:4] <- c("Tissue", "Diet")
pathCounts$Diet <- as.factor(pathCounts$Diet)
levels(pathCounts$Diet) <- c("TRF", "EXF")

write.xlsx(pathCounts, file = "figures/F6b.xlsx")

##################################################################################
#serum waffle plot
colorsWaffle <- c("Amino Acid" = '#1b9e77',
                  "Carbohydrate" = '#d95f02',
                  "Energy" = '#e7298a', 
                  "Nucleotide" = 'steelblue4',
                  "Other" = "gray", 
                  "Peptide" = '#e6ab02' , 
                  "Lipid" = '#7570b3')

dataWaffle <- NULL
enrich <- NULL
for(tiss in c("muscle", "serum")){
  U <- paste0(tiss,".U")
  R <- paste0(tiss,".R")
  
  temp <- rbind(cbind(resSig[[U]],cond = U), cbind(resSig[[R]],cond = R))
  
  temp$diet <- strsplit2(temp$cond, "\\.")[,2]
  temp$diet <- as.factor(temp$diet)
  levels(temp$diet) <- c("TRF", "EXF")
  temp$time <- sapply(temp$Acrophase, nearestTime)
  
  tempAnnot <- dataAnnotation[[tiss]]
  tempAnnot[tempAnnot$SUPERPATHWAY %in% c("Xenobiotics", "Partially Characterized Molecules", "Cofactors and Vitamins"),]$SUPERPATHWAY <- "Other"
  # temp[pathway %in% c("Xenobiotics", "Partially Characterized Molecules", "Cofactors and Vitamins"), pathway] <- "Other"
  
  enirchTime <- compareCluster(feature ~ diet + time, data = temp,
                               fun = "enricher", minGSSize = 5,
                               TERM2GENE = data.frame(term = tempAnnot$SUPERPATHWAY, 
                                                      feature = tempAnnot$Label),
                               universe = dataAnnotation[[tiss]]$feature)
  
  enrich[[tiss]] <- enirchTime@compareClusterResult
  
  temp <- cbind(temp[,names(table(SUPERPATHWAY)), by = c("diet", "time")], temp[,table(SUPERPATHWAY), by = c("diet", "time")][,"V1"])
  colnames(temp) <- c("diet", "time", "pathway", "value")
  temp$group <- paste0(temp$diet, "-", temp$time)
  temp[pathway %in% c("Xenobiotics", "Partially Characterized Molecules", "Cofactors and Vitamins")]$pathway <- "Other"
  temp[pathway == "Other", value:=sum(value), by = "group"]
  temp <- temp[!duplicated(paste0(group, pathway)), ]
  
  temp$pathway <- as.factor(temp$pathway)
  temp$pathway <- gdata::reorder.factor(temp$pathway, new.order = names(colorsWaffle))
  
  temp$group <- as.factor(temp$group)
  temp$group <- gdata::reorder.factor(temp$group, new.order = paste0(rep(c("EXF-", "TRF-"), 6),
                                                                     sort(rep(seq(0,20,4), 2))))
  
  dataWaffle[[tiss]] <- temp
}

facetRenamer <- function(variable,value){
  return(gsub("-.*", "", value))
}

labeler <- function(x) labels(x * 5)

table(dataWaffle$muscle$group)
dataWaffle$muscle <- rbind(dataWaffle$muscle, dataWaffle$muscle[1,])
dataWaffle$muscle[nrow(dataWaffle$muscle), "group"] <- "TRF-8"
dataWaffle$muscle[nrow(dataWaffle$muscle), "value"] <- 1
dataWaffle$muscle[nrow(dataWaffle$muscle), "pathway"] <- "stupid_ass_hack" #as the name implies, stupid ass hack to get a blank facet
dataWaffle$muscle[nrow(dataWaffle$muscle), "diet"] <- "TRF"

colorsWaffle <- c(colorsWaffle, "stupid_ass_hack" = "white" )
plotMuscle <- NULL
plotSerum <- NULL

for(i in c("EXF", "TRF")){
  plotSerum[[i]] <- ggplot(dataWaffle$serum[order(pathway)][diet == i], aes(fill = pathway, values = value)) +
    geom_waffle(color = "white", size = .25, n_rows = 5, flip = TRUE) +
    facet_wrap(~group, nrow = 1, strip.position = "bottom", drop = T, labeller = facetRenamer)  +
    scale_x_discrete() + 
    scale_fill_manual(values = colorsWaffle, breaks = names(colorsWaffle), limits = names(colorsWaffle)) +
    scale_y_continuous(labels = seq(0,150,50),
                       breaks = seq(0, 30, 10),
                       limits = c(0, 30),
                       expand = c(0,0)) +
    coord_equal() +
    ylab("Counts") +
    theme_minimal() +
    guides(fill = guide_legend(reverse = F)) +
    theme(legend.position = "none",
          #plot.margin = unit(c(0,0,0,0), "cm"),
          # strip.text = element_text(size = 12),
          panel.grid = element_blank(), 
          axis.line.y = element_line(),
          axis.ticks.y = element_line(),
          axis.text = element_text(color = "black")) 
  
  plotMuscle[[i]] <- ggplot(dataWaffle$muscle[order(pathway)][diet==i], aes(fill = pathway, values = value)) +
    geom_waffle(color = "white", size = .25, n_rows = 5, flip = TRUE) +
    facet_wrap(~group, nrow = 1, strip.position = "bottom", drop = T, labeller = facetRenamer)  +
    scale_x_discrete() + 
    scale_alpha_manual() +
    scale_fill_manual(values = colorsWaffle, breaks = names(colorsWaffle), limits = names(colorsWaffle)) +
    scale_y_continuous(labels = c(0,25),
                       breaks = c(0,6),
                       limits = c(0,6),
                       expand = c(0,0)) +
    coord_equal() +
    ylab("Counts") +
    theme_minimal() +
    guides(fill = guide_legend(reverse = F)) +
    theme(#plot.margin = unit(c(0,0,0,0), "cm"),
      legend.position = "none",
      # strip.text = element_text(size = 12),
      panel.grid = element_blank(), 
      axis.line.y = element_line(),
      axis.ticks.y = element_line(),
      axis.text = element_text(color = "black")) 
  
}

###############
#assemble final

theme_waffle_muscle <- theme(plot.margin = unit(c(0,0,0,0), "cm"),
                             axis.line.y = element_line(),
                             strip.text.x = element_blank(),
                             axis.ticks = element_line(),
                             strip.text = element_blank(), 
                             axis.title = element_blank(),
                             legend.position = "none")
theme_waffle_serum <- theme(axis.line.y = element_line(),
                            strip.text.x = element_blank(),
                            axis.ticks = element_line(),
                            strip.text = element_blank(), 
                            axis.title = element_blank(),
                            plot.margin = unit(c(-2,0,0,0), "cm"),
                            legend.position = "none")

plotMuscle$EXF + theme_waffle_muscle +
  # ggsave("figures/waffle.muscle.EXF.pdf", height = 5, width = 6)
  plotMuscle$TRF + theme_waffle_muscle +
  # ggsave("figures/waffle.muscle.TRF.pdf", height = 5, width = 6)
  plotSerum$EXF + theme_waffle_serum +
  # ggsave("figures/waffle.serum.EXF.pdf", height = 5, width = 6)
  plotSerum$TRF + theme_waffle_serum + plot_layout(ncol = 2, widths = c(.1,.1,.1,.1))
# ggsave("figures/waffle.serum.TRF.pdf", height = 5, width = 6)

ggsave("figures/final.waffle.pdf", height = 3, width = 5)

#############
#export final for publication

dataWaffle <- lapply(dataWaffle, function(x){
  x <- x[,-5]
  x$time <- x$time + 1
  colnames(x)[4] <- "counts"
  x
})

write.xlsx(dataWaffle, "figures/F6c.xlsx")

enrich <- lapply(enrich, function(x){
  x$time <- as.numeric(x$time) + 1
  x[,c(-1,-4,-10)]
})

write.xlsx(enrich, "figures/F6cEnrich.xlsx")

##################################################################################
#core circadian

coreCirc <- c("Clock", "ARNTL", "Npas2", "Per1", "Per2", "Per3", "Cry1", "Cry2", "Nr1d1", "Nr1d2", "Dbp", "RORA")
coreCirc <- toupper(coreCirc)
coreCirc <- dataAnnotation$genes$Symbol[dataAnnotation$genes$Symbol %in% coreCirc] 
names(coreCirc) <- dataAnnotation$genes$genes[dataAnnotation$genes$Symbol %in% coreCirc] 

coreCirc["ENSG00000126368"] <- "REVERBalpha"
coreCirc["ENSG00000174738"] <- "REVERBbeta"

pointsPlot <- dataRaw$genes
pointsPlot$feature <- dataAnnotation$genes$genes
pointsPlot <- pointsPlot[feature %in% names(coreCirc)]
pointsPlot$feature == names(coreCirc)
pointsPlot$feature <- coreCirc

resPlot <- res[c("genes.U", "genes.R")]
resPlot$genes.U$feature <- dataAnnotation$genes$Symbol
resPlot$genes.R$feature <- dataAnnotation$genes$Symbol
resPlot$genes.U[feature == "NR1D1","feature"] <- temp2$genes.R[feature == "NR1D1","feature"] <- "REVERBalpha"
resPlot$genes.U[feature == "NR1D2","feature"] <- temp2$genes.R[feature == "NR1D2","feature"] <- "REVERBbeta"

#reorder datacurves to match dataAnnotation
curvesPlot <- dataCurves$genes[,names(coreCirc), with =F]
colnames(curvesPlot) <- coreCirc
curvesPlot$group <- dataCurves$genes$group

plotCircadian(raw = pointsPlot, 
              resU = resPlot$genes.U, 
              resR = resPlot$genes.R, 
              curves = curvesPlot, 
              features = coreCirc) +
  ylab("Expression (counts per million)") +
  theme(plot.margin=unit(c(.5,1,.5,.5,.5),"cm"))

#Cairo fails to export the triangles..
#ggsave("figures/circadian.core.clock.pdf", device = cairo_pdf())
#export manually...

#############
#export for publication

pointsPlot <- cbind(feature = pointsPlot$feature, pointsPlot[,-ncol(pointsPlot), with = F])
curvesPlot <- cbind(strsplit2(curvesPlot$group, "_"), curvesPlot[,-ncol(curvesPlot), with = F])
colnames(curvesPlot)[1:2] <- c("time", "diet")

write.xlsx(list(rawDataPoints = pointsPlot,
                fittedDataLine = curvesPlot),
           file = "figures/F5a.xlsx")

##################################################################################
#cortisol and stuff
curvesPlot <- dataCurves$serum[,c("corticosterone", "cortisol", "cortisone"), with =F]
colnames(curvesPlot) <- c("corticosterone", "cortisol", "cortisone")
curvesPlot$group <- dataCurves$serum$group

pointsPlot <- dataRaw$serum[feature %in% c("corticosterone", "cortisol", "cortisone"),]

plotCircadian(raw = pointsPlot, 
              resU = res$serum.U, 
              resR = res$serum.R, 
              curves = curvesPlot, 
              features = c("corticosterone", "cortisol", "cortisone"),
              log = T) +
  ylab("Metabolite content (arbitrary values)")

#############
#export for publicatoins

pointsPlot <- cbind(feature = pointsPlot$feature, pointsPlot[,-ncol(pointsPlot), with = F])
curvesPlot <- cbind(strsplit2(curvesPlot$group, "_"), curvesPlot[,-ncol(curvesPlot), with = F])
colnames(curvesPlot)[1:2] <- c("time", "diet")

write.xlsx(list(rawDataPoints = pointsPlot,
                fittedDataLine = curvesPlot),
           file = "figures/F5b.xlsx")


##################################################################################
#Plot the significantly different

resDiff$genes <- merge(resDiff$genes, dataAnnotation$genes[,c("genes", "Symbol")], by.x = "features", by.y = "genes")

diffSerum <- lapply(c("fdr_Acr", "fdr_amp", "fdr_MESOR"), function(x){
  resDiff$serum$features[resDiff$serum[[x]]<0.05 & resDiff$serum$period == 24]
})
names(diffSerum) <- c("Acrophase", "Amplitude", "MESOR")
upset(fromList(diffSerum))

x <- venn.diagram(diffSerum,
                  fill = c("gray80"),
                  col = "black",
                  cex = 7, 
                  reverse = T,
                  fontfamily = "sans", 
                  # direct.area = T,
                  # category.names = c("",""),
                  rotation.degree = 0,
                  cat.fontfamily = "sans", 
                  filename = NULL, )
pdf("figures/venn.params.serum.diff.pdf")
grid.draw(x)
dev.off()

#############
#export data for publication

max(unlist(lapply(diffSerum, length)))
out <- Reduce(cbind,lapply(diffSerum, filler, total = 41))
colnames(out) <- c("Acrophase", "Amplitude", "MESOR")

write.xlsx(out, file = "figures/F7b.xlsx")

###############
#heatmap of only acrophase

plotHeatmapData <- dataHeatmap(resU = res$serum.U, 
                               resR = res$serum.R, 
                               dataU = dataRaw$serum.U, 
                               dataR = dataRaw$serum.R,
                               targets = diffSerum$Acrophase)

order <- plotHeatmapData[diet == "TRF", which(max(scale)==scale), by = c("feature")][order(V1), feature]

acro <- plotHeatmap(plotHeatmapData, order = rev(order), cutline = c(0)) +
  theme(axis.text.y = element_text(size=11),
        legend.title.align = .5,
        plot.margin = unit(c(0,1,0,0), "cm"))

colorsWaffle <- c("Amino Acid" = '#1b9e77',
                  "Carbohydrate" = '#d95f02',
                  "Energy" = '#e7298a', 
                  "Nucleotide" = 'steelblue4',
                  "Other" = "gray", 
                  "Peptide" = '#e6ab02' , 
                  "Lipid" = '#7570b3')

plotHeatmapData[SUPERPATHWAY %in% c("Partially Characterized Molecules", "Xenobiotics"), "SUPERPATHWAY"] <- "Other"
path <- ggplot(plotHeatmapData, aes(x = feature, y= 1, fill = SUPERPATHWAY)) + 
  geom_tile()+ coord_flip() + 
  scale_x_discrete(limits = rev(order)) + 
  scale_fill_manual(values = colorsWaffle, name = "")

acro + 
  path + theme_void() + theme(plot.margin = unit(c(0,0,0,-1), "cm"), legend.position = "right") + 
  plot_layout(widths = c(7,1), ncol = 2, nrow = 1)

ggsave("figures/heatmap.differential.acrophase.serum.pdf", width = 8, height = 8)

#############
#export data for publication

out <- plotHeatmapData
out <- out[,c(-3)]

write.xlsx(out, file = "figures/F7c.xlsx")

###############
#scatterplot x y amplitudes 
resDiffPlot <- resDiff$serum[fdr_amp<0.05 & period == 24,]
resDiffPlot <- data.frame(EXF = resSig$serum.R[match(resDiffPlot$features, feature), amplitude], 
                          TRF = resSig$serum.U[match(resDiffPlot$features, feature), amplitude],
                          feature = resSig$serum.U[match(resDiffPlot$features, feature), feature])

resDiffPlot <- as.data.table(merge(resDiffPlot, 
                                   dataAnnotation$serum[,c("feature", "SUPERPATHWAY")], 
                                   by = "feature"))

resDiffPlot[SUPERPATHWAY %in% c("Partially Characterized Molecules", "Xenobiotics"), "SUPERPATHWAY"] <- "Other"
table(resDiffPlot$SUPERPATHWAY)

resDiffPlot <- list(amplitude = resDiffPlot)

amp <- ggplot(resDiffPlot$amplitude, aes( x = EXF, y = TRF, group = feature, color = SUPERPATHWAY)) + 
  geom_point(size = 1.5) +
  scale_color_manual(values = colorsWaffle, name = "") +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(10^4, 10^8)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(10^4, 10^8)) +
  geom_abline(intercept = 0, slope = 1) + 
  theme_linedraw() +
  theme(axis.title = element_text(size = 20), 
        axis.text = element_text(size = 12), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "gray")
  ) +
  labs(x="Amplitude EXF", y = "Amplitude TRF") 
amp

resDiffPlot$mesor <- resDiff$serum[fdr_MESOR<0.05 & period == 24,]
resDiffPlot$mesor <- data.frame(EXF = resSig$serum.R[match(resDiffPlot$mesor$features, feature), mesor], 
                                TRF = resSig$serum.U[match(resDiffPlot$mesor$features, feature), mesor],
                                feature = resSig$serum.U[match(resDiffPlot$mesor$features, feature), feature])

resDiffPlot$mesor <- as.data.table(merge(resDiffPlot$mesor, 
                                         dataAnnotation$serum[,c("feature", "SUPERPATHWAY")], 
                                         by = "feature"))
resDiffPlot$mesor[SUPERPATHWAY %in% c("Partially Characterized Molecules", "Xenobiotics"), "SUPERPATHWAY"] <- "Other"
table(resDiffPlot$mesor$SUPERPATHWAY)

mesor <- ggplot(resDiffPlot$mesor, aes( x = EXF, y = TRF, group = feature, color = SUPERPATHWAY)) + geom_point() + 
  geom_point(size = 1.5) +
  scale_color_manual(values = colorsWaffle, name = "") +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(10^4.7, 10^8)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(10^4.7, 10^8)) +
  geom_abline(intercept = 0, slope = 1) +
  theme_linedraw() +
  theme(axis.title = element_text(size = 20), 
        axis.text = element_text(size = 12), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "gray")) +
  labs(x="MESOR EXF", y = "MESOR TRF")
mesor

amp + mesor + plot_layout(ncol = 1) & theme(legend.position = "none")
ggsave("figures/scatter.amplitude.mesor.differences.pdf", height = 5, width = 3)

#############
#export data for plublication

write.xlsx(resDiffPlot, file = "figures/F7d.xlsx")

###############
#plot the differences in genes

diffGenes <- resDiff$genes[fdr_amp<0.05|fdr_Acr<0.05|fdr_MESOR<0.05,]$features
diffGenes <- dataAnnotation$genes[dataAnnotation$genes$feature %in% diffGenes,]
temp <- diffGenes$feature
diffGenes <- diffGenes$Symbol
names(diffGenes) <- temp

pointsPlot <- dataRaw$genes
pointsPlot$feature <- dataAnnotation$genes$Symbol
pointsPlot <- pointsPlot[feature %in% diffGenes]

statsPlot <- res
statsPlot$genes.U$feature <- dataAnnotation$genes$Symbol
statsPlot$genes.R$feature <- dataAnnotation$genes$Symbol

#reorder datacurves to match dataAnnotation
curvesPlot <- dataCurves$genes[,names(diffGenes), with =F]
colnames(curvesPlot) <- diffGenes
curvesPlot$group <- dataCurves$genes$group

plotCircadian(raw = pointsPlot, 
              resU = statsPlot$genes.U, 
              resR = statsPlot$genes.R, 
              curves = curvesPlot, 
              features = diffGenes) +
  facet_wrap(~feature, scale = "free_y", drop = T, nrow = 2) + 
  ylab("Expression (counts per million)") +
  theme(plot.margin=unit(c(.5,1,.5,.5,.5),"cm"))

#############
#export data for publicatoin

pointsPlot <- cbind(feature = pointsPlot$feature, pointsPlot[,-ncol(pointsPlot), with = F])
curvesPlot <- cbind(strsplit2(curvesPlot$group, "_"), curvesPlot[,-ncol(curvesPlot), with = F])
colnames(curvesPlot)[1:2] <- c("time", "diet")

write.xlsx(list(rawDataPoints = pointsPlot,
                fittedDataLine = curvesPlot),
           file = "F7e.xlsx")

##################################################################################
#chrononogram

lims <- list(genes = .3001, muscle = .5001, serum = .5)
exclude <- list(genes = c("10%"), serum = c("10%"), muscle = c("10%", "20%", "30%"))

plotsChrono <- NULL
dataChrono <- NULL

for(tiss in c("muscle", "serum", "genes")){
  U <- paste0(tiss,".U")
  R <- paste0(tiss,".R")
  
  temp <- rbind(cbind(resSig[[U]], diet = "EXF"), 
                cbind(resSig[[R]], diet = "TRF"))
  temp$phase <- sapply(temp$Acrophase, nearestTime)
  
  temp <- cbind(temp[, table(phase), by = "diet"], temp[, names(table(phase)), by = "diet"][,"V1"])
  colnames(temp) <- c("diet", "counts", "time")
  
  temp[,`:=`(value = counts/sum(counts), time = time), by = "diet"]
  
  temp$time <- factor(temp$time, levels = seq(0, 20, 4)) 
  levels(temp$time) <- c("07:00", "11:00", "15:00", "19:00", "23:00", "03:00")
  
  dataChrono[[tiss]] <- temp
  
  plotsChrono[[tiss]] <- plotChrono(data = temp, 
                                    limits = c(0,lims[[tiss]]),
                                    labelsExclude = exclude[[tiss]])
}

plotsChrono$genes + plotsChrono$muscle + plotsChrono$serum + 
  plot_layout(ncol = 3, nrow=1) & theme(plot.margin = unit(c(0,0.5,1,1), "cm"))

# ggsave("figures/final.chrono.pdf",width = 12, height = 7, device = cairo_pdf)
#this doesnt work???!
#save manualy as a jpg...

#############
#export data for publication

dataChrono <- lapply(dataChrono, function(x){
  colnames(x)[4] <- "proportion"
  x
})

write.xlsx(dataChrono, file = "figures/F3.chrono.xlsx")

##################################################################################
#amplitude and mesor distributions

amplitude <- NULL
mesor <- NULL

for(tiss in c("muscle", "serum", "genes")){
  U <- paste0(tiss,".U")
  R <- paste0(tiss,".R")
  temp <- data.table(value = resSig[[R]]$amplitude, diet = "TRF")
  temp <- rbind(temp,
                data.table(value = resSig[[U]]$amplitude, diet = "EXF"))
  temp$value <- log10(temp$value)
  amplitude[[tiss]] <- temp
  
  temp <- data.table(value = resSig[[R]]$mesor, diet = "TRF")
  temp <- rbind(temp,
                data.table(value = resSig[[U]]$mesor, diet = "EXF"))
  temp$value <- log10(temp$value)
  mesor[[tiss]] <- temp
  
}

limsAmp <- lapply(amplitude, function(x){
  max(c(density(na.omit(x$value)[x$diet == "TRF"])$y*1.1, density(na.omit(x$value[x$diet == "EXF"]))$y*1.1))
})
limsMesor <- lapply(mesor, function(x){
  max(density(na.omit(x$value))$y)*1.1
})

plotDensityAmp <- NULL
for(tiss in c("muscle", "serum", "genes")){
  
  plotDensityAmp[[tiss]] <- ggplot(amplitude[[tiss]], 
                                   aes(x = value, color = diet, group = diet)) +
    geom_line(stat="density", size = 2, adjust = 1)+
    scale_y_continuous(expand=c(0,0), limits = c(0, limsAmp[[tiss]])) +
    scale_x_continuous(expand=c(0,0)) +
    scale_color_manual(values=c("black","red1")) + 
    xlab("Amplitude (log10)") + 
    ylab("Density") +
    theme_bw() + 
    theme(
      plot.title = element_text(vjust = 4),
      axis.title.x = element_text(size=25, color="black"),
      axis.title.y = element_text(size=25, color="black"),
      axis.text.x = element_text(size=22, color="black"),
      axis.text.y = element_text(size=22, color="black"),
      panel.border = element_blank(), 
      panel.grid = element_blank(),
      axis.line = element_line(),
      legend.position="none") 
  
  
}
plotDensityAmp$genes
plotDensityAmp$muscle
plotDensityAmp$serum

#test diff
lapply(amplitude, function(x) ks.test(x[diet == "EXF", value], x[diet == "TRF", value]))

plotDensityMESOR <- NULL
for(tiss in c("muscle", "serum", "genes")){

  plotDensityMESOR[[tiss]] <- ggplot(mesor[[tiss]], 
              aes(x = value, color = diet)) +
    geom_line(stat="density", size = 2, adjust = 1)+
    scale_y_continuous(expand=c(0,0), limits = c(0, limsMesor[[tiss]])) +
    scale_x_continuous(expand=c(0,0)) +
    scale_color_manual(values=c("black","red1")) + 
    xlab("MESOR (log10)") + 
    ylab("Density") +
    theme_bw() + 
    theme(
      plot.title = element_text(vjust = 4),
      axis.title.x = element_text(size=25, color="black"),
      axis.title.y = element_text(size=25, color="black"),
      axis.text.x = element_text(size=22, color="black"),
      axis.text.y = element_text(size=22, color="black"),
      
      panel.border = element_blank(), 
      panel.grid = element_blank(),
      axis.line = element_line(),
      legend.position="none") 
  
}
plotDensityMESOR$genes
plotDensityMESOR$muscle
plotDensityMESOR$serum

lapply(mesor, function(x) ks.test(x[diet == "EXF", value], x[diet == "TRF", value]))

plotDensityAmp$genes + plotDensityAmp$muscle + ylab("") + plotDensityAmp$serum + ylab("")+
  plotDensityMESOR$genes + plotDensityMESOR$muscle + ylab("") + plotDensityMESOR$serum + ylab("") +
  plot_layout(ncol = 3, nrow = 2)
ggsave("figures/final.density.params.pdf",width = 12, height = 7)

#############
#export data for publication

write.xlsx(amplitude, file = "figures/F3.amp.xlsx")
write.xlsx(mesor, file = "figures/F3.mes.xlsx")

##############################################################################################################################
#enrichment of unique TRF and EXF features

U <- table(dataAnnotation$muscle[dataAnnotation$muscle$feature %in% resSig$muscle.U[!feature %in% intersect(resSig$muscle.U$feature,
                                                                                                            resSig$muscle.R$feature), 
                                                                                    feature],"SUPERPATHWAY"])
R <- table(dataAnnotation$muscle[dataAnnotation$muscle$feature %in% resSig$muscle.R[!feature %in% intersect(resSig$muscle.U$feature,
                                                                                                            resSig$muscle.R$feature), 
                                                                                    feature],"SUPERPATHWAY"])

x <- venn.diagram(list(TRF = resSig$muscle.R$feature,
                       EXF = resSig$muscle.U$feature),
                  fill = c("red", "black"),
                  col = "transparent",
                  cex = c(7,.1,7), reverse = T,
                  fontfamily = "sans", 
                  category.names = c("",""),
                  rotation.degree = 0,
                  cat.fontfamily = "sans", 
                  filename = NULL, )


pdf("figures/barplot.counts.right.join.metabolites.muscle.EXF.pdf")
barplot(-U, horiz = T, xlab = "Counts of metabolites",axes = F, width = 1, cex.names = .9)
axis(1, at = -seq(0,20, 5), labels = seq(0,20, 5))
dev.off()

pdf("figures/barplot.counts.left.join.metabolites.muscle.TRF.pdf")
barplot(R, horiz = T, xlab = "Counts of metabolites",axes = F, width = 1, cex.names = .9)
axis(1, at = seq(0,20, 5), labels = seq(0,20, 5))
dev.off()

pdf("figures/venn.enrichment.muscle.pdf")
grid.draw(x)
dev.off()

#############
#export data for publication
inters <- intersect(resSig$muscle.U$feature, resSig$muscle.R$feature)

write.xlsx(list(EXF = resSig$muscle.U[!feature %in% inters, c("feature", "SUPERPATHWAY")],
                TRF = resSig$muscle.R[!feature %in% inters, c("feature", "SUPERPATHWAY")]),
           file = "figures/S3b.xlsx")

##############################################################################################################################
#supplementary table 1, and 2

load("dataPublication/circadian.metabolites.Rdata")

resSig$genes.U$Symbol <- dataAnnotation$genes$Symbol[match(resSig$genes.U$feature, dataAnnotation$genes$genes)]
resSig$genes.R$Symbol <- dataAnnotation$genes$Symbol[match(resSig$genes.R$feature, dataAnnotation$genes$genes)]

lapply(resSig, colnames)

names(resSig) == names(resSigHalf)

temp <- lapply(names(resSig), function(x){
  y <- rbind(resSig[[x]][,c("feature", "pVal", "pValAdj", "period", "mesor", "amplitude", "Acrophase")],
        resSigHalf[[x]][,c("feature", "pVal", "pValAdj", "period", "mesor", "amplitude", "Acrophase")])
  colnames(y) <- c("Feature", "RAIN.pVal", "RAIN.FDR", "RAIN.period", "MESOR", "Amplitude", "Acrophase")
  y
})

names(temp) <- names(resSig)

temp$genes.U$Symbol <- dataAnnotation$genes$Symbol[match(temp$genes.U$Feature, dataAnnotation$genes$genes)]
temp$genes.R$Symbol <- dataAnnotation$genes$Symbol[match(temp$genes.R$Feature, dataAnnotation$genes$genes)]

names(temp) <- gsub("\\.R", ".TRF", names(temp))
names(temp) <- gsub("\\.U", ".EXF", names(temp))

openxlsx::write.xlsx(temp, file = "dataPublication/supplementary.table.1.xlsx")

coreCirc <- c("Clock", "ARNTL", "Npas2", "Per1", "Per2", "Per3", "Cry1", "Cry2", "Nr1d1", "Nr1d2", "Dbp", "RORA")
coreCirc <- toupper(coreCirc)
coreCirc <- dataAnnotation$genes$Symbol[dataAnnotation$genes$Symbol %in% coreCirc] 
names(coreCirc) <- dataAnnotation$genes$genes[dataAnnotation$genes$Symbol %in% coreCirc] 

temp2 <- resDiff$genes[features %in% names(coreCirc)]
temp2$features == names(coreCirc)[c(-7,-9)]
temp2$Symbol <- coreCirc[c(-7,-9)]
temp2 <- temp2[,c("Symbol", "MESOR_pVal", "Amp_pVal", "Acr_pVal", "fdr_MESOR", "fdr_amp", "fdr_Acr")]
colnames(temp2) <- c("Symbol", "MESOR.pVal", "Amplitude.pVal", "Acrophase.pval", "MESOR.FDR", "Amplitude.FDR", "Acrophase.FDR")

openxlsx::write.xlsx(temp2, file = "dataPublication/supplementary.table.2.xlsx")
