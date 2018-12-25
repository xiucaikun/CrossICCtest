library(ConsensusClusterPlus)
library(survival)
library(survminer)
library(ggsci)
library(sva)
library(CrossICC)
library(RColorBrewer)
library(tidyverse)

# library(ggbiplot)

#### prepare function for ploting survival curve ####
consensusSurvival <- function(SurvivalData){
    SurvivalData <<- SurvivalData
    gene_surv <- Surv(SurvivalData$Time.days./30, SurvivalData$Event)~SurvivalData$Cluster
    gene_survfit <- survfit(gene_surv)
    gene_survdiff <- tryCatch(survdiff(gene_surv), error = function(e) return(NA))
    
    plotcorlor <- pal_d3(alpha = 1)(10)
    fit <- survfit(Surv(SurvivalData$Time.days./30, SurvivalData$Event)~SurvivalData$Cluster, data = SurvivalData)
    surv_plot <- ggsurvplot(fit, main = "", xlab = "Time(Months)",
                            conf.int = FALSE, pval = TRUE, legend = c(0.8, 0.75),
                            legend.title = "Clusters",
                            risk.table.fontsize = 5,
                            risk.table.height = 0.30,
                            break.x.by = 48,
                            risk.table = TRUE, risk.table.y.text.col = TRUE, pval.size = 6,
                            font.main = c(16, "plain", "black"),
                            font.x = c(14, "plain", "black"),font.y = c(16, "plain", "black"),
                            font.tickslab = c(16, "plain", "black"), palette = plotcorlor,
                            ggtheme = theme(axis.text.x = element_text(angle=0, hjust=0.5, vjust=0.8, size=14, face="plain"),
                                            axis.text.y = element_text(angle=0, hjust=0.5, vjust=0.5, size=14, face="plain"),
                                            axis.title.x = element_text(angle=0, size=14, face="plain"),
                                            legend.title = element_text(angle=0, hjust=0.8, vjust=0.8, size=16, face="plain"),
                                            legend.text = element_text(angle=0, hjust=0.8, vjust=0.8, size=16, face="plain")
                            )+theme_bw()+theme_classic())
    surv_plot$plot <- surv_plot$plot+ theme(legend.text = element_text(size = 16, color = "black", face = "plain"),
                                            axis.title.x = element_text(color = "black", hjust=0.5, vjust=0.8, size=20, face="plain"),
                                            axis.title.y = element_text(angle=90, hjust=0.5, vjust=1, size=14, face="plain"),
                                            legend.title = element_text(angle=0, hjust=0.8, vjust=0.8, size=18, face="plain"))
    surv_plot$table <- surv_plot$table + theme(axis.title.y = element_blank(),
                                               axis.text.x = element_text(color = "black", hjust=0.5, vjust=0.8, size=18, face="plain"))
    remove(SurvivalData)
    return(surv_plot)
}


#### prepare function for rbind ####
srbind <- function(m1, m2){
    sharerow <- intersect(rownames(m1), rownames(m2))
    newmatrix <- cbind(m1[sharerow,], m2[sharerow,])
    return(newmatrix)
}

### prepare function for plot CrossICC result heatmap ###
plot_expression_heatmap_with_cluster <- function(df,sample.cluster, genes,cluster_row=FALSE,showRowname=FALSE){
    plot.matrix <- df
    samplename <- colnames(plot.matrix)
    annotation.list <- sample.cluster[samplename]
    annotation.list <- sort(annotation.list)
    plot.matrix <- plot.matrix[,names(annotation.list)]
    annotation.frame <- data.frame(cluster=as.factor(annotation.list))
    rownames(annotation.frame) <- names(annotation.list)
    #heatmap colors
    colorlength <- 3
    if(length(unique(annotation.frame[,1]))>3){
        colorlength <- length(unique(annotation.frame[,1]))
        color.list<-brewer.pal(colorlength, "Set2")
    }else{
        colorlength <- length(unique(annotation.frame[,1]))
        color.list<-brewer.pal(3, "Set2")[1:colorlength]
    }
    
    names(color.list) <- unique(annotation.frame[,1])
    #plot heatmap
    pheatmap::pheatmap(plot.matrix[genes,],
                       scale = 'none',
                       border_color = NA,
                       cluster_cols = FALSE,
                       cluster_rows = cluster_row,
                       annotation_col = annotation.frame,
                       show_rownames = showRowname,
                       show_colnames = FALSE,
                       colorRampPalette(c("blue", "white", "red"))(100))
    
}



setwd("/data3/liuzk/CrossICCtest/data/")
testfiles <- c("GPL6883_GSE45725.exp.RData",  "GPL6098_GSE22219.exp.RData",
               "GPL13667_GSE65095.exp.RData", "GPL13667_GSE53031.exp.RData",
               "GPL6244_GSE37751.exp.RData",  "GPL570_GSE31448.exp.RData",
               "BRCA_TCGA.RData")

#### get a intergrated expression matrix for consensusclusterPlus ####
# prepare datalist for CrossICC #
# get sample and batch information #


datalist <- vector("list", length = 7)
for(i in 1:length(testfiles)){
    load(testfiles[i])
    GEO_matrix0 <- codingExp
    datalist[[i]] <- as.matrix(codingExp)
    GEO_matrix0_mad <- apply(GEO_matrix0, 1, mad)
    GEO_matrix0 <-  GEO_matrix0[GEO_matrix0_mad > 0.1, ]
    
    if(i == 1){GEO_matrix <- GEO_matrix0
    batchlist <- as.character(rep(i, ncol(GEO_matrix0)))
    patientsID <- colnames(GEO_matrix0)
    }else{
        GEO_matrix <- srbind(GEO_matrix, GEO_matrix0)
        batchlist <- c(batchlist, rep(i, ncol(GEO_matrix0)))
        patientsID <- c(patientsID, colnames(GEO_matrix0))
    }
}

patientsID_batch <- data.frame(samples = patientsID,
                               batch = batchlist,
                               row.names = patientsID,
                               stringsAsFactors = FALSE)


#### Prepare Survival Data ####

brcaGEOTCGA_file7 <- read.delim("/data3/liuzk/CrossICCtest/data/BreastCancer.GEO-TCGA-Survival.txt", sep = "\t", stringsAsFactors = FALSE)
row.names(brcaGEOTCGA_file7) <- brcaGEOTCGA_file7$GSM

###### test not remove batch effect result ####
# GEO_matrix[GEO_matrix == "-Inf"] <- NA
# GEO_matrix[GEO_matrix == "Inf"] <- NA
# GEO_matrix_imputed <- impute.knn(as.matrix(GEO_matrix), k = 10, rowmax = 0.5, colmax = 0.8)
# GEO_matrix_i <- as.data.frame(GEO_matrix_imputed$data)

NotRemoveBatch <- ConsensusClusterPlus(as.matrix(GEO_matrix),
                                       maxK = 10,
                                       reps = 100,
                                       pItem = 0.8,
                                       pFeature = 1,
                                       plot = 'png',
                                       title = "/data3/liuzk/CrossICCtest/plot/NotRemoveBatch",
                                       distance = "pearson",
                                       clusterAlg = "hc")

# When k = 10, show best CDF result, but to banlance the relationship between #
# efficacy and number of each group, we choice the result of k = 6 to show.   #
NotRemoveBatchResultClass <- data.frame(group = NotRemoveBatch[[6]]$consensusClass) %>% rownames_to_column(var = "sample")

### finish survival analysis for not removig batch effect  ####
BRCAsurvival_NRB <- merge(x = NotRemoveBatchResultClass, y = brcaGEOTCGA_file7, by.x = "sample", by.y = "GSM") %>%
    mutate(Cluster = paste("Cluster", group, sep = "")) %>%
    filter(!Cluster %in% c("Cluster2", "Cluster6"))

BRCAsurvival_NRB.Plot <- consensusSurvival(SurvivalData = BRCAsurvival_NRB)
pdf(file = "/data3/liuzk/CrossICCtest/plot/BreastSurvival-NotRemoveBatch.pdf", width = 6.5, height = 6)
print(BRCAsurvival_NRB.Plot, newpage = FALSE)
dev.off()


#### PCA analysis for not remove batch effect ####
notRemoveBatchPCA <- prcomp(t(GEO_matrix), scale. = TRUE)
NRB_matrix <- notRemoveBatchPCA$x
NRBpcaPlotData <- data.frame(NRB_matrix[, c(1:2)],
                             Batch = batchlist,
                             ConsensusGroup = as.character(NotRemoveBatchResultClass$group),
                             stringsAsFactors = FALSE)

NotRemoveBatchPCAplot <- ggplot()+theme_bw()+
                        geom_point(data = NRBpcaPlotData,
                                   aes(x = PC1, y = PC2, color = Batch, shape = ConsensusGroup), size =2)+
                        scale_colour_nejm()+
                        xlab("PC1(76.2% explained var.)")+ylab("PC2(7.1% explained var.)")+
                        theme(axis.text.x = element_text(angle=0, hjust=0.5, vjust=0.8, size=14, face="plain"),
                              axis.text.y = element_text(angle=0, hjust=0.8, vjust=0.5, size=14, face="plain"),
                              axis.title.y = element_text(angle=90, vjust=0.8, size=14, face="plain"),
                              axis.title.x = element_text(angle=0, size=14, face="plain"),
                              legend.title = element_text(angle=0, hjust=0.8, vjust=0.8, size=13, face="plain"),
                              legend.text = element_text(angle=0, hjust=0.8, vjust=0.8, size=12.5, face="plain"),
                              legend.position = c(0.8, 0.6),
                              legend.background = element_blank(),
                              legend.key = element_blank())

ggsave(NotRemoveBatchPCAplot, filename = "/data3/liuzk/CrossICCtest/plot/NotRemoveBatchPCAplot.pdf", width = 7, height = 7)
# ggbiplot(notRemoveBatchPCA, obs.scale = 1, var.scale = 1,
#          groups = as.factor(batchlist), ellipse = TRUE,var.axes = FALSE)


##### test for Removing batch effect ####
batch <- as.factor(batchlist)
RemoveBatchMatrix <- ComBat(as.matrix(GEO_matrix), batch)

RemoveBatch <- ConsensusClusterPlus(as.matrix(RemoveBatchMatrix),
                                    maxK = 10,
                                    reps = 100,
                                    pItem = 0.8,
                                    pFeature = 1,
                                    plot = 'png',
                                    title = "/data3/liuzk/CrossICCtest/plot/RemoveBatch",
                                    distance = "pearson",
                                    clusterAlg="hc")
# When k = 2, show the best CDF value #
RemoveBatchResultClass <- data.frame(group = RemoveBatch[[2]]$consensusClass) %>% rownames_to_column(var = "sample")

### finish survival analysis for  removig batch effect  ####
BRCAsurvival_RB <- merge(x = RemoveBatchResultClass, y = brcaGEOTCGA_file7, by.x = "sample", by.y = "GSM") %>%
    mutate(Cluster = paste("Cluster", group, sep = ""))

BRCAsurvival_RB.Plot <- consensusSurvival(SurvivalData = BRCAsurvival_RB)
pdf(file = "/data3/liuzk/CrossICCtest/plot/BreastSurvival-RemoveBatcPlot.pdf", width = 6.5, height = 6)
print(BRCAsurvival_RB.Plot, newpage = FALSE)
dev.off()


#### PCA analysis for removing batch effect ####
RemoveBatchPCA <- prcomp(t(RemoveBatchMatrix), scale. = TRUE)
RB_matrix <- RemoveBatchPCA$x
RBpcaPlotData <- data.frame(RB_matrix[, c(1:2)],
                            Batch = batchlist,
                            ConsensusGroup = as.character(RemoveBatchResultClass$group),
                            stringsAsFactors = FALSE)

RemoveBatchPCAplot <- ggplot()+
                        theme_bw()+
                        geom_point(data = RBpcaPlotData,
                                   aes(x = PC1, y = PC2, color = Batch, shape = ConsensusGroup), size =2)+
                        scale_colour_nejm()+
                        xlab("PC1(8.3% explained var.)")+ylab("PC2(6.9% explained var.)")+
                        theme(axis.text.x = element_text(angle=0, hjust=0.5, vjust=0.8, size=14, face="plain"),
                              axis.text.y = element_text(angle=0, hjust=0.8, vjust=0.5, size=14, face="plain"),
                              axis.title.y = element_text(angle=90, vjust=0.8, size=14, face="plain"),
                              axis.title.x = element_text(angle=0, size=14, face="plain"),
                              legend.title = element_text(angle=0, hjust=0.8, vjust=0.8, size=13, face="plain"),
                              legend.text = element_text(angle=0, hjust=0.8, vjust=0.8, size=12.5, face="plain"),
                              legend.position = c(0.85, 0.7),
                              legend.background = element_blank(),
                              legend.key = element_blank())

ggsave(RemoveBatchPCAplot, filename = "/data3/liuzk/CrossICCtest/plot/RemoveBatchPCAplot.pdf", width = 7, height = 7)
# ggbiplot(RemoveBatchPCA, obs.scale = 1, var.scale = 1,
#          groups = as.factor(batchlist), ellipse = TRUE,var.axes = FALSE)

### Visualization of Principal component analysis ####
## finish Principal component analysis for merged dataset
importancedata <- summary(RemoveBatchPCA)$importance
importancedata <- t(importancedata)[1:10,]
importancedata <- data.frame(importancedata, comp = rownames(importancedata), Sets = "Total")
importancedata$Proportion.of.Variance = importancedata$Proportion.of.Variance*100

importancePlot <- ggplot(importancedata, aes(x = factor(importancedata$comp, levels = paste("PC", 1:10, sep = "")),
                                             y = Proportion.of.Variance,
                                             group = Sets, color = Sets))+
    theme_classic()+geom_line(show.legend = FALSE)+geom_point(size = 3, show.legend = FALSE)+
    xlab("Components")+ ylab("Proportion of Variance (%)")+
    theme(axis.text = element_text(angle = 0, size = 14, face = "plain"),
          axis.title = element_text(angle = 0, size = 14, face = "plain"))+
    scale_color_jco()

importancePlot
ggsave(importancePlot, filename = "/data3/liuzk/CrossICCtest/plot/importancePlot.pdf", width = 7, height = 3)

## finish Principal component analysis for each dataset
datatestImpTotal <- data.frame()
for(i in 1:7){
    datatest <- datalist[[i]]
    removeGene <- rownames(datatest)[which(apply(datatest, 1, var) == 0)]
    datatest <- datatest[!(rownames(datatest) %in% removeGene), ]
    datatestPCA <- prcomp(t(as.matrix(datatest)), scale. = TRUE)
    datatestImp <- summary(datatestPCA)$importance
    datatestImp <- t(datatestImp)[1:10,]
    datatestImp <- data.frame(datatestImp, comp = rownames(datatestImp), Sets = paste("Dataset", i, sep = ""))
    datatestImpTotal <- rbind(datatestImp, datatestImpTotal)
}

datatestImpTotal$Proportion.of.Variance <- datatestImpTotal$Proportion.of.Variance*100
setlevels <- c("Total", paste("Dataset", 1:10, sep = ""))
All2show <- rbind(importancedata, datatestImpTotal)
datatestImpPlot <- ggplot(All2show, aes(x = factor(All2show$comp, levels = paste("PC", 1:10, sep = "")),
                                        y = Proportion.of.Variance,
                                        group = factor(All2show$Sets, levels = setlevels),
                                        color = factor(All2show$Sets, levels = setlevels)))+
    theme_classic()+geom_line()+geom_point(size = 3)+
    xlab("Components")+ ylab("Proportion of Variance (%)")+labs(color = "Sets")+
    theme(axis.text = element_text(color = "black", angle = 0, size = 14, face = "plain"),
          axis.title = element_text(angle = 0, size = 14, face = "plain"),
          legend.title = element_text(angle = 0, size = 13, face = "plain"),
          legend.text = element_text(angle = 0, size = 13, face = "plain"))+
    scale_color_jco()

datatestImpPlot
ggsave(datatestImpPlot, filename = "/data3/liuzk/CrossICCtest/plot/datatestImpPlot.pdf", width = 7, height = 4)



################################
##### Analysis by CrossICC #####
################################

CrossICCtestRes <- CrossICC(datalist,
                            max.iter = 100, use.shiny = TRUE, cross = "cluster",
                            fdr.cutoff = 0.1, ebayes.cutoff = 0.1,
                            filter.cutoff = 0.1, n.platform = 4,
                            output.dir = "/data3/liuzk/")
CrossICCtestRes <- readRDS("/data3/liuzk/CrossICC.object.rds")
# runShinyCrossICC()

# get cluster result of CrossICC #
CrossICCclusterResult <- data.frame(cluster = CrossICC::summary.CrossICC(CrossICCtestRes)$clusters, stringsAsFactors = FALSE) %>%
                         rownames_to_column(var = "sampleid") %>%
                         dplyr::filter(cluster != "") %>%
                         mutate(Cluster = paste("Cluster", cluster, sep = ""))

CrossICCclusterResult <- merge(x = CrossICCclusterResult, y = brcaGEOTCGA_file7, by.x = "sampleid", by.y = "GSM") %>%
                         dplyr::filter(!Cluster %in% c("Cluster8", "Cluster3", "Cluster9"))


### plot CrossICC result survival curve ###
breast_survivalGEOTCGA_CrossICC_plot <- consensusSurvival(SurvivalData = CrossICCclusterResult)
pdf(file = "/data3/liuzk/CrossICCtest/plot/BreastSurvival-CrossICCplot.pdf")
print(breast_survivalGEOTCGA_CrossICC_plot, newpage = FALSE)
dev.off()

### PCA analysis of CrossICC result ####
CrossICCMatrixTotal <- GEO_matrix[, CrossICCclusterResult$sampleid]
CrossICCpca <- prcomp(t(CrossICCMatrixTotal), scale. = TRUE)

CrossICC.PCA.PlotData <- CrossICCpca$x %>% as.data.frame() %>% .[, c(1:3)]

patientsIDbatchCrossICC <- patientsID_batch[rownames(CrossICC.PCA.PlotData), ]

CrossICC.PCA.PlotData <- mutate(CrossICC.PCA.PlotData,
                                Batch = as.character(patientsIDbatchCrossICC$batch),
                                CrossICC.Cluster = CrossICCclusterResult$Cluster)

CrossICC_pca_plot <- ggplot()+theme_bw()+
    geom_point(data = CrossICC.PCA.PlotData,
               aes(x = PC1, y = PC2, color = Batch, shape = CrossICC.Cluster), size =2)+
    scale_colour_nejm()+
    xlab("PC1(76.2% explained var.)")+ylab("PC2(7.1% explained var.)")+
    theme(axis.text.x = element_text(angle=0, hjust=0.5, vjust=0.8, size=14, face="plain"),
          axis.text.y = element_text(angle=0, hjust=0.8, vjust=0.5, size=14, face="plain"),
          axis.title.y = element_text(angle=90, vjust=0.8, size=14, face="plain"),
          axis.title.x = element_text(angle=0, size=14, face="plain"),
          legend.title = element_text(angle=0, hjust=0.8, vjust=0.8, size=13, face="plain"),
          legend.text = element_text(angle=0, hjust=0.8, vjust=0.8, size=12.5, face="plain"),
          legend.position = c(0.8, 0.7),
          legend.background = element_blank(),
          legend.key = element_blank())

CrossICC_pca_plot
ggsave(CrossICC_pca_plot, filename = "/data3/liuzk/CrossICCtest/plot/CrossICC-PCAplot.pdf", width = 7, height = 7)

# ggbiplot(CrossICCpca, obs.scale = 1, var.scale = 1,
#       groups = as.factor(CrossICCclusterResult$cluster), ellipse = TRUE,var.axes = FALSE)


#### plot heatmap for each matrix by signature genes ####
for(k in 1:7){
    plot.matrix<-as.data.frame(CrossICCtestRes$platforms[[k]])
    if(class(CrossICCtestRes$clusters$clusters)=="list"){
        cluster.table<-CrossICCtestRes$clusters$clusters[[k]]
    }else{
        cluster.table<-CrossICCtestRes$clusters$clusters
    }
    gsig<-CrossICCtestRes$gene.order[[k]]
    pdf(paste("/data3/liuzk/CrossICCtest/plot/CrossICC-Dataset", k, "-heatmap.pdf", sep = ""))
    plot_expression_heatmap_with_cluster(plot.matrix,cluster.table,gsig, cluster_row = FALSE, showRowname = FALSE)
    dev.off()
}









