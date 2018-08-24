library(ConsensusClusterPlus)
library(survival)
library(survminer)
library(ggsci)
library(sva)
library(CrossICC)
library(RColorBrewer)
# library(ggbiplot)

#### prepare function for ploting survival curve ####
consensusSurvival = function(SurvivalData){
  SurvivalData <<- SurvivalData
  gene_surv <- Surv(SurvivalData$Time.days./30, SurvivalData$Event)~SurvivalData$Cluster
  gene_survfit <- survfit(gene_surv)
  gene_survdiff <- tryCatch(survdiff(gene_surv), error = function(e) return(NA))
  gene_p <- pchisq(gene_survdiff$chisq, length(gene_survdiff$n) - 1, lower.tail = FALSE)

  plotcorlor = pal_d3(alpha = 1)(10)
  fit<- survfit(Surv(SurvivalData$Time.days./30, SurvivalData$Event)~SurvivalData$Cluster, data = SurvivalData)
  surv_plot = ggsurvplot(fit, main = "READ", xlab = "Time(Months)", #ylab = "Progression-free survival (%)",
                         conf.int = F, pval = TRUE, legend = c(0.8, 0.75), #scale_y_continuous(expand = c(0,0)),
                         legend.title = "Cluster of Breast Cancer",
                         risk.table.fontsize = 5,
                         risk.table = TRUE, risk.table.y.text.col = TRUE, pval.size = 6,
                         font.main = c(16, "plain", "black"),#
                         font.x = c(14, "plain", "black"),font.y = c(14, "plain", "black"),
                         font.tickslab = c(14, "plain", "black"), palette = plotcorlor,
                         ggtheme = theme(axis.text.x = element_text(angle=0, hjust=0.5, vjust=0.8, size=14, face="plain"),
                                         axis.text.y = element_text(angle=0, hjust=0.8, vjust=0.5, size=14, face="plain"),
                                         axis.title.y = element_text(angle=90, hjust=0.8, vjust=0.8, size=14, face="plain"),
                                         axis.title.x = element_text(angle=0, size=14, face="plain"),
                                         legend.title = element_text(angle=0, hjust=0.8, vjust=0.8, size=13, face="plain"),
                                         legend.text = element_text(angle=0, hjust=0.8, vjust=0.8, size=12.5, face="plain")
                         )+theme_bw()+theme_classic())
  return(surv_plot)
}


#### prepare function for rbind ####
srbind = function(m1, m2){
  sharerow = intersect(rownames(m1), rownames(m2))
  newmatrix = cbind(m1[sharerow,], m2[sharerow,])
  return(newmatrix)
}

### prepare function for plot CrossICC result heatmap ###
plot_expression_heatmap_with_cluster<-function(df,sample.cluster, genes,cluster_row=FALSE,showRowname=FALSE){
  plot.matrix<-df
  samplename<-colnames(plot.matrix)

  annotation.list<-sample.cluster[samplename]

  annotation.list<-sort(annotation.list)
  plot.matrix<-plot.matrix[,names(annotation.list)]
  annotation.frame<-data.frame(cluster=as.factor(annotation.list))
  rownames(annotation.frame)<-names(annotation.list)
  #heatmap colors
  colorlength <- 3
  if(length(unique(annotation.frame[,1]))>3){
    colorlength <- length(unique(annotation.frame[,1]))
    color.list<-brewer.pal(colorlength, "Set2")
  }else{
    colorlength <- length(unique(annotation.frame[,1]))
    color.list<-brewer.pal(3, "Set2")[1:colorlength]
  }

  names(color.list)<-unique(annotation.frame[,1])
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
testfiles = c(
"GPL6883_GSE45725.exp.RData",  "GPL6098_GSE22219.exp.RData",  "GPL13667_GSE65095.exp.RData", "GPL13667_GSE53031.exp.RData",
"GPL6244_GSE37751.exp.RData",  "GPL570_GSE31448.exp.RData",   "BRCA_TCGA.RData")

#### get a intergrated expression matrix for consensusclusterPlus ####
# prepare datalist for CrossICC #
# get sample and batch information #


datalist = vector("list", length = 7)
for(i in 1:length(testfiles)){
  load(testfiles[i])
  GEO_matrix0 = codingExp
  datalist[[i]] =codingExp
  GEO_matrix0_mad = apply(GEO_matrix0, 1, mad)
  GEO_matrix0 = GEO_matrix0[which(as.numeric(GEO_matrix0_mad) > 0.1),]
  if(i == 1){GEO_matrix = GEO_matrix0
            batchlist = as.character(rep(i, ncol(GEO_matrix0)))
            patientsID = colnames(GEO_matrix0)
  }else{
    GEO_matrix = srbind(GEO_matrix, GEO_matrix0)
    batchlist = c(batchlist, rep(i, ncol(GEO_matrix0)))
    patientsID = c(patientsID, colnames(GEO_matrix0))
  }
}

patientsID_batch = data.frame(samples = patientsID, batch = batchlist, stringsAsFactors = F)
rownames(patientsID_batch) = patientsID_batch$samples


#### Prepare Survival Data ####

brcaGEOTCGA_file7 = read.delim("/data3/liuzk/CrossICCtest/data/BreastCancer.GEO-TCGA-Survival.txt", sep = "\t", stringsAsFactors = F)
row.names(brcaGEOTCGA_file7) = brcaGEOTCGA_file7$GSM

###### test not remove batch effect result ####
# GEO_matrix[GEO_matrix == "-Inf"] = NA
# GEO_matrix[GEO_matrix == "Inf"] = NA
# GEO_matrix_imputed = impute.knn(as.matrix(GEO_matrix), k = 10, rowmax = 0.5, colmax = 0.8)
# GEO_matrix_i = as.data.frame(GEO_matrix_imputed$data)

NotRemoveBatch = ConsensusClusterPlus(as.matrix(GEO_matrix),
                                      maxK=10,
                                      reps=100,
                                      pItem=0.8,
                                      pFeature=1,
                                      plot='png',
                                      title="/data3/liuzk/CrossICCtest/plot/NotRemoveBatch",
                                      distance="pearson",
                                      clusterAlg="hc")
# When k = 10, show best CDF result, but to banlance the relationship between #
# efficacy and number of each group, we choice the result of k = 6 to show.   #
NotRemoveBatchResultClass = data.frame(group = NotRemoveBatch[[6]]$consensusClass)

### finish survival analysis for not removig batch effect  ####
breast_survivalGEOTCGA_NRB = brcaGEOTCGA_file7[intersect(brcaGEOTCGA_file7$GSM, rownames(NotRemoveBatchResultClass)), ]
NotRemoveBatch_4survival = NotRemoveBatchResultClass[rownames(breast_survivalGEOTCGA_NRB), , drop = FALSE]
breast_survivalGEOTCGA_NRB$Cluster = paste("Cluster", NotRemoveBatch_4survival$group, sep = "")

BreastSurvivalNotRemoveBatcPlot = consensusSurvival(SurvivalData = breast_survivalGEOTCGA_NRB)
pdf(file = "/data3/liuzk/CrossICCtest/plot/BreastSurvival-NotRemoveBatcPlot.pdf")
print(BreastSurvivalNotRemoveBatcPlot, newpage = FALSE)
dev.off()


#### PCA analysis for not remove batch effect ####
notRemoveBatchPCA = prcomp(t(GEO_matrix), scale. = TRUE)
NRB_matrix = notRemoveBatchPCA$x
NRBpcaPlotData = data.frame(NRB_matrix[, c(1:2)],
                               Batch = batchlist,
                               ConsensusGroup = as.character(NotRemoveBatchResultClass$group),
                               stringsAsFactors = F)

NotRemoveBatchPCAplot = ggplot()+theme_bw()+
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
#          groups = as.factor(batchlist), ellipse = TRUE,var.axes = F)


##### test for Removing batch effect ####
batch=as.factor(batchlist)
RemoveBatchMatrix = ComBat(as.matrix(GEO_matrix), batch)

RemoveBatch = ConsensusClusterPlus(as.matrix(RemoveBatchMatrix),
                                    maxK=10,
                                    reps=100,
                                    pItem=0.8,
                                    pFeature=1,
                                    plot='png',
                                    title="/data3/liuzk/CrossICCtest/plot/RemoveBatch",
                                    distance="pearson",
                                    clusterAlg="hc")
# When k = 2, show the best CDF value #
RemoveBatchResultClass = data.frame(group = RemoveBatch[[2]]$consensusClass)

### finish survival analysis for  removig batch effect  ####
breast_survivalGEOTCGA_RB = brcaGEOTCGA_file7[intersect(brcaGEOTCGA_file7$GSM, rownames(RemoveBatchResultClass)), ]
RemoveBatchResult_surv = RemoveBatchResultClass[rownames(breast_survivalGEOTCGA_RB), , drop = FALSE]
breast_survivalGEOTCGA_RB$Cluster = paste("Cluster", RemoveBatchResult_surv$group, sep = "")

breast_survivalGEOTCGA_RB_plot = consensusSurvival(SurvivalData = breast_survivalGEOTCGA_RB)
pdf(file = "/data3/liuzk/CrossICCtest/plot/BreastSurvival-RemoveBatcPlot.pdf")
print(breast_survivalGEOTCGA_RB_plot, newpage = FALSE)
dev.off()

#### PCA analysis for removing batch effect ####
RemoveBatch_pca = prcomp(t(RemoveBatchMatrix), scale. = TRUE)
RB_matrix =  RemoveBatch_pca$x
RBpcaPlotData = data.frame(RB_matrix[, c(1:2)],
                                 Batch = batchlist,
                                 ConsensusGroup = as.character(RemoveBatchResultClass$group),
                                 stringsAsFactors = F)

RemoveBatchPCAplot = ggplot()+theme_bw()+
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
# ggbiplot(RemoveBatch_pca, obs.scale = 1, var.scale = 1,
#          groups = as.factor(batchlist), ellipse = TRUE,var.axes = F)


##### Analysis by CrossICC ####
CrossICCtestRes <- CrossICC(datalist, max.iter = 100, use.shiny = TRUE, cross = "cluster",
                                fdr.cutoff = 0.1, ebayes.cutoff = 0.1, filter.cutoff = 0.1)
CrossICCtestRes <- readRDS("~/CrossICC.object.rds")
# runShinyCrossICC()

# get cluster result of CrossICC #
CrossICCclusterResult = data.frame(cluster = CrossICC::summary.CrossICC(CrossICCtestRes)$clusters, stringsAsFactors = F)
CrossICCclusterResult = CrossICCclusterResult[CrossICCclusterResult$cluster != "", , drop = FALSE]
CrossICCclusterResult$cluster = paste("Cluster", CrossICCclusterResult$cluster, sep = "")

### plot CrossICC result survival curve ###
breastSurvivalCrossICC = brcaGEOTCGA_file7[rownames(CrossICCclusterResult), ]
CrossICCclusterResult = CrossICCclusterResult[rownames(breastSurvivalCrossICC), , drop = FALSE]
breastSurvivalCrossICC$Cluster = CrossICCclusterResult$cluster

breast_survivalGEOTCGA_CrossICC_plot = consensusSurvival(SurvivalData = breastSurvivalCrossICC)
pdf(file = "/data3/liuzk/CrossICCtest/plot/BreastSurvival-CrossICCplot.pdf")
print(breast_survivalGEOTCGA_CrossICC_plot, newpage = FALSE)
dev.off()

### PCA analysis of CrossICC result ####
CrossICCMatrixTotal = GEO_matrix[, rownames(CrossICCclusterResult)]
CrossICCpca = prcomp(t(CrossICCMatrixTotal), scale. = TRUE)

CrossICCMatrix =  CrossICCpca$x
CrossICC.PCA.PlotData = as.data.frame(CrossICCMatrix[, c(1:3)])
patientsIDbatchCrossICC = patientsID_batch[rownames(CrossICC.PCA.PlotData), ]

CrossICC.PCA.PlotData$Batch = as.character(patientsIDbatchCrossICC$batch)
CrossICC.PCA.PlotData$CrossICC.Cluster = CrossICCclusterResult$cluster

CrossICC_pca_plot = ggplot()+theme_bw()+
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

ggsave(CrossICC_pca_plot, filename = "/data3/liuzk/CrossICCtest/plot/CrossICC-PCAplot.pdf", width = 7, height = 7)

# ggbiplot(CrossICCpca, obs.scale = 1, var.scale = 1,
#       groups = as.factor(CrossICCclusterResult$cluster), ellipse = TRUE,var.axes = F)


#### plot heatmap for each matrix by signature genes ####
for(k in 1:length(unique(CrossICCclusterResult$cluster))){
  plot.matrix<-as.data.frame(CrossICCtestRes$platforms[[k]])
  if(class(CrossICCtestRes$clusters$clusters)=="list"){
    cluster.table<-CrossICCtestRes$clusters$clusters[[k]]
  }else{
    cluster.table<-CrossICCtestRes$clusters$clusters
  }
  gsig<-CrossICCtestRes$gene.order[[k]]
  pdf(paste("/data3/liuzk/CrossICCtest/plot/CrossICC-Cluster", k, "-heatmap.pdf", sep = ""))
  plot_expression_heatmap_with_cluster(plot.matrix,cluster.table,gsig,cluster_row = FALSE,showRowname = FALSE)
  dev.off()
}


