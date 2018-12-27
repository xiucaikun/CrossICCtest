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
