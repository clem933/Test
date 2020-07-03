rm(list = ls())
install.packages("Seurat")
if(F){
BiocManager::install("DESeq2")
BiocManager::install("edgeR")
}
rm(list=ls())
options(stringsAsFactors = F)
b<-read.csv(file='D:csv',header =T,row.names = 1 )
                    
exprSet = b
suppressMessages(library(DESeq2))
group_list=c('young','young','young','young','young','young','young','ageing','ageing','ageing','ageing','ageing','ageing','ageing')
(colData <- data.frame(row.names = colnames(exprSet),
                       group_list=group_list))
dds <- DESeqDataSetFromMatrix(countData = exprSet,
                            colData = colData,
                            design = ~ group_list)
dds<-DESeq(dds)
res <- results(dds,
               contrast = c("group_list","ageing","young"))
resOrdered <- res[order(res$padj),]
head(resOrdered)
DEG=as.data.frame(resOrdered)
DEG1 = na.omit(DEG)

nrDEG = DEG
library(pheatmap)
choose_gene = head(rownames(nrDEG),50)
choose_matrix = exprSet[choose_gene,]
choose_matrix = t(scale(t(choose_matrix)))
pheatmap(choose_matrix,fontsize_col=8, fontsize_row=8, cluster_cols=F,filename = 'DEG_heatmap.png')

logFC_cutoff<-with(DEG,mean(abs(log2FoldChange))+2*sd(abs(log2FoldChange)))
DEG$change = as.factor(ifelse(DEG$pvalue < 0.05 & abs(DEG$log2FoldChange)>logFC_cutoff,
                              ifelse(DEG$log2FoldChange>logFC_cutoff,'UP','DOWN'),'NOT'))
this_tile <- paste0('Cutoff for logFC is',round(logFC_cutoff,3),
                   '\nThe number of up gene is', nrow(DEG[DEG$change =='up',]),
                   '\nThe number of down gene is', nrow(DEG[DEG$change =='down',]))
library(ggplot2)
g = ggplot(data=DEG,
           aes(x=log2FoldChange, y=-log10(pvalue),
               color=change))+
        geom_point(alpha=0.4, size=1.75)+
        theme_set(theme_set(theme_bw(base_size=20)))+
        xlab("log2 fold change")+ylab("-log10 p-value")+
        ggtitle( this_tile )+theme(plot.title=element_text(size=15,hjust=0.5))+
        scale_colour_manual(values=c('blue','black','red'))
print(g)
ggsave(g,filename = 'volcano_.png')

write.csv(data.frame(DEG1),"D:/Tools/R/Rdata/DEG1.csv")

png("qc_dispersions.png",1000,1000,pointsize = 20)
 plotDispEsts(dds, main="Dispersion plot")
 dev.off()
 
