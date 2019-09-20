library(cummeRbund)
cuff <- readCufflinks("/dev/datasets/FairWind/_results/Fatima/cuffdiff")
gene_diff_data <- diffData(genes(cuff))
typeof(gene_diff_data)

write.table(as.data.frame(gene_diff_data),file="/dev/datasets/FairWind/_results/Fatima/genes_diff.csv", quote=F,sep="\t",row.names=F)


#capture.output(summary(gene_diff_data), file = "/dev/datasets/FairWind/genes2.csv")

#lapply(gene_diff_data, function(x) write.table( data.frame(x), '/dev/datasets/FairWind/genes.csv'  , append= T, sep=',' ))



head(gene_diff_data)
diffGeneIDs <- getSig(cuff,level="genes",alpha=0.05)
diffGenes<-getGenes(cuff,diffGeneIDs)
a=sort(diffGenes)

