#Split Toxoplasma genes into groups depending on their expression level

#Load data (downloaded from ToxoDB): Transcriptome RNA-seq data Lorenzi 2011
exprs<-read.table("~/Documents/Kim_lab/GSEA/Expression level/Percentile_RNAseq_data.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
exprssplit<-split(exprs, cut(exprs$Expression_percentile, c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100), include.lowest=TRUE))
lapply(names(exprssplit), function(x) {
  write.table(as.data.frame(exprssplit[[x]]$Gene_ID), file = paste(x, ".txt", collapse = "\t"),
              append = FALSE, row.names = FALSE, quote=FALSE, col.names=FALSE)
})

#Load data (downloaded from ToxoDB): Microarray Roos
exprs2<-read.table("Documents/Kim_lab/Ubiquitination project/Pilot/GSEA/V10/Percentile_tachy_microarray_Roos.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
exprssplit2<-split(exprs2, cut(exprs2$Expression_percentile, c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100), include.lowest=TRUE))
lapply(names(exprssplit2), function(x) {
  write.table(as.data.frame(exprssplit2[[x]]$Gene_ID), file = paste(x, "_2.txt", collapse = "\t"),
              append = FALSE, row.names = FALSE, quote=FALSE, col.names=FALSE)
})