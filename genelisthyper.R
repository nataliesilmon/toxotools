#R Script for Gene Set Enrichment Testing				17 March 2016
#Natalie Silmon de Monerri - Kim Lab

#Read in genome and specific gene list from experiment
#Files must be: tab delimited lists of gene accession numbers. 
#Gene lists and gene sets should all be UNIQUE

#BEFORE RUNNING SCRIPT
#read in genes from genome transpose to get a vector
#genomet<-c(t(read.table("TgondiiV26_allgenes.txt", sep="\t", stringsAsFactors=FALSE, header=FALSE)))

#genelist<-read.table("Documents/Kim_lab/Ubiquitination project/Pilot/GSEA/V10/PTM_genelists/Ubi_All_ku80.txt", sep="\t", stringsAsFactors=FALSE, header=FALSE)
#transpose
#genelistt<-c(t(genelist))

#SCRIPT
#source mydhyper function: runs a hypergeometric statistical test comparing observed overlap between a list of genes (genelist) and a predefined gene list (geneset - see GSEA paper by Croken et al 2014 PLoS One for gene sets example), with the background number of genes. Output: P-value of enrichment.
mydhyper<-function(genelistt, genesett, genomet){
  myX<-length(which(genelistt %in% genesett))
  myM<-length(genesett)
  myN <- length(genomet) - length(genesett)
  myK <- length(genelistt)
  return(dhyper(x=myX,m=myM, n=myN, k=myK))
}

#source hyperRandom function: runs a hypergeometric statistical test which randomly samples all genes in the TGME49 genome for a set of genes the same length as the gene SET being tested. This is used to test for the likelihood of random enrichment. Output: p-value of enrichment.
dhyperRandom <- function(genelistt, genesett, genomet){
  myRandomGS <- sample(genomet,size=length(genesett) )
  myX <- length(which(genelistt %in% myRandomGS))
  myM <- length(myRandomGS)
  myN <- length(genomet) - length(myM)
  myK <- length(genelistt)
  return(dhyper(x=myX, m=myM, n=myN, k=myK))
}

#Source genelisthyper function: automates the hypergeometric testing to test one list of genes against ALL predefined gene sets in a defined folder (directory) and incorporates the random enrichment testing - sampling 1000 times. The output is a dataframe containing: P-values of enrichment, Adjusted P-values (adjusted for multiple hypothesis testing by Bonferroni method), Length of Gene Set, Overlap between gene set and gene list, and fold enrichment.
#Inputs: gene list to be tested, directory of gene sets, name of output file
genelisthyper<-function(genelist, directory, output, outputPDF, plottitle){
  #genelist as a tab delimited list; directory as a path to a folder, output as a path to a file with the desired name
  #library(stringr)
  output<-c(output)
  direc<-list.files(path=directory, full.names=TRUE)
  id<-1:length(direc)
  pvalue<-c()
  hypervalues<-c()
  randomvalues<-c()
  lengths<-c()
  overlap<-c()
  foldenrich<-c()
  foldenrichrandom<-c()
  names<-c()
  meanenrich<-c()
  print("Examining datasets…hang on…")
  for(i in id){
    geneset<-read.table(direc[i], sep="\t", stringsAsFactors=FALSE, header=FALSE)
    genesett<-c(t(geneset))
    hypervalues<-rbind(hypervalues, mydhyper(genelistt, genesett, genomet))
    lengths<-rbind(lengths, length(genesett))
    overlap<-rbind(overlap, length(which(genelistt %in% genesett)))
    #names<-rbind(names, str_c(genelistt[which(genelistt %in% genesett)], collapse=","))
    names<-rbind(names, paste(intersect(genesett, genelistt), collapse=","))
    foldenrich<-rbind(foldenrich, ((length(which(genelistt %in% genesett))/length(genesett))/(length(genesett)/length(genomet))))
    for(i in 1:1000){
      pvalue[i] <- dhyperRandom(genelistt, genesett, genomet)
    }
    randomvalues<-c(randomvalues, mean(pvalue))
    for(i in 1:1000){
      myRandomGS<-sample(genomet,size=length(genesett) )
      foldenrichrandom<-rbind(foldenrichrandom, ((length(which(myRandomGS %in% genesett))/length(genesett))/(length(genesett)/length(genomet))))
    }
    meanenrich<-rbind(meanenrich, mean(foldenrichrandom))
    
  }
  adjusthyper<-p.adjust(hypervalues, method="bonferroni")
  randomvalues<-p.adjust(randomvalues, method="bonferroni")	
  #all_genes<-cbind(hypervalues, adjusthyper, randomvalues, lengths, overlap, foldenrich, meanenrich, names)
  all_genes<-data.frame(hypervalues, adjusthyper, randomvalues, lengths, overlap, foldenrich, meanenrich, names)
  rownames(all_genes)<-list.files(directory)
  rownames(all_genes)<-gsub(".txt", "", rownames(all_genes))
  colnames(all_genes)<-c("Enrichment_Pvalue", "Bonferroni_Pvalue", "Random_Pvalue", "Gene_Set_Length", "Overlap", "Fold_Enrichment", "Random_Fold_Enrichment", "Genes")
  all_genes<-as.data.frame(all_genes)
  print("Done! Check Output Files.")
  #write output to file, with the name specified at the beginning of the function.
  write.table(all_genes, file=output, sep="\t", row.names=TRUE, quote=FALSE)

#export as a pdf
barplotvalues<-(-log((all_genes$Bonferroni_Pvalue/all_genes$Random_Pvalue),2))
barplotvalues[abs(barplotvalues)==Inf]<-0
pdf(outputPDF)	
	barplot(barplotvalues, names=rownames(all_genes), las=2, col="deepskyblue", main=plottitle, ylab="-log(adjusted pvalue)", xlab="Gene set", cex.names=0.5)
	abline(h=4.32, col="green", lty=2)
	legend("topright", col="green", legend="pvalue=0.05", lty=2, bty = "n") 
dev.off()
x<-list(barplotvalues, all_genes)
names(x)<-c("barplotvalues", "gseaoutput")
return(x)
}

#Example usage:
#source("genelisthyper.R")
#genelisthyper(genelistt, "GeneSetDirectory", "Analysisoutputname", "Plotoutputname", "Plottitle")
