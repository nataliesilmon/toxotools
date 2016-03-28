#Script for converting GMT files into individual text files - individual gene lists

#Read in gmt file
gmt<-data.frame(t(read.table("Tgondii_GO_Terms.gmt", sep="\t", stringsAsFactors=FALSE, row.names=1, na.strings =c("	na", ""))))

#Loop through gmt file - remove NA values 
for (i in 1:ncol(gmt)){
	a<-data.frame(gmt[,i])
	a<-a[!is.na(a)]
	myfile<-file.path("GO_terms", paste0(colnames(gmt)[i], ".txt"))
	write.table(a, file=myfile, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, append=FALSE)
}


KEGG<-data.frame(t(read.table("Tgondii_kegg_pathways.curated.gmt", sep="\t", stringsAsFactors=FALSE, row.names=1, na.strings =c("	na", ""))))

#export as text files with column name as file name
sapply(gmt, write.table, sep="\t". file=colnames

for (i in 1:ncol(KEGG)){
	a<-data.frame(KEGG[,i])
	a<-a[!is.na(a)]
	myfile<-file.path("KEGG", paste0(colnames(KEGG)[i], ".txt"))
	write.table(a, file=myfile, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, append=FALSE)
}


#import LAMP pathways
LAMP<-read.table("Toxoplasma_gondii_llamp_annotation_release2.txt", sep="\t", stringsAsFactors=FALSE, header=TRUE)
LAMPS<-split(LAMP, f=LAMP$Metabolic.pathway)
LAMPIDs<-sapply(LAMPS, function(x) x$T..gondii.gene.id)
for (i in 1:38){
	a<-data.frame(LAMPIDs[i])
	myfile<-file.path("LAMP", paste0(names(a), ".txt"))
	write.table(a, file=myfile, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, append=FALSE)
}

LAMPreact<-read.table("Toxoplasma_gondii_llamp_reactions_list_release2.txt", sep="\t", stringsAsFactors=FALSE, header=TRUE)
LAMPreact<-split(LAMPreact, f=LAMPreact$Pathway.name)
LAMPIDs<-sapply(LAMPreact, function(x) x$Gene.Ids)
for (i in 1:51){
	a<-data.frame(LAMPIDs[i])
	myfile<-file.path("LAMP_react", paste0(names(a), ".txt"))
	write.table(a, file=myfile, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, append=FALSE)
}

