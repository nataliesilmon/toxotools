library(ChIPpeakAnno)
library(rtracklayer)
gff<-read.table("~/Documents/Kim_lab/Computational analysis/ChIP-seq/Annotation files for ChIPpeakAnno/ToxoDB-26_TgondiiME49_annotation.gff.txt", sep="\t", stringsAsFactors=FALSE, header=FALSE)
message("file 1 loaded")
colnames(gff)<-c("chromosome", "source", "feature", "start", "end", "score", "strand", "frame", "group")
#Use only genes - not mRNA, 5' or 3' UTR or exons as well
gff<-gff[gff$feature=="gene", ]
#Replace values in strand, score and frame column with Ranged Data "friendly" values - see GFF2RangedData usage
gff$strand[gff$strand=="+"]<-1
gff$strand[gff$strand=="-"]<--1
gff$score[gff$score=="."]<-NA
gff$frame[gff$frame=="."]<-NA
gffRD<-GFF2RangedData(gff)


#GFF2RangedData messes up file - changes order, doesn't carry over IDs.
#So need to get back IDs. By ordering the RD object and original file to then be able to extract IDs in correct order
#Reordering objects based on chromosome and start position
gffRDo<-gffRD[order(space(gffRD),start(gffRD)),]
gffso<-gff[order(gff$chromosome, gff$start),]
#Extract IDs from column 9 from annotation file
IDs<-strsplit(gffso$group, split=";")
IDs<-sapply(IDs, function(X){X[1]}) #extracts 1st part of the element of each vector of #the list
IDs<-strsplit(IDs, split="=")
IDs<-sapply(IDs, function(X){X[2]})
#head(IDs)
#Make IDs row names
rownames(gffRDo)<-IDs #If crashes/has difficulty, close X11 and open browsers!!
#head(gffRDo)
message("gene RD object created")


#TSS annotation
TSS<-read.table("~/Documents/Kim_lab/Computational analysis/ChIP-seq/Annotation files for ChIPpeakAnno/TgondiiTSS_V11.bed" , sep="\t", stringsAsFactors=FALSE)
message("file 2 loaded")

TSS$V1<-gsub("TGME49_chr", "TGME49_", TSS$V1)
TSSRD<-RangedData(ranges=IRanges(start=TSS[,2], end=TSS[,3]), space=TSS[,1])
#for some reason these files have duplicated rows - remove
TSSRD<-TSSRD[!duplicated(TSSRD),]
message("TSS RD object created")


#OR
TSS2<-read.table("tss_out.txt" , sep="\t", stringsAsFactors=FALSE)
message("file 3 loaded")

TSS2$V2<-gsub("TGME49_chr", "TGME49_", TSS2$V2)
TSS2<-TSS2[!is.na(TSS2$V6),]
TSSRD2<-RangedData(ranges=IRanges(start=TSS2[,8], end=TSS2[,7]), space=TSS2[,2])
#for some reason these files have duplicated rows - remove
TSSRD2<-TSSRD2[!duplicated(TSSRD2),]
message("TSS RD object 2 created")


###Annotation by promoter
promo<-read.table("~/Documents/Kim_lab/Computational analysis/ChIP-seq/Annotation files for ChIPpeakAnno/ToxDB-26.0_TgondiiME49-promo.gff", sep="\t", stringsAsFactors=FALSE)
message("file 4 loaded")
colnames(promo)<-c("chromosome", "source", "feature", "start", "end", "score", "strand", "frame", "group")
#For some reason promoter file contains duplicates - remove
promo<-promo[!duplicated(promo$group),]
#Replace values in strand, score and frame column with Ranged Data "friendly" values - see GFF2RangedData usage
#Move gene to feature type
promo$strand[promo$strand=="+"]<-1
promo$strand[promo$strand=="-"]<--1
promo$score[promo$score=="."]<-NA
promo$frame[promo$frame=="."]<-NA
promoRD<-GFF2RangedData(promo)
#loses Gene IDs!
promoRDo<-promoRD[order(space(promoRD),start(promoRD)),]
promoso<-promo[order(promo$chromosome, promo$start),]
#Extract IDs from column 9 from annotation file
IDs<-promoso$group

#Make IDs row names
rownames(promoRDo)<-IDs #If crashes/has difficulty, close X11 and open browsers!!

message("Promoter RD object created")
