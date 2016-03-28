# toxotools
Tools for analysis of Toxoplasma gondii experimental datasets

Toxotools is a selection of R scripts for interrogation of a variety of Toxoplasma gondii datasets. 

1. Next generation sequencing data analysis
getanno.R
getanno.R is an R script that loads ToxoDB genome annotation data for Toxoplasma gondii and converts this to a RangedData object for analysis for downstream analysis with e.g. ChIPpeakAnno. User can specify which genome version to be used.

moreanno.R
moreanno.R is an R script that loads genome annotation data 


2. Gene lists
genelists.R
genelists.R is an R script that loads gene lists defined by Croken et al (2014) BMC Genomics as well as novel gene lists into R for downstream analysis. Currently available gene lists are:
- Cell cycle G1 
- Cell cycle SM
- Stage
- Localisation
- GO Terms
- KEGG pathways
- LAMP pathways
- Posttranslational modification proteomes

2. Functional gene enrichment analysis
genelisthyper.R 
genelisthyper.R performs enrichment analysis using user-specified gene lists. 

