# toxotools
Tools for analysis of *Toxoplasma gondii* experimental datasets

Development: Kim Laboratory, Albert Einstein College of Medicine, Bronx, NY, USA


*Toxoplasma gondii* is a protozoan parasite. Since it is not a widely used model organism, many tools often used for e.g. human or mouse datasets can not easily be applied to this organism and custom analysis methods need to be developed. toxotools is a selection of R scripts designed to aid in the analysis of *Toxoplasma gondii* proteomic, next generation sequencing or array data. 

1. Next generation sequencing data analysis

**getanno.R**

**getanno.R** loads ToxoDB genome annotation data for *Toxoplasma gondii* and converts this to a RangedData object for downstream analysis with e.g. ChIPpeakAnno. User can specify which genome version to be used.

**moreanno.R**

**moreanno.R** loads genome annotations (that are not available on ToxoDB) into RangedData objects for downstream analysis with e.g. ChIPpeakAnno. Annotations are in version 26.


- Gene lists

 **genelists.R**

**genelists.R** loads gene sets defined in Croken et al (2014) BMC Genomics as well as novel gene sets into R for downstream analysis. Currently available gene sets are:

Cell cycle G1, Cell cycle SM, Stage, Localisation, GO Terms, KEGG pathways, LAMP pathways, Posttranslational modification proteomes

- Functional gene enrichment analysis

**genelisthyper.R** 

**genelisthyper.R** performs enrichment analysis against predefined gene sets (specified by **genelists.R**)


- Conversion of versions
ChIP-seq data
ChIP-chip data

