# toxotools
**Tools for analysis of *Toxoplasma gondii* experimental datasets**

Development: Natalie Silmon de Monerri, Kim Laboratory, Department of Medicine, Albert Einstein College of Medicine, Bronx, NY, USA

*Toxoplasma gondii* is a protozoan parasite. Since it is a non-model organism, many tools often used for e.g. human or mouse datasets can not easily be applied to this organism and custom methods for annotation and analysis  are required. toxotools is a selection of R scripts and annotation files designed to aid in the analysis of *Toxoplasma gondii* proteomic, next generation sequencing or array data. 

This is an ongoing project and is being expanded to include *Plasmodium falciparum* scripts.

1. Next generation sequencing data analysis

**getanno.R**
Loads ToxoDB genome annotation data for *Toxoplasma gondii* and converts this to a RangedData object for downstream analysis with e.g. ChIPpeakAnno. User can specify which genome version to be used.

**moreanno.R**
Loads genome annotations (that are not available on ToxoDB) into RangedData objects for downstream analysis with e.g. ChIPpeakAnno. Annotations are in version 26.

2. Gene lists for download and loading into R.

**genelists.R**
Loads Toxoplasma gondii gene sets defined in Croken *et al* (2014) BMC Genomics as well as novel gene sets into R for downstream analysis. 
Currently available gene sets are located in the file 'GENE SETS 2016 TXT.zip' and are in .txt format. **GMT files for use in GSEA analysis will be available soon.
Cell cycle G1, Cell cycle SM, Stage, Localisation, GO Terms, KEGG pathways, LAMP pathways, Posttranslational modification proteomes

3. Functional gene enrichment analysis
**genelisthyper.R** 
Performs enrichment analysis against predefined gene sets (specified by **genelists.R**)

4. Genome annotation
**annotation.zip**: Set of BED files containing annotations of *Toxoplasma gondii* genome. Currently available with version 26 ToxoDB IDs.
- Introns
- Exons
- Intergenic regions
- Promoters (defined as 1kb upstream of 5'UTR)
- 5'UTR (defined from T. gondii RH strain RNA-seq)
- 3'UTR (defined from T. gondii RH strain RNA-seq)
- Genes active in tachyzoites 
- Genes inactive in tachyzoites
- Promoters of genes active in tachyzoites
- Promoters of genes inactive in tachyzoites
- GFF file containing annotation data for genes with different expression levels in tachyzoites split by 5 (1,2,3,4,5 (low to high))

5. Conversion of *T. gondii* ME49 annotations between versions
Available chain files for use with UCSC LiftOver tool:
tgme49_6.1_To_tgme49_26.1.over.chain - converts annotation version 6.1 to version 26.1
tgme49_9.0_To_tgme49_26.0.over.chain - converts annotation version 9.0 to version 26.0
tgme49_11.0_To_tgme49_26.0.over.chain - converts annotation version 11.0 to 26.0

6. Indexes for read alignments
Derived from genome version 26 (ToxoDB)
- Bowtie
- Bowtie2
- bwa




