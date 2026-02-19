# cryptic_choice_anubis

This project contains code for analyses and figures related to:

Mechanisms of cryptic female choice in the vaginal tract of a primate species. 
Rachel M. Petersen, Lee (Emily) M. Nonnamaker, Jaclyn Anderson, Christina Bergey, Christian Roos, Amanda D. Melin, and James P. Higham

Specifically, the scripts provided here are:

  1. RNAseq_data_processing.sh : trim and map reads, build and intersect with kracken metagenomic database, generate gene counts matrix

  2. DE_phase.Rmd : filter samples, normalize counts, and run differential expression analysis to test for differences in vaginal gene expression across ovarian cycle phases (in the absense of mating) and create figures in Fig 1

  3. DE_postcop.Rmd : filter samples, normalize counts, and run differential expression analyses to test for differences in vaginal gene expression in noncopulatory versus postcopulatory samples and create figures in Fig 2

  4. DE_genotype.Rmd : filter samples, normalize counts, and run differential expression analyses to test for an interactive effect between postcopulatory status and 10 different aspects of genetic make-up and create figures in Fig 3

  5. pH analyses.Rmd : run linear mixed effects models to test whether 1) ovarian cycle phase, 2) copulation, or 3) male genotype (following copulation) are associated with vaginal pH and create figures in Fig 4

  6. Semen_gene_analyses.Rmd: identifies genes overexpressed in semen samples and creates "semengenevector.txt" used to filter genes in postcopulatory DE analyses

Data required as input for each script are hosted on Zenodo (10.5281/zenodo.14976902). The raw RNA-seq data (FASTQ files) are available on NCBI SRA at accession PRJNA1232174.

Note that specific paths in the scripts will not work on your computer and you will need to change them accordingly.

Please contact me at rpetersen42@gmail.com with any questions.
