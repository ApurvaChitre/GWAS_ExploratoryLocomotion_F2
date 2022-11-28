# Prerequisites

• This code was run on a Linux based high-performance computing cluster
(HPC system) with a PBS Resource Manager. On this system we had access
to four 2 x 12-core Intel Haswell processors with 64GB of main memory.  
• Each standalone script of the pipeline can be adapted for other linux/
HPC systems with adequate adjustments such as inputting the paths for
the locations of files/ programs and configuring appropriate system
specific scheduler directives.  
• The data required to run the pipeline can be found at \[UCSD library
link\]

# GRM

*Software : GCTA version 1.26.0*  
*Command :*

``` bash

INFILE='/path/to/geno/u01_huda_akil_genotypes'
OUTFILE='/path/to/output/u01_huda_akil_genotypes'


gcta \
--bfile  $INFILE \ 
--thread-num 10 \
--make-grm-bin \
--autosome \
--autosome-num 20 \
--out $OUTFILE 
```

# SNP heritability *h<sup>2</sup>*

*Software : GCTA version 1.26.0*  
*Command :*

``` bash
GENO='/path/to/geno/u01_huda_akil_genotypes'
PHENO='/path/to/pheno/pheno.txt'
GRM='/path/to/grm/u01_huda_akil_genotypes'
OUTFILE='/path/to/out/snp_h2'



gcta \
--bfile  $GENO \ 
--thread-num 5 \
--grm $GRM \
--reml \
--pheno  $PHENO \
--out $OUTFILE 
```

# GWAS

*Software : GCTA version 1.26.0*  
*Command :*

``` bash
GENO='/path/to/geno/u01_huda_akil_genotypes'
PHENO='/path/to/pheno/pheno.txt'
OUTFILE='/path/to/out/gwas_out'



gcta \
--bfile  $GENO \ 
--thread-num 5 \
--mlma-loco \
--pheno  $PHENO \
--out $OUTFILE 
```

# Calling QTLs

To identify QTLs, we scanned each chromosome to determine if there was
at least one SNP that exceeded the permutation-derived threshold, which
was supported by a second SNP within 0.5 Mb that had a p-value that was
within 2 - log10(p) units of the index SNP. Other QTLs on the same
chromosome were tested to ensure that they were conditionally
independent of the first.

Script name: **calling_qtls.R**

# Conditional analysis

To establish conditional independence, we used the top SNP from the
first QTL as a covariate and performed a second GWAS of the chromosome
in question. If the resulting GWAS had an additional SNP with a p-value
that exceeded our permutation-derived threshold, it was considered to be
a second, independent locus. This process was repeated (including all
previously significant SNPs as covariates), until no more QTLs were
detected on a given chromosome.

Script names: **conditional_analysis.R conditional_analysis.sh**

# Annotation of QTL intervals

• Intervals for each identified QTL were determined by identifying all
markers that had a high correlation with the peak marker \>0.6
(r<sup>2</sup> = 0.6). **PLINK-1.90** was used to calculate the LD
r<sup>2</sup> values.  
• Genes within the QTL intervals were queried from the Rat Genome
Database (RGD).  
• Variant annotation for all SNPs within a QTL interval using SNPEff.

Script names: **conditional_analysis.R conditional_analysis.sh**  
**annotating_qtl_intervals.R** **variant_annotation_snpeff.R**
