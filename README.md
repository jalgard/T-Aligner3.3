# T-Aligner3.3
Legacy version of T-Aligner, U-insertion/deletion analysis toolkit

## Description

T-Aligner was specially written for RNA editing analysis in mitochondrial genomes of kinetoplastids using NGS data as input. 

### Typical pipeline is:

1. Map edited reads on the maxicircle reference contig using 'alignlib' tool
2. Define edited regions, prepare flanked references for each potential cryptogene
3. Align reads on each cryptogene separately with 'alignlib' and reconstruct transcriptome with 'findorfs'
4. Generate per-cryptogene T-Aligner dotmatrices using 'coverage.py' plotter
5. Determine canonical mRNA edited sequence using combination of 'blastp', coverage and T-Aligner dotmatrix plots
6. Study differential editing (if there are several repicates of RNAseq prepared for two experimental conditions) with DEST-pipeline (based on EdgeR)
7. Find gRNAs using 'findorfs' tool and verify hits with smallRNA sequencing data (if possible)
8. Build gRNA:mRNA editing maps (various tools depending on the situation)

## Legacy versions

We develop and maintain T-Aligner constantly. Most relevant version of T-Aligner will be accessible at https://github.com/jalgard/T-Aligner.
Please use currently maintained version. 

First version of T-Aligner (v. 1-2) was using Qt library to plot dot matrices and other stuff, but the code was written for Qt4.x, which 
is currently not available. After v. 3 we completely switched to python's matplotlib for plotting.
