# rna-devo-analysis
Methods for differential gene expression and co-expression network analysis of bulk RNA-seq data from developmental samples

THIS WORK IS IN PROGRESS

## Overview

Many methods have been created to analyze RNA-seq data for differential gene expression. Many of these methods are relatively strightforward for comparing gene expression between a 
wild-type and a modified organism, a healthy and diseased varient, or an untreated and treated system. But, what if you're working with RNA-seq samples through multiple stages 
of an organism's development? Developmental biology can be complex, and there are many factors to take into account when analyzing developmental data. I wanted to put together a guide
for bulk RNA-seq analysis taking into account best practices for developmental data. 

For the past few years, I have been researching and learning methods for analyzing multiomics datasets. The most recent dataset I worked with consisted of bulk RNA-seq and ATAC-seq tissue samples through multiple stages of development. 
This guide will focus on the RNA-seq part of the analysis pipeline. 

<img width="740" height="214" alt="multiomics data analysis workflow" src="https://github.com/user-attachments/assets/34c506a3-b08e-491b-967d-265ae3e1377c" />

## Read Counts 

The first step in RNA-seq analysis is to determine read counts for genes. There are many packages that can be used to count how many reads map to a specific genomic feature, in the case of RNA-seq, the features we're interested in are genes. Here, I chose to use `HTseq` to determine read counts for RNA-seq. Take a look at the [HTseq documentation](https://htseq.readthedocs.io/en/latest/htseqcount.html) for careful consideration of cases when reads might overlap or align to more than one gene. 

To determine read count, you will need the final cleaned `.bam` files that you generated using a processing pipeline. You can use my processing pipeline from my [multiomics-pipeline](https://github.com/oefitch/multiomics-pipeline) repository. You will also need a `.gtf` file that contains the features of interested, in the case of RNA-seq this will be the gene annotation `.gtf` file from the genome you aligned to. 

Make sure to download the most [recent version](https://htseq.readthedocs.io/en/latest/install.html) of `HTseq`

```
htseq-count -f bam $FINAL.bam REFERENCE_GENOME.gtf -c $COUNTS.csv
```

## Batch Effects



## Differential Analysis

### Normalization
In the case of developmental data , we're often not comparing upregulated and downregulated genes within one sample, but upregulated and downregulated genes between samples. 
