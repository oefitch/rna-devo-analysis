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

After running this, you should have a `.csv` file with your read counts for your sample. Yay! 

## Batch Effects

After identifying read counts, now we are going to consider batch effects. Often times in genomics research, we will send off our sampels for sequencing in batches. Maybe you'll run a few just to test, and then you'll send the rest when you have the funds. Samples being sequenced through different methods, potenitally on different machines, this can lead to batch effects that artificially make your samples seem more different or more similar than they actually are. You can minimize these effects using software tools designed to adjust for batch effects. 

I used `ComBat-seq` to minimize batch effects. I also assigned some additonal covariates because when considering developmental samples, there may be multiple biological variables in your samples. In my case, I was dealing with samples from different tissues at different time points, so I wanted to account for that in my model. Take a look at the `covar_mod` parameter in the `ComBat-seq` [documentation](https://github.com/zhangyuqing/ComBat-seq).

Use the [combatseq.R](https://github.com/oefitch/rna-devo-analysis/blob/main/combatseq.R) script in this repository to account for batch effects and make an adjusted counts `.csv' file that you can use for further analysis.

Now that you have your adjusted counts, you can down perform downstream analysis and make some pretty pictures! Let's get started. 

## Differential Analysis: Identifying Differentially Expressed Genes (DEGs)

### Normalization
In the case of developmental data, we're often not comparing upregulated and downregulated genes within one sample, but upregulated and downregulated genes between samples. Because of this, we need to pick the right type of normalization. 

There are several ways to normalize gene counts, two of the popular ways are to calculate transcripts per million (TPM) or trimmed mean of m-values (TMM), but there are plenty others. When choosing the normalization method to use, it's important to consider what types of comparisons you'll be making. With developmental data, it is likely that you will be comparing gene expression between samples, not just within one sample. It is common to make the mistake of choosng the wrong normalization method when comparing between samples, methods like TPM and RPKM are not designed for compairng between samples [(Zhao et al., 2020)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7373998/). This means that we will need to choose a normalization method that normalizes across samples so that samples can be compared. TMM is a good method for this. The normalization methods that are built in to edgeR and DEseq are also acceptable methods.


