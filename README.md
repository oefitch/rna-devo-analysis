# rna-devo-analysis
Methods for differential gene expression and co-expression network analysis of bulk RNA-seq data from developmental samples

THIS WORK IS IN PROGRESS

## Overview

Many methods have been created to analyze RNA-seq data for differentail gene expression. Many of these methods are relatively strightforward for comparing gene expression between a 
wild-type and a modified organism, a healthy and diseased varient, or an untreated and treated system. But, what if you're working with with RNA-seq samples through multiple stages 
of an organisms development? Developmental biology can be complex, and there are many factors to take into account when analyzing developmental data. I wanted to put together a guide
for bulk RNA-seq analysis taking into account best practices for developmental data. 

For the past 6 years I have been working on analyzing a multiomics dataset consisting of bulk RNA-seq and ATAC-seq tissue samples through multiple stages of development. 
This guide will focus on the RNA-seq analysis part of the analysis pipeline. 

<img width="740" height="214" alt="Image" src="https://github.com/user-attachments/assets/34c506a3-b08e-491b-967d-265ae3e1377c" />

## Read Counts 

## Batch Effects

## Differential Analysis

### Normalization
In the case of developmental data , we're often not comparing upregulated and downregulated genes within one sample, but upregulated and downregulated genes between samples. 
