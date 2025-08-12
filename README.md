# Introduction to imaging-based spatial transcriptomics analysis

Biological systems are complex and heterogeneous, involving a diverse suite of cells with unique functions, interacting and communicating to produce complex behaviour at the tissue-level. Understanding such complexity requires detailed measurements that are now made possible through emerging sub-cellular resolution spatial molecular technologies. These technologies allow us to measure the activity of 100s of genes across >100,000 cells from a biological tissue. In the context of cancer, these technologies allow us to study how cancer cells interact with their environment to survive. To unleash the power of these technologies, we need appropriate bioinformatics analysis to generate biological insight.

## Overview

In this workshop, we will analyse a sub-sampled spatial transcriptomics dataset to demonstrate quality control, normalisation, cell typing, spatial domain identification, and domain-specific functional analysis. We will use spatially aware computational methods to perform analyses where such methods are available. At the end of the workshop, attendees will be equipped with the computational tools and data structures used to analyse spatial transcriptomics datasets. They will also be able to analyse their own spatial transcriptomics datasets and decipher complex behaviour in their biological systems of interest.

## Pre-requisites 

This workshop will be relevant to anyone interested in analysing data from emerging high-resolution imaging-based spatial molecular technologies such as 10x Xenium, NanoString CosMx, and Vizgen MERSCOPE.

Attendees will require a laptop with internet access and should have familiarity with R, some familiarity with Bioconductor. Some basic knowledge on spatial transcriptomics technologies and the standard steps needed to gain insight from the data they generate is desirable.

## Time outline

### Session 1 - An introduction to spatial analysis

| Activity                                                             | Time |
|----------------------------------------------------------------------|------|
| Introduction & setup                                                 | 30m  |
| Part 1. Preprocess spatial `omics data (QC, Normalisation, SVG, PCA) | 45m  |
| Part 2. Infer spatial biology (cell typing, spatial domains)         | 45m  |
| Part 3. Interpret spatial biology (DE, GSEA)                         | 45m  |
| Q & A                                                                | 15m  |

### Session 2 - Exploring cell states using PhiSpace

| Activity                                                             | Time |
|----------------------------------------------------------------------|------|
| Introduction & setup                                                 | 30m  |
| Part 1. How to run PhiSpace: the basics.                             | 45m  |
| Part 2. Identify spatial niches/doma                                 | 45m  |
| Part 3. Learn biology from multiple tissue sections                  | 45m  |
| Q & A                                                                | 15m  |

### Session 3 - Deciphering spatial marker genes using JazzPanda

| Activity                                                             | Time |
|----------------------------------------------------------------------|------|
| Introduction & setup                                                 | 30m  |
| Part 1. Marker gene concepts in spatial transcriptomics              | 45m  |
| Part 2. One-sample scenario                                          | 45m  |
| Part 3. Multiple-sample scenario                                     | 45m  |
| Q & A                                                                | 15m  |


## Workshop goals and objectives

### Learning goals

 - Preprocess spatial transcriptomics datasets.
 - Infer cell types, spatial domains, and biological function from spatial transcriptomics datasets.
 - Understand the importance of visualisation in bioinformatics and computational biology.

### Learning objectives

 - Exectue a standard spatial transcritpomics analysis pipeline.
 - Interpret spatial biology from spatial transriptomics datasets.

## Workshop package installation 

### Guide

This is necessary in order to reproduce the code shown in the workshop. 
The workshop is designed for R `4.5.0` and Bioconductor `3.21`, and can be installed as instructed below.

### Via GitHub

Alternatively, you could install the workshop using the commands below in R `>4.5.0` and BioConductor `>3.21`.

```
install.packages('remotes')

# Install workshop package
remotes::install_github("bhuvad/CosMxSpatialAnalysisWorkshop")

# Alternatively, if choose to build vignettes:
remotes::install_github("bhuvad/CosMxSpatialAnalysisWorkshop", build_vignettes = TRUE)
library(CosMxSpatialAnalysisWorkshop)
browseVignettes("CosMxSpatialAnalysisWorkshop")
```
