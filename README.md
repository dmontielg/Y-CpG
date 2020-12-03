
# Male-specific age estimation based on Y-chromosomal DNA methylation


#### Athina Vidaki, Diego Montiel González, Benjamin Planterose Jiménez, Manfred Kayser
Department of Genetic Identification, Erasmus MC University Medical Center Rotterdam, Rotterdam, The Netherlands.

## Requirements

    Operating system: tested on Ubuntu 18.04LTS
    R: tested on R version 3.6.1 (2019-07-05) -- "Action of the toes"
    RAM requirements: Especially for data preprocessing, and data normalization use at least 120 GB of RAM.


## Structure
    
    Datasets: Scripts employed in the preprocessing of 450K data whose raw IDAT data are available in the Gene [Expression Omnibus database [GEO](https://www.ncbi.nlm.nih.gov/geo/):
    
    [GSE128235](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE128235), [GSE100386](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100386), [GSE125105](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125105), [GSE61496](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE61496), [GSE87571](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87571), and [GSE115278](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115278).

    *quality_control.R: quality control assessment of probes/cpg-sites, samples and sex prediction

    *normalization.R: normalization pipeline for all raw IDATs. !! Warning, using all 1057 samples requires approximately ~160GB RAM to store matrix transformation and tested with 40 CPUs.

    *train.R: model training using Support vector Machines with Radial Kernel and eps-regression technique

    *predict.R: Script for predicting age using pre-normalized beta values

    *plots.R: Plots generated as seen in the paper. Including violing plots, histograms and scatter-plots

    *annotation.R: Employed in the functional annotation of evCpGs.

    *probes_correlation.R: age correlation among all Y-Cpg probes as in the paper

### data/ folder contains the following: 
    *qc/ list of probes used for data preprocessing/normalization for train + validation and test set
    *annotation/ annotation and correlation files to be used on Integrative Genomics Viewer (IGV)
    *feature_selection/ list of CpG sites based on IQR (> 1.0) and Stepwise-Forward feature selection
    *normalized/ contains normalized methylation beta values BMIQ + ENmix for horvath and Y-chromosome


Please contact me at d.montielgonzalez@erasmusmc.nl for any questions or issues concerning the scripts.

### References and Supporting Information
A. Vidaki *et al* (**2020**). Male-specific age estimation based on Y-chromosomal DNA methylation. *Aging*



