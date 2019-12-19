# RNA-velociraptors

Alice Eddershaw<sup>1</sup>, Ester Paolocci<sup>1</sup>, Ashwin Jainarayanan<sup>1</sup>, Jakke Neiro<sup>12</sup>, Szymon Stodolak<sup>1</sup>

1. Doctoral Training Centre for Interdisciplinary Bioscience DTP, University of Oxford, United Kingdom
2. Department of Mathematics and Statistics, Faculty of Science, University of Helsinki, Finland

Bioinformatics Hackathon 

20th December 2019

# Summary
Glioblastoma is a devastating central nervous system tumour arising from glial cells that causes the death of about 250,000 people each year. In our Hackathon project, we aimed to investigate whether RNA velocity can be used to understand the differences between tumour layers in glioblastoma and help us identify the origin of tumour cells. Furthermore, we aimed to see whether RNA velocity could be used to predict tumour migration and to identify relevant molecular markers for the progression of the disease. Our presentation of the project can be found [here](RNAVelocity_Hackathon.pdf).


## 1. Installing WSL on Windows 10
Before even starting, we we were posed with a computational problem: the PCs provided by the Doctoral Training Centre (DTC) did not have storage capacities required for storing downloaded RNAseq data and perform computationally heavy tasks. Therefore, we downloaded the Windows Subsystem for Linux (WSL) on our own personal laptops to be able to perform our project. The documentation for how we set up our data science environment can be seen [here](WSLInstallation.md)

## 2. Downloading raw data, STAR indexing and alignment
However, a normal PC does not enough computational capacities for creating genome indices and aligning several data sets of single-cell data, so we ended up with performing this step on computing clusters (Preliminary data set of 4 cells on ARC, Oxford and 100 cells on Taito, CSC, Kajaani). The raw data (SRA indices [here](FINAL_TEST_SET.txt)) was downloaded with SRA tools and the alignment of the reads performed with STAR. Full documentation of this step (100 cells analyzed done on Taito) can be seen [here](STAR.md).

## 3. Counting unspliced/spliced transcripts with velocyto CLI
The bam files generated by STAR alignment were analyzed with the velocyto CLI tool to count the proportion of unspliced/spliced transcripts for each gene. The full documentation is to be seen [here](velocytoCLI.md).

## 4. velocyto.R


