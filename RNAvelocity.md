---
title: "RNA velocity"
author: "Our group"
date: "19 December 2019"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Downloading data from SRA

We chose 100 good cells and saved the indices as FINAL_SET_TEST, which we transferred to the taito-shell. 

```{bash, eval=FALSE}
scp FINAL_SET_TEST neiroja1@taito-shell.csc.fi:/wrk/neiroja1/human_genome/FINAL_SET_TEST
ssh neiroja1@taito-shell.csc.fi
cd $WRKDIR
mkdir human_genome
cd human_genome
```

We downloaded the compiled sra-tookit

```{bash, eval=FALSE}
wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.0/sratoolkit.2.10.0-ubuntu64.tar.gz
tar -zxvf sratoolkit.2.10.0-ubuntu64.tar.gz
mv sratoolkit.2.10.0-ubuntu64 sratoolkit
mkdir FINAL_SRA
cd FINAL_SRA
while read p; do ../SRA/sratoolkit/bin/fastq-dump --split-files $p; done < ../FINAL_TEST_SET
```

## Building the index

We donwloaded the human genome and the corresponding annotation file

```{bash, eval=FALSE}
wget ftp://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_rm.primary_assembly.fa.gz
tar -zxvf Homo_sapiens.GRCh.38.dna_rm.primary_assembly.fa.gz

wget ftp://ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh38.98.gtf.gz
tar -zxvf Homo_sapiens.GRCh38.98.gtf.gz
```

Then we created the file star_index_final.sh and composed the following sbatch script.

```{bash, eval=FALSE}
#!/bin/bash -l
# author: neiroja1
#SBATCH -J StarIndex
#SBATCH --constraint="snb|hsw"
#SBATCH -o StarIndex_out_%j.txt
#SBATCH -e StarIndex_err_%j.txt
#SBATCH -p hugemem
#SBATCH -n 1
#SBATCH --cpus-per-task=32
#SBATCH -t 12:00:00
#SBATCH --mem=300000 ## 32GB memory


genomeGtf="/wrk/neiroja1/human_genome/Homo_sapiens.GRCh38.98.gtf"
if [ ! -f "$genomeGtf" ]
then
        echo "gtf not found"
fi

genomeFasta="/wrk/neiroja1/human_genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
if [ ! -f "$genomeFasta" ]
then
        echo "genome not found"
fi

mkdir GENOME_Indices_final

STAR --runMode genomeGenerate --genomeDir GENOME_Indices_final \
--genomeFastaFiles $genomeFasta \
--sjdbGTFfile $genomeGtf \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 65 --limitGenomeGenerateRAM 160263683456
```

The file was executed

```{bash, eval=FALSE}
sbatch star_align_final
```

