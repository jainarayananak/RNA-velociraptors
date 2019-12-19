# Downloading raw data, STAR indexing and alignment

## 2.1 Downloading data from SRA

We chose 100 good cells and saved the indices as FINAL_SET_TEST, which we transferred to the Taito cluster (Jakke Neiro's account). 

```{bash, eval=FALSE}
scp FINAL_SET_TEST neiroja1@taito-shell.csc.fi:/wrk/neiroja1/human_genome/FINAL_SET_TEST
ssh neiroja1@taito-shell.csc.fi
cd $WRKDIR
mkdir human_genome
cd human_genome
```

We downloaded the compiled sra-toolkit (2.10.0) and used it to download the sequence files from the SRA (Sequence read archive). In order to use the sra-toolkit, the compatible version perl was activated by loading the cluster's biokit.   

```{bash, eval=FALSE}
wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.0/sratoolkit.2.10.0-ubuntu64.tar.gz
tar -zxvf sratoolkit.2.10.0-ubuntu64.tar.gz
mv sratoolkit.2.10.0-ubuntu64 sratoolkit
mkdir FINAL_SRA
cd FINAL_SRA
module load biokit
while read p; do ../SRA/sratoolkit/bin/fastq-dump --split-files $p; done < ../FINAL_TEST_SET
```

## 2.2 Building the index

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

## 2.3 Alignment to the genome
Thhe reads were aligned to the genome with STAR (version 2.2.1), which was pre-installed on the cluster. The alignment was performed as long-running screen process as no SLURM quota was left at this point. 

First a connction to taito.csc.fi was opened, and then a screen job was initiated

```{bash, eval=TRUE}
screen -R
```

In the screen session, a taito-shell session was launched

```{bash, eval=TRUE}
sinteractive
```

Then the STAR alignment was initiated using 5 threads

```{bash, eval=TRUE}
while read p; do STAR --genomeDir GENOME_Indices_final --readFilesIn FINAL_SRA/"$p"_1.fastq FINAL_SRA/"$p"_2.fastq --outFileNamePrefix FINAL_align/star_final_$p --outSAMtype BAM SortedByCoordinate --runThreadN 5; done < FINAL_TEST_SET 
```