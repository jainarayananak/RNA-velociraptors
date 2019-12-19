# 3. Counting unspliced/spliced transcripts with velocyto CLI

## 3.1 Downloading and configuring Anaconda

Before we could use the velocyto Anaconda, we had to download the newest version of Anaconda with python newer than 3.6 on the Taito cluster. The Taito cluster is being deprecated and therefore the cluster has been poorly maintained during 2019, meaning that the most recent version of python had not been pre-installed on the computer. 

Anaconda was downloaded and installed into the root directory:
```{bash}
cd ~
wget https://repo.anaconda.com/archive/Anaconda3-2019.10-Linux-x86_64.sh
bash Anaconda3-2019.10-Linux-x86_64.sh
nano condaBash.sh
```

The following script was filled into the file condaBash.sh:
```{bash, eval=FALSE}
export CONDA_PATH=~/anaconda3 export PATH=$CONDA_PATH:bin:$PATH export CPATH=$CONDA_PATH:include:$CPATH export LIBRARY_PATH=$CONDA_PATH:lib:$LIBRARY_PATH export LD_LIBRARY_PATH=$CONDA_PATH:lib:$LD_LIBRARY_PATH
```

The conda environment was activated:
```{bash}
source condaBash.sh
```

## 3.2 Installing velocyto
In the activated conda environment, the necessary dependencies of velocyto were installed (the pip installer is more robust for pysam and therefore preferable):

```{bash}
conda install numpy scipy cython numba matplotlib scikit-learn h5py click
pip install pysam
```

Then the velocyto package was downloaded with pip:

```{bash}
pip install velocyto
```

## 3.3 Counting and unspliced/spliced transcripts
The velocyto CLI was performed with as a long-running screen process as in [2.3](STAR.md). In other words, a connction to taito.csc.fi was opened, and then a screen job was initiated

```{bash, eval=TRUE}
screen -R
```

In the screen session, a taito-shell session was launched

```{bash, eval=TRUE}
sinteractive
```

Then the velocyto CLI tool was launched

```{bash, eval=TRUE}
while read p; do STAR --genomeDir GENOME_Indices_final --readFilesIn FINAL_SRA/"$p"_1.fastq FINAL_SRA/"$p"_2.fastq --outFileNamePrefix FINAL_align/star_final_$p --outSAMtype BAM SortedByCoordinate --runThreadN 5; done < FINAL_TEST_SET 
```