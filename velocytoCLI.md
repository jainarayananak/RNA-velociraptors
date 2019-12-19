# 3. Counting unspliced/spliced transcripts with velocyto CLI

## 3.1 Downloading and configuring Anaconda and velocyto CLI

Before we could use the velocyto Anaconda, we had to download the newest version of Anaconda with python newer than 3.6 on the Taito cluster. The Taito cluster is being deprecated and therefore the cluster has been poorly maintained during 2019, meaning that the most recent version of python had not been pre-installed on the computer. 

Anaconda was downloaded and installed into the root directory:
```{bash}
cd ~
wget https://repo.anaconda.com/archive/Anaconda3-2019.10-Linux-x86_64.sh
bash Anaconda3-2019.10-Linux-x86_64.sh
nano condaBash.sh
```

The following script was filled into the file condaBash.sh:
```{bash}
export CONDA_PATH=~/anaconda3 export PATH=$CONDA_PATH:bin:$PATH export CPATH=$CONDA_PATH:include:$CPATH export LIBRARY_PATH=$CONDA_PATH:lib:$LIBRARY_PATH export LD_LIBRARY_PATH=$CONDA_PATH:lib:$LD_LIBRARY_PATH
```

The conda environment was activated:
```{bash}
source condaBash.sh
```

