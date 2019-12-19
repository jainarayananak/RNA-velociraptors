# WSL installation

## velocyto.R package

velocyto.R (velox + κύτος, quick cell) is a package for the analysis of expression dynamics in single cell RNA seq data. In particular, it enables estimations of RNA velocities of single cells by distinguishing unspliced and spliced mRNAs in standard single-cell RNA sequencing protocols. For more details on this R package see <http://velocyto.org/> or <https://github.com/velocyto-team/velocyto.R>.

## prerequisite One - linux kernel 

velocyto.R requires unix-flavoured operating system with good storage capacities to enable processing large RNA sequencing data. However, DTC computers do not have enough storage capacities (10GB or more is recommened). In addition, their general computing performance is only moderate. One convenient solution is to install linux kernel on a Windows device with better specifications. Note, everything that follows in this markdown was done with a lenovo s740-14iil laptop (intel i7-1065G7, 16 GB DDR4 RAM, 1 TB SSD).

- Update Windows 10 using Windows Insider Program to Windows Build 18917 or higher
- Launch PowerShell (run as administrator) and enable Windows Subsystem for Linux (WSL)
```{bash, eval = FALSE}
Enable-WindowsOptionalFeature -Online -FeatureName Microsoft-Windows-Subsystem-Linux 
```
- Download and install Ubuntu 18.04 LTS distro from Microsoft Store <https://www.microsoft.com/en-gb/p/ubuntu-1804-lts/9n9tngvndl3q?activetab=pivot:overviewtab>
- Enable 'Virtual Machine Platform'
```{bash, eval = FALSE}
dism.exe /online /enable-feature /featurename:Microsoft-Windows-Subsystem-Linux /all /norestart
dism.exe /online /enable-feature /featurename:VirtualMachinePlatform /all /norestart
```
- Set a distro to be backed by WSL version 2 using PowerShell (run as administrator)
```{bash, eval = FALSE}
wsl --set-version Ubuntu 2
```
- Verify what version of WSL Ubuntu is using by the following command (only available in Windows Build 18917 or higher)
```{bash, eval = FALSE}
wsl -l -v
```

Source: <https://docs.microsoft.com/en-us/windows/wsl/wsl2-install> 

## prerequisite Two - R base
- Launch Ubuntu 18.04 terminal to get started 
- Install the packages necessary to add a new repository over HTTPS
```{bash, eval = FALSE}
sudo apt install apt-transport-https software-properties-common
```
- Enable the CRAN repository and add the CRAN GPG key to your system using the following commands
```{bash, eval = FALSE}
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'
```
- Now that the apt repository is added, update the packages list and install the R package by typing
```{bash, eval = FALSE}
sudo apt update
sudo apt install r-base
```
- To verify that the installation was successful run the following command which will print the R version
```{bash, eval = FALSE}
R --version
```

Source: <https://linuxize.com/post/how-to-install-r-on-ubuntu-18-04/>

## prerequisite Three - RStudio Server

RStudio Server provides a browser based interface to a version of R running on a remote Linux server.

- again in the Ubuntu terminal, install gdebi-core which is a simple tool to install deb files
```{bash, eval = FALSE}
sudo apt-get install gdebi-core
```
- download Rstudio Server
```{bash, eval = FALSE}
wget https://download2.rstudio.org/server/bionic/amd64/rstudio-server-1.2.5019-amd64.deb
```
- install Rstudio Server which is a deb file
```{bash, eval = FALSE}
sudo gdebi rstudio-server-1.2.5019-amd64.deb
```
- access Rstudio Server on a web browser by typing following in the address bar
```{bash, eval = FALSE}
http://<server-ip>:8787
```
- to find your local server ip <server-ip> , in the Ubuntu terminal type
```{bash, eval = FALSE}
hostname -I
```

Source: <https://rstudio.com/products/rstudio/download-server/debian-ubuntu/>

## prerequisite Four - devtools
The devtools package is useful for installing R packages directly from GitHub. Later on, velocyto.R will be downloaded in this way.

- install system dependencies for devtools
```{bash, eval = FALSE}
sudo apt install build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev
```
- install the devtools package by opening R in the Ubuntu terminal. Note installing this package takes a few minutes
```{bash, eval = FALSE}
sudo -i R
install.packages('devtools')
```

Source: <https://www.digitalocean.com/community/tutorials/how-to-install-r-packages-using-devtools-on-ubuntu-18-04#conclusion>

## prerequisite Five - system dependencies for velocyto.R  
- Run update command to update package repositories and get latest package information
```{bash, eval = FALSE}
sudo apt-get update -y
```
- Run the install command with -y flag to quickly install the packages and dependencies
```{bash, eval = FALSE}
sudo apt-get install -y libhdf5-dev
sudo apt-get install -y libboost-all-dev
sudo apt-get install -y gcc-multilib
```
Sources:

- velocyto dependecies  <http://velocyto.org/>
- openmp library <http://manueldeveloper.blogspot.com/2012/04/how-to-install-openmp-in-ubuntu-linux.html>
- libhdf5-dev package <https://zoomadmin.com/HowToInstall/UbuntuPackage/libhdf5-dev>

## prerequisite Six - pcaMethods R package
This dependency comes from the Bioconductor. It is not in CRAN depository so it has to be installed manually.
- Launch R in the Ubuntu terminal
- install BiocManager package which helps to install and manage packages from the Bioconductor 
```{r, eval = FALSE}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
```
- install the pcaMethods package
```{r, eval = FALSE}
BiocManager::install("pcaMethods")
```

Source <http://www.bioconductor.org/packages/release/bioc/html/pcaMethods.html>

## install the velocyto.R package
- Launch R in the Ubuntu terminal
- install the velocyto.R package using devtools
```{r, eval = FALSE}
devtools::install("velocyto.R")
```

Source <http://velocyto.org/>