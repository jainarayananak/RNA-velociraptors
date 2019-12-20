# 4. velocyto.R
This markdown attempts to follow the guidelines on the RNA velocity analysis using `velocyto.R` package. For more information visit [velocyto.org](http://velocyto.org/).

## 4.1 Load libraries
To get started, necessary libraries are loaded. Make sure you have them and their dependencies already installed.
```{r}
library(velocyto.R)
library(pagoda2)
```

## 4.2 Read in cell cluster assignment
Here the table with information on each cell is loaded. In addition, named vector `cell.colors` is created to cluster cells later on in the analysis.
```{r}
SraRunInfo <- read_csv('SraRunInfoTable.csv') # Read the table
SraRunInfoTableNeo <- filter(SraRunInfo, Cell_type == 'Neoplastic') # Filter for neoplastic cells only
SraRunInfoTableNeo <- mutate(SraRunInfoTableNeo, ColourTissue = ifelse(tissue == 'Tumor', '#00CD00B3', '#FF00FFB3')) # Assign distinct colours to tumour core and tumour periphery
cell.colors <- SraRunInfoTableNeo$ColourTissue # extract the colour column as a vector
names(cell.colors) <- SraRunInfoTableNeo$Run # name each element based on the run accession
```

## 4.3 Read in count matrices from the loom files
Those loom files contain counts of spliced, unspliced, and spanning RNA reads.  
```{r}
# Read in a loom file
ldat <- read.loom.matrices("My77cells.loom")

## 4.4 Assign run accessions for each column in ldat to represent each cell
My7cells_info <- read.csv('My7cells_info.csv') # loom fille which was created first
My43cells_info <- read.csv('My43cells_info.csv') # created second
My27cells_info <- read.csv('My27cells_info.csv') # created third

colnames(ldat$spliced) <-  c(levels(My7cells_info$Run), levels(My43cells_info$Run), levels(My27cells_info$Run))
                             
colnames(ldat$unspliced) <- colnames(ldat$spliced)
colnames(ldat$ambiguous) <- colnames(ldat$spliced)
colnames(ldat$spanning) <- colnames(ldat$spliced)
```

Set up expression matrices.
```{r}
## 4.5 Exonic read (spliced) expression matrix
emat <- ldat$spliced

## 4.6 Normalise
emat <- emat[,colSums(emat)>=1e3]

## 4.7 Remove dupliate row names (genes)
emat <- emat[!duplicated(rownames(emat)),]

## 4.8 Filter out all non-neoplastic cells
all_types <- tibble(Run = colnames(ldat$spliced))
Neoplastic_info <- inner_join(SraRunInfoTableNeo, all_types, by="Run")
emat <- emat[, Neoplastic_info$Run]

```

## 4.9 Reduction of dimension analysis
Create pagoda2 objects and adjust variances.
```{r}
r <- Pagoda2$new(emat,modelType='plain',trim=10,log.scale=T)
r$adjustVariance(plot=T,do.par=T,gam.k=10)
```

Run basic analysis steps to generate cell embedding and clustering.
```{r}
r$calculatePcaReduction(nPcs=100,n.odgenes=3e3,maxit=5000)
r$makeKnnGraph(k=5,type='PCA',center=T,distance='cosine')
r$getKnnClusters(method=multilevel.community,type='PCA',name='multilevel')
r$getEmbedding(type='PCA',embeddingType='tSNE',perplexity=5,verbose=F)
```

Plot embedding, labeling clusters (left) and expression of a gene of interest.
```{r}
par(mfrow=c(1,2))
r$plotEmbedding(type='PCA',embeddingType='tSNE',show.legend=F,mark.clusters=F,min.group.size=2,shuffle.colors=F,mark.cluster.cex=1,alpha=1,main='cell clusters')
r$plotEmbedding(type='PCA',embeddingType='tSNE',colors=r$counts[,"GAPDH"], alpha = 1, main='GAPDH')
```

### Velocity estimation

Restrict cells that passed p2 filter
```{r}
emat <- ldat$spliced
nmat <- ldat$unspliced
emat <- emat[,rownames(r$counts)]
nmat <- nmat[,rownames(r$counts)]
smat <- ldat$spanning
smat <- smat[,rownames(r$counts)]
```

Take cluster labels and tSNE embeddings
```{r}
cluster.label <- r$clusters$PCA[[1]]
#cell.colors <- pagoda2:::fac2col(cluster.label)
emb <- r$embeddings$PCA$tSNE
```

In addition to clustering and the t-SNE embedding, from the p2 processing we will also take a cell-cell distance, which will be better than the default whole-transcriptome correlation distance that velocyto.R would normally use.
```{r}
cell.dist <- as.dist(1-armaCor(t(r$reductions$PCA)))
```

Filter genes based on the minimum average expresion magnitude (in at least one of the clusters), output total number of resulting valid genes
```{r}
emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.5)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.05)
smat <- filter.genes.by.cluster.expression(smat,cluster.label,min.max.cluster.average = 0.5)
length(intersect(rownames(emat),rownames(emat)))
```
### Estimate RNA velocity
```{r}
fit.quantile <- 0.05
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=5,cell.dist=cell.dist,fit.quantile=fit.quantile)
```

Visualize velocity on the t-SNE embedding, using velocity vector fields
```{r}
vel <- rvel.cd; arrow.scale=18; cell.alpha=0.7; cell.cex=1; fig.height=4; fig.width=4.5;
show.velocity.on.embedding.cor(emb,vel,n=300,scale='sqrt',cell.colors=ac(cell.colors,alpha=1),cex=cell.cex,arrow.scale=arrow.scale,arrow.lwd=1)
```

Visualize a fit for a particular gene.
```{r}
gene <- "FGFR3"
gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells = 2,kGenes=1,fit.quantile=fit.quantile,cell.emb=emb,cell.colors=ac(cell.colors, alpha=1), cell.dist=cell.dist,show.gene=gene,old.fit=rvel.cd,do.par=T)
```

Alternatively visualise RNA velocity by projecting observed and extrapolated cells onto PCs:
```{r}
pca.velocity.plot(vel,nPcs=2,plot.cols=2,cell.colors=ac(cell.colors,alpha=0.7),cex=1.2,pcount=0.1,pc.multipliers=c(1,-1))
```




