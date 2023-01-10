#### Computes SSE between matrix (dat) and its rank-K approximation. 
reconstructionError <- function(dat, this.svd, k) {
  if (k==1) {
    approx <- this.svd$u[,1:k, drop='F']%*%(this.svd$d[1:k])%*%t(this.svd$v[,1:k, drop='F'])
  } else {
    approx <- this.svd$u[,1:k]%*%diag(this.svd$d[1:k])%*%t(this.svd$v[,1:k])
  }
  return(sum((dat-approx)^2))
}


library(Seurat)
library(countsplit)
library(ggplot2)
library(patchwork)
library(mclust)


#### The data from 10X genomics was previously downloaded and included in the countsplit R package. 
data(pbmc.counts, package="countsplit")
#### This is necessary to avoid later errors in preprocessing. 
rownames(pbmc.counts) <- sapply(rownames(pbmc.counts), function(u) stringr::str_replace_all(u, "_","-"))


#### Reproduce the analysis from the Seurat Guided Clustering tutorial. 
pbmc <- CreateSeuratObject(counts = pbmc.counts, min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc)
pbmc  <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc,features = all.genes)
all_var_genes <- VariableFeatures(pbmc)
X_preproc <- GetAssayData(pbmc, slot="scale.data")[all_var_genes,]

#### The two equivalent formations of the error using the "naive" method. 
full.svd <- svd(t(as.matrix(X_preproc)))
res.naive <- sapply(1:20, function(u) reconstructionError(t(X_preproc), full.svd,u))
res.naive2 <- sapply(1:20, function(u) var(full.svd$u[,u]*full.svd$d[u]))


#### Redoing the analysis via data thinning. 
set.seed(1)
split <- countsplit(pbmc.counts, epsilon=0.5)
Xtrain <- split$train
Xtest <- split$test

## Set up the training set object and run the suggested Seurat normalization. 
pbmc.train <- CreateSeuratObject(counts = Xtrain, min.cells = 3, min.features = 200)
pbmc.train[["percent.mt"]] <- PercentageFeatureSet(pbmc.train, pattern = "^MT-")
pbmc.train <- subset(pbmc.train, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc.train <- NormalizeData(pbmc.train)
pbmc.train <- FindVariableFeatures(pbmc.train, selection.method = "vst", nfeatures = 2000)
all.genes.train <- rownames(pbmc.train)
pbmc.train <- ScaleData(pbmc.train,features = all.genes.train)
var_genes <- VariableFeatures(pbmc.train)
Xtrain_preproc <- GetAssayData(pbmc.train, slot="scale.data")[var_genes,]

### Set up the test set object.  We need to make sure that the rows/columns match those of the 
## training set. We do not do additional filtering on the test set cells. 
rows <- rownames(pbmc.train)
cols <- colnames(pbmc.train)
Xtestsubset <- Xtest[rows,cols]
pbmc.test <- CreateSeuratObject(counts = Xtestsubset, min.cells = 0, min.features = 0)
pbmc.test <- NormalizeData(pbmc.test)
pbmc.test <- ScaleData(pbmc.test,features = all.genes.train)
Xtest_preproc <- GetAssayData(pbmc.test, slot="scale.data")[var_genes,]

### Recontruction set error between TRAINING SET svd and TEST SET counts. 
train.svd <- svd(t(as.matrix(Xtrain_preproc)))
SSEs_cs <- sapply(1:20, function(u) reconstructionError(t(Xtest_preproc), train.svd,u))


p11 <- ggplot(data=NULL, aes(x=1:20, y=res.naive))+geom_point()+theme_bw()+xlab("# of PCs (K)")+ylab("Sum of squared errors (SSE)")+
  ggtitle(expression("(b) SSE, computed on"~tilde(Y))) 

p13 <- ggplot(data=NULL, aes(x=1:20, y=sapply(1:20, function(u) sd(full.svd$u[,u]*full.svd$d[u]))))+geom_point()+theme_bw()+xlab("# of PCs (K)")+ylab("Standard deviation")+
  ggtitle(expression("(a) SD of PCs, computed on"~tilde(Y)))



p21 <- ggplot(data=NULL, aes(x=1:20, y=SSEs_cs))+geom_point()+theme_bw()+xlab("# of PCs (K)")+ylab("Sum of squared errors (SSE)")+
  ggtitle("(c) SSE, data thinning")#+scale_y_log10()

p13+p11+
  p21+
  plot_layout(nrow=1, byrow=TRUE)
  #plot_annotation(tag_level="a")

ggsave("~/Dropbox/Generalized Count Splitting/Figures/seuratplot.png", width=10, height=4)










