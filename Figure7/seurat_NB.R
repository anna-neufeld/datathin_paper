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
library(Matrix)


#### The data from 10X genomics was previously downloaded and included in the countsplit R package. 
data(pbmc.counts, package="countsplit.tutorials")
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

#### Redoing the analysis via Poisson data thinning. 
set.seed(1)
split <- countsplit(pbmc.counts, epsilon=0.5)
Xtrain <- split[[1]]
Xtest <- split[[2]]

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


#### Redoing the analysis via NB data thinning. Need to estimate the overdispersions first. 
set.seed(1)
temp.ests <- sctransform::vst(pbmc.counts, n_genes=10000)$model_pars[,1]
### If it didn't get an estimate, assume it is Poisson
overdisp.ests <- rep(Inf, NROW(pbmc.counts))
names(overdisp.ests) <- rownames(pbmc.counts)
overdisp.ests[names(temp.ests)] <- temp.ests


split.nb <- countsplit(t(pbmc.counts), epsilon=0.5, overdisps = overdisp.ests)
Xtrain.nb <- split.nb[[1]]
Xtest.nb <- split.nb[[2]]

## Set up the training set object and run the suggested Seurat normalization. 
pbmc.train.nb <- CreateSeuratObject(counts = Xtrain.nb, min.cells = 3, min.features = 200)
pbmc.train.nb[["percent.mt"]] <- PercentageFeatureSet(pbmc.train.nb, pattern = "^MT-")
pbmc.train.nb <- subset(pbmc.train.nb, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc.train.nb <- NormalizeData(pbmc.train.nb)
pbmc.train.nb <- FindVariableFeatures(pbmc.train.nb, selection.method = "vst", nfeatures = 2000)
all.genes.train.nb <- rownames(pbmc.train.nb)
pbmc.train.nb <- ScaleData(pbmc.train.nb,features = all.genes.train.nb)
var_genes <- VariableFeatures(pbmc.train.nb)
Xtrain.nb_preproc <- GetAssayData(pbmc.train.nb, slot="scale.data")[var_genes,]

### Set up the test set object.  We need to make sure that the rows/columns match those of the 
## train.nbing set. We do not do additional filtering on the test set cells. 
rows <- rownames(pbmc.train.nb)
cols <- colnames(pbmc.train.nb)
Xtest.nbsubset <- Xtest.nb[rows,cols]
pbmc.test.nb <- CreateSeuratObject(counts = Xtest.nbsubset, min.cells = 0, min.features = 0)
pbmc.test.nb <- NormalizeData(pbmc.test.nb)
pbmc.test.nb <- ScaleData(pbmc.test.nb,features = all.genes.train.nb)
Xtest.nb_preproc <- GetAssayData(pbmc.test.nb, slot="scale.data")[var_genes,]

### Recontruction set error between TRAINING SET svd and TEST SET counts. 
train.svd.nb <- svd(t(as.matrix(Xtrain.nb_preproc)))
SSEs_cs.nb <- sapply(1:20, function(u) reconstructionError(t(Xtest.nb_preproc), train.svd.nb,u))


p11 <- ggplot(data=NULL, aes(x=1:20, y=res.naive))+geom_point()+theme_bw()+
  xlab("Number of principal components (K)")+
  ylab("Sum of squared errors")+
  ggtitle(expression("(b) Sum of squared errors,"), expression("computed on"~tilde(Y))) 

p13 <- ggplot(data=NULL, aes(x=1:20, y=sapply(1:20, function(u) sd(full.svd$u[,u]*full.svd$d[u]))))+geom_point()+theme_bw()+
  xlab("Principal component (K)")+ylab("Standard deviation")+
  ggtitle(expression("(a) Standard deviation of principal"), 
          expression("components, computed on"~tilde(Y))) 



p21 <- ggplot(data=NULL, aes(x=1:20, y=SSEs_cs))+geom_point()+theme_bw()+
  xlab("Number of principal components (K)")+ylab("Sum of squared errors")+
  ggtitle("(c) Sum of squared errors,", "computed with Poisson thinning")#+scale_y_log10()


p22 <- ggplot(data=NULL, aes(x=1:20, y=SSEs_cs.nb))+geom_point()+theme_bw()+
  xlab("Number of principal components (K)")+
  ylab("Sum of squared errors")+
  ggtitle("(d) Sum of squared errors,", "computed with Neg. bin. thinning")#+scale_y_log10()

p13+p11+
  p21+p22+
  plot_layout(nrow=1, byrow=TRUE) &
  theme(
    plot.title=element_text(size=12),
    plot.subtitle=element_text(size=12)
  )
  #plot_annotation(tag_level="a")

ggsave("~/Dropbox/Generalized Count Splitting/JMLR-resubmit-feb-2024/Figures/seuratplotNEW.png", width=16, height=5)










