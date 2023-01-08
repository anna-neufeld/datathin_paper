reconstructionError <- function(dat, this.svd, k) {
  if (k==1) {
    approx <- this.svd$u[,1:k, drop='F']%*%(this.svd$d[1:k])%*%t(this.svd$v[,1:k, drop='F'])
  } else {
    approx <- this.svd$u[,1:k]%*%diag(this.svd$d[1:k])%*%t(this.svd$v[,1:k])
  }
  return(mean((dat-approx)^2))
}


#### Setup that is shared by both. 
library(Seurat)
library(countsplit)
library(ggplot2)
library(patchwork)
library(mclust)

data(pbmc.counts, package="countsplit")
rownames(pbmc.counts) <- sapply(rownames(pbmc.counts), function(u) stringr::str_replace_all(u, "_","-"))


#### NAIVE METHOD
pbmc <- CreateSeuratObject(counts = pbmc.counts, min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc)
pbmc  <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc,features = all.genes)
all_var_genes <- VariableFeatures(pbmc)
X_preproc <- GetAssayData(pbmc, slot="scale.data")[all_var_genes,]

full.svd <- svd(t(as.matrix(X_preproc)))
res.naive <- sapply(1:20, function(u) reconstructionError(t(X_preproc), full.svd,u))
res.naive2 <- sapply(1:20, function(u) sd(full.svd$u[,u]*full.svd$d[u]))


p1 <- ggplot(data=NULL, aes(x=1:20, y=res.naive2))+geom_point()+theme_bw()+xlab("PC")+ylab("Standard Deviation")+
  ggtitle("My SD Elbow Plot", "(Naive Method)")
p2 <- ggplot(data=NULL, aes(x=1:20, y=res.naive))+geom_point()+theme_bw()+xlab("PC")+ylab("MSE")+
  ggtitle("My MSE Plot", "(Naive Method")
p1+p2


#### COUNT SPLIT
set.seed(1)
split <- countsplit(pbmc.counts, epsilon=0.1)
Xtrain <- split$train
Xtest <- split$test

#Xtrain_subset <- Xtrain[all_var_genes, colnames(X_preproc)]
#Xtest_subset <- Xtest[all_var_genes, colnames(X_preproc)]

#XtrainNorm <- apply(Xtrain, 2, function(u) log1p(u/rowTotals*10000))
#XtrainScale <- apply(XtrainNorm, 2, function(u) (u - mean(u)) / (sd(u)))
#svdTrain <- svd(XtrainScale)

#XtestNorm <- apply(Xtest, 2, function(u) log1p(u/rowTotals*10000))
#XtestScale <- apply(XtestNorm, 2, function(u) (u-mean(u))/sd(u))

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

### Set up the test set object. 
rows <- rownames(pbmc.train)
cols <- colnames(pbmc.train)
Xtestsubset <- Xtest[rows,cols]
pbmc.test <- CreateSeuratObject(counts = Xtestsubset, min.cells = 0, min.features = 0)
### This divides by TEST SET totals and I sort of wish it divided by TRAINING SET totals.
### We can discuss I guess. 
### This is gonna get wonky for multiple folds. 
pbmc.test <- NormalizeData(pbmc.test)
pbmc.test <- ScaleData(pbmc.test,features = all.genes.train)
Xtest_preproc <- GetAssayData(pbmc.test, slot="scale.data")[var_genes,]


### This is so slow lol. 
train.svd <- svd(t(as.matrix(Xtrain_preproc)))
### This step is where we could eventually SCALE if needed. 
train.svd$d <- 1/9*train.svd$d
res_cs <- sapply(1:20, function(u) reconstructionError(t(Xtest_preproc), train.svd,u))

plot(res_cs)


#### COUNT SPLIT MANY FOLDS. 
### multiple fold splitting.
### I really want to use fewer genes lol.
### But its cheating to select the highly variable genes FIRST. 
### Bummer because this code is gonna be SO SLOW. 
set.seed(2)
folds <- 8
Xsubset_counts <- t(as.matrix(pbmc.counts[all_var_genes,]))
Xfolds <- apply(Xsubset_counts, 1:2, function(u) rmultinom(1, u, rep(1/folds,folds)))
res_cs_10fold <- matrix(NA, nrow=folds, ncol=20)

rowTotals <- rowSums(Xsubset_counts)

for (fold in 1:folds) {
  print(fold)
  Xtrain <- Xsubset_counts-Xfolds[fold,,]
  Xtest <- Xfolds[fold,,]
  
  pbmc.train <- CreateSeuratObject(counts = Xtrain, min.cells = 3/8, min.features = 200/8)
  pbmc.train[["percent.mt"]] <- PercentageFeatureSet(pbmc.train, pattern = "^MT-")
  pbmc.train <- subset(pbmc.train, subset = nFeature_RNA > 200/8 & nFeature_RNA < 2500/8 & percent.mt < 5/8)
  pbmc.train <- NormalizeData(pbmc.train)
  pbmc.train <- FindVariableFeatures(pbmc.train, selection.method = "vst", nfeatures = 2000)
  all.genes.train <- rownames(pbmc.train)
  pbmc.train <- ScaleData(pbmc.train,features = all.genes.train)
  var_genes <- VariableFeatures(pbmc.train)
  Xtrain_preproc <- GetAssayData(pbmc.train, slot="scale.data")[var_genes,]
  
  ### Set up the test set object. 
  rows <- rownames(pbmc.train)
  cols <- colnames(pbmc.train)
  Xtestsubset <- Xtest[rows,cols]
  pbmc.test <- CreateSeuratObject(counts = Xtestsubset, min.cells = 0, min.features = 0)
  ### This divides by TEST SET totals and I sort of wish it divided by TRAINING SET totals.
  ### We can discuss I guess. 
  ### This is gonna get wonky for multiple folds. 
  pbmc.test <- NormalizeData(pbmc.test)
  pbmc.test <- ScaleData(pbmc.test,features = all.genes.train)
  Xtest_preproc <- GetAssayData(pbmc.test, slot="scale.data")[var_genes,]
  

  #XtrainNorm <- apply(Xtrain, 2, function(u) log1p(u/rowTotals*10000))
  #XtrainScale <- apply(XtrainNorm, 2, function(u) (u - mean(u)) / (sd(u)))
  svdTrain <- svd(t(Xtrain_preproc))
  
  #svdTest <- svd(Xtest_preproc)
  
  #XtestNorm <- apply(Xtest, 2, function(u) log1p(u/rowTotals*10000))
  #XtestScale <- apply(XtestNorm, 2, function(u) (u-mean(u))/sd(u))
  #### Shit lol how does this work if we scaled??? 
  #svdTrain$d <- (1)/(folds-1)*(svdTrain$d)
  res_cs_10fold[fold,] <- sapply(1:20, function(u) reconstructionError(t(Xtest_preproc), svdTrain,u))
}

cs_10fold_res <- colSums(res_cs_10fold)
resCS <- colMeans(res_cs)
ggplot(data=NULL)+
  geom_point(aes(x=1:20, y=(res-min(res))/(max(res)-min(res)), col="Naive"))+
  geom_line(aes(x=1:20, y=(res-min(res))/(max(res)-min(res)), col="Naive"))+
  geom_point(aes(x=1:20, y=(res_cs-min(res_cs))/(max(res_cs)-min(res_cs)), col="Data thinning, train/test"))+
  geom_line(aes(x=1:20, y=(res_cs-min(res_cs))/(max(res_cs)-min(res_cs)), col="Data thinning, train/test"))+
  geom_point(aes(x=1:20, y=(cs_10fold_res-min(cs_10fold_res))/(max(cs_10fold_res)-min(cs_10fold_res)), col="Data thinning, 8 folds"))+
  geom_line(aes(x=1:20, y=(cs_10fold_res-min(cs_10fold_res))/(max(cs_10fold_res)-min(cs_10fold_res)), col="Data thinning, 8 folds"))+
  xlab("Number of PCs")+ylab("Normalized Reconstruction Error")+labs(col="Method")+
  theme_bw()
  
plot(colSums(res_cs_10fold))



