reconstructionError <- function(dat, this.svd, k) {
  if (k==1) {
    approx <- this.svd$u[,1:k, drop='F']%*%(this.svd$d[1:k])%*%t(this.svd$v[,1:k, drop='F'])
  } else {
    approx <- this.svd$u[,1:k]%*%diag(this.svd$d[1:k])%*%t(this.svd$v[,1:k])
  }
  return(sum((dat-approx)^2))
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
res.naive2 <- sapply(1:20, function(u) var(full.svd$u[,u]*full.svd$d[u]))



#### COUNT SPLIT
set.seed(1)
split <- countsplit(pbmc.counts, epsilon=0.8)
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

train.svd <- svd(t(as.matrix(Xtrain_preproc)))
### This step is where we could eventually SCALE if needed. 
#train.svd$d <- 1/9*train.svd$d
res_cs <- sapply(1:20, function(u) reconstructionError(t(Xtest_preproc), train.svd,u))

p11 <- ggplot(data=NULL, aes(x=1:20, y=sapply(1:20, function(u) reconstructionError(t(X_preproc), full.svd,u))))+geom_point()+theme_bw()+xlab("PC")+ylab("MSE")+
  ggtitle("Naive method: SSE") 
#p12 <- ggplot(data=NULL, aes(x=1:20, y=cumsum(full.svd$d[1:20]^2)), full.svd,u)+geom_point()+theme_bw()+xlab("PC")+ylab("MSE")+
#  ggtitle("Sum D^2", "(Naive Method: Direct)") 
p13 <- ggplot(data=NULL, aes(x=1:20, y=sapply(1:20, function(u) sd(full.svd$u[,u]*full.svd$d[u]))))+geom_point()+theme_bw()+xlab("PC")+ylab("Standard Deviation")+
  ggtitle("Naive method: SD of PCs")
#p14 <- ggplot(data=NULL, aes(x=1:20, y=(full.svd$d[1:20]^2)))+geom_point()+theme_bw()+xlab("PC")+ylab("Standard Deviation")+
#  ggtitle("D^2", "(Naive Method: Direct)")+
#  geom_hline(yintercept=0, col="red")

#SSEs <- sapply(1:20, function(u) reconstructionError(t(X_preproc), full.svd,u))
#p31 <- ggplot(data=NULL, aes(x=1:20, y=SSEs))+geom_point()+theme_bw()+xlab("PC")+ylab("MSE")+
 # ggtitle("SSE ", "(Naive Method: Direct)") 
#Cumulative_Ds <- sum(X_preproc^2)-SSEs
#p32 <- ggplot(data=NULL, aes(x=1:20, y=Cumulative_Ds))+geom_point()+theme_bw()+xlab("PC")+ylab("Standard Deviation")+
#  ggtitle("Sum D^2", "(Naive Method: from SSE)")
#Ds <- c(Cumulative_Ds[1], sapply(2:20, function(u) Cumulative_Ds[u]- Cumulative_Ds[u-1]))
#p33 <- ggplot(data=NULL, aes(x=1:20, y=sqrt(Ds/sum(X_preproc^2))))+geom_point()+theme_bw()+xlab("PC")+ylab("MSE")+
 # ggtitle("SD of PCs", "(Naive Method: from SSE)")+
#  geom_hline(yintercept=0, col="red")
#p34 <- ggplot(data=NULL, aes(x=1:20, y=Ds))+geom_point()+theme_bw()+xlab("PC")+ylab("MSE")+
#  ggtitle("D^2 ", "(Naive Method: from SSE)")+
#  geom_hline(yintercept=0, col="red")

#p11+p12+p13+p14+p31+p32+p33+p34+plot_layout(nrow=2, byrow=TRUE)

SSEs_cs <- sapply(1:20, function(u) reconstructionError(t(Xtest_preproc), train.svd,u))

#Cumulative_Ds_cs<- sum(Xtest_preproc^2)-SSEs_cs

#Ds_cs <- c(Cumulative_Ds_cs[1], sapply(2:20, function(u) Cumulative_Ds_cs[u]- Cumulative_Ds_cs[u-1]))
p21 <- ggplot(data=NULL, aes(x=1:20, y=SSEs_cs))+geom_point()+theme_bw()+xlab("PC")+ylab("MSE")+
  ggtitle("Data Thinning: SSE")
#p22 <- ggplot(data=NULL, aes(x=1:20, y=Cumulative_Ds_cs))+
#                geom_point()+theme_bw()+xlab("PC")+ylab("MSE")+
#  ggtitle("Sum D^2 kinda", "(Count split: 50/50. From SSE")
  
  
#p23 <- ggplot(data=NULL, aes(x=1:20, y=sqrt(Ds_cs*sum(Xtest_preproc^2))))+geom_point()+theme_bw()+xlab("PC")+ylab("MSE")+
#  ggtitle("SD of PCs", "(Count split: 50/50. From SSE)")+
#  geom_hline(yintercept=0, col="red")
#p24 <- ggplot(data=NULL, aes(x=1:20, y=Ds_cs))+geom_point()+theme_bw()+xlab("PC")+ylab("MSE")+
#  ggtitle("D^2", "(Count split: 50/50. From SSE)")+
#  geom_hline(yintercept=0, col="red")

p13+p11+
  p21+
  plot_layout(nrow=1, byrow=TRUE)
ggsave("figures/seuratplot.png")










