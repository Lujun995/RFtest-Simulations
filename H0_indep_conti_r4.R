#H0_i_c: type I error, with independent Xi, continuous outcome
#sample size=50

#initialization
setwd("/hpc/home/lz197/RFomnibus/")#to modify according to computer
rm(list=ls())
source("randomForestTest_parallel_omnibus2_updated9.R")
source("getDescendants.R")
library(GUniFrac)
library(ranger)
library(ape)
library(vegan)
library(parallel)

n <- 50 #number of observations
iter <- 1000 #iterations, number of simulations

method=c("wRF","uwRF", "Omnibus.RF")
beta0= 10
pv.mat <- array(NA, 
                dim=c(iter, length(method)),
                dimnames = list(NULL, method)
                )#iter, method; indexed by h,i

load("/hpc/home/lz197/RFomnibus/adenoma_Li.RData")

#to generate the sample distribution and simulate a trait
set.seed(seed=23456)
data.obj$otu.tab->otu.tab
otu.tab<-t(otu.tab)
otu.tab[rowSums(otu.tab)>=20000, ]->data #sequencing depth >=20000
Rarefy(data)$otu.tab->data
data[,colSums(data)!=0]->data #439*2100 after rarefication
data.obj$tree->tree.rooted

set.seed(seed=2345610)
#for different methods, using the same x.index, x and phylogenetic tree
#for different iter, using different x and trees
x.index<-matrix(NA, nrow=n, ncol=iter)#sample size, iter; indexed by g, h
for(h in 1:iter) x.index[,h]<-sample(1:nrow(data),n)

#for different samples, iter, using different y
#for different methods, using the same y
y<-array(NA, dim=c(n, iter),
         dimnames = list(NULL, NULL))#sample size, iter; indexed by g,h
z <- array(NA, dim=c(n, iter),
         dimnames = list(NULL, NULL))#sample size, iter; indexed by g,h

for(h in 1:iter){
  z[,h] <- scale(rnorm(n)) #independent covariate z1
  #generate y
  temp <- rnorm(n, sd=3)
  y[ ,h] <- beta0 + z[,h] + temp    #to simulate a trait
}

#on computer cluster
cores = 10 #to modify according to computer

for(h in 1:iter){
  if(h %% round(iter/100) == 0) cat(".")
  
  #generate x and tree which are determined by the x.index
  x <- data[x.index[,h],]
  x <- x[,colSums(x)>=1]  #50*1506
  droptips<-tree.rooted$tip.label[!(tree.rooted$tip.label %in% colnames(x))]
  tree <- drop.tip(phy = tree.rooted, tip = droptips) #a real tree
  if(sum(!(colnames(x) %in% tree$tip.label))==0) 
    x<-x[,tree$tip.label] else stop("don't match") #order x according to the tree
	
  meta<- data.frame(y=y[,h],z=z[,h])
  
  pv1 <- randomForestTest_parallel_omnibus2(comm = x , meta.data = meta, tree=tree, response.var="y", 
                  adjust.vars = "z", perm.no = 999, n.cores=cores, test.type = "Training",method = "o",
				  prediction.type = 'Regression')
  pv1$p.value.perm[c("weighted", "unweighted","Omnibus")]->pv.mat[h,c("wRF", "uwRF","Omnibus.RF")]
  
  if(h %% round(iter/100) == 0) save.image(file="/hpc/home/lz197/RFomnibus/H0_i_c_r4.RData")
}
save.image(file="/hpc/home/lz197/RFomnibus/H0_i_c_r4.RData")