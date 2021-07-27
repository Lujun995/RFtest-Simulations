#Scenario 3: different linking functions
#beta= 0, 0.25, 0.5, 0.75, 1
#signal density=15%
#four methods= RFtest, MiRKAT, OMiAT, aMiSPU
#two signal types = non-phylogenetically clustered OTUs and phylogenetically clustered OTUs
#linking function = log2(x+1) and presence/absence

#initialization
setwd("/hpc/home/lz197/RFomnibus/")#to modify according to computer
rm(list=ls())
source("randomForestTest_parallel_omnibus2.R")
source("getDescendants.R")
library(GUniFrac)
library(ranger)
library(ape)
library(vegan)
library(parallel)
library(ecodist)

n <- 50 #number of observations
iter <- 1000 #iterations, number of simulations

method=c("wRF","uwRF", "Omnibus.RF", "Omnibus.MiRKAT", "OMiAT", "MiHC", "aMiSPU")
beta0= 10
beta = c(0, 0.25, 0.5, 0.75, 1)
signal<-c("non-phylo","phylo")
Density=0.15
linking.method=c("log2","P/A")
linking <- function(x, method=c("log2","P/A")) 
  switch(method, 'log2'=log2(x+1), 'P/A' = (1*(x!=0)) )
pv.mat <- array(NA, 
                dim=c(iter, length(method), length(beta), length(signal), length(linking.method)),
                dimnames = list(NULL, method, as.character(beta),signal, linking.method)
                )#iter, method, beta, signal; indexed by h,i,j,k,l

load("/hpc/home/lz197/RFomnibus/adenoma_Li.RData")

#to generate the sample distribution and simulate a trait
set.seed(seed=23456)
data.obj$otu.tab->otu.tab
otu.tab<-t(otu.tab)
otu.tab[rowSums(otu.tab)>=20000, ]->data #sequencing depth >=20000
Rarefy(data)$otu.tab->data
data[,colSums(data)!=0]->data #439*2100 after rarefication
data.obj$tree->tree.rooted
#to locate the lineage A comprising the following tips
lineageA<-getDescendants.tip(tree.rooted,node = 4190) #a abundant lineage
lineageA<-tree.rooted$tip.label[lineageA] #having ~15% OTUs in x, ~21.1% abundance of total
#length(lineageA)/length(tree.rooted$tip.label)
# droptips<-tree.rooted$tip.label[!(tree.rooted$tip.label %in% colnames(data))]
# setdiff(lineageA,droptips)->tempA
# sum(colSums(data[,tempA]))/sum(colSums(data))

set.seed(seed=33456) #not equal to s3_b
#for different method, beta, signal, linking functions using the same x.index, x and phylogenetic tree
#for different iter, using different x and trees
x.index<-matrix(NA, nrow=n, ncol=iter)#sample size, iter; indexed by g, h
for(h in 1:iter) x.index[,h]<-sample(1:nrow(data),n)

#for different samples, iter, beta, signal, Density, using different y
#for different methods, using the same y
y<-array(NA, 
         dim=c(n, iter, length(beta), length(signal), length(linking.method)),
         dimnames = list(NULL, NULL, as.character(beta),signal, linking.method)
)#sample size, iter, beta, signal, linking.method; indexed by g,h,j,k,l

for(h in 1:iter){
  #generate x which are determined by the x.index
  #the x.index is determined by h
  x <- data[x.index[,h],]
  x <- x[,colSums(x)>=1]  #50*1506
  droptips<-tree.rooted$tip.label[!(tree.rooted$tip.label %in% colnames(x))]
  #generate y
  #y varies according to different samples, iter, beta, signal, linking.method
  for(k in 1:length(signal)){ #for each iter, this loop is excuted twice, however, each "if" is excuted only once 
    if(signal[k]=="non-phylo"){
      OTU.list<-sample(colnames(x)) #randomize names of x for each iter
	  signal.tips<- OTU.list[1:round(ncol(x)*Density)]
    }else if(signal[k]=="phylo"){
	  if(abs(Density-0.15)<1e-5)
		setdiff(lineageA,droptips)->signal.tips else #the tips in the new tree in lineage A in each iter
	  stop("Density is not 0.15")
    }else stop("signal is not non-phylo or phylo")
	
	for(l in 1:length(linking.method)){
	#for each linking.method, signal tips are the same
	  x.trans<-linking(x, linking.method[l])
	  signal.strength <- rowSums(x.trans[,signal.tips]) #signal.strength varies according to linking functions
	  q=beta0+matrix(beta, ncol=1) %*% matrix(scale(signal.strength), nrow=1)
	  #for different beta, the scaled signal strength are the same
	  rownames(q)<-as.character(beta);colnames(q)<-names(signal.strength)
	  temp<-rnorm(n)
	  y[ ,h, ,k,l] <- apply(X=q, MARGIN=1, FUN="+", temp) #to simulate a trait
	}
  }
}



save.image(file="/hpc/home/lz197/RFomnibus/s3_c.RData")


#on computer cluster
cores = 37 #to modify according to computer

#set.seed(seed=3456122)
for(h in 1:iter){
  if(h %% round(iter/100) == 0) cat(".")
  #generate x and tree which are determined by the x.index
  #the x.index is determined by h
  x <- data[x.index[,h],]
  x <- x[,colSums(x)>=1]  #50*1506
  droptips<-tree.rooted$tip.label[!(tree.rooted$tip.label %in% colnames(x))]
  tree <- drop.tip(phy = tree.rooted, tip = droptips) #a real tree
  if(sum(!(colnames(x) %in% tree$tip.label))==0) 
    x<-x[,tree$tip.label] else stop("don't match") #order x according to the tree

#   #for MiRKAT
#   unifracs = GUniFrac(x, tree=tree, alpha = c(1))$unifracs
#   Ds = list(w = unifracs[,,"d_1"], uw = unifracs[,,"d_UW"],
# 			BC= as.matrix(vegdist(x, method="bray")))
#   Ks = lapply(Ds, FUN = function(d) MiRKAT::D2K(d))
  
  for(k in 1:length(signal)){
	  for(l in 1:length(linking.method)){
		for(j in 1:length(beta)){
			#method[i]=c("wRF","uwRF", "Omnibus.RF", "Omnibus.MiRKAT", "OMiAT", "MiHC", "aMiSPU")
			#obtain y=y[ ,h,j,k,l]
			#RF omnibus test
			pv1 <- randomForestTest_parallel_omnibus2(comm = x , meta.data = data.frame(y=y[ ,h,j,k,l]), tree=tree,
                  response.var="y", perm.no = 999, n.cores=cores, test.type = "Training", method = "o", 
				  prediction.type = 'Regression')
			pv1$p.value.perm[c("weighted", "unweighted","Omnibus")]->pv.mat[h,c("wRF", "uwRF","Omnibus.RF"),j,k,l]			
			# #MiRKAT
			# pv2 <- MiRKAT::MiRKAT(y[ ,h,j,k,l], Ks = Ks, out_type="C",nperm=999)
			# pv2$omnibus_p -> pv.mat[h,c("Omnibus.MiRKAT"),j,k,l]
			# #OMiAT
			# pv3 <- OMiAT::OMiAT(Y=y[ ,h,j,k,l], otu.tab=x, tree=tree, model="continous",n.perm=999)
			# pv3$OMiAT.pvalue -> pv.mat[h,c("OMiAT"),j,k,l]
			# #aMiSPU
			# pv4 <- MiSPU::MiSPU(y = y[ ,h,j,k,l], X= as.matrix(x), tree=tree, model="gaussian", n.perm = 999)
			# pv4$aMiSPU$pvalue -> pv.mat[h,c("aMiSPU"),j,k,l]
			# #MiHC
			# pv5 <- MiHC::MiHC(y=y[ ,h,j,k,l], otu.tab=x, tree = tree, model="gaussian", n.perm=999)
			# pv5$ada.pvs["MiHC"] -> pv.mat[h,c("MiHC"),j,k,l]
		}
	  }
  }
  if(h %% round(iter/100) == 0) save.image(file="/hpc/home/lz197/RFomnibus/s3_c.RData")
}
save.image(file="/hpc/home/lz197/RFomnibus/s3_c.RData")