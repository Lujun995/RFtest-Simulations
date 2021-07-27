#Scenario 1: three testing methods
#beta=0, 0.75, 1.5, 2.25, 3
#n=50
#signal density=5%, 15%
#three methods=RFtest, MiRKAT, OMiAT
#three signal types=rare OTUs, abundant OTUs and phylogenetically clustered OTUs

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
beta0= 0
beta = c(0, 0.75, 1.5, 2.25, 3)
signal<-c("non-phylo","phylo")
Density=c(0.05, 0.15)
pv.mat <- array(NA, 
                dim=c(iter, length(method), length(beta), length(signal), length(Density)),
                dimnames = list(NULL, method, as.character(beta),signal, as.character(Density))
                )#iter, method, beta, signal, Density; indexed by h,i,j,k,l

load("/hpc/home/lz197/RFomnibus/adenoma_Li.RData")

#to generate the sample distribution and simulate a trait
set.seed(seed=23456)
data.obj$otu.tab->otu.tab
otu.tab<-t(otu.tab)
otu.tab[rowSums(otu.tab)>=20000, ]->data #sequencing depth >=20000
Rarefy(data)$otu.tab->data
data[,colSums(data)!=0]->data #439*2100 after rarefication
data.obj$tree->tree.rooted
colMeans(data)-> OTU.mean
#to locate the lineage (lineage A, C) comprising the following tips
lineageA<-getDescendants.tip(tree.rooted,node = 4190) #a abundant lineage
lineageA<-tree.rooted$tip.label[lineageA] #having ~15% OTUs in x, ~21.1% abundance of total
#length(lineageA)/length(tree.rooted$tip.label)
# droptips<-tree.rooted$tip.label[!(tree.rooted$tip.label %in% colnames(data))]
# setdiff(lineageA,droptips)->tempA
# sum(colSums(data[,tempA]))/sum(colSums(data))
lineageC<-getDescendants.tip(tree.rooted,node = 4300) #lineageC is a lineage under lineageA
lineageC<-tree.rooted$tip.label[lineageC] #having ~5% OTUs in x, ~10.9% abundance of total
# droptips<-tree.rooted$tip.label[!(tree.rooted$tip.label %in% colnames(data))]
# setdiff(lineageC,droptips)->tempC
# sum(colSums(data[,tempC]))/sum(colSums(data))
# length(lineageC)/length(tree.rooted$tip.label)

set.seed(seed=23456)
#for different method, beta, signal, Density, using the same x.index, x and phylogenetic tree
#for different iter, using different x and trees
x.index<-matrix(NA, nrow=n, ncol=iter)#sample size, iter; indexed by g, h
for(h in 1:iter) x.index[,h]<-sample(1:nrow(data),n)

#for different samples, iter, beta, signal, Density, using different y
#for different methods, using the same y
y<-array(NA, 
         dim=c(n, iter, length(beta), length(signal), length(Density)),
         dimnames = list(NULL, NULL, as.character(beta),signal, as.character(Density))
)#sample size, iter, beta, signal, Density; indexed by g,h,j,k,l

for(h in 1:iter){
  #generate x which are determined by the x.index
  #the x.index is determined by h
  x <- data[x.index[,h],]
  x <- x[,colSums(x)>=1]  #50*1506
  droptips<-tree.rooted$tip.label[!(tree.rooted$tip.label %in% colnames(x))]
  #generate y
  #y varies according to different samples, iter, beta, signal
  for(k in 1:length(signal)){
	signal.tips<-vector(mode="list", length=length(Density))
	names(signal.tips)<-paste("Density=",c(0.05,0.15),sep="")
    if(signal[k]=="non-phylo"){
      OTU.list<-sample(colnames(x)) #randomize names of x for each iter
	  for(l in 1:length(Density)){
		#for different Density, use the same order of OTU.list
		#therefore, 1:round(ncol(x)*0.15) > 1:round(ncol(x)*0.05)
		signal.tips[[l]]<- OTU.list[1:round(ncol(x)*Density[l])]
	  }
    }else if(signal[k]=="phylo"){
	  for(l in 1:length(Density)){
		if(abs(Density[l]-0.05)<1e-5) 
			setdiff(lineageC,droptips)->signal.tips[[l]] else #the tips in the new tree in lineage C in each simulation
		if(abs(Density[l]-0.15)<1e-5)
			setdiff(lineageA,droptips)->signal.tips[[l]] else
		stop("Density is not 0.05 or 0.15")
	  }
    }else stop("signal is not non-phylo or phylo")
	for(l in 1:length(Density)){
	  if(signal[k]=="non-phylo") #Zij/mean(Zj) to avoid a single or a few OTUs dominating the total effect.
		signal.strength<-rowSums(t(t(x[,signal.tips[[l]]]) / colMeans(x[,signal.tips[[l]]]))) else # matrix/vector element-wise
	  if(signal[k]=="phylo")
	    signal.strength<-rowSums(x[,signal.tips[[l]]])
	  q=beta0+matrix(beta, ncol=1) %*% matrix(scale(signal.strength), nrow=1)#for each beta and sample, signal tips are the same
      rownames(q)<-as.character(beta);colnames(q)<-names(signal.strength)
      p=exp(q)/(1+exp(q)) #inverse logit function
      y[ ,h, ,k,l] <- apply(X=p, MARGIN=1, FUN=rbinom, n=n, size=1) #to simulate a trait
	  #NOTE on rbinom: rbinom(n, size, prob), if n < length(prob), it will generate the first nth binomial random variable 
	  #according to p; if n > length(prob), it will recycle prob; if n = length(prob), it will give the binomial sample according
	  #to the probability given by prob. Therefore, as the row of p is the probability of one sample under a specific beta,
	  #this apply(X=p, MARGIN=1, FUN=rbinom, n=n, size=1) should generate the desired binomial random sample (with different
      #probability).
	  #temp<-sapply(1:1e4, function(i) rbinom(n=5, size=1, prob=c(0.1,0.2,0.4,0.6,0.8)))
	  #apply(temp,MARGIN=1,sum)/1e4
	}
  }
}


# tempA<-tempB<-tempC<-temp2<-rep(NA,iter);for(h in 1:iter) {
#   x <- data[x.index[,h],]
#   x <- x[,colSums(x)>=1]  #50*1506
#   droptips<-tree.rooted$tip.label[!(tree.rooted$tip.label %in% colnames(x))]
#   length(setdiff(lineageA,droptips))/ncol(x)->tempA[h]
#   length(setdiff(lineageB,droptips))/ncol(x)->tempB[h]
#   length(setdiff(lineageC,droptips))/ncol(x)->tempC[h]
#   length(union(setdiff(lineageA,droptips),setdiff(lineageB,droptips)))/ncol(x)->temp2[h]
# }
# mean(tempA)
# mean(tempB)
# mean(tempC)
# mean(temp2)

save.image(file="/hpc/home/lz197/RFomnibus/s1_b_run3.RData")




#on computer cluster
cores = 37 #to modify according to computer

#set.seed(seed=234561)
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
	
  # #For MiRKAT
  # unifracs = GUniFrac(x, tree=tree, alpha = c(1))$unifracs
  # Ds = list(w = unifracs[,,"d_1"], uw = unifracs[,,"d_UW"], 
			# BC= as.matrix(vegdist(x, method="bray")))
  # Ks = lapply(Ds, FUN = function(d) MiRKAT::D2K(d))
  
  for(k in 1:length(signal)){
	for(l in 1:length(Density)){
		for(j in 1:length(beta)){   #for(j in 1:length(beta)){ #beta=0, 0.75, 1.5, 2.25, 3
			#method[i]=c("wRF","uwRF", "Omnibus.RF", "Omnibus.MiRKAT", "OMiAT", "MiHC", "aMiSPU")
			#obtain y=y[ ,h,j,k,l]
			#RF omnibus test
			pv1 <- randomForestTest_parallel_omnibus2(comm = x , meta.data = data.frame(y=factor(y[ ,h,j,k,l])), tree=tree, 
                  response.var="y", perm.no = 999, n.cores=cores, test.type = "Training",method = "o")
			pv1$p.value.perm[c("weighted", "unweighted","Omnibus")]->pv.mat[h,c("wRF", "uwRF","Omnibus.RF"),j,k,l]
			# #MiRKAT
			# pv2 <- MiRKAT::MiRKAT(y[ ,h,j,k,l], Ks = Ks, out_type="D",nperm=999)
			# pv2$omnibus_p -> pv.mat[h,c("Omnibus.MiRKAT"),j,k,l]
			# #OMiAT
			# pv3 <- OMiAT::OMiAT(Y=y[ ,h,j,k,l], otu.tab=x, tree=tree, model="binomial",n.perm=999)
			# pv3$OMiAT.pvalue -> pv.mat[h,c("OMiAT"),j,k,l]
			# #aMiSPU
			# pv4 <- MiSPU::MiSPU(y = y[ ,h,j,k,l], X= as.matrix(x), tree=tree, model="binomial", n.perm = 999)
			# pv4$aMiSPU$pvalue -> pv.mat[h,c("aMiSPU"),j,k,l]
			# #MiHC
			# pv5 <- MiHC::MiHC(y=y[ ,h,j,k,l], otu.tab=x, tree = tree, model="binomial", n.perm=999)
			# pv5$ada.pvs["MiHC"] -> pv.mat[h,c("MiHC"),j,k,l]
		}
	}
  }
  if(h %% round(iter/100) == 0) save.image(file="/hpc/home/lz197/RFomnibus/s1_b_run3.RData")
}
save.image(file="/hpc/home/lz197/RFomnibus/s1_b_run3.RData")