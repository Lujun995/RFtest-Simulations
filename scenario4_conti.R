#Scenario 4: different lineages and different linking function, continuous outcomes
#beta = 0.35
#n = 50
#signal density = 10~15%
#five methods = RFtest, MiRKAT, OMiAT, aMiSPU
#signal type = phylogenetically clustered OTUs
#seven different lineages

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
beta = 0.35
signal<-"phylo"
Density= 0.15
lineage.nodes = as.character(c(4190,2514,3590,2519,2734,3399,3171))
linking.method=c("identity","log2","P/A")
linking <- function(x, method=c('identity','log2','P/A')) 
  switch(method, 'identity'=x, 'log2'=log2(x+1), 'P/A' = (1*(x!=0)) )
pv.mat <- array(NA, 
                dim=c(iter, length(method), length(linking.method), length(lineage.nodes)),
                dimnames = list(NULL, method, linking.method, lineage.nodes)
                )#iter, method, linking, lineag.nodes; indexed by h,i,k,l

load("/hpc/home/lz197/RFomnibus/adenoma_Li.RData")

#to generate the sample distribution and simulate a trait
set.seed(seed=23456)
data.obj$otu.tab->otu.tab
otu.tab<-t(otu.tab)
otu.tab[rowSums(otu.tab)>=20000, ]->data #sequencing depth >=20000
Rarefy(data)$otu.tab->data
data[,colSums(data)!=0]->data #439*2100 after rarefication
data.obj$tree->tree.rooted
#to locate the lineages under the nodes in the lineage.nodes
lineage.list<-vector(length=length(lineage.nodes), mode="list")
names(lineage.list)<-lineage.nodes
for(l in 1:length(lineage.nodes)){
  temp<-getDescendants.tip(tree.rooted,node = as.numeric(lineage.nodes[l])) #get the tip number of these nodes
  lineage.list[[l]]<-tree.rooted$tip.label[temp]
}

# #the total abundance of OTUs in these lineages
# droptips<-tree.rooted$tip.label[!(tree.rooted$tip.label %in% colnames(data))]
# for(l in 1:length(lineage.list)) 
#   cat(round(sum(colSums(data[,setdiff(lineage.list[[l]],droptips)]))/sum(colSums(data))*100,2)," ")
# #ranging from 1% to 40%, comprising more than 80% of total abundance

# #check if there is an intersection between lineages
# temp=rep(NA,length=choose(length(lineage.list), 2))
# l=1
# for(l1 in 1:(length(lineage.list)-1)){
  # for(l2 in (l1+1):length(lineage.list)) {
    # temp[l]<-length(intersect(lineage.list[[l1]],lineage.list[[l2]]))
	# l=l+1
  # }
# }
# if(sum(temp!=0)!=0) 
  # stop("There is some common OTUs between two lineages") #no intersection between every pair of lineages


set.seed(seed=4456)
#for different method, beta, signal, Density, using the same x.index, x and phylogenetic tree
#for different iter, using different x and trees
x.index<-matrix(NA, nrow=n, ncol=iter)#sample size, iter; indexed by g, h
for(h in 1:iter) x.index[,h]<-sample(1:nrow(data),n)

#for different samples, iter, beta, signal, Density, using different y
#for different methods, using the same y
y<-array(NA, 
         dim=c(n, iter, length(linking.method), length(lineage.nodes)),
         dimnames = list(NULL, NULL, linking.method, lineage.nodes)
)#sample size, iter, linking.method, lineage.nodes; indexed by g,h,k,l

for(h in 1:iter){
  #generate x which are determined by the x.index
  #the x.index is determined by h
  x <- data[x.index[,h],]
  x <- x[,colSums(x)>=1]  #50*1506
  droptips<-tree.rooted$tip.label[!(tree.rooted$tip.label %in% colnames(x))]
  
  #generate y
  #y varies according to different samples, iter, beta, signal
  #signal.tip, determined by lineage.list and droptips, is the same for each linking function
  signal.tips<-lapply(lineage.list, FUN= setdiff, y=droptips) 
  names(signal.tips)<-paste("Node ",lineage.nodes,sep="")
  for(k in 1: length(linking.method)){
    x.trans<-linking(x, linking.method[k])
	  signal.strength <- #signal.strength varies according to different linking functions
	    sapply(1:length(signal.tips),  function(l) rowSums(x.trans[,signal.tips[[l]]])) #a matrix of g-by-l
	  colnames(signal.strength) <- names(signal.tips)
      q = beta0+ beta * apply(signal.strength, MARGIN=2, FUN= scale) #a matrix of g-by-l
	  rownames(q)<-rownames(signal.strength);colnames(q)<-names(signal.tips)
	  temp<-rnorm(n)
	  y[ ,h,k, ]<- apply(X=q, MARGIN=2, FUN="+", temp) #to simulate a trait
  }
}

#the proportion of OTUs in each lineage
# temp<-matrix(NA, nrow=iter, ncol=length(lineage.list))
# colnames(temp) <- names(lineage.list)
# for(h in 1:iter) {
#   x <- data[x.index[,h],]
#   x <- x[,colSums(x)>=1]  #50*1506
#   droptips<-tree.rooted$tip.label[!(tree.rooted$tip.label %in% colnames(x))]
#   for(l in 1: length(lineage.list))
# 	  length(setdiff(lineage.list[[l]],droptips))/ncol(x)->temp[h,l]
# }
# apply(temp,MARGIN=2,mean) #comprising 5%~20% of OTUs in x, a total of 80%

save.image(file="/hpc/home/lz197/RFomnibus/s4_c.RData")




#on computer cluster
cores = 37 #to modify according to computer

#set.seed(seed=4561233)
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
  
  for(k in 1){ #for(k in 1:length(linking.method)){
	for(l in 1:length(lineage.nodes)){
		#method[i]=c("wRF","uwRF", "Omnibus.RF", "Omnibus.MiRKAT", "OMiAT", "MiHC", "aMiSPU")
		#y: sample size, iter, linking.method, lineage.nodes; indexed by g,h,k,l
		#obtain y=y[ ,h,k,l]
		#pv.mat:iter, method, linking, lineag.nodes; indexed by h,i,k,l
		#RF omnibus test
		pv1 <- randomForestTest_parallel_omnibus2(comm = x , meta.data = data.frame(y=y[ ,h,k,l]), tree=tree, 
               response.var="y", perm.no = 999, n.cores=cores, test.type = "Training", method = "o",
               prediction.type = 'Regression')
		pv1$p.value.perm[c("weighted", "unweighted","Omnibus")]->pv.mat[h,c("wRF", "uwRF","Omnibus.RF"),k,l]
		# #MiRKAT
		# pv2 <- MiRKAT::MiRKAT(y[ ,h,k,l], Ks = Ks, out_type="C",nperm=999)
		# pv2$omnibus_p -> pv.mat[h,c("Omnibus.MiRKAT"),k,l]
		# #OMiAT
		# pv3 <- OMiAT::OMiAT(Y=y[ ,h,k,l], otu.tab=x, tree=tree, model="continous",n.perm=999)
		# pv3$OMiAT.pvalue -> pv.mat[h,c("OMiAT"),k,l]
		# #aMiSPU
		# pv4 <- MiSPU::MiSPU(y = y[ ,h,k,l], X= as.matrix(x), tree=tree, model="gaussian", n.perm = 999)
		# pv4$aMiSPU$pvalue -> pv.mat[h,c("aMiSPU"),k,l]
		# #MiHC
		# pv5 <- MiHC::MiHC(y=y[ ,h,k,l], otu.tab=x, tree = tree, model="gaussian", n.perm=999)
		# pv5$ada.pvs["MiHC"] -> pv.mat[h,c("MiHC"),k,l]
	}
  }
  if(h %% round(iter/100) == 0) save.image(file="/hpc/home/lz197/RFomnibus/s4_c.RData")
}
save.image(file="/hpc/home/lz197/RFomnibus/s4_c.RData")