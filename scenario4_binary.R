#Source code for the scenario 4 (clear): interaction
#beta=5
#n=50, 100, 150, 200, 250
#signal density= 13% * 15%
#three methods=RFtest, MiRKAT, OMiAT, aMiSPU
#two signal types=non-clustered random signal and phylogenetically clustered signal

#initialization
setwd("/hpc/home/lz197/RFomnibus/")#to modify according to computer
rm(list=ls())
source("randomForestTest.R")
source("getDescendants.R")
library(GUniFrac)
library(ranger)
library(ape)
library(vegan)
library(parallel)
library(ecodist)

n = c(50, 100, 150, 200, 250) #number of observations
iter <- 1000 #iterations, number of simulations

method=c("wRF","uwRF", "Omnibus.RF", "Omnibus.MiRKAT", "OMiAT", "aMiSPU")
beta0= 0
beta1 = 5 #for interaction
beta2 = 1.5 #for synergistic
signal<-c("non-phylo","phylo")
Density = 0.15 #the second set of signaling OTUs
Switch = 0.13 #the first set set of signaling OTUs

inter = c("interaction", "synergistic")
pv.mat <- array(NA, 
                dim=c(iter, length(method), length(n), length(signal), length(inter)),
                dimnames = list(NULL, method, as.character(n),signal, inter)
                )#iter, method, n, signal, inter; indexed by h,i,j,k,l. g=1:250 or j=1:5

load("/hpc/home/lz197/RFomnibus/adenoma_Li.RData")

#to generate the sampling distribution and simulate a trait
set.seed(seed=23456)
data.obj$otu.tab->otu.tab
otu.tab<-t(otu.tab)
otu.tab[rowSums(otu.tab)>=20000, ]->data #sequencing depth >=20000
Rarefy(data)$otu.tab->data
data[,colSums(data)!=0]->data #439*2100 after rarefication
data.obj$tree->tree.rooted
#to locate the lineage (lineage A, S) comprising the following tips
lineageA<-getDescendants.tip(tree.rooted,node = 4190) #a abundant lineage
lineageA<-tree.rooted$tip.label[lineageA] #having ~15% OTUs in x, ~21.1% abundance of total
#length(lineageA)/length(tree.rooted$tip.label)
# droptips<-tree.rooted$tip.label[!(tree.rooted$tip.label %in% colnames(data))]
# setdiff(lineageA,droptips)->tempA
# sum(colSums(data[,tempA]))/sum(colSums(data))
lineageS<-getDescendants.tip(tree.rooted,node = 2514) #switch lineage, disjoint with lineageA
lineageS<-tree.rooted$tip.label[lineageS] #having ~13% OTUs in x, ~12% abundance of total
# intersect(lineageS, lineageA)
# droptips<-tree.rooted$tip.label[!(tree.rooted$tip.label %in% colnames(data))]
# setdiff(lineageS,droptips)->tempS
# sum(colSums(data[,tempS]))/sum(colSums(data)) #the average abundance of lineageS in the population
# length(lineageS)/length(tree.rooted$tip.label) #the proportion of OTUs in S to population
# #not proportion of OTUs in lineageS to x because rare OTU may be singled out when resampled from the population

set.seed(seed=6)
#for different method, sample size(n), signal, inter, using the same x.index, x and phylogenetic tree
#for different iter, using different x and trees
x.index<-matrix(NA, nrow=n[length(n)], ncol=iter)#sample size, iter; indexed by g, h.
#Using the largest n to generate the sampling distribution. For a smaller n, use the first 1:n[j] observations
#Therefore, 1:50 is in 1:100, and 1:100 is in 1:150 ...
for(h in 1:iter) x.index[,h]<-sample(1:nrow(data),n[length(n)])

# tempS<-integer(iter)
# for(h in 1:iter) {
  # x <- data[x.index[,h],]
  # x <- x[,colSums(x)>=1]  #50*1506
  # droptips<-tree.rooted$tip.label[!(tree.rooted$tip.label %in% colnames(x))]
  # length(setdiff(lineageS,droptips))/ncol(x)->tempS[h]
# }
# mean(tempS) #having ~2% OTUs in x


#for different samples, iter, signal, inter, using different y
#for different methods, using the same y. For different sample sizes, use the first several records.
#in a sample with a size other than 250, the first, for exmaple, 50 will have a valid record while the 
#following numbers would be NAs. 
y<-array(NA, 
         dim=c(n[length(n)], iter, length(n), length(signal), length(inter)),
         dimnames = list(NULL, NULL, as.character(n), signal, inter)
)#sample size, iter, n, signal, inter; indexed by g,h,j,k,l. g=1:250 or j=1:5

for(h in 1:iter){
  #generate x which are determined by the x.index
  #the x.index is determined by h
  x.all <- data[x.index[,h],] #retrieve all (250) x.index to generate signal.tips
  x.all <- x.all[,colSums(x.all)>=1]  #250*2000
  #droptips<-tree.rooted$tip.label[!(tree.rooted$tip.label %in% colnames(x))]
  #generate y
  #y varies according to different samples, iter, n, beta, signal, inter
  for(k in 1:length(signal)){
    signal.tips <- vector(mode = "list", length=2)
    names(signal.tips) <- c("switch", "density")
    if(signal[k]=="non-phylo"){
      OTU.list<-sample(colnames(x.all)) #randomize names of x for each iter
      signal.tips[[1]] <- OTU.list[ 1                             :  round(ncol(x.all)*Switch)                         ]
      signal.tips[[2]] <- OTU.list[ (round(ncol(x.all)*Switch)+1) : (round(ncol(x.all)*Switch)+round(ncol(x.all)*Density)) ]
      #the first 2% and the following 15%
    }else if(signal[k]=="phylo"){
      if(abs(Switch-0.13)<1e-5)  signal.tips[[1]] <- lineageS else stop("Switch is not 0.13")
      if(abs(Density-0.15)<1e-5) signal.tips[[2]] <- lineageA else stop("Density is not 0.15")
    }else stop("signal is not non-phylo or phylo")
    for(j in 1:length(n)){
      x <- data[x.index[1:n[j],h],] #retrieve the first 50 or 100 records x.index to generate y
      x <- x[,colSums(x)>=1] #50*1511
      droptips<-tree.rooted$tip.label[!(tree.rooted$tip.label %in% colnames(x))]
      lapply(signal.tips, setdiff,  y = droptips) -> signal.tips.filtered
      
      signal.strength <- vector(mode = "list", length=2)
      names(signal.strength) <- c("switch", "density")
      for(ii in 1:2){
        if(signal[k]=="non-phylo"){
          signal.strength[[ii]] <- 
            rowSums(t(t(x[,signal.tips.filtered[[ii]]]) / colMeans(x[,signal.tips.filtered[[ii]]])))
          # matrix / vector element-wise by col
          # Zij/mean(Zj) to avoid a single or a few OTUs dominating the total effect.
        }else if(signal[k]=="phylo")
          signal.strength[[ii]]<-rowSums(x[,signal.tips.filtered[[ii]]])
      }
      for(l in 1:length(inter)){
        if(inter[l]=="interaction"){
          q=beta0 + beta1 * scale(signal.strength[[1]])*scale(signal.strength[[2]]) 
          # for each beta and sample, signal tips are the same
          # center the signal strength at 0 to reduce the marginal effect
          # because signal.strength is not symmetrically distributed, q here may be still correlated with signal.strength?
          # cor(x=q[2,],y=signal.strength[[2]],method="s")
          p=exp(q)/(1+exp(q)) #inverse logit function
          y[1:n[j],h,j,k,l] <- rbinom(n=n[j], size= 1, prob = p[,1])
          #NOTE on rbinom: rbinom(n, size, prob), if n < length(prob), it will generate the first nth binomial random variable 
          #according to p; if n > length(prob), it will recycle prob; if n = length(prob), it will give the binomial sample according
          #to the probability given by prob. Therefore, as the row of p is the probability of one sample under a specific beta,
          #this apply(X=p, MARGIN=1, FUN=rbinom, n=n, size=1) should generate the desired binomial random sample (with different
          #probability).
          #temp<-sapply(1:1e4, function(i) rbinom(n=5, size=1, prob=c(0.1,0.2,0.4,0.6,0.8)))
          #apply(temp,MARGIN=1,sum)/1e4
        }else if(inter[l]=="synergistic"){
          lapply(signal.strength, beta0=beta0, beta2=beta2, FUN = function(.x, beta0, beta2){
            q = beta0 + beta2 * scale(.x)
            p = exp(q)/(1+exp(q)) #inverse logit function
            y = rbinom(n=length(.x), size= 1, prob = p[,1])
          }) -> temp
          y[1:n[j],h,j,k,l] <- temp[[1]]*temp[[2]]
        }
      }
    }
  }
}


#on computer cluster
cores =10 #to modify according to computer

#set.seed(seed=612345)
for(h in 1:iter){
  if(h %% round(iter/100) == 0) cat(".")
  for(j in 1:length(n)){
    #generate x and tree which are determined by the x.index and j (sample size)
    #the x.index is determined by h
    x <- data[x.index[1:n[j],h],] #the same x used to generate signal.tips.filtered, signal strength and y
    x <- x[,colSums(x)>=1]
    droptips<-tree.rooted$tip.label[!(tree.rooted$tip.label %in% colnames(x))]
    tree <- drop.tip(phy = tree.rooted, tip = droptips) #a real tree
    if(sum(!(colnames(x) %in% tree$tip.label))==0) 
      x<-x[,tree$tip.label] else stop("don't match") #order x according to the tree
    
    #For MiRKAT
    unifracs = GUniFrac(x, tree=tree, alpha = c(1))$unifracs
    Ds = list(w = unifracs[,,"d_1"], uw = unifracs[,,"d_UW"], 
              BC= as.matrix(vegdist(x, method="bray")))
    Ks = lapply(Ds, FUN = function(d) MiRKAT::D2K(d))
    
    for(k in 1:length(signal)){
      for(l in 1){ #interaction only
        #method[i]=c("wRF","uwRF", "Omnibus.RF", "Omnibus.MiRKAT", "OMiAT", "MiHC", "aMiSPU")
        #obtain y=y[1:n[j],h,j,k,l] #sample size, iter, n, signal, inter; indexed by g,h,j,k,l
		    #the first 1:n[j] ys should be non-NAs
        #pv.mat = iter, method, n, signal, inter; indexed by h,i,j,k,l. g=1:250 or j=1:5
        #RF omnibus test
        pv1 <- randomForestTest(comm = x , meta.data = data.frame(y=factor(y[1:n[j],h,j,k,l])), tree=tree, 
                                                  response.var="y", perm.no = 999, n.cores=cores, test.type = "OOB",method = "w",presence.in=1)
        pv1$p.value.perm->pv.mat[h,c("wRF"),j,k,l]
        #MiRKAT
        pv2 <- MiRKAT::MiRKAT(y[1:n[j],h,j,k,l], Ks = Ks, out_type="D",nperm=999)
        pv2$omnibus_p -> pv.mat[h,c("Omnibus.MiRKAT"),j,k,l]
        #OMiAT
        pv3 <- OMiAT::OMiAT(Y=y[1:n[j],h,j,k,l], otu.tab=x, tree=tree, model="binomial",n.perm=999)
        pv3$OMiAT.pvalue -> pv.mat[h,c("OMiAT"),j,k,l]
        #aMiSPU
        pv4 <- MiSPU::MiSPU(y = y[1:n[j],h,j,k,l], X= as.matrix(x), tree=tree, model="binomial", n.perm = 999)
        pv4$aMiSPU$pvalue -> pv.mat[h,c("aMiSPU"),j,k,l]
      }
    }
  }
  if(h %% round(iter/100) == 0) save.image(file="/hpc/home/lz197/RFomnibus/s4_b.RData")
}
save.image(file="/hpc/home/lz197/RFomnibus/s4_b.RData")