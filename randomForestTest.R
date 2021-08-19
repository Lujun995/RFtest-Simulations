###############################################################################
# 2021/07/21 randomForestTest_parallel_omnibus2_updated9.R
# Author: Lujun Zhang, Jun Chen
###############################################################################
require(ranger)
require(ape)
require(vegan)
require(parallel)

# Calculate the abundance of the internal nodes
calculate.cum <- function (otu.tab, tree) {
  #normalization (relative abundance)
  if (!is.rooted(tree)) stop("Rooted phylogenetic tree required!")
  otu.tab <- as.matrix(otu.tab)
  row.sum <- rowSums(otu.tab)
  otu.tab <- otu.tab / row.sum
  n <- nrow(otu.tab)
  
  #exception handling
  if (sum(!(colnames(otu.tab) %in% tree$tip.label)) != 0) {
    stop("The OTU table contains unknown OTUs! OTU names\n\t\t\t\t\tin the OTU table and the tree should match!")
  }
  absent <- tree$tip.label[!(tree$tip.label %in% colnames(otu.tab))]
  if (length(absent) != 0) {
    tree <- drop.tip(tree, absent)
    warning("The tree has more OTU than the OTU table!")
  }
  
  #calculate the accumulative abundance of nodes in the phylo tree
  #the accumulative abundance = sum of the abundance of all the child tips
  tip.label <- tree$tip.label
  otu.tab <- otu.tab[, tip.label]
  ntip <- length(tip.label)
  nbr <- nrow(tree$edge) #number of edges
  edge <- tree$edge #tree$edge[,1] where the edge is from, tree$edge[,2] where the edge is to
  edge2 <- edge[, 2]
  br.len <- tree$edge.length
  cum <- matrix(0, nbr, n) #nbr edges X n OTUs
  for (i in 1:ntip) {
    tip.loc <- which(edge2 == i) #the node connecting tip/leaf/OTU i
    cum[tip.loc, ] <- cum[tip.loc, ] + otu.tab[, i] #add the abundance of tip i to the node
    node <- edge[tip.loc, 1] 
    node.loc <- which(edge2 == node) #the node connects the node which connects tip i
    while (length(node.loc)) {
      cum[node.loc, ] <- cum[node.loc, ] + otu.tab[, i] #add the abundance of tip i to the node
      node <- edge[node.loc, 1]
      node.loc <- which(edge2 == node) #the next node which connects the current node
    }
  }
  return(cum)
}

# Internal function
# Need to revise!  Prediction with probability
calculateError <- function (rf.obj, data, response.var, prediction.type, test.type, w = NULL) {
  yy.o <- data[, response.var]
  yy.b <- rf.obj$predictions
  
  if (test.type == 'OOB') {
    # Not clear how the prediction error is calcuated for multi-class proabilities and survival
    stat <- rf.obj$prediction.error
  }
  
  if (test.type == 'Training') {
    if (prediction.type %in% 'Classification') {
      yy.o <- model.matrix(~.- 1, data.frame(yy.o))
      yy.t <- predict(rf.obj, data = data)$prediction 
      
      if (is.factor(yy.t)) {
        yy.t <- model.matrix(~.- 1, data.frame(yy.t))
      }
      # Brier's score
      stat <- mean(rowSums((yy.o - yy.t)^2))
    }
    
    if (prediction.type %in% 'Regression') {
      yy.t <- predict(rf.obj, data = data)$prediction 
      stat <- mean((yy.o - yy.t)^2)
    }
  }
  
  if (test.type %in% c('0.632', '0.632+')) {
    if (test.type == '0.632') {
      w <- 0.632
    }
    
    if (prediction.type %in% 'Classification') {
      yy.o <- model.matrix(~.- 1, data.frame(yy.o))
      yy.t <- predict(rf.obj, data = data)$prediction 
      
      if (is.factor(yy.t)) {
        yy.t <- model.matrix(~.- 1, data.frame(yy.t))
      }
      if (is.factor(yy.b)) {
        yy.b <- model.matrix(~.- 1, data.frame(yy.b))
      }
      
      # Brier's score
      err.t <- mean(rowSums((yy.o - yy.t)^2))
      err.b <- mean(rowSums((yy.o - yy.b)^2))
      
      if (test.type == '0.632+' & is.null(w)) {
        err.n <- sum(sapply(1:ncol(yy.o), function (i) {
          mean(outer(yy.t[, i], yy.o[, i], function (x, y) (x - y)^2))
        }))
        
        ratio <- (err.b - err.t) / (err.n - err.t)
        if (ratio < 0) ratio <- 0
        if (ratio > 1) ratio <- 1
        w <- 0.632 / (1 - 0.368 * ratio)
      }
      
      stat <- (1 - w) * err.t + w * err.b
    }
    
    if (prediction.type %in% 'Regression') {
      yy.t <- predict(rf.obj, data = data)$prediction 
      err.t <- mean((yy.o - yy.t)^2)
      err.b <- mean((yy.o - yy.b)^2)
      
      if (test.type == '0.632+' & is.null(w)) {
        err.n <- mean(outer(yy.t, yy.o, function (x, y) (x - y)^2))
        ratio <- (err.b - err.t) / (err.n - err.t)
        if (ratio < 0) ratio <- 0
        if (ratio > 1) ratio <- 1
        w <- 0.632 / (1 - 0.368 * ratio)
      }
      
      stat <- (1 - w) * err.t + w * err.b
    }
  }
  
  
  return(list(stat = stat, w = w))
}

# This is the major function
randomForestTest <- 
  function (comm, meta.data, response.var, status.var = NULL, adjust.vars = NULL, 
            prediction.type = c('Classification', 'Regression', 'Survival'),
            method = c('weighted','omnibus','unweighted'), tree = NULL, probability = TRUE,
            presence.in = 2,  test.type = c('OOB', 'Training', '0.632', '0.632+'),
            perm.no = 999, pearson.fit = FALSE,  trace = FALSE, n.cores=detectCores()-1,
			reps = 3, ...) {
  call <- match.call()
  
  prediction.type <- match.arg(prediction.type)
  method <- match.arg(method)
  test.type <- match.arg(test.type)
  
  # Future revision to allow formula, data.frame, etc.
  # Survial prediction to be implemented later
  # Multi-class prediction to be further checked
  
  meta.data <- data.frame(meta.data)
  if (prediction.type %in% 'Classification') {
    meta.data[, response.var] <- factor(meta.data[, response.var])
  }
  
  if (prediction.type %in% c('Regression', 'Survival')) {
    meta.data[, response.var] <- as.numeric(meta.data[, response.var])
    probability <- FALSE
  }
  
  if(!is.null(adjust.vars)){ #fitting the model outside the loop to save time
    fml = as.formula(paste(response.var, paste(adjust.vars,collapse = "+"), sep = "~"))
    if(prediction.type %in% c('Regression')) {
      model =  lm(formula = fml, data = meta.data)
    } else 
      if(prediction.type %in% 'Classification') {
        model = glm(formula = fml, data = meta.data, family = "binomial")
      }
  }
  
  drop.otus <- NULL
  if(!is.null(presence.in)) {
    if(presence.in < 1) 
      drop.otus.number <- which( colSums(1*(comm !=0)) < as.numeric(presence.in)*nrow(comm) ) else
        drop.otus.number <- which( colSums(1*(comm !=0)) < as.numeric(presence.in) )
    if(length(drop.otus.number) != 0) {
      drop.otus <- names(drop.otus.number)
      comm<-comm[,-drop.otus.number]
    }
  }
  
  # if (two.stage) {# to select important variables, is not currently supported
  #   rf.obj <- ranger(dependent.variable.name = response.var, data = data, importance = 'permutation',  probability = probability,
  #                    always.split.variables = adjust.vars, status.variable.name = status.var, ...)
  #   feature.sel <- c(response.var, union(names(sort(rf.obj$variable.importance, decr = TRUE))[1: feature.no], adjust.vars))
  # } else {
  #   feature.sel <- colnames(data)
  # }
  # if((method %in% c('omnibus','phylogenetic')) & (is.null(tree))) 
    # stop('Please specify a phylogenetic tree.')
  if(!probability & test.type != "OOB"){
    warning("test.type is set to OOB when probability is FALSE")
    test.type = "OOB"
  }
  if(method=='omnibus' & test.type=="0.632+") 
    stop('Omnibus RFtest can not support 0.632+ currently.')
  if(is.null(tree)) cum.mat<-comm else 
    if (class(tree) == "phylo") cum.mat <- (t(calculate.cum(comm, tree))) else
	warning("Tree is not of class \"phylo\"!")
  

  if(method=="omnibus"){
    stat.o<-rep(0,2)
    names(stat.o)<-c("weighted", "unweighted")
	  cum.mat.pa <- 1*(cum.mat != 0) # cum.mat and cum.mat.pa will also be used to calculate stat.p
    for(i in 1:reps){# take average of reps runs to estimate stat.o
      #weighted
      data <- data.frame(cum.mat, meta.data[, c(response.var, adjust.vars), drop = FALSE])
      feature.sel <- colnames(data)
      rf.obj <- ranger(dependent.variable.name = response.var, data = data[, feature.sel, drop = FALSE], 
                       importance = 'none',  probability = probability,
                       always.split.variables = adjust.vars, status.variable.name = status.var, ...)
      stat.o["weighted"]<-stat.o["weighted"]+calculateError(rf.obj, data, response.var, prediction.type, test.type)$stat
      #unweighted
      data <- data.frame(cum.mat.pa, meta.data[, c(response.var, adjust.vars), drop = FALSE])
      feature.sel <- colnames(data)
      rf.obj <- ranger(dependent.variable.name = response.var, data = data[, feature.sel, drop = FALSE], 
                       importance = 'none',  probability = probability,
                       always.split.variables = adjust.vars, status.variable.name = status.var, ...)
	    temp<-calculateError(rf.obj, data, response.var, prediction.type, test.type)
      stat.o["unweighted"]<-stat.o["unweighted"]+temp$stat
          
      w<-temp$w #for omnibus RFtest, the w should be 0.632 or NULL, otherwise an error will be thrown above.
    }
  }else if(method=="weighted"){
    stat.o<-0
    names(stat.o)<-"weighted"
    for(i in 1:reps){
	  data <- data.frame(cum.mat, meta.data[, c(response.var, adjust.vars), drop = FALSE])
      feature.sel <- colnames(data)
      rf.obj <- ranger(dependent.variable.name = response.var, data = data[, feature.sel, drop = FALSE], 
                       importance = 'none',  probability = probability,
                       always.split.variables = adjust.vars, status.variable.name = status.var, ...)
	  temp<-calculateError(rf.obj, data, response.var, prediction.type, test.type)
	  stat.o = stat.o+temp$stat
	  w <- temp$w
    }
  }else if(method=="unweighted"){
    stat.o<-0
    names(stat.o)<-"unweighted"
    cum.mat.pa <- 1*(cum.mat != 0)
    for(i in 1:reps){
	  data <- data.frame(cum.mat.pa, meta.data[, c(response.var, adjust.vars), drop = FALSE])
      feature.sel <- colnames(data)
      rf.obj <- ranger(dependent.variable.name = response.var, data = data[, feature.sel, drop = FALSE], 
                       importance = 'none',  probability = probability,
                       always.split.variables = adjust.vars, status.variable.name = status.var, ...)
	  temp<-calculateError(rf.obj, data, response.var, prediction.type, test.type)
      stat.o = stat.o+temp$stat
      w<-temp$w
    }
  }
  stat.o<-stat.o/reps

  
  
  if (trace)  cat('Begin permutation ...\n')
  if(perm.no >= 1){
    cl<-makeCluster(n.cores)
    clusterExport(cl, c("ranger","calculateError", "calculate.cum")) 
    clusterEvalQ(cl, {library(ape)})
    #check if I need to export variables such as "meta", "comm" .... 
    #In this version, it seems that variables in this namespace (the function namespace) will
    #autonomously be exported into the cluster. However, these variables should be evaluated 
    #once but not just passed via the arguments.
    #check if this is needed:
    #clusterExport(cl, varlist = ls(), envir = environment()) #export all the variables in the function namespace
    stat.p <- parSapply(cl,1:perm.no, function(i) {
      #    if (trace) { #doesn't work in parallel version
      #      if (i %% 100 == 0) cat('.') 
      #    }
      # permutate the covariates, simulating the response var
      meta.perm<-meta.data
      if(!is.null(adjust.vars)){
        nr = nrow(meta.data)
        #meta.perm <- meta.perm[sample.int(nr), ] 
        if(prediction.type %in% c('Regression')) {
          meta.perm[,response.var] <- 
            predict(object = model, newdata = meta.perm) + sample(residuals(model)) #fitted values and permutated residuals
        } else if(prediction.type %in% 'Classification'){
          pr <- predict(object = model, newdata = meta.perm , type= "response") #the expected probability for meta.perm
          temp <- rbinom(n=nr, size = 1, prob = pr)
          names(temp) <- names(pr)
          temp <- factor(temp, levels = c(0,1)) #make sure 1 will be the second level (success)
          if(length(levels(meta.data[,response.var])) != length(levels(temp))) 
            stop("Currently, RFtest only accepts a dichotomous response variable!")
          levels(temp) <- levels(meta.data[,response.var])
          meta.perm[,response.var] <- temp
        }
      } else
        meta.perm[,response.var]<-sample(meta.data[, response.var])
      
      # if (two.stage) { #not currently supported
      #   rf.obj <- ranger(dependent.variable.name = response.var, data = data.perm, importance = 'permutation',  probability = probability,
      #                    always.split.variables = adjust.vars, status.variable.name = status.var, ...)
      #   feature.sel <- c(response.var, union(names(sort(rf.obj$variable.importance, decr = TRUE))[1: feature.no], adjust.vars))
      # } else {
      #  feature.sel <- colnames(data.perm)
      #}
      if(method=="omnibus"){
        stat.p<-rep(NA,2)
        names(stat.p)<-c("weighted","unweighted")
        #weighted
        data <- data.frame(cum.mat, meta.perm[, c(response.var, adjust.vars), drop = FALSE])
        feature.sel <- colnames(data)
        rf.obj <- ranger(dependent.variable.name = response.var, data = data[, feature.sel, drop = FALSE], 
                         importance = 'none',  probability = probability,
                         always.split.variables = adjust.vars, status.variable.name = status.var, ...)
        stat.p["weighted"]<-calculateError(rf.obj, data, response.var, prediction.type, test.type, w)$stat
        #unweighted
        data <- data.frame(cum.mat.pa, meta.perm[, c(response.var, adjust.vars), drop = FALSE])
        feature.sel <- colnames(data)
        rf.obj <- ranger(dependent.variable.name = response.var, data = data[, feature.sel, drop = FALSE], 
                         importance = 'none',  probability = probability,
                         always.split.variables = adjust.vars, status.variable.name = status.var, ...)
        stat.p["unweighted"]<-calculateError(rf.obj, data, response.var, prediction.type, test.type, w)$stat
      }else if(method=="weighted"){
        data <- data.frame(cum.mat, meta.perm[, c(response.var, adjust.vars), drop = FALSE])
        feature.sel <- colnames(data)
        rf.obj <- ranger(dependent.variable.name = response.var, data = data[, feature.sel, drop = FALSE], 
                         importance = 'none',  probability = probability,
                         always.split.variables = adjust.vars, status.variable.name = status.var, ...)		
        stat.p<-calculateError(rf.obj, data, response.var, prediction.type, test.type, w)$stat
        names(stat.p)<-"weighted"
      }else if(method=="unweighted"){		
        data <- data.frame(cum.mat.pa, meta.perm[, c(response.var, adjust.vars), drop = FALSE])
        feature.sel <- colnames(data)
        rf.obj <- ranger(dependent.variable.name = response.var, data = data[, feature.sel, drop = FALSE], 
                         importance = 'none',  probability = probability,
                         always.split.variables = adjust.vars, status.variable.name = status.var, ...)
        stat.p<-calculateError(rf.obj, data, response.var, prediction.type, test.type, w)$stat
        names(stat.p)<-"unweighted"
      }
      return(stat.p)
    }
    )
    stopCluster(cl)
    #  cat('\n')
    
    if(method=="omnibus"){
      p.value.perm<-rep(NA,3)
      names(p.value.perm)<-c("weighted", "unweighted", "Omnibus")
      cbind(stat.o,stat.p)->stat
      t(apply(stat, MARGIN=1, rank))->stat.r #the p.values
      stat.r[,1]/(perm.no+1)->p.value.perm[1:(length(p.value.perm)-1)]
      apply(stat.r, MARGIN=2, min) -> stat.r.m #the best p.value each run, also the sample distribution of Omnibus
      p.value.perm["Omnibus"] <- sum(stat.r.m <= stat.r.m[1])/(perm.no+1)
    }else p.value.perm <- (sum(stat.p <= stat.o) + 1) / (perm.no + 1)
    
    
    
    
    if (pearson.fit) {
      moments <- empMoments(stat.p)
      pearson.para <- pearsonFitM(moments = moments)
      p.value.pearson <- ppearson(stat.o, pearson.para)
      #		hist(stat.p, prob=TRUE)
      #		x <- seq(0.20, 0.32, len=1000)
      #		lines(x, dpearson(x, moments = empMoments(stat.p), lty = 2, col = 'red'))
    } else {
      p.value.pearson <- NA
      pearson.para <- NA
    }
  } else{
    p.value.perm <- NA
    p.value.pearson <- NA
    pearson.para <- NA
    stat.p <- NA
  }
  
  return(list(call = call, method = method, p.value.perm = p.value.perm, p.value.pearson = p.value.pearson, 
              pearson.para = pearson.para, stat.o = stat.o, w = w, drop.otus=drop.otus, stat.p=stat.p))
}
