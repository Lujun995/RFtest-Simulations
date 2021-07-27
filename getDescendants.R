getDescendants.node<-function(tree,node,curr=NULL){
  if(is.null(curr)) curr<-vector()
  daughters<-tree$edge[which(tree$edge[,1]%in%node),2] #if node<=Ntip, it will return null.
  curr<-c(curr,daughters[which(daughters>Ntip(tree))])
  if(length(daughters)>0) curr<-getDescendants.node(tree,daughters,curr)
  return(curr)
}

getDescendants<-function(tree,node,curr=NULL){
  if(is.null(curr)) curr<-vector()
  daughters<-tree$edge[which(tree$edge[,1]%in%node),2] #if node<=Ntip, it will return null.
  curr<-c(curr,daughters)
  if(length(daughters)>0) curr<-getDescendants(tree,daughters,curr)
  return(curr)
}

getDescendants.tip<-function(tree,node,curr=NULL){
  if(is.null(curr)) curr<-vector()
  daughters<-tree$edge[which(tree$edge[,1]%in%node),2] #if node<=Ntip, it will return null.
  curr<-c(curr,daughters[which(daughters<=Ntip(tree))])
  if(length(daughters)>0) curr<-getDescendants.tip(tree,daughters,curr)
  return(curr)
}
