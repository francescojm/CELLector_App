library(shiny)
library(shinythemes)
#library(stringr)
library(plotrix)

#polar.plot(), pie3D()
#library(devtools)
#library(RColorBrewer)
library(CELLector)
library(collapsibleTree)
#library(pander)

#install_github("sjp/grImport2")
library(grImport2)

library(data.tree)
library(igraph)
library(sunburstR)

## Loading Primary Tumour Binary Event Matrices
data(CELLector.PrimTum.BEMs)
## Loading Cell Lines's Binary Event Matrices
data(CELLector.CellLine.BEMs)

data(CELLector.CFEs)
data(CELLector.Pathway_CFEs)

## Deriving available TCGA labels
TCGALabels<-names(CELLector.PrimTum.BEMs)

tumours<-CELLector.PrimTum.BEMs$COREAD
features<-CELLector.CFEs
pathways<-names(CELLector.Pathway_CFEs)

CELLector_App.sunBurstFormat<-function(table_tree){
  
  table_tree<-data.frame(lapply(table_tree[,1:11], as.character), stringsAsFactors=FALSE)
  table_tree$Left.Child.Index[table_tree$Left.Child.Index==0]<- -1
  table_tree$Right.Child.Index[table_tree$Right.Child.Index==0]<- -1
  
  table_tree<-rbind(c(0,'TOTAL','TOTAL','root','-1',table_tree$CurrentTotal[1],table_tree$CurrentTotal[1],1,1,0),table_tree)
  table_tree$Type[2]<-'Left.Child'
  
  table_tree$Idx<-as.numeric(table_tree$Idx)+1
  table_tree$Parent.Idx<-as.numeric(table_tree$Parent.Idx)+1
  table_tree$Left.Child.Index<-as.numeric(table_tree$Left.Child.Index)+1
  table_tree$Right.Child.Index<-as.numeric(table_tree$Right.Child.Index)+1
  
  leaves<-which(table_tree$Left.Child.Index==0 & table_tree$Right.Child.Index==0)
  
  stable_tree<-table_tree
  
  for (i in 1:length(leaves)){
    
    CurrentTotal<-as.numeric(table_tree$CurrentTotal[leaves[i]])-as.numeric(table_tree$AbsSupport[leaves[i]])
    
    
    
    stable_tree<-rbind(stable_tree,c(nrow(stable_tree)+1,'Others','Others','Right.Child',leaves[i],CurrentTotal,CurrentTotal,
                                     1,CurrentTotal/as.numeric(stable_tree$CurrentTotal[1]),0,0))
  }
  
  
  nnodes<-nrow(stable_tree)
  edgeList<-NULL
  
  for (i in 1:nnodes){
    
    currentNode<-i
    
    
    if (stable_tree$Type[currentNode]=='Left.Child'){
      edgeList<-c(edgeList,as.numeric(c(stable_tree$Parent.Idx[currentNode],stable_tree$Idx[currentNode])))
    #  print(c(stable_tree$Parent.Idx[currentNode],stable_tree$Idx[currentNode]))
    }else{
      startingNode<-currentNode
      while (stable_tree$Type[startingNode]=='Right.Child'){
        startingNode<-as.numeric(stable_tree$Parent.Idx[startingNode])
        
      }
      if(currentNode!=1){
        edgeList<-c(edgeList,as.numeric(c(stable_tree$Parent.Idx[startingNode],currentNode)))
        #print(c(stable_tree$Parent.Idx[startingNode],currentNode))
      }
    }
  }
  
  G<-make_graph(edgeList)
  leaves<-V(G)[which(degree(G,v=V(G),'out')==0)]
  paths<-all_simple_paths(G,1,leaves)
  nleaves<-length(paths)
  
  chainP<-vector()
  npat<-vector()
  for (i in 1:nleaves){
    
    currentId<-as.numeric(paths[[i]][length(paths[[i]])])
    
    
    chainP[i]<-paste(paste(as.numeric(paths[[i]])-1,stable_tree$ItemsDecoded[as.numeric(as.character(paths[[i]]))]),
                     collapse='-')
    
    npat[i]<-as.numeric(stable_tree$AbsSupport[currentId])
  }
  
  sequences<-data.frame(V1=chainP,V2=npat,stringsAsFactors = FALSE)
  
  
  # labes<-unique(unlist(strsplit(sequences$V1,"-")))
  # 
  # n <- length(labes)
  # qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  # col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  # 
  # # define specific colors
  # colors <- list(
  #   domain=unique(unlist(strsplit(sequences$V1,"-"))),
  #   range=sample(col_vector, n)
  # )
  # 
  # 
  # colors$range[colors$domain=='Others']<-'white'
  # 
  
  return(sequences)
  
}
CELLector_App.complementarPieChart<-function(Tree,NavTab,nodeIdx){
  
  supports<-vector()
  supports[1]<-NavTab$GlobalSupport[[nodeIdx]]
  
  flag<-2
  names(supports)<-paste('SubT.',nodeIdx,sep='')
  
  NIDX<-nodeIdx
  
  while(NavTab$Right.Child.Index[nodeIdx]>0){
    
    nodeIdx<-NavTab$Right.Child.Index[nodeIdx]
    supports[flag]<-NavTab$GlobalSupport[[nodeIdx]]
    names(supports)[flag]<-paste('SubT.',nodeIdx,sep='')
    flag<-flag+1
    NIDX<-c(NIDX,nodeIdx)
    
  }
  
  tmpCol<-Get(Traverse(Tree,traversal = 'level'),'Colors')
  nn<-names(tmpCol)
  nn<-str_split(nn,' ')
  nn<-as.numeric(unlist(lapply(nn,function(x){x[[1]][1]})))
  
  id<-match(NIDX,nn)
  COLORS<-tmpCol[id]
  
  supports<-c(100*supports,100-100*sum(supports))
  names(supports)[flag]<-'Others'  
  
  return(list(supports=supports,COLORS=COLORS))
}

