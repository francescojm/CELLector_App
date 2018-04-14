library(shiny)
library(shinythemes)
library(beeswarm)
library(plotrix)

#polar.plot(), pie3D()
library(devtools)
#library(RColorBrewer)
#install_github("najha/CELLector")
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

CELLector.PrimTum.BEMs<-CELLector.PrimTum.BEMs[2:length(CELLector.PrimTum.BEMs)]
## Loading Cell Lines's Binary Event Matrices
data(CELLector.CellLine.BEMs)

data(CELLector.CFEs.CNAid_decode)

CNAid_decode<-data.frame(Id=as.character(CELLector.CFEs.CNAid_decode$CNA_Identifier),
                   CancerType=as.character(CELLector.CFEs.CNAid_decode$CancerType),
                   Gain_Loss=as.character(CELLector.CFEs.CNAid_decode$Recurrent),
                   Locus=as.character(CELLector.CFEs.CNAid_decode$locus),
                   n.Genes=as.numeric(CELLector.CFEs.CNAid_decode$nGenes),
                   Genes=as.character(CELLector.CFEs.CNAid_decode$ContainedGenes),stringsAsFactors = FALSE)
CNAid_decode$Gain_Loss[CNAid_decode$Gain_Loss=='Amplification']<-'Gain'
CNAid_decode$Gain_Loss[CNAid_decode$Gain_Loss=='Deletion']<-'Loss'


for (i in 1:nrow(CNAid_decode)){
   
   tmp<-CNAid_decode$Genes[i]
   
   tmp<-unlist(str_split(tmp,','))
   if(length(tmp)>5){
     
     ntmp<-intersect(tmp,c(CELLector.HCCancerDrivers,'ERBB2','MYC','PTEN'))
     ntmp<-sort(ntmp)[1:5]
     tmp<-paste(c(sort(c(ntmp,tmp[1:(5-length(ntmp))])),'...'),collapse=', ')
   }else{
     tmp<-paste(tmp,collapse=', ')
   }
   
   CNAid_decode$Genes[i]<-tmp
   
 }


data(CELLector.CFEs)
data(CELLector.Pathway_CFEs)
data(CELLector.MSIstatus)

## Deriving available TCGA labels
TCGALabels<-names(CELLector.PrimTum.BEMs)

tumours<-CELLector.PrimTum.BEMs$COREAD
colnames(tumours)<-paste(colnames(tumours),'_',1:ncol(tumours),sep='')

features<-CELLector.CFEs
pathways<-names(CELLector.Pathway_CFEs)


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

