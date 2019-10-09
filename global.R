library(shiny)

library(shinyBS)
library(shinythemes)
library(beeswarm)
library(plotrix)
library(xfun)

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
options(warn=-1)
## Loading Primary Tumour Binary Event Matrices
data(CELLector.PrimTum.BEMs)
data(CELLector.PrimTum.BEMs_v2)

OrPrimTumBEMs<-CELLector.PrimTum.BEMs
OrPrimTumBEMs_v2<-CELLector.PrimTum.BEMs_v2

CELLector.PrimTum.BEMs<-OrPrimTumBEMs[2:length(OrPrimTumBEMs)]

## Loading Cell Lines's Binary Event Matrices
data(CELLector.CellLine.BEMs)
data(CELLector.CellLine.BEMs_v2)

OrCellLineBEMs<-CELLector.CellLine.BEMs
OrCellLineBEMs_v2<-CELLector.CellLine.BEMs_v2

data(CELLector.CFEs.CNAid_decode)
data(CELLector.CFEs.HMSid_decode)

CNAid_decode<-data.frame(Id=as.character(CELLector.CFEs.CNAid_decode$CNA_Identifier),
                   CancerType=as.character(CELLector.CFEs.CNAid_decode$CancerType),
                   Gain_Loss=as.character(CELLector.CFEs.CNAid_decode$Recurrent),
                   Locus=as.character(CELLector.CFEs.CNAid_decode$locus),
                   n.Genes=as.numeric(CELLector.CFEs.CNAid_decode$nGenes),
                   Genes=as.character(CELLector.CFEs.CNAid_decode$ContainedGenes),stringsAsFactors = FALSE)
CNAid_decode$Gain_Loss[CNAid_decode$Gain_Loss=='Amplification']<-'Gain'
CNAid_decode$Gain_Loss[CNAid_decode$Gain_Loss=='Deletion']<-'Loss'

HMSid_decode<-data.frame(Id=as.character(CELLector.CFEs.HMSid_decode$hms_id),
                         CancerType=as.character(CELLector.CFEs.HMSid_decode$Cancer.Types),
                         GenomicCoords=as.character(CELLector.CFEs.HMSid_decode$Genomic.Coordinates),
                         DownStream.Genes=as.character(CELLector.CFEs.HMSid_decode$GN))

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
data(CELLector.PrimTumVarCatalog)


## Deriving available TCGA labels
TCGALabels<-names(CELLector.PrimTum.BEMs)

tumours<-CELLector.PrimTum.BEMs$COREAD
colnames(tumours)<-paste(colnames(tumours),'_',1:ncol(tumours),sep='')

features<-CELLector.CFEs
pathways<-names(CELLector.Pathway_CFEs)

## Downloading Cell Model Passports annotations
CMPs_model_annotations<-CELLector.CMPs_getModelAnnotation()
CMPs_model_annotations$cancer_type_detail<-
  str_sub(CMPs_model_annotations$cancer_type_detail,3,end = str_length(CMPs_model_annotations$cancer_type_detail)-3)
CMPs_driverGenes<-CELLector.CMPs_getDriverGenes()

#Iorios_driverGenes<-
#CMPs_variants<-CELLector.CMPs_getVariants()

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
CELLector_App.checkBEMformats<-function(primTumBEMs,cellLinBEMS){
  
  errFlag<- 0
  
  if (!is.list(primTumBEMs) | !is.list(cellLinBEMs)){
    errFlag<-1
    errMessage<-paste('Wrong BEM formats!\nA named list of binary matrices (with TCGA cancer type labels as names). The entries of each of these matrices indicate the status (Present/Absent) of each CFE (one per row) across primary tumors samples (one per column).')
  }
  
}

CELLector_App.current_Model_ids<-function(input,annotation,genomicData=TRUE){
  ids<-which(
          CMPs_model_annotations$tissue==input$CMP_selectTissue & 
          is.element(CMPs_model_annotations$cancer_type,input$CMP_selectCancerType) &
          is.element(CMPs_model_annotations$cancer_type_detail,input$CMP_selectCancerType_detailed) &
          is.element(CMPs_model_annotations$sample_site,input$CMP_selectSample_site) &
          (input$CMP_exclude_organoids * (CMPs_model_annotations$model_type!='Organoid') |
             !input$CMP_exclude_organoids * rep(1,nrow(CMPs_model_annotations))) &
          (input$CMP_human_samples_only * (CMPs_model_annotations$species=="Homo Sapiens") |
             !input$CMP_human_samples_only * rep(1,nrow(CMPs_model_annotations))) &
          (CMPs_model_annotations$gender==input$CMP_gender | 
             (input$CMP_gender=='All (including Unknown)')*rep(1,nrow(CMPs_model_annotations))) &
          (CMPs_model_annotations$msi_status==input$CMP_msi_status | 
             (input$CMP_msi_status=='All (including NA)')*rep(1,nrow(CMPs_model_annotations))) &
          ((!input$CMP_based_on_mut_burden*rep(1,nrow(CMPs_model_annotations))) |
             (input$CMP_based_on_mut_burden * (round(CMPs_model_annotations$mutational_burden) >= input$CMP_mutBurdend_slide[1] &
                                                 round(CMPs_model_annotations$mutational_burden) <= input$CMP_mutBurdend_slide[2])
             )
          ) &
          ((!input$CMP_based_on_ploidy*rep(1,nrow(CMPs_model_annotations))) |
             (input$CMP_based_on_ploidy * (round(CMPs_model_annotations$ploidy) >= input$CMP_ploidy_slide[1] &
                                             round(CMPs_model_annotations$ploidy) <= input$CMP_ploidy_slide[2])
             )
          ) &
          ((!input$CMP_age_at_sampling*rep(1,nrow(CMPs_model_annotations))) |
             (input$CMP_age_at_sampling * (round(CMPs_model_annotations$age_at_sampling) >= input$CMP_age_at_sampling_slide[1] &
                                             round(CMPs_model_annotations$age_at_sampling) <= input$CMP_age_at_sampling_slide[2])
             )
          ) &
          ((!input$CMP_based_on_etnicity*rep(1,nrow(CMPs_model_annotations))) |
             (input$CMP_based_on_etnicity * !is.element(CMPs_model_annotations$ethnicity,input$CMP_etnicity))
          )
  )
  

  
  if(genomicData){
    ids<-ids[which(CMPs_model_annotations$mutation_data[ids]=='True')]
  }
  
  return(ids)
}


