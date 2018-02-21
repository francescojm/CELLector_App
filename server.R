
#testing the new repo

# Define server logic required to draw a histogram
#' Title
#'
#' @param input 
#' @param output 
#'
#' @return
#' @export
#'
#' @examples
server <- function(input, output) {
  
  output$cnaDecodeTable <- renderDataTable(CNAid_decode[CNAid_decode$CancerType==input$selectCancerType,],
                                     options = list(
                                       pageLength = 10,autoWidth = FALSE)
                                     )
     
  
  
  output$str <- renderPrint({
   if(!is.null(NT$data)){
     
     if (input$subSet!=''){
       if(!input$checkboxNegation){
         if (is.element(input$subSet,rownames(TUMOURS$data))){
          nn <- length(which(TUMOURS$data[input$subSet,]>0))
          }else{
            nn <- 0
            }
         if (is.element(input$subSet,colnames(CELLlineData$data))){
           nc <- length(which(CELLlineData$data[,input$subSet]>0))
           CL <- rownames(CELLlineData$data)[which(CELLlineData$data[,input$subSet]>0)]
           }else{
             nc <- 0
             CL<-NULL
             }
         }else{
           if (is.element(input$subSet,rownames(TUMOURS$data))){
             nn <- length(which(TUMOURS$data[input$subSet,]==0))
           }else{
             nn <- ncol(TUMOURS$data)
             
           }
           if (is.element(input$subSet,colnames(CELLlineData$data))){
             nc <- length(which(CELLlineData$data[,input$subSet]==0))
             CL <- rownames(CELLlineData$data)[which(CELLlineData$data[,input$subSet]==0)]
           }else{
             nc <- nrow(CELLlineData$data)
             CL<-rownames(CELLlineData$data)
           }
         }
       }else{
         nn <- ncol(TUMOURS$data)
         nc <- nrow(CELLlineData$data)
         CL<-rownames(CELLlineData$data)
       }
     
      if(input$whatToInclude2=='Microsatellite stable'){
        CTE<-rownames(CELLlineData$data)[match(names(which(CELLector.MSIstatus=="MSI-H")),COSMICids$data)]
      }else{
        if(input$whatToInclude2=='Microsatellite instable'){
           CTE<-rownames(CELLlineData$data)[match(names(which(CELLector.MSIstatus!="MSI-H")),COSMICids$data)]
        }else{
           CTE<-NULL
        } 
      }
     
      nc<-length(setdiff(CL,CTE))
     
      pander::pander(paste(nc,input$selectCancerType,'cell lines and ',
                           nn,'patients considered in this session [Re-build CELLector Search Space to Make any change to the criteria below effective]'))
      
    }else{
      pander::pander('Build CELLector Search Space to START')
   }
  })
  
  TUMOURS <- reactiveValues(data = NULL)
  FEATURES <- reactiveValues(data = NULL)
  CTYPE <- reactiveValues(data = NULL)
  NT <- reactiveValues(data = NULL)
  
  PATIENTcoords<-reactiveValues(data = NULL)
  
  RULES <- reactiveValues(data = NULL)
  SELECTEDNODE <- reactiveValues(data = NULL)
  SIGNATURES <- reactiveValues(data = NULL)
  encodedSIGNATURES <- reactiveValues(data = NULL)
  CELLlineData <- reactiveValues(data = NULL)
  COSMICids<-reactiveValues(data = NULL)
  
  
  STATUS <- reactiveValues(data = NULL)
  SunBurstSequences <- reactiveValues(data = NULL)
  
  
  output$DownSearchSpace <- downloadHandler(
    filename = function(){
      paste("CELLector_SearchSpace_", input$selectCancerType, "_", Sys.Date(), ".tsv", sep = "")
    },
    content = function(file) {
      if(length(NT$data)>0){
        toSave<-cbind(NT$data$navTable[,c(1,3:9)],SIGNATURES$data)
        
        colnames(toSave)<-c('Node Index',
                            'Genomic Alterations',
                            'Node Type (Left = Refinement, Right = Complement)',
                            'Parent Node',
                            'N. Represented Patients (Globally)',
                            'N. Patients in the considered sub-population',
                            'Relative Support %',
                            'Global Support %',
                            'Genomic Signature'
        )
        
        nsig<-length(encodedSIGNATURES$data)
        
        FTE<-NULL
        if (input$subSet!='' & !input$checkboxNegation){
          if (is.element(input$subSet,colnames(CELLlineData$data))){
            FTE<-rownames(CELLlineData$data)[which(CELLlineData$data[,input$subSet]==0)]
          }else{
            FTE<-rownames(CELLlineData$data)
          }    
        }
        if (input$subSet!='' & input$checkboxNegation){
          if (is.element(input$subSet,colnames(CELLlineData$data))){
            FTE<-rownames(CELLlineData$data)[which(CELLlineData$data[,input$subSet]==1)]
          }else{
            FTE<-NULL
          }
        }
        
        if(input$whatToInclude2=='Microsatellite stable'){
          CTE<-rownames(CELLlineData$data)[match(names(which(CELLector.MSIstatus=="MSI-H")),COSMICids$data)]
        }else{
          if(input$whatToInclude2=='Microsatellite instable'){
            CTE<-rownames(CELLlineData$data)[match(names(which(CELLector.MSIstatus!="MSI-H")),COSMICids$data)]
          }else{
            CTE<-NULL
          } 
        }
        
        MODELS<-vector()
        suppressWarnings(
          for (cc in 1:nsig){
            
            solved<-CELLector.solveFormula(encodedSIGNATURES$data[[cc]],dataset = CELLlineData$data,To_beExcluded = union(CTE,FTE))    
            
            MODELS[cc]<-paste(sort(solved$PS),collapse=', ')
          }
        )
        toSave<-cbind(toSave,MODELS)
        colnames(toSave)[ncol(toSave)]<-'Representative Cell Line'
        
        visit<-CELLector.selectionVisit(NT$data$navTable)
        
        TOS<-toSave[visit,]
        
        colnames(TOS)[1]<-'Tumour SubType Index'
        colnames(TOS)[2]<-'Node Genomic Alteration'
        
        TOS<-TOS[,c(1,9,2,3,4,5,6,7,8,10)]
        write.table(TOS, file, sep='\t',row.names = FALSE,quote=FALSE)  
        
      }
      
    }
  )
  
  output$CELLect <- downloadHandler(
    filename = function(){
      paste("CELLected_CellLines_", input$selectCancerType, "_", Sys.Date(), ".tsv", sep = "")
    },
    content = function(file) {
      if(length(NT$data)>0){
        toSave<-cbind(NT$data$navTable[,c(1,3:9)],SIGNATURES$data)
        
        nsig<-length(encodedSIGNATURES$data)
        
        MODELS<-vector()
        
        FTE<-NULL
        if (input$subSet!='' & !input$checkboxNegation){
          if (is.element(input$subSet,colnames(CELLlineData$data))){
            FTE<-rownames(CELLlineData$data)[which(CELLlineData$data[,input$subSet]==0)]
          }else{
            FTE<-rownames(CELLlineData$data)
          }    
        }
        if (input$subSet!='' & input$checkboxNegation){
          if (is.element(input$subSet,colnames(CELLlineData$data))){
            FTE<-rownames(CELLlineData$data)[which(CELLlineData$data[,input$subSet]==1)]
          }else{
            FTE<-NULL
          }
        }
        
        if(input$whatToInclude2=='Microsatellite stable'){
          CTE<-rownames(CELLlineData$data)[match(names(which(CELLector.MSIstatus=="MSI-H")),COSMICids$data)]
        }else{
          if(input$whatToInclude2=='Microsatellite instable'){
            CTE<-rownames(CELLlineData$data)[match(names(which(CELLector.MSIstatus!="MSI-H")),COSMICids$data)]
          }else{
            CTE<-NULL
          } 
        }
        
        
        
        suppressWarnings(
          for (cc in 1:nsig){
            solved<-CELLector.solveFormula(encodedSIGNATURES$data[[cc]],dataset = CELLlineData$data,To_beExcluded = union(CTE,FTE))    
            MODELS[cc]<-paste(sort(solved$PS),collapse=', ')
          }
        )
        
        toSave<-cbind(toSave,MODELS)
        colnames(toSave)[ncol(toSave)]<-'Representative Cell Line'
        
        visit<-CELLector.selectionVisit(NT$data$navTable)
        
        ncellLines<-input$N.CellLines
        sortedModels<-MODELS[visit]
        
        NODEidx<-visit[which(sortedModels!='')]
        sortedModels<-sortedModels[which(sortedModels!='')]
        modelMat<-CELLector.buildModelMatrix(sortedModels)
        
        res<-CELLector.makeSelection(modelMat,ncellLines)
        res$modelAccounted<-NODEidx[res$modelAccounted]
        
        colnames(res)<-c('Tumour SubType Index','Representative Cell Line')
        write.table(res, file,sep='\t', row.names = FALSE, quote=FALSE)    
      }
    }
  )
  
  observeEvent(input$action, {
    CELLlineData$data<-CELLector.CellLine.BEMs[[input$selectCancerType]]
    r<-CELLlineData$data[,2]
    COSMICids$data<-CELLlineData$data[,1]
    CELLlineData$data<-CELLlineData$data[,3:ncol(CELLlineData$data)]
    rownames(CELLlineData$data)<-r
    SELECTEDNODE$data <- NULL
    TUMOURS$data <- CELLector.PrimTum.BEMs[[input$selectCancerType]]
    colnames(TUMOURS$data)<-paste(colnames(TUMOURS$data),'_',1:ncol(TUMOURS$data),sep='')
    FEATURES$data <- rownames(TUMOURS$data)
    
    PATIENTcoords$data<-sample(ncol(TUMOURS$data))
    names(PATIENTcoords$data)<-colnames(TUMOURS$data)
    
    
    #.............. 
    progress <- shiny::Progress$new(style='old')
    progress$set(message = "Loading primary tumour data and building CELLector Search Space... Please Wait", value = 0)
    #.............. 
    ctype <- input$selectCancerType
    CTYPE$data <- ctype

    minLen <- input$minSetSize
    minGlobSupp <- input$minGlobalSupport/100
    pathway <- input$pathFocus
    #..............
    mutonly <- FALSE
    cnaonly <- FALSE

    if(input$whatToInclude=='Mutations in high confidence cancer genes'){
      mutonly <- TRUE
    }
    if(input$whatToInclude=='Recurrently CN altered chromosomal segments'){
      cnaonly <- TRUE
    }
    #..............
    #if (input$whereToNeglect=='While building search space' | input$whereToNeglect=='Always'){
    #  toRemove <- input$toExclude
    #}else{ toRemove <- NULL
    #}

    #..............
    progress$set(message = "Loading primary tumour data and building CELLector Search Space... Please Wait", value = 0.5)
       NT$data<- CELLector.Build_Search_Space(ctumours = t(TUMOURS$data),
                                           cancerType = input$selectCancerType,
                                           minlen=input$minSetSize,
                                           mutOnly = (input$whatToInclude=='Mutations in high confidence cancer genes'),
                                           cnaOnly = (input$whatToInclude=='Recurrently CN altered chromosomal segments'),
                                           minGlobSupp = input$minGlobalSupport/100,
                                           FeatureToExclude = input$toExclude,
                                           pathwayFocused = input$pathFocus,
                                           pathway_CFEs = CELLector.Pathway_CFEs,
                                           cnaIdMap = CELLector.CFEs.CNAid_mapping,
                                           cnaIdDecode = CELLector.CFEs.CNAid_decode,
                                           cdg = CELLector.HCCancerDrivers,
                                           subCohortDefinition=input$subSet,
                                           NegativeDefinition=input$checkboxNegation)
    if(nrow(NT$data$navTable)>0){

      S <- CELLector.createAllSignatures(NavTab = NT$data$navTable)

      SIGNATURES$data <- S$S
      encodedSIGNATURES$data <- S$ES

      CLD <- CELLector.CellLine.BEMs[[input$selectCancerType]]
      CIDS <- CLD$COSMIC_identifier
      rn <- CLD[,2]
      CLD <- as.matrix(CLD[,3:ncol(CLD)])
      rownames(CLD) <- rn

    }

    progress$set(message = "Done!", value = 1)
    progress$close()

    #cnaLookUp$data <- CELLector.cna_look_up(input$cnaID , input$selectCancerType)
    
    #CLG_feature$data <- CELLector.get_cell_line_genomic_features(input$CLname, input$selectCancerType)
    
    
  })

  output$plot <- renderCollapsibleTree({
    
    if(length(NT$data)>0){
      
      if (nrow(NT$data$navTable)>1){
        
        CC <- colors(distinct = TRUE)
        CC <- CC[setdiff(1:length(CC),grep('gray',CC))]
        CC <- rgb(t(col2rgb(CC)),maxColorValue = 255)
        
        COLORSbyLev <- CC[sample(length(CC))][1:NT$data$TreeRoot$totalCount]
        
        RelatesToFatherAs <- rep('-',NT$data$TreeRoot$totalCount)
        RelatesToFatherAs[which(Get(Traverse(NT$data$TreeRoot,traversal = 'level'),attribute = 'NodeType')=='Right.Child')]<-'Complement'
        RelatesToFatherAs[which(Get(Traverse(NT$data$TreeRoot,traversal = 'level'),attribute = 'NodeType')=='Left.Child')]<-'Refinement'
        
        NT$data$TreeRoot$Set(Colors=COLORSbyLev,traversal = 'level')
        NT$data$TreeRoot$Set(RelatesToFatherAs=RelatesToFatherAs,traversal = 'level')
        
        SunBurstSequences$data<-CELLector_App.sunBurstFormat(NT$data$navTable)
        
        collapsibleTree(NT$data$TreeRoot,
                        fill = 'Colors',
                        inputId = 'searchSpace',
                        tooltip = TRUE,
                        attribute = 'RelatesToFatherAs')
        
      } 
    }
  })
  
  output$sunburst <- renderSunburst({
    #invalidateLater(1000, session)
    if(length(NT$data)>0){
      sequences <- SunBurstSequences$data
      
      tmpCol <- Get(Traverse(NT$data$TreeRoot,traversal = 'level'),'Colors')
      ttmp<-tmpCol
      
      names(ttmp)<-NULL
      
      nvoid<-grep('Others',unique(unlist(strsplit(sequences$V1,'-'))),value = TRUE)
      
      stpes<-nvoid
      
      colors <- list(
        domain=c('0 TOTAL',names(tmpCol),stpes),
        range=c('black',ttmp,rep('white',length(stpes)))
      )
      
      #add_shiny(
      #          )
      


      htmlwidgets::onRender(
        sunburst(sequences,breadcrumb = list(w = 400),percent = FALSE,count = FALSE,colors=colors,

                 explanation = "function(d) {     var ssr = d.data.name
                 if (!ssr.match(/Others/gi)){
                 return ssr
                 }
    }"),
        "
        function(el,x){
        d3.select(el).select('.sunburst-sidebar').remove()
        }
        "
        )
    }
  })
  
  output$NodeDetails<-renderTable({
    
    if(length(input$searchSpace)>0){
      nodeLabel<-input$searchSpace[[length(input$searchSpace)]]
      
      nodeIdx<-str_split(nodeLabel,' ')
      itemSet<-nodeIdx[[1]][2]
      nodeIdx<-nodeIdx[[1]][1]
      
      itt<-as.numeric(nodeIdx)
      
      if (itt <= length(SIGNATURES$data)){
        SELECTEDNODE$data<-itt
        data.frame(`Subpopulation_Index`=nodeIdx,
                   Molecular_Signature=SIGNATURES$data[[itt]])
      }else{
        SELECTEDNODE$data<-1
        data.frame(`Subpopulation_Index`=format(1,digits=1),
                   Molecular_Signature=SIGNATURES$data[[1]])
      }
    }else{
      if(length(SIGNATURES$data)>0){
        SELECTEDNODE$data<-1
        data.frame(`Subpopulation_Index`=format(1,digits=1),
                   Molecular_Signature=SIGNATURES$data[[1]]) 
      }else{
        SELECTEDNODE$data<-NULL
        return(NULL)
      }
    }
  })
  
  
  output$CellLineDetails<-renderTable({
    
    if(is.null(SELECTEDNODE$data)){
      if (!is.null(NT$data)){
        SELECTEDNODE$data<-1
      }
    }
    
    if(length(SELECTEDNODE$data)>0){
      FTE<-NULL
      if (input$subSet!='' & !input$checkboxNegation){
        if (is.element(input$subSet,colnames(CELLlineData$data))){
            FTE<-rownames(CELLlineData$data)[which(CELLlineData$data[,input$subSet]==0)]
        }else{
            FTE<-rownames(CELLlineData$data)
        }    
      }
      if (input$subSet!='' & input$checkboxNegation){
        if (is.element(input$subSet,colnames(CELLlineData$data))){
          FTE<-rownames(CELLlineData$data)[which(CELLlineData$data[,input$subSet]==1)]
        }else{
          FTE<-NULL
        }
      }
      
      if(input$whatToInclude2=='Microsatellite stable'){
          CTE<-rownames(CELLlineData$data)[match(names(which(CELLector.MSIstatus=="MSI-H")),COSMICids$data)]
      }else{
        if(input$whatToInclude2=='Microsatellite instable'){
          CTE<-rownames(CELLlineData$data)[match(names(which(CELLector.MSIstatus!="MSI-H")),COSMICids$data)]
        }else{
          CTE<-NULL
        } 
      }
      
      tmp <- CELLector.solveFormula(encodedSIGNATURES$data[[SELECTEDNODE$data]],dataset = CELLlineData$data,
                                    To_beExcluded = union(CTE,FTE))
      
      if(length(tmp)>0){
        N<-tmp$N
        Ids<-tmp$PS
      }else{
        N<-0
        Ids<-NULL
      }
      data.frame(`Representative_Models`=N, Model_ID=paste(Ids,collapse = ', '))
    } 
  })
  
  output$GlobalPieChart<-renderPlot({
    
    if(length(NT$data$TreeRoot)>0 & length(SELECTEDNODE$data)>0){
      
      currentGlobalSupport <- NT$data$navTable$GlobalSupport[SELECTEDNODE$data]
      supportingPatients <- NT$data$navTable$positivePoints[SELECTEDNODE$data]
      
      if (input$subSet!=''){
        
        if(!input$checkboxNegation){
          nn <- length(which(TUMOURS$data[input$subSet,]>0))
          patients <- rep(0,nn)
          names(patients) <- names(which(TUMOURS$data[input$subSet,]>0))
        }else{
          nn <- length(which(TUMOURS$data[input$subSet,]==0))
          patients <- rep(0,nn)
          names(patients) <- names(which(TUMOURS$data[input$subSet,]==0))
        }
      }else{
        nn <- ncol(TUMOURS$data)
        patients <- rep(0,nn)
        names(patients) <- colnames(TUMOURS$data)
        
      }
      
      
      
      supportingPatients <- unlist(str_split(supportingPatients,','))
      patients[supportingPatients] <- 1
      
      tmpCol <- Get(Traverse(NT$data$TreeRoot,traversal = 'level'),'Colors')
      
      non <- names(tmpCol)
      non <- str_split(non,' ')
      non <- unlist(lapply(non,function(x){x[[1]][1]}))
      
      id <- which(non==SELECTEDNODE$data)
      
      CCOL<-tmpCol[id]
      
      if(is.na(CCOL[1])){
        CCOL<-'darkgray'
      }
      
      CCCOL<-rep('lightgray',length(patients))
      CCCOL[which(patients>0)]<-CCOL
      
      
      
      plotCoords<-beeswarm(PATIENTcoords$data[names(patients)],method = 'hex',xaxt='n',yaxt='n',pwcol = CCCOL,pch=16,cex=1.5,
               main=paste("SubType ", SELECTEDNODE$data,' (', format(100*currentGlobalSupport, digits=3),'% of ', nn, ' patients)',
                          sep=''))
      
      plot(plotCoords$x,plotCoords$y,xaxt='n',yaxt='n',col = CCCOL,pch=16,cex=1.5,frame.plot = FALSE,xlab='',ylab='',
           main=paste("SubType ", SELECTEDNODE$data,' (', format(100*currentGlobalSupport, digits=3),'% of ', nn, ' patients)',
                                                                                             sep=''))
      
    } 
  })
  
  output$ComplementPieChart<-renderPlot({
     
      if(length(NT$data$TreeRoot)>0 & length(SELECTEDNODE$data)>0){
        
        res<-CELLector_App.complementarPieChart(NT$data$TreeRoot,NT$data$navTable,SELECTEDNODE$data)
        
        
        COLORS<-res$COLORS
        res<-res$supports
        
        
        if(is.na(COLORS[1])){
          COLORS<-'darkgray'
        }
        
        if(is.element(100,res)){res<-100}
        pie3D(res,explode = 0,labels =names(res),labelcex = 0.8,
              col=c(COLORS,NA),
              radius = 1,las=2,main = paste('Selected subpopulation \nwith respect to the entire cohort (total ',
                                            NT$data$navTable$CurrentTotal[SELECTEDNODE$data],')',sep=''))
      }
    }
    )
  
  
}
