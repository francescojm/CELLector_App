
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
server <- function(input, output, session) {
  
  
  
  output$cnaDecodeTable <- renderDataTable(CNAid_decode[CNAid_decode$CancerType==input$selectCancerType,],
                                     options = list(
                                       pageLength = 5,autoWidth = FALSE)
                                     )
  output$hmsDecodeTable <- renderDataTable(HMSid_decode[HMSid_decode$CancerType==input$selectCancerType,],
                                           options = list(
                                             pageLength = 5,autoWidth = FALSE)
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
        nc<-length(setdiff(CL,CTE))
         
      }else{
        if(input$whatToInclude2=='Microsatellite instable'){
           CTE<-rownames(CELLlineData$data)[match(names(which(CELLector.MSIstatus!="MSI-H")),COSMICids$data)]
           nc<-length(setdiff(CL,CTE))
        }else{
           CTE<-NULL
        } 
      }
     
     pander::pander(paste(nc,'cell lines and ',
                          nn,'patients considered in this session\n\n[Re-build CELLector Search Space to Make any change to the criteria below effective]'))
     
    }else{
      pander::pander('Build CELLector Search Space to START')
   }
  })
  
  
  output$str_UD_CL_BEM_STATUS <- renderPrint({
    if(length(UD_CL_BEM$data)==0){
      pander::pander('No in-vitro models BEM built yet.')
      
      }else{
      
        tmp<-as.matrix(UD_CL_BEM$data[,3:ncol(UD_CL_BEM$data)])
        pander::pander(paste(paste(UD_CL_BEM_TISSUE$data,': ',
                                   paste(UD_CL_BEM_CTYPE$data,collapse=', '),
                                   '\nBEM built and ready to be saved.\n',
                                   dim(tmp)[1],' in-vitro models x ',
                                   dim(tmp)[2],' mutated genes',sep='')))
        
      }})
  
  output$str_UD_TCGA_BEM_STATUS <- renderPrint({
    if(length(UD_TUM_BEM$data)==0){
      pander::pander('No tumours BEM built yet.')
    }else{
      
      pander::pander(paste(paste(UD_TUM_BEM_CTYPE$data,': ',
                                 '\nBEM built and ready to be saved.\n',
                                 dim(UD_TUM_BEM$data)[2],' tumours x ',
                                 dim(UD_TUM_BEM$data)[1],' mutated genes',sep='')))
      
    }})
  
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
  
  
  UserDefinedPrimTum<-reactiveValues(data = NULL)
  UserDefinedCellLines<-reactiveValues(data = NULL)
  
  UD_CL_BEM<-reactiveValues(data = NULL)
  UD_CL_BEM_TISSUE<-reactiveValues(data = NULL)
  UD_CL_BEM_CTYPE<-reactiveValues(data = NULL)
  
  UD_TUM_BEM<-reactiveValues(data = NULL)
  UD_TUM_BEM_CTYPE<-reactiveValues(data = NULL)
  
  STATUS <- reactiveValues(data = NULL)
  
  
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
         res<-CELLector.makeSelection(modelMat = CELLector.buildModelMatrix(encodedSIGNATURES$data,
                                                                            CELLlineData$data,
                                                                            NT$data$navTable),
                                      input$N.CellLines,
                                      NT$data$navTable)
        )
         
        write.table(res, file,sep='\t', row.names = FALSE, quote=FALSE)    
      }
    }
  )
  
  observeEvent(input$action, {
  
    if(!input$UDgenomic){
      if(length(input$useMeth)>0){
        CELLector.PrimTum.BEMs<-OrPrimTumBEMs_v2[2:length(OrPrimTumBEMs)]
        CELLlineData$data<-OrCellLineBEMs_v2[[input$selectCancerType]]
      }else{
        CELLector.PrimTum.BEMs<-OrPrimTumBEMs[2:length(OrPrimTumBEMs)]
        CELLlineData$data<-OrCellLineBEMs[[input$selectCancerType]]
      }
      TUMOURS$data <- CELLector.PrimTum.BEMs[[input$selectCancerType]]
      TUMOURS$data <- CELLector.unicizeSamples(TUMOURS$data)
    
      toInclude<-CELLlineData$data$COSMIC_identifier
      
      if(input$whatToInclude2=='Microsatellite stable'){
        toInclude<-names(which(CELLector.MSIstatus!="MSI-H"))
      }
      
      if(input$whatToInclude2=='Microsatellite instable'){
        toInclude<-names(which(CELLector.MSIstatus=="MSI-H"))
      }
      
      toInclude<-intersect(toInclude,CELLlineData$data$COSMIC_identifier)
      
      tmp<-setdiff(input$whatToInclude2,'All')
      
      if(length(toInclude)<2){
        if(length(toInclude)==0){
          
          
          message=paste('Searching Space not built: no',tmp,
                        input$selectCancerType,'cell lines available')
        }else{
          message=paste('Searching Space not built: only 1',tmp,
                        input$selectCancerType,'cell line available:',
                        CELLlineData$data$CellLine[
                          match(as.character(toInclude),CELLlineData$data$COSMIC_identifier)])
          
        }
        session$sendCustomMessage(type = 'testmessage',
                                  message = message)
        return()
      }
      
      CELLlineData$data<-CELLlineData$data[match(toInclude,CELLlineData$data$COSMIC_identifier),]
      
      SELECTEDNODE$data <- NULL
      FEATURES$data <- rownames(TUMOURS$data)
      
      PATIENTcoords$data<-sample(ncol(TUMOURS$data))
      names(PATIENTcoords$data)<-colnames(TUMOURS$data)
      
      
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
      
      
      pathways<-input$pathFocus
      
      if(sum(is.element(pathways,names(CELLector.Pathway_CFEs)))!=length(pathways)){
        session$sendCustomMessage(type = 'testmessage',
                                  message = 'Searching Space not built: Invalid pathway identifier inserted')
        return()
      }
      
      progress <- shiny::Progress$new(style='notification')
      
      progress$set(message = "Loading primary tumour data and building CELLector Search Space... Please Wait", value = 0.5)
      
      
      NT$data<-CELLector.Build_Search_Space(ctumours = t(TUMOURS$data),
                                            minGlobSupp = input$minGlobalSupport/100,
                                            minlen=input$minSetSize,
                                            cancerType = input$selectCancerType,
                                            pathwayFocused = input$pathFocus,
                                            pathway_CFEs = CELLector.Pathway_CFEs,
                                            cnaIdMap = CELLector.CFEs.CNAid_mapping,
                                            cnaIdDecode = CELLector.CFEs.CNAid_decode,
                                            hmsIdDecode = CELLector.CFEs.HMSid_decode,
                                            includeHMS = length(input$useMeth)>0,
                                            cdg = CELLector.HCCancerDrivers,
                                            subCohortDefinition=input$subSet,
                                            NegativeDefinition=input$checkboxNegation,
                                            mutOnly = (input$whatToInclude=='Mutations in high confidence cancer genes'),
                                            cnaOnly = (input$whatToInclude=='Recurrently CN altered chromosomal segments'),
                                            FeatureToExclude=input$toExclude)
      
      
      if(nrow(NT$data$navTable)>0){
        
        S <- CELLector.createAllSignatures(NavTab = NT$data$navTable)
        
        SIGNATURES$data <- S$S
        encodedSIGNATURES$data <- S$ES
        
        CLD <- CELLector.CellLine.BEMs[[input$selectCancerType]]
        CIDS <- CLD$COSMIC_identifier
        rn <- CLD[,2]
        CLD <- as.matrix(CLD[,3:ncol(CLD)])
        rownames(CLD) <- rn
        
        
        
        progress$set(message = "Done!", value = 1)
        progress$close()
        }
      }else{
        if (length(UserDefinedPrimTum$data)==0 | length(UserDefinedCellLines$data)==0){
          showModal(modalDialog(
            title = "No Data found!",
            "Please upload and validate BEMs respectively for Primary Tumours and Cell Lines (or other in-vitro models) first."
          ))
          }else{
            CELLector.PrimTum.BEMs<-UserDefinedPrimTum$data
            CELLlineData$data<-UserDefinedCellLines$data
            TUMOURS$data <- CELLector.PrimTum.BEMs
            TUMOURS$data <- CELLector.unicizeSamples(TUMOURS$data)
             
            SELECTEDNODE$data <- NULL
            FEATURES$data <- rownames(TUMOURS$data)
             
            PATIENTcoords$data<-sample(ncol(TUMOURS$data))
            names(PATIENTcoords$data)<-colnames(TUMOURS$data)
             
            minLen <- input$minSetSize
            minGlobSupp <- input$minGlobalSupport/100
            
             
            progress <- shiny::Progress$new(style='notification')
             
            progress$set(message = "Building CELLector Search Space... Please Wait", value = 0.5)
             
             
            NT$data<-CELLector.Build_Search_Space(ctumours = t(TUMOURS$data),
                                                  minGlobSupp = input$minGlobalSupport/100,
                                                  minlen=input$minSetSize,
                                                  cancerType = 'User Defined',
                                                  cnaIdDecode = CELLector.CFEs.CNAid_decode,
                                                  hmsIdDecode = CELLector.CFEs.HMSid_decode,
                                                  includeHMS = FALSE,
                                                  cdg = colnames(TUMOURS$data),
                                                  FeatureToExclude=input$toExclude)
             
             
             if(nrow(NT$data$navTable)>0){
               
               S <- CELLector.createAllSignatures(NavTab = NT$data$navTable)
               
               SIGNATURES$data <- S$S
               encodedSIGNATURES$data <- S$ES
               
               CLD <- CELLector.CellLine.BEMs[[input$selectCancerType]]
               CIDS <- CLD$COSMIC_identifier
               rn <- CLD[,2]
               CLD <- as.matrix(CLD[,3:ncol(CLD)])
               rownames(CLD) <- rn
               
               
               
              
             progress$set(message = "Done!", value = 1)
             progress$close()
            }
          }
      }
    })

  
  observeEvent(input$changeColors, {
    if(length(NT$data)>0){
      NT$data<-CELLector.changeSScolors(NT$data)
    }
  })
  
  observeEvent(input$ValidateUpdateBEMs, {
    
    inFile_primTum <- input$ud_tumourBEMs
    inFile_cellLin <- input$ud_cellLineBEMs
    
    ErrFlag<-0

    
    if(length(inFile_primTum)==0 | length(inFile_cellLin)==0){
      ErrFlag<-1
      showModal(modalDialog(
        title = "Warning!",
        "Please select two 2 RData object files, containing Binary genomic Event Matrices for primary tumours and cell lines (respectively)."
      ))
    }else{
      if(tolower(file_ext(inFile_primTum$name))!='rdata' | tolower(file_ext(inFile_cellLin$name))!='rdata'){
        ErrFlag<-1
        showModal(modalDialog(
          title = "Warning!",
          "Wrong file format. Please select two RData object files"
          ))  
      }
      if(inFile_primTum$name==inFile_cellLin$name){
        ErrFlag<-1
        showModal(modalDialog(
          title = "Warning!",
          "Please select two distict files"
        ))  
      }
      
    }
    
    if(!ErrFlag){
      oldWS<-ls()
      load(inFile_primTum$datapath)
      primTum_object<-res
      rm(res)
      
      load(inFile_cellLin$datapath)
      cellLin_object<-res
      rm(res)
            
      if(!is.matrix(primTum_object)){
        ErrFlag<-1
        showModal(modalDialog(
          title = "Warning!",
          "The Binary genomic Event Matrices for primary Tumours shoud be a binary event matrix (BEM) modeling a cohort of cancer patients. With cancer functional events (CFEs) on the columns and sample identifers on the rows. See the documentation entry for the CELLector.PrimTum.BEMs object of the CELLector R package for further details"
        ))  
        }
      
      if(!ErrFlag){
        UserDefinedPrimTum$data<-primTum_object
        UserDefinedCellLines$data<-cellLin_object

        showModal(modalDialog(
          title = "Binary event matrix file names and format look good",
          "Build the CELLector search space to start"
        ))  
      }else{
        updateCheckboxInput(session,'UDgenomic', value = FALSE)
      }
    }else{
      updateCheckboxInput(session,'UDgenomic', value = FALSE)
    }
  })
  
  output$plot <- renderCollapsibleTree({
    
    if(length(NT$data)>0){
      
      if (nrow(NT$data$navTable)>1){
        
        CELLector.visualiseSearchingSpace(searchSpace = NT$data,CLdata=CELLlineData$data)
        
        } 
    }
  })
  
  output$sunburst <- renderSunburst({
    if(length(NT$data)>0){
    
      CELLector.visualiseSearchingSpace_sunBurst(NT$data)
        
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
         data.frame(`Patient SubType`=nodeIdx,
                    `Underlying Signature`=SIGNATURES$data[[itt]])
       }else{
         SELECTEDNODE$data<-1
         data.frame(`Patient SubType`=format(1,digits=1),
                    `Underlying Signature`=SIGNATURES$data[[1]])
       }
    }else{
       if(length(SIGNATURES$data)>0){
         SELECTEDNODE$data<-1
         data.frame(`Patient SubType`=format(1,digits=1),
                    `Underlying Signature`=SIGNATURES$data[[1]]) 
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
      data.frame(`n. Representative Cell Lines`=N, `Names`=paste(Ids,collapse = ', '))
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
        
        if(NT$data$navTable$Type[SELECTEDNODE$data]=='Left.Child')
        {cp<-paste('SubType',NT$data$navTable$Parent.Idx[SELECTEDNODE$data])}
        else{
          if(NT$data$navTable$Type[SELECTEDNODE$data]=='Right.Child'){
            cp<-paste('Complement of SubType',NT$data$navTable$Parent.Idx[SELECTEDNODE$data])  
          }else{
            cp<-paste('whole cohort')
          }
          
        }
        
        pie3D(res,explode = 0,labels =names(res),labelcex = 0.8,
              col=c(COLORS,NA),
              radius = 1,las=2,main = paste('Derived from: \n',cp, '\n(n. Patients = ',
                                            NT$data$navTable$CurrentTotal[SELECTEDNODE$data],')',sep=''))
      }
    }
    )
  
  output$score <- downloadHandler(
    filename = function(){
      paste("CELLector_scores_", input$selectCancerType, "_", Sys.Date(), ".tsv", sep = "")
    },
    content = function(file) {
      if(length(NT$data)>0){
        res<-CELLector.Score(NT$data$navTable,CELLlineData$data,input$scoreAlpha) 
        
        write.table(res, file,sep='\t', row.names = FALSE, quote=FALSE)    
      }
    }
  )
  
  output$subTypeMap <- downloadHandler(
    filename = function(){
      paste("CELLector_subTypeMap_", input$selectCancerType, "_", Sys.Date(), ".tsv", sep = "")
    },
    content = function(file) {
      if(length(NT$data)>0){
        #res<-CELLector.Score(NT$data$navTable,CELLlineData$data,input$scoreAlpha) 
        
        ### take all the signatures from the searching space
        Signatures <- CELLector.createAllSignatures(NT$data$navTable)
        
        ### mapping colorectal cancer cell lines onto the CELLector searching space
        ModelMat<-CELLector.buildModelMatrix(Signatures$ES,CELLlineData$data,NT$data$navTable)
        
        rownames(ModelMat)<-Signatures$S[rownames(ModelMat)]
        
        ModelMat<-cbind(rownames(ModelMat),ModelMat)
        colnames(ModelMat)[1]<-"Signature"
        
        write.table(ModelMat, file,sep='\t', row.names = FALSE, quote=FALSE)    
      }
    }
  )
  
  
  output$CMP_selectCancerType_uiOutput <- renderUI({
    choices_ = sort(unique(CMPs_model_annotations$cancer_type[which(CMPs_model_annotations$tissue==input$CMP_selectTissue)]))
    selectInput("CMP_selectCancerType", "Cancer Type:", choices = choices_,selected = choices_,multiple = TRUE)
  })
  
  output$CMP_selectCancerType_details_uiOutput <- renderUI({
    choices_ = sort(unique(CMPs_model_annotations$cancer_type_detail[which(CMPs_model_annotations$tissue==input$CMP_selectTissue & 
                                                                      is.element(CMPs_model_annotations$cancer_type,input$CMP_selectCancerType))]))
    selectInput("CMP_selectCancerType_detailed", "Cancer Type Details:", choices = choices_,selected = choices_,multiple = TRUE)
  })
  
  output$CMP_selectSample_site_uiOutput <- renderUI({
    choices_ = sort(unique(CMPs_model_annotations$sample_site[which(CMPs_model_annotations$tissue==input$CMP_selectTissue & 
                                                                 is.element(CMPs_model_annotations$cancer_type,input$CMP_selectCancerType) &
                                                                   is.element(CMPs_model_annotations$cancer_type_detail,input$CMP_selectCancerType_detailed))]))
    selectInput("CMP_selectSample_site", "Sample site:", choices = choices_,selected = choices_,multiple = TRUE)
  })
  
  output$CMP_mutBurdend_slide_uiOutput <- renderUI({
      ids = which(CMPs_model_annotations$tissue==input$CMP_selectTissue & 
                    is.element(CMPs_model_annotations$cancer_type,input$CMP_selectCancerType) &
                    is.element(CMPs_model_annotations$cancer_type_detail,input$CMP_selectCancerType_detailed) &
                    is.element(CMPs_model_annotations$sample_site,input$CMP_selectSample_site))
      
      range=range(CMPs_model_annotations$mutational_burden[ids],na.rm = TRUE)
      sliderInput("CMP_mutBurdend_slide", "Mutation burden range (n.Mut x Mb):",
                  min = round(range[1]), max = round(range[2]),
                  value = round(range))
  })
  
  output$CMP_ploidy_slide_uiOutput <- renderUI({
    ids = which(CMPs_model_annotations$tissue==input$CMP_selectTissue & 
                  is.element(CMPs_model_annotations$cancer_type,input$CMP_selectCancerType) &
                  is.element(CMPs_model_annotations$cancer_type_detail,input$CMP_selectCancerType_detailed) &
                  is.element(CMPs_model_annotations$sample_site,input$CMP_selectSample_site))
    
    range=range(CMPs_model_annotations$ploidy[ids],na.rm = TRUE)
    sliderInput("CMP_ploidy_slide", "Ploidy range:",
                min = round(range[1]), max = round(range[2]),
                value = round(range))
  })
  
  output$CMP_age_at_sampling_slide_uiOutput <- renderUI({
    ids = which(CMPs_model_annotations$tissue==input$CMP_selectTissue & 
                  is.element(CMPs_model_annotations$cancer_type,input$CMP_selectCancerType) &
                  is.element(CMPs_model_annotations$cancer_type_detail,input$CMP_selectCancerType_detailed) &
                  is.element(CMPs_model_annotations$sample_site,input$CMP_selectSample_site))
    
    range=range(CMPs_model_annotations$age_at_sampling[ids],na.rm = TRUE)
    sliderInput("CMP_age_at_sampling_slide", "Age at sampling range:",
                min = round(range[1]), max = round(range[2]),
                value = round(range))
  })
  
  output$CMP_etnicity_uiOutput<-renderUI({
    ids = which(CMPs_model_annotations$tissue==input$CMP_selectTissue & 
                  is.element(CMPs_model_annotations$cancer_type,input$CMP_selectCancerType) &
                  is.element(CMPs_model_annotations$cancer_type_detail,input$CMP_selectCancerType_detailed) &
                  is.element(CMPs_model_annotations$sample_site,input$CMP_selectSample_site))
    
    ranges=unique(CMPs_model_annotations$ethnicity[ids])
    
    if(length(ranges)==1){
      h6(paste('Only one ethnic group value present (',ranges,')',sep=''))
    }else{
      selectInput("CMP_etnicity", paste(length(ranges),' ethnic group values avaible, select what to exclude',sep=''), choices = ranges,selected = NULL,multiple = TRUE)  
    }
  })
  
  
  output$CMP_N_cell_lines<-renderUI({
    
    ids = CELLector_App.current_Model_ids(input,CMPs_model_annotations,TRUE)
    ids_all = CELLector_App.current_Model_ids(input,CMPs_model_annotations,FALSE)
  
    h4(paste(length(ids),' in-vitro models with genomic data available (out of ',
             length(ids_all),') based on current filter settings.', sep=''))
  })
  
  observeEvent(input$CL_BEM_generation, {
    if(input$USE_cellModelPassports=='Use Variants Catalogue from Cell Model Passports (CMPs)'){
      ids = CELLector_App.current_Model_ids(input,CMPs_model_annotations,TRUE)

      if(length(ids)>0){
        if(length(input$CMP_selectCancerType)>0){
          if(length(input$CMP_selectSample_site)>0){
            if(input$CMP_genes=='Iorio et al. 2016 drivers'){
              data(CELLector.HCCancerDrivers)
              genesToConsider<-CELLector.HCCancerDrivers
            }
            if(input$CMP_genes=='CMPs drivers'){
              genesToConsider<-CELLector.CMPs_getDriverGenes()
            }
            if(input$CMP_genes=='All'){
              genesToConsider<-NULL
            }
            if(input$CMP_genes=='User defined list'){
              fn<-input$CMP_ud_geneList
              genesToConsider<-unlist(read.table(fn$datapath,stringsAsFactors = FALSE))
              names(genesToConsider)<-NULL
            }
            
            if(input$CMP_variants=='Iorio et al. 2016 variants (COSMIC filtered)'){
              data(CELLector.RecfiltVariants)
              variantsToConsider<-CELLector.RecfiltVariants
            }
            if(input$CMP_variants=='All'){
              variantsToConsider<-NULL
            }
            if(input$CMP_variants=='User defined list'){
              fn<-input$CMP_ud_variants
              variantsToConsider<-unlist(read.table(fn$datapath,stringsAsFactors = FALSE,sep='\t'))
              names(variantsToConsider)<-NULL
            }
            
            progress <- shiny::Progress$new(style='notification')
            
            progress$set(message = "Building Genomic Binary Event matrix for In-vitro models... Please Wait", value = 0.5)
            
            GENDER<-input$CMP_gender
            if (GENDER == 'all (including Unknown'){
              GENDER <- NULL
            }
            
        
            
            if(input$CMP_gender=="All (including Unknown)"){
              gender<-NULL
            }else{
              gender<-input$CMP_gender
            }
            
            
            if(input$CMP_msi_status=="All (including NA)"){
              msi<-NULL
            }else{
              msi<-input$CMP_msi_status
            }
            
            if(input$CMP_age_at_sampling){
              age<-input$CMP_age_at_sampling_slide
            }else{
              age<-NULL
            }
            
            
            
            if(input$CMP_based_on_mut_burden){
               mutBurd<-input$CMP_mutBurdend_slide
             }else{
               mutBurd<-NULL
             }
            
            if(input$CMP_based_on_mut_burden){
              mutBurd<-input$CMP_mutBurdend_slide
            }else{
              mutBurd<-NULL
            }
            
            if(input$CMP_based_on_ploidy){
              ploidy<-input$CMP_ploidy_slide
            }else{
              ploidy<-NULL
            }
            
            if(input$CMP_based_on_etnicity){
              etno<-input$CMP_etnicity
            }else{
              etno<-NULL
            }
            
            BEM<-CELLector.CELLline_buildBEM(
              Tissue = input$CMP_selectTissue,
              Cancer_Type = input$CMP_selectCancerType,
              
              excludeOrganoids = input$CMP_exclude_organoids,
              humanonly = input$CMP_human_samples_only,
              Cancer_Type_details = input$CMP_selectCancerType_detailed,
              sample_site = input$CMP_selectSample_site,
              gender_select = gender,
              msi_status_select = msi,
              age_at_sampling = age,
              mutational_burden_th = mutBurd,
              ploidy_th = ploidy,
              ethnicity_to_exclude = etno,
              
              GenesToConsider = genesToConsider,
              VariantsToConsider = variantsToConsider
              )
            
            progress$set(message = "Done!", value = 1)
            progress$close()
            
            tmp<-as.matrix(BEM[,3:ncol(BEM)])
            
            showModal(modalDialog(paste("Binary Event Matrix created: ",
                                        dim(tmp)[1],
                                        ' in-vitro models x ',
                                        dim(tmp)[2],
                                        ' mutated genes.\nDensity: ',
                                        format(100*sum(c(tmp))/prod(dim(tmp)),digits = 3),'%',
                                        sep='')))
            
            UD_CL_BEM$data<-BEM
            UD_CL_BEM_TISSUE$data<-input$CMP_selectTissue
            UD_CL_BEM_CTYPE$data<-input$CMP_selectCancerType
          }
        }
        else{showModal(modalDialog(
          title = "Warning!",
          paste("Select a sample site first.")
        ))}
      }else{
        showModal(modalDialog(
          title = "Warning!",
          paste("0 in-vitro models available based on current filter settings.")
        ))
      }
    }
    else{
      if(length(input$CMP_ud_variant_catalogue)==0){
        showModal(modalDialog('Upload variant catalogue file first.'))
      }else{
        
        progress <- shiny::Progress$new(style='notification')
        
        progress$set(message = "Building Genomic Binary Event matrix for In-vitro models... Please Wait", value = 0.5)
        varCat<-read.table(input$CMP_ud_variant_catalogue$datapath,sep='\t',header=TRUE,stringsAsFactors = FALSE)
        
        if(length(intersect(c("model_name","model_id","gene_symbol"),colnames(varCat)))!=3){
          showModal(modalDialog('Incorrect file format: model_name and/or model_id and/or gene_symbol column not found!'))
          progress$close()
        }else{
           BEM<-CELLector.CELLline_buildBEM(varCat = varCat)  
           progress$set(message = "Done!", value = 1)
        
           progress$close()
        
           tmp<-as.matrix(BEM[,3:ncol(BEM)])
           showModal(modalDialog(paste("Binary Event Matrix created: ",
                                       dim(tmp)[1],
                                       ' in-vitro models x ',
                                       dim(tmp)[2],
                                       ' mutated genes.\nDensity: ',
                                       format(100*sum(c(tmp))/prod(dim(tmp)),digits = 3),'%',
                                       sep='')))
           UD_CL_BEM$data<-BEM
           UD_CL_BEM_TISSUE$data<-"User Defined"
           UD_CL_BEM_CTYPE$data<-NULL
        }
      }
    }
  })
  
  output$SAVE_CL_BEM <- 
    downloadHandler(
      filename = function(){
        if(input$CL_BEM_AS_WHAT=='R object'){paste("in-vitro_BEM_",Sys.Date(), ".RData", sep = "")}
        else{paste("in-vitro_BEM_",Sys.Date(), ".tsv", sep = "")}
        },
      content = function(file){
        res<-UD_CL_BEM$data
        #tmp<-as.matrix(res[,3:ncol(res)])
        #rownames(tmp)<-res[,2]
        #colnames(tmp)<-colnames(res)[3:ncol(res)]
        #tmp<-t(tmp)
        if(input$CL_BEM_AS_WHAT=='R object'){save(res,file = file)}
        else{write.table(res,file=file,quote=FALSE,sep='\t',row.names=FALSE)}
        }
      )
  
  observeEvent(input$TGCA_BEM_generation, {
    if(input$USE_whatTumData=='Use curated TCGA data from Iorio et al. 2016'){
      if(input$TCGA_genes=='Iorio et al. 2016 drivers'){
        data(CELLector.HCCancerDrivers)
        genesToConsider<-CELLector.HCCancerDrivers
      }
      if(input$TCGA_genes=='CMPs drivers'){
        genesToConsider<-CELLector.CMPs_getDriverGenes()
      }
      if(input$TCGA_genes=='All'){
        genesToConsider<-NULL
      }
      if(input$TCGA_genes=='User defined list'){
        fn<-input$TCGA_ud_geneList
        genesToConsider<-unlist(read.table(fn$datapath,stringsAsFactors = FALSE))
        names(genesToConsider)<-NULL
      }
      if(input$TCGA_variants=='Iorio et al. 2016 variants (COSMIC filtered)'){
        data(CELLector.RecfiltVariants)
        variantsToConsider<-CELLector.RecfiltVariants
      }
      if(input$TCGA_variants=='All'){
        variantsToConsider<-NULL
      }
      if(input$TCGA_variants=='User defined list'){
        fn<-input$TCGA_ud_variants
        variantsToConsider<-unlist(read.table(fn$datapath,stringsAsFactors = FALSE,sep='\t'))
        names(variantsToConsider)<-NULL
      }
         
         progress <- shiny::Progress$new(style='notification')
         
         progress$set(message = "Building Genomic Binary Event matrix for Tumours... Please Wait", value = 0.5)
         
         BEM<-CELLector.Tumours_buildBEM(
           Cancer_Type = input$TCGA_selectCancerType,
           GenesToConsider = genesToConsider,
           VariantsToConsider = variantsToConsider)
         progress$set(message = "Done!", value = 1)
         progress$close()
         
         tmp<-BEM
         
         showModal(modalDialog(paste("Binary Event Matrix created: ",
                                     dim(tmp)[2],
                                     ' Tumours x ',
                                     dim(tmp)[1],
                                     ' mutated genes.\nDensity: ',
                                     format(100*sum(c(tmp))/prod(dim(tmp)),digits = 3),'%',
                                     sep='')))
         
         UD_TUM_BEM$data<-BEM
         UD_TUM_BEM_CTYPE$data<-input$TCGA_selectCancerType
      }else{
       if(length(input$TCGA_ud_variant_catalogue)==0){
         showModal(modalDialog('Upload variant catalogue file first.'))
       }else{
         progress <- shiny::Progress$new(style='notification')
         progress$set(message = "Building Genomic Binary Event matrix for Tumours... Please Wait", value = 0.5)
         
         varCat<-read.table(input$TCGA_ud_variant_catalogue$datapath,sep='\t',header=TRUE,stringsAsFactors = FALSE)
         
         if(length(intersect(c("SAMPLE","gene_symbol"),colnames(varCat)))!=3){
           showModal(modalDialog('Incorrect file format: SAMPLE and/or gene_symbol column not found!'))
           progress$close()
         }else{
           BEM<-CELLector.Tumours_buildBEM(varCat = varCat)  
           progress$set(message = "Done!", value = 1)
           
           progress$close()
           
           tmp<-as.matrix(BEM[,3:ncol(BEM)])
           showModal(modalDialog(paste("Binary Event Matrix created: ",
                                       dim(tmp)[1],
                                       ' tumour samples x ',
                                       dim(tmp)[2],
                                       ' mutated genes.\nDensity: ',
                                       format(100*sum(c(tmp))/prod(dim(tmp)),digits = 3),'%',
                                       sep='')))
           UD_TUM_BEM$data<-BEM
           UD_TUM_BEM_CTYPE$data<-NULL
         }
       }
     }
    })
  
  output$SAVE_TCGA_BEM <- 
    downloadHandler(
      filename = function(){
        if(input$TGCA_BEM_AS_WHAT=='R object'){paste("tumours_BEM_",Sys.Date(), ".RData", sep = "")}
        else{paste("tumours_BEM_",Sys.Date(), ".tsv", sep = "")}
      },
      content = function(file){
        res<-UD_TUM_BEM$data
        if(input$TGCA_BEM_AS_WHAT=='R object'){save(res,file = file)}
        else{write.table(res,file=file,quote=FALSE,sep='\t',row.names=TRUE)}
      }
    )
  
}
