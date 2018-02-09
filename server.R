
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
  
  output$str <- renderPrint({
   if(!is.null(NT$data)){
      CELLlineData$data<-CELLector.CellLine.BEMs[[input$selectCancerType]]
      TUMOURS$data<-CELLector.PrimTum.BEMs[[input$selectCancerType]]
      pander::pander(paste(dim(CELLlineData$data)[1],input$selectCancerType,'cell lines and ',
                           dim(TUMOURS$data)[1],'patients considered in this session'))
      
    }else{
      pander::pander('Build CELLector Search Space to START')
   }
  })
  
  TUMOURS <- reactiveValues(data = NULL)
  FEATURES <- reactiveValues(data = NULL)
  CTYPE <- reactiveValues(data = NULL)
  NT <- reactiveValues(data = NULL)
  
  RULES <- reactiveValues(data = NULL)
  SELECTEDNODE <- reactiveValues(data = NULL)
  SIGNATURES <- reactiveValues(data = NULL)
  encodedSIGNATURES <- reactiveValues(data = NULL)
  CELLlineData <- reactiveValues(data = NULL)
  STATUS <- reactiveValues(data = NULL)
  SunBurstSequences <- reactiveValues(data = NULL)
  
  
  
  observeEvent(input$action, {
    
    SELECTEDNODE$data <- NULL
    TUMOURS$data <- CELLector.PrimTum.BEMs[[input$selectCancerType]]
    FEATURES$data <- colnames(TUMOURS$data)
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
    if (input$whereToNeglect=='While building search space' | input$whereToNeglect=='Always'){
      toRemove <- input$toExclude
    }else{ toRemove <- NULL
    }

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

      #..............
      if (input$whereToNeglect=='While selecting cell lines' | input$whereToNeglect=='Always'){
        toRemove <- input$toExclude

        if(length(toRemove)>1){
          CLD <- CLD[which(rowSums(CLD[,toRemove])==0),]
        }else{
          if (length(toRemove)>0){
            CLD <- CLD[which(CLD[,toRemove]==0),]
          }
        }
      }
      #..............
      if (input$whatToInclude2!='Everything'){

        if (input$whatToInclude2=='Microsatellite stable'){
          CLD <- CLD[which(MSI_STATUS[as.character(CIDS)]==0),]
        } else{
          CLD <- CLD[which(MSI_STATUS[as.character(CIDS)]==1),]
        }
      }
      #..............

      CELLlineData$data <- CLD
    }
    #..............

    progress$set(message = "Done!", value = 1)
    progress$close()
    #..............
    
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
      
      tmp <- CELLector.solveFormula(encodedSIGNATURES$data[[SELECTEDNODE$data]],dataset = CELLlineData$data)
      
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
      
      nn <- nrow(TUMOURS$data)
      patients <- rep(0,nn)
      names(patients) <- rownames(TUMOURS$data)
      supportingPatients <- unlist(str_split(supportingPatients,','))
      patients[supportingPatients] <- 1
      
      tmpCol <- Get(Traverse(NT$data$TreeRoot,traversal = 'level'),'Colors')
      non <- names(tmpCol)
      non <- str_split(non,' ')
      non <- unlist(lapply(non,function(x){x[[1]][1]}))
      
      id <- which(non==SELECTEDNODE$data)
      
      polar.plot(patients, 1:nn, main=paste("SubType ", SELECTEDNODE$data, ' (', format(100*currentGlobalSupport, digits=3),'% of ', nrow(TUMOURS$data), ' patients)',sep=''), lwd=3, line.col=tmpCol[id],rp.type = 'r', labels=NULL,show.grid=FALSE, show.radial.grid=FALSE, show.grid.labels=FALSE)
      
    } else{
      if(!is.null(NT$data$navTable)){
        
        currentGlobalSupport <- NT$data$navTable$GlobalSupport[1]
        supportingPatients <- NT$data$navTable$positivePoints[1]
        
        nn<-nrow(TUMOURS$data)
        patients<-rep(0,nn)
        names(patients)<-rownames(TUMOURS$data)
        supportingPatients<-unlist(str_split(supportingPatients,','))
        patients[supportingPatients]<-1
        
        polar.plot(patients,1:nn,
                   main=paste("SubType ", 1, ' (', format(100*currentGlobalSupport, digits=3),'% of ',
                              nrow(TUMOURS$data),' patients)',sep=''),lwd=3,
                   line.col='blue',rp.type = 'r',
                   labels=NULL,show.grid=FALSE,
                   show.radial.grid=FALSE,show.grid.labels=FALSE)
      }
    }
  })
  
  output$ComplementPieChart<-renderPlot({
    
    if(length(NT$data$TreeRoot)>0 & length(SELECTEDNODE$data)>0){
      res<-CELLector_App.complementarPieChart(NT$data$TreeRoot,NT$data$navTable,SELECTEDNODE$data)
      COLORS<-res$COLORS
      res<-res$supports
      
      pie3D(res,explode = 0,labels =names(res),labelcex = 0.8,
            col=c(COLORS,NA),
            radius = 1,las=2,main = paste('Selected subpopulation \nwith respect to the entire cohort (total ',
                                          NT$data$navTable$CurrentTotal[SELECTEDNODE$data],')',sep=''))
    }
  })
  
  
}
