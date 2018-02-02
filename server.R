
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
  
  
  observeEvent(input$action, {
    
    SELECTEDNODE$data <- NULL
    TUMOURS$data <- CELLector.PrimTum.BEMs[[input$selectCancerType]]
    FEATURES$data <- colnames(TUMOURS$data)
    if(length(input$subSet)>0 & input$subSet!=''){
       TUMOURS$data <- TUMOURS$data[which(TUMOURS$data[,input$subSet]),]
       rownames(TUMOURS$data) <- as.character(1:nrow(TUMOURS$data))
    }
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
    # #.............. 
    
    #### other parameters to be added.
    #NT$data <- CELLector.Build_Search_Space(
    #                                    ctumours = TUMOURS$data,
    #                                    cancerType = input$selectCancerType,
    #                                    minlen = minLen,
    #                                    minGlobSupp = minGlobSupp/100,
    #                                    pathwayFocused = pathway,
    #                                    mutOnly = mutonly,
    #                                    cnaOnly = cnaonly,
    #                                    subCohortDefinition = 
    #                                    FeatureToExclude = toRemove)
    # #..............
    # if(nrow(NT$data$navTable)>0){
    #   
    #   S <- createAllSignatures(NavTab = NT$data$navTable)
    #   
    #   SIGNATURES$data <- S$S
    #   encodedSIGNATURES$data <- S$ES
    #   
    #   CLD <- CellLine.list[[which(TCGALabels$CancerType_TCGALabel==input$selectCancerType)]]
    #   CIDS <- CLD$COSMIC_identifier
    #   rn <- CLD[,2]
    #   CLD <- as.matrix(CLD[,3:ncol(CLD)])
    #   rownames(CLD) <- rn
    #   
    #   #..............
    #   if (input$whereToNeglect=='While selecting cell lines' | input$whereToNeglect=='Always'){
    #     toRemove <- input$toExclude
    #     
    #     if(length(toRemove)>1){
    #       CLD <- CLD[which(rowSums(CLD[,toRemove])==0),]    
    #     }else{
    #       if (length(toRemove)>0){
    #         CLD <- CLD[which(CLD[,toRemove]==0),]    
    #       }
    #     }
    #   }
    #   #.............. 
    #   if (input$whatToInclude2!='Everything'){
    #     
    #     if (input$whatToInclude2=='Microsatellite stable'){
    #       CLD <- CLD[which(MSI_STATUS[as.character(CIDS)]==0),]
    #     } else{
    #       CLD <- CLD[which(MSI_STATUS[as.character(CIDS)]==1),]
    #     } 
    #   }
    #   #..............
    #   
    #   CELLlineData$data <- CLD
    # }
    # #.............. 
    # 
    # progress$set(message = "Done!", value = 1)
    # progress$close()
    #..............
    
    # cnaLookUp$data <- cna_look_up(input$cnaID , input$selectCancerType)
    
    # CLG_feature$data <- get_cell_line_genomic_features(input$CLname, input$selectCancerType)
  })

}
