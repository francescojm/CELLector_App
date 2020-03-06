
ui <- fluidPage(
  tags$head(tags$script(src = "message-handler.js")),
  navbarPage(theme = shinytheme("yeti"), # other options: "cosmo","yeti","lumen"
             title = "CELLector",
             tabPanel(title = "Select Cell Lines",
                      fluidRow(
                        column(4,
                               a(img(src="cti_ot_primary_logo_blk_hr.jpg", height = 57, width = 156, align = "middle"),
                                 href='https://www.opentargets.org/',target="_blank")
                               ),
                        column(4,
                               a(img(src="Wellcome_Sanger_Institute_Logo_Landscape_Digital_RGB_Monotone_Black.png",
                                     height = 58, width = 169, align = "middle"),
                                  href='http://www.sanger.ac.uk/',target="_blank")
                               ),
                        column(4,
                               a(img(src="EMBL_EBI_Logo_black.jpg",
                                   height = 46, width = 155, align = "middle"),
                                 href='https://www.ebi.ac.uk/',target="_blank")
                               )
                      ),
                      br(),
                      titlePanel(div(HTML("Genomics Guided Selection of Cancer <em>in vitro</em> Models"))),
                      
                      hr(),
                      h3('v1.0.0 (beta)'),
                      p(h5(a("Tutorial", href="https://www.dropbox.com/s/djaaj2b33hqv4w1/Supplemental_Information.pdf?dl=1" ,
                             target="_blank"))),
                      p(h5(a("Code", href="https://github.com/francescojm/CELLector_App" ,
                             target="_blank"))),
                      p(h5(a("CELLector R Package", href="https://github.com/francescojm/CELLector" ,
                                                             target="_blank"))),
                      p(h5(a("CELLector R Package interactive vignette", href="https://rpubs.com/francescojm/CELLector" ,
                                                              target="_blank"))),
                      verbatimTextOutput("str"),
                      fluidRow(
                        # - - - - - - - - - - - - - - - - - - - - -
                        column(3,
                               wellPanel(
                                 img(src="cellcultures.jpg", height = 100, width = 240, align = "middle"),
                                 br(), br(),
                                 wellPanel(
                                   checkboxInput('UDgenomic', 'User Defined Binary genomic Event Matrices (BEMs)', FALSE),
                                   conditionalPanel(condition = 'input.UDgenomic',
                                                   
                                                     fileInput("ud_tumourBEMs", "Select Primary Tumour BEMs file:",
                                                            multiple = FALSE,
                                                            accept = c(".RData")),
                                                  
                                                     fileInput("ud_cellLineBEMs", "Select Cell Line BEMs file:",
                                                            multiple = TRUE,
                                                            accept = c(".RData")),
                                                     actionButton('ValidateUpdateBEMs','Validate and update BEMs')
                                                   ),
                                   conditionalPanel(condition = '!input.UDgenomic',
                                                    selectInput("selectCancerType", label = h5("Cancer Type:"),
                                                                choices = sort(as.character(TCGALabels)), selected = 'COREAD'),
                                                    p(h6("Underlying data available at the", a("GDSC1000 data portal.", href="http://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Home.html" ,
                                                                                               target="_blank")))
                                   )
                                 ),
                                 actionButton("action", label = tags$strong(em("Build Search Space"))),
                                 br(), br(),
                                 downloadButton("DownSearchSpace", label = "Download Search Space"),
                                 br(), br()
                               ),
                               
                               fluidRow(column(12, 
                                                wellPanel(
                                                  h5("Representative Cell Line Selection"),
                                                  
                                                  numericInput("N.CellLines", h5(strong(em("N. of Cell Lines to Select:"))), 10, min = 1, max = 100),
                                                  downloadButton("CELLect", label = "CELLect Cell Lines"),
                                                  
                                                  br(), 
                                              
                                                  
                                                  downloadButton("score", "Score cell lines"),
                                                  downloadButton("subTypeMap", "Cell lines SubTypes Map"),
                                                      
                                                  sliderInput("scoreAlpha",min = 0,max=1,
                                                              strong(paste("\nSignature Length weight (= 1 - n. Patients score weight)")),
                                                              value = 0.75,
                                                              step = 0.05,round = TRUE),
                                                  
                                                  p(h6("Please cite:",
                                                       a("Najgebauer et al. 2018 - CELLector: Genomics Guided Selection of Cancer in vitro Models", href = "https://www.biorxiv.org/content/early/2018/03/03/275032" , target="_blank")))
                                                )) )
                        ),
                        
                        # - - - - - - - - - - - - - - - - - - - - -
                        column(9,
                                wellPanel(
                                  h5("Primary Tumours: Subtyping Criteria"),
                                  fluidRow(
                                    column(5,
                                           conditionalPanel("!input.UDgenomic",
                                             radioButtons('whatToInclude', strong('Cancer Functional Events (CFEs) to consider:'),
                                                          choices = c('Mutations in high confidence cancer genes',
                                                                      'Recurrently CN altered chromosomal segments',
                                                                      'Both'),
                                                          selected = 'Both',
                                                          inline = FALSE),
                                             checkboxGroupInput("useMeth",label=NULL,
                                                                choices = list("Include Methylation data" = 1),
                                                                selected = NULL)
                                           )
                                    ), 
                                    column(7,
                                           fluidRow(
                                             column(6,sliderInput("minSetSize",
                                                                  strong("Alteration set size:"), 1, 5, value = 1,step = 1,round = TRUE)),
                                             column(6,
                                                    sliderInput("minGlobalSupport",
                                                                strong("Global support (%):"), 1, 50, value = 5,step = 1,round = TRUE)
                                             )
                                           )
                                    )
                                  )
                                )
                        ),
                        # - - - - - - - - - - - - - - - - - - - - -
                        column(3,conditionalPanel("!input.UDgenomic",
                                wellPanel(
                                  p(h5("Supervised Search Space Construction")),
                                  br(),
                                  # - - - - - - - - - - - - - - - - - - - - -
                                  fluidRow(
                                    column(6,
                                           selectizeInput(
                                             'subSet',
                                             label = "1. Define subcohort based on the status of an individual CFE:",
                                             choices = c('', features),
                                             selected = '',
                                             options = list(create = TRUE, maxItems = 1)
                                           ),

                                           checkboxInput("checkboxNegation", label = "wild-type", value = FALSE),
                                           
                                          selectizeInput(
                                            'pathFocus',
                                            label = "2. Focus on CFEs in cancer pathways (max 3):",
                                            choices = c('',pathways),
                                            selected="",
                                            options = list(create = TRUE, maxItems = 3)
                                          ),

                                           radioButtons('whatToInclude2', '3. Consider only cell lines that are:',
                                                        choices = c('Microsatellite stable',
                                                                    'Microsatellite instable',
                                                                    'All'),
                                                        selected = 'All',
                                                        inline = FALSE)
                                          )
                                    )
                                
                                  )
                        
                                )
                        ),
                        
                        # - - - - - - - - - - - - - - - - - - - - -
                         column(6,conditionalPanel("!input.UDgenomic",
                                checkboxInput('showCNAdecode', 'Show CNA id decoding table', TRUE),
                                conditionalPanel(condition = 'input.showCNAdecode',
                                
                                  wellPanel(
                                    p(h5("Id decoding for recurrently CN altered chromosomal segments (values in Id column can be used in box 1)")),
                                    dataTableOutput('cnaDecodeTable'),
                                     p(h6("Full decoding table available at the", a("GDSC1000 data portal",
                                                                                 href = "http://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources//Data/suppData/TableS2D.xlsx" , target="_blank")))
      
                                    )
                                  
                                  ),
                                checkboxInput('showHMSdecode', 'Show HyperMeth. id decoding table', FALSE),
                                conditionalPanel(condition = 'input.showHMSdecode',
                                wellPanel(
                                  p(h5("Id decoding for informative CpG island hyper-methylations (values in Id column can be used in box 1)")),
                                  dataTableOutput('hmsDecodeTable'),
                                  p(h6("Full decoding table available at the", a("GDSC1000 data portal",
                                                                                 href = "http://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources//Data/suppData/TableS2D.xlsx" , target="_blank")))
                                  )
                                
                                )
                         
                                )
                         )
                      ),
                      # # - - - - - - - - - - - - - - - - - - - - -
                      actionButton('changeColors','Use Different Color Scheme'),
                      fluidRow(
                        column(9,
                               wellPanel(
                               collapsibleTreeOutput("plot"))
                               ),
                        column(3,
                               sunburstOutput("sunburst"))
                      ),
                      
                      fluidRow(   
                         column(6,
                                tableOutput('NodeDetails'),
                                tableOutput('CellLineDetails')),
                         column(3,plotOutput('GlobalPieChart')),
                         column(3,plotOutput('ComplementPieChart'))
                      )
             ),
             tabPanel(title = "BEM builder",
                      fluidRow(
                        column(4,
                               a(img(src="cti_ot_primary_logo_blk_hr.jpg", height = 57, width = 156, align = "middle"),
                                 href='https://www.opentargets.org/',target="_blank")
                        ),
                        column(4,
                               a(img(src="Wellcome_Sanger_Institute_Logo_Landscape_Digital_RGB_Monotone_Black.png",
                                     height = 58, width = 169, align = "middle"),
                                 href='http://www.sanger.ac.uk/',target="_blank")
                        ),
                        column(4,
                               a(img(src="EMBL_EBI_Logo_black.jpg",
                                     height = 46, width = 155, align = "middle"),
                                 href='https://www.ebi.ac.uk/',target="_blank")
                        )
                      ),
                      br(),
                      titlePanel(div(HTML("Genomic Binary Event Matrix (BEM) builder"))),
                      hr(),
                      h3('v1.0.0 (beta)'),
                      column(6,
                             wellPanel(
                               h4(strong("Primary tumours")),
                               radioButtons('USE_whatTumData',"Variant catalogue to consider:",
                                            choices = c('Use curated TCGA data from Iorio et al. 2016',
                                                        'Upload Variants'),
                                            selected = 'Use curated TCGA data from Iorio et al. 2016'),
                               conditionalPanel(condition = "input.USE_whatTumData == 'Use curated TCGA data from Iorio et al. 2016'",
                                                wellPanel(
                                                  selectInput("TCGA_selectCancerType", label = h5("Cancer Type:"),
                                                              choices = sort(unique(CELLector.PrimTumVarCatalog$Cancer.Type)),
                                                              selected = sort(unique(CELLector.PrimTumVarCatalog$Cancer.Type))[1]),
                                            
                                                  radioButtons("TCGA_genes", 'Genes to consider:',
                                                               choices = c('All','Iorio et al. 2016 drivers','CMPs drivers','User defined list'),
                                                               selected = 'Iorio et al. 2016 drivers',
                                                               inline = TRUE),
                                                   
                                                  conditionalPanel(condition = "input.TCGA_genes == 'User defined list'",
                                                                    fileInput("TCGA_ud_geneList",
                                                                              "Upload a plain .txt file with Gene identifiers (one per line):",
                                                                              multiple = FALSE,
                                                                              accept = c(".txt"))),
                                                   
                                                  radioButtons("TCGA_variants", 'Variants to consider:',
                                                                choices = c('All','Iorio et al. 2016 variants (COSMIC filtered)','User defined list'),
                                                                selected = 'Iorio et al. 2016 variants (COSMIC filtered)',
                                                                inline = TRUE),
                                                   
                                                  conditionalPanel(condition = "input.TCGA_variants == 'User defined list'",
                                                                    fileInput("TCGA_ud_variants",
                                                                              "Upload a plain .txt file with tab separated [Gene identifier - variant] (one pair per line):",
                                                                              multiple = FALSE,
                                                                              accept = c(".txt")))
                                                )
                               ),
                             conditionalPanel(condition = "input.USE_whatTumData !='Use curated TCGA data from Iorio et al. 2016'",
                                              wellPanel(
                                                fileInput("TCGA_ud_variant_catalogue",
                                                          "Upload the variant catalogue as a plain tab separated .txt file:",
                                                          multiple = FALSE,
                                                          accept = c(".txt"))
                                              )
                             ),
                             fluidRow(
                               column(4,actionButton('TGCA_BEM_generation','Make new BEM')),
                               column(4,radioButtons('TGCA_BEM_AS_WHAT', strong('BEM file format:'),
                                                     choices = c('R object','.tsv'),
                                                     selected = 'R object',
                                                     inline = TRUE)),
                               column(4,downloadButton('SAVE_TCGA_BEM','Save BEM'))
                             ),
                             verbatimTextOutput("str_UD_TCGA_BEM_STATUS")
                             )
                      ),
                      column(6,
                             wellPanel(
                               h4(strong("in-vitro models")),
                               radioButtons('USE_cellModelPassports',"Variant catalogue to consider:",
                                            choices = c('Use Variants Catalogue from Cell Model Passports (CMPs)',
                                                        'Upload Variants'),
                                            selected = 'Use Variants Catalogue from Cell Model Passports (CMPs)'),
                               
                               conditionalPanel(condition = "input.USE_cellModelPassports == 'Use Variants Catalogue from Cell Model Passports (CMPs)'",
                                                wellPanel(
                                                  
                                                  fluidRow(
                                                    column(3,
                                                           checkboxInput('CMP_exclude_organoids', 'Exclude organoids', FALSE),
                                                           checkboxInput('CMP_human_samples_only', 'Human derived only', TRUE),
                                                           wellPanel(
                                                           checkboxInput('CMP_age_at_sampling', 'Filter based on age at sampling', TRUE),
                                                           conditionalPanel(condition = "input.CMP_age_at_sampling",
                                                                            uiOutput("CMP_age_at_sampling_slide_uiOutput")
                                                           ))
                                                           ),
                                                    column(3,radioButtons('CMP_gender',label = 'Gender:',choices = c('Male','Female','All (including Unknown)'),
                                                                          selected = 'All (including Unknown)'),
                                                           wellPanel(
                                                           checkboxInput('CMP_based_on_etnicity',label='Filter based on etnicity'),
                                                           conditionalPanel(condition = "input.CMP_based_on_etnicity",
                                                                            uiOutput("CMP_etnicity_uiOutput")
                                                           ))
                                                           ),
                                                    column(3,radioButtons('CMP_msi_status',label = 'MSI status:',
                                                                          choices = c('MSS','MSI','All (including NA)'),
                                                                          selected = 'All (including NA)')),
                                                    column(3,
                                                           wellPanel(
                                                           checkboxInput('CMP_based_on_mut_burden',label='Filter based on mutation burden',value = TRUE),
                                                           conditionalPanel(condition = "input.CMP_based_on_mut_burden",
                                                                            uiOutput("CMP_mutBurdend_slide_uiOutput")
                                                                            )),
                                                           wellPanel(
                                                           checkboxInput('CMP_based_on_ploidy',label='Filter based on ploidy',value = TRUE),
                                                           conditionalPanel(condition = "input.CMP_based_on_ploidy",
                                                                            uiOutput("CMP_ploidy_slide_uiOutput")
                                                                            ))
                                                          )
                                                    
                                                  ),
                                                  
                                                  uiOutput("CMP_N_cell_lines"),
                                                  selectInput("CMP_selectTissue", label = h5("Tissue:"),
                                                              choices = sort(unique(CMPs_model_annotations$tissue)),
                                                              selected = sort(unique(CMPs_model_annotations$tissue))[1]),
                                                  uiOutput("CMP_selectCancerType_uiOutput"),
                                                  uiOutput("CMP_selectCancerType_details_uiOutput"),
                                                  uiOutput("CMP_selectSample_site_uiOutput"),
                                                  
                                                  radioButtons("CMP_genes", 'Genes to consider:',
                                                               choices = c('All','Iorio et al. 2016 drivers','CMPs drivers','User defined list'),
                                                               selected = 'Iorio et al. 2016 drivers',
                                                               inline = TRUE),
                                                  
                                                  conditionalPanel(condition = "input.CMP_genes == 'User defined list'",
                                                                   fileInput("CMP_ud_geneList",
                                                                             "Upload a plain .txt file with Gene identifiers (one per line):",
                                                                             multiple = FALSE,
                                                                             accept = c(".txt"))),
                                                  
                                                  radioButtons("CMP_variants", 'Variants to consider:',
                                                               choices = c('All','Iorio et al. 2016 variants (COSMIC filtered)','User defined list'),
                                                               selected = 'Iorio et al. 2016 variants (COSMIC filtered)',
                                                               inline = TRUE),
                                                  
                                                  conditionalPanel(condition = "input.CMP_variants == 'User defined list'",
                                                                   fileInput("CMP_ud_variants",
                                                                             "Upload a plain .txt file with tab separated [Gene identifier - variant] (one pair per line):",
                                                                             multiple = FALSE,
                                                                             accept = c(".txt")))
                                                  )
                                                ),
                               conditionalPanel(condition = "input.USE_cellModelPassports != 'Use Variants Catalogue from Cell Model Passports (CMPs)'",
                                                wellPanel(
                                                  fileInput("CMP_ud_variant_catalogue",
                                                            "Upload the variant catalogue as a plain tab separated .txt file:",
                                                            multiple = FALSE,
                                                            accept = c(".txt"))
                                                          )
                                                ),
                              
                               fluidRow(
                                 column(4,actionButton('CL_BEM_generation','Make new BEM')),
                                 column(4,radioButtons('CL_BEM_AS_WHAT', strong('BEM file format:'),
                                                                     choices = c('R object','.tsv'),
                                                                     selected = 'R object',
                                                                     inline = TRUE)),
                                 column(4,downloadButton('SAVE_CL_BEM','Save BEM'))
                               ),  
                               
                            
                               verbatimTextOutput("str_UD_CL_BEM_STATUS")
                              
                               )
                      
                             )
                      )
  
             )

  )

