
ui <- fluidPage(
  tags$head(tags$script(src = "message-handler.js")),
  navbarPage(theme = shinytheme("yeti"), # other options: "cosmo","yeti","lumen"
             title = "CELLector",
             tabPanel(title = "Select",
                      
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
                      verbatimTextOutput("str"),
                      fluidRow(
                        # - - - - - - - - - - - - - - - - - - - - -
                        column(3,
                               wellPanel(
                                 img(src="cellcultures.jpg", height = 100, width = 240, align = "middle"),
                                 br(), br(),
                                 selectInput("selectCancerType", label = h5("Cancer Type:"),
                                             choices = sort(as.character(TCGALabels)), selected = 'COREAD'),
                                 actionButton("action", label = tags$strong(em("Build Search Space"))),
                                 br(), br(),
                                 downloadButton("DownSearchSpace", label = "Download Search Space"),
                                 br(), br(),
                                 p(h6("Tutorial available", a("here.", href="https://drive.google.com/open?id=1vTachWWefdS8sWFunBDh6ssIfwzhivJ6" ,
                                                              target="_blank"))),
                                 p(h6("Underlying data available at the", a("GDSC1000 data portal.", href="http://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Home.html" ,
                                                              target="_blank")))
                                 
                                 
                               ),
                               
                                
                               fluidRow(column(12, 
                                                wellPanel(
                                                  
                                                  h5("Representative Cell Line Selection"),
                                                  
                                                  numericInput("N.CellLines", h5(strong(em("N. of Cell Lines to Select:"))), 10, min = 1, max = 100),
                                                  downloadButton("CELLect", label = "CELLect Cell Lines"),
                                                  br(), #br(),
                                                  
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
                                           radioButtons('whatToInclude', strong('Cancer Functional Events (CFEs) to consider:'),
                                                        choices = c('Mutations in high confidence cancer genes',
                                                                    'Recurrently CN altered chromosomal segments',
                                                                    'Both'),
                                                        selected = 'Both',
                                                        inline = FALSE)
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
                        column(3,
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
                        ),
                        
                        # - - - - - - - - - - - - - - - - - - - - -
                         column(6, 
                                  wellPanel(
                                    p(h5("Id decoding for recurrently CN altered chromosomal segments (values in Id column should be used in box 1.)")),
                                    dataTableOutput('cnaDecodeTable'),
                                     p(h6("Full decoding table available at the", a("GDSC1000 data portal",
                                                                                 href = "http://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources//Data/suppData/TableS2D.xlsx" , target="_blank")))
                                  ))
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
             )
  )
)

