
ui <- fluidPage(
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
                                 actionButton("action", label = tags$strong(em("Build Search Space"))), #style="color: #fff"
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
                                                  
                                                  p(h6("Please cite:", a("xxxx", href = "https://github.com/francescojm/CELLector" , target="_blank")))
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
                                    # h3("Recurrent Copy Number Alterations"),
                                    
                                    # selectInput("cnaID", label = h4("Look up:"), 
                                    #              choices = c('', cnaidlist), 
                                    #             selected = ''),
                                    
                                    # p("XxXXXXXxxXxxidahfklshgishgishglskhgsytwhflalt"),
                                    # 
                                    # hr(),
                                    # fluidRow(column(12, verbatimTextOutput("value")))
                                    p(h6("Full decoding table available at the", a("GDSC1000 data portal",
                                                                                 href = "http://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources//Data/suppData/TableS2D.xlsx" , target="_blank")))
                                  ))
                      ),
                      # # - - - - - - - - - - - - - - - - - - - - -
                      
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
             
             
             # - - - - - - - - - - - - - - - - - - - - -
             # tabPanel(title = "Explore",
             #          titlePanel(em("Genomics Guided Selection of Cancer in vitro Models")),
             #          hr(),
             #          
             #          p(h5("Primary tumour and cell line genomic datasets are described in", a("Iorio et al., 2016", href="http://www.cell.com/cell/fulltext/S0092-8674(16)30746-2" , target="_blank"), "and relevant data were obtained from the accompanied", a("web-portal.", href="http://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Home.html" , target="_blank"))),
             #          br(),
             #          
             #          fluidPage(
             #            sidebarPanel(
             #              h4(em("Primary Tumour Genomic Features")),
             #              
             #              img(src="tumoursA.pdf", height = 130, width = 280, align = "middle"),
             #              br(), br(),
             #              
             #              # Not implemented
             #              # Input: Choose dataset ----
             #              selectInput("selectCancerType", label = em(h5("Select Cancer Type:")),
             #                          choices = sort(as.character(TCGALabels[,1])), selected = 'COREAD'),
             #              radioButtons("whatToInclude", "Feature type:",
             #                           choices = c("Mutations", "CNAs", "Both")),
             #              # Button
             #              downloadButton("downloadData", "Download")
             #              
             #            ),
             #            
             #            # - - - - - - - - - - - - - - - - - - - - -
             #            sidebarPanel(
             #              h4(em("List Examined Cell Lines")),
             #              
             #              img(src="celllines.pdf", height = 130, width = 260, align = "middle"),
             #              br(), br(),
             #              
             #              # Not implemented
             #              # Input: Choose dataset ----
             #              selectInput("selectCancerType", label = em(h5("Select Cancer Type:")),
             #                          choices = c('',sort(as.character(TCGALabels[,1]))), 
             #                          selected = ''),
             #              # Button
             #              actionButton("goButton", "Go"),
             #              br(), br(),
             #              
             #              verbatimTextOutput("nText")
             #            ),
             #            # - - - - - - - - - - - - - - - - - - - - -
             #            
             #            
             #            
             #            
             #            
             #            
             #            # - - - - - - - - - - - - - - - - - - - - -
             #            sidebarPanel(
             #              h4(em("Explore Cell Line Genomic Features")),
             #              
             #              img(src="cellecta-stable-cell-line-engineering.jpg", height = 130, width = 280, align = "middle"),
             #              br(), br(),
             #              
             #              # Not implemented
             #              # Input: Choose dataset ----
             #              textInput("CLname", em(h5("Enter cell line of interest:")), placeholder="Cell Line name..."),
             #              # Button
             #              actionButton("go", "Go"),
             #              br(), br(),
             #              
             #              verbatimTextOutput("value2")
             #            )
             #          )
             # ),
             # - - - - - - - - - - - - - - - - - - - - -
             # tabPanel(title = "Help",
             #          titlePanel(em("Genomics Guided Selection of Cancer in vitro Models")),
             #          hr(), br(),
             #          
             #          p("This web application provides a suite of tools that enable researchers to make appropriate, informed choices about", em("in vitro"), "model inclusion/exclusion in retrospective analyses and future studies.", em("CELLector"), "assists in the selection of the best-representative cancer", em("in vitro"), "models based on their molecular similarity to their primary disease of origin and enables identification of disease subtypes currently lacking appropriate", em("in vitro"), "models."),
             #          br(),
             #          p("The", em("CELLector algorithm"), "is described fully in our manuscript", a("CELLector: Genomics Guided Selection of Cancer in vitro Models", href="https://github.com/HNaj/CELLector" , target="_blank"), icon("newspaper-o")),
             #          p("A tutorial demonstrating the full functionality of this web application is available", a("here.", href="https://github.com/HNaj/CELLector" , target="_blank"), icon("github-square")),
             #          br(),
             #          p("The", em("CELLector algorithm"), "and interactive visualisation tools are implemented in R and available as", a("an open-source R package from GitHub.", href="https://github.com/HNaj/CELLector" , target="_blank"), icon("github-square")),
             #          p("These pages were built using R Shiny, and the source code is freely available", a("here.", href="https://github.com/HNaj/CELLector" , target="_blank"), icon("github-square")),
             #          br(),
             #          # p("The app is released under a ........   License"),
             #          # br(),
             #          # p("Please cite:", a("CELLector R package", href = "https://github.com/HNaj/CELLector" , target="_blank"),icon("github-square"), "(only temp example)"),
             #          # p("Please report any issues using", a("the issue tracker.", href="https://github.com/HNaj/CELLector" , target="_blank"), icon("github-square")),
             #          p("Please report any issues to: email", icon("envelope")),
             #          # div("Please report any issues to: email", icon("envelope"), style = "color:red"),
             #          p(h5("CELLector v0.1"))
             # )
  )
)

