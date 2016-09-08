print("Starting structureXploR ...")
library(shiny) 
library(shinydashboard)
library(ape)
library(rjson)
library(jsonlite)
library(pvclust)
library(colorspace)
library(DT)
library(cluster)
library(plyr)
library(shinyjs)
library(devtools)
library(rCharts)
library(Biostrings)
library(BiocGenerics)
library("shinyBS")

header <- dashboardHeader(title = tagList(tags$em(strong(HTML('StructureXpl<i class="fa fa-compass"></i>R')))),titleWidth = 200, disable = FALSE)
#header <- dashboardHeader( title = span(tagList(icon("calendar"), "Example")))

#Sidebar elements for the search visualizations

sidebar <- dashboardSidebar(
  width = 200,
 # tags$em(h4(align="center",strong(HTML('StructureXpl<i class="fa fa-compass"></i>R')))),
  sidebarMenu(id="menu",
    menuItem(strong("Home"),tabName = "about", icon = icon("home")),
    menuItem(strong("Prepare"), tabName = "prepare", icon = icon("gears")),
    menuItem(strong("Explore"),tabName = "explore", icon = icon("compass")),


     div(id="clustConfig",style="display: none;",# better than using conditional panel ( panel is shown and after it is hidden)
       hr(),
       h5(align="center",icon("sliders"),strong("Clust. config.")),
       selectInput("hcExploreParam", label = tags$h6("Hierarchical clustering"), 
                   choices = list("UPGMA" = 0, "WARD" = 1,
                                  "SLINK" = 2,"CLINK" = 3), selected = 0)
        ,
       selectInput("setnbcluster", label = tags$h6("Nb. of clusters"), 
                   choices = list("Opt. clusters"=0,"2" = 2, "3" = 3, "4" = 4,"5"=5, "6" = 6, "7" = 7, "8" = 8,"9"=9,
                                  "10"=10,"11"=11,"12"=12,"13"=13,"14"=14,"15"=15,"16"=16,"17"=17,"18"=18,"19"=19,"20"=20,
                                  "21"=21,"22"=22,"23"=23,"24"=24,"25"=25,"26"=26,"27"=27,"28"=28,"29"=29,"30"=30)
                   , selected = 0)
       ,
       conditionalPanel(
         condition = "input.distChoiceParam == 1",
         selectInput("setnmotifs", label = tags$h6("Top rep. regions"), 
                     choices = list("1" = 1,"2" = 2, "3" = 3, "4" = 4, "5" = 5)
                     , selected = 1)
       )
     )
    , 
   tags$head(
              tags$link(rel = "stylesheet", type = "text/css", href = "custom.css"),
              tags$link(rel = "stylesheet", type = "text/css", href = "phylotree.css"),
              tags$link(rel = "stylesheet", type = "text/css", href = "introjs.css"),
              tags$link(rel = "stylesheet", type = "text/css", href = "label.css"),
              
              tags$script(src="fornac.js"),
              tags$script(src="d3.v3.min.js"),
              tags$script(src="underscore-min.js"),
              tags$script(src="phylotree.js"),
              tags$script(src="intro.js"),
              tags$script(src="custom.js"),
              
              tags$script(src="g_a.js")
              

      )
    
  ) # /sidebarMenu
) # /dashboardSidebar

#Body elements for the search visualizations.
body <- dashboardBody(

  useShinyjs(),
  fluidRow(column(12,
  

    tabItems(
    tabItem(tabName = "about",

            div(id="mainPanelStructurXploR",style='background-color:white;'
                ,br(),br(),br(),br(),
                
            fluidRow(
                  
              div(id="step11",
              tags$h1(strong(HTML('StructureXpl<i class="fa fa-compass"></i>R')),align="center",style = "font-size:50pt;"),br(),
              #p("StructureXpl",icon("gears","fa-stack-1x"),"R"),
              em(tags$h1("An interactive plateform to explore structural features of RNA secondary structures",align="center")),br()
              ),

            tags$h2(strong("Get started with your secondary structures"),align="center"),
            fluidRow(
              column(12,align="center",
                     h4(actionButton("buttonPrepare1","Prepare",icon("gears"),class="btn-primary",style='font-size:250%;color:white'),
                     actionButton("buttonExplore","Explore",icon("compass"),class="btn-success",style='font-size:250%;color:white;margin-left: 50px'))
              ))
              ),
            br(),br(),br(),br()
            ),

            div(id="step12",
            tags$h2(strong("Find, assess and explore clustering"),align="center"),
            br(),
            
                box(width = 12,status = "warning",
                    fluidRow( 
                              #column(4,align="center",h3(icon("share-alt","fa-rotate-180 fa-2x"),br(),"Find clusters of related structures")),
                              #column(4,align="center", h3(icon("check fa-2x"),br(),"Assess cluster quality")),
                              #column(4,align="center", h3(icon("sliders fa-2x"),br(),"Explore clustering configurations")
                              column(4,align="center",h3(HTML('<i class="fa fa-share-alt" style="color: #619CFF;"></i>'),"Find clusters of related structures" ),
                                     img(src='img/clusters.png',width="60%", height="60%", align = "center",style="box-shadow: 1px 1px 1px 1px  DarkGrey;	border-radius: 10px;")
                              ),
                              column(4,align="center",h3(HTML('<i class="fa fa-check" style="color: green;"></i>'),"Assess clustering quality" ),
                                     img(src='img/qual.png',width="50%", height="50%", align = "center",style="box-shadow: 1px 1px 1px 1px  grey;border-style: solid;	border-radius: 10px;"),
                                     img(src='img/qual2.png',width="50%", height="50%", align = "center",style="box-shadow: 1px 1px 1px 1px  grey;border-style: solid;	border-radius: 10px;"),
                                     img(src='img/qual3.png',width="50%", height="50%", align = "center",style="box-shadow: 1px 1px 1px 1px  grey;border-style: solid;	border-radius: 10px;")
                                     
                              ),
                              column(4,align="center",h3(HTML('<i class="fa fa-sliders" style="color: grey;"></i>'),"Explore clustering configurations" ),
                                     img(src='img/clustConf.png',width="40%", height="40%", align = "center",style="box-shadow: 1px 1px 1px 1px  grey;border-style: solid;	border-radius: 10px;")
                              )                      
                             ),
                    br(),br()
                )
                
                ,
                tags$h2(strong("Explore cluster features"),align="center"),
                br(),
                box(width = 12, status = "warning",
                    fluidRow( 
                      #column(4,align="center",h3(icon("sitemap fa-2x" ),br(),"Explore cluster hierarchy" )),
                      #column(4,align="center",h3(icon("dot-circle-o fa-2x" ),br(),"Identify representative structures" )),
                      # column(4,align="center",h3(icon("area-chart fa-2x" ),br(),"Visualize representative regions"), 
                      #        img(src='img/struct.png',width="50%", height="50%", align = "center"),
                      #        )
                      
                      column(4,align="center",h3(HTML('<i class="fa fa-sitemap" style="color: #F8766D;"></i>'),"Explore cluster hierarchy" ),
                             img(src='img/tree.png',width="100%", height="130%", align = "center",style="box-shadow: 1px 1px 1px 1px  grey;border-style: solid;border-radius: 10px;")
                             ),
                      column(4,
                             align="center",h3(icon("square"),"Identify representative structures" ),
                             img(src='img/rep.png',width="60%", height="60%", align = "center",style="box-shadow: 1px 1px 1px 1px  grey;border-style: solid;border-radius: 10px;")
                             ),
                      column(4, 
                             align="center",h3(HTML('<i class="fa fa-circle" style="color: #36C3E0;"></i>'),"Visualize representative regions"),
                             img(src='img/struct.png',width="50%", height="50%", align = "center",style="box-shadow: 1px 1px 1px 1px  grey;border-style: solid;border-radius: 10px;")
                              )                      
                  ),
                  br(),
                  em(h4(align="center","and many more features..." ))
                )
            
                                           
            ),
            hr(),  
            tags$h2(strong("About ", HTML('StructureXpl<i class="fa fa-compass"></i>R')),align="center"),
            tags$h3(align="center",
                    HTML('StructureXpl<i class="fa fa-compass"></i>R (v0.3)'),
                    "is maintained by Jean-Pierre Glouzon ",HTML('  <a href="https://github.com/jpsglouzon" target="_blank"><i class="fa fa-github"></i></a>'),br(),
                    "Suggestions or issues can be reported ",HTML('  <a href="https://github.com/jpsglouzon/structurexplor/issues" target="_blank">here</a>'),br(),br(),
                  
                    HTML(' <a href="http://info.usherbrooke.ca/Prospectus" target="_blank">Prospectus </a> ')," & ",
                    HTML(' <a href="http://jpperreaultlab.recherche.usherbrooke.ca/" target="_blank">Jean-Pierre Perreault Lab.</a> '),br(),
                    "GPL licence",br(),"2016"
                    )
          
      ),
    tabItem(tabName = "prepare",
            fluidRow(
            column(6,
              box(width=12,title = tagList(shiny::icon("gears"), "Prepare your secondary structures"), status='primary', solidHeader = TRUE,
                        'Prepare your data for exploration:', 
                        br(),
                        tags$ol(
                          tags$li(icon("reply-all","fa-rotate-180"),strong('Select your file of RNA scondary structures'),' : ' ,strong('click') ,' on',strong(' browse '),' button'),
                          tags$li(strong('Adjust '),' (or use default) ', icon("share-alt","fa-rotate-180"),strong(' structural comparison and clustering with bootstrap parameters') ,'(boxes below). '),
                          tags$li('Press',strong(icon("play-circle"),'to compute structural features'), ' button and explore','.')
                        ),
                        'Alternatively' ,strong(icon("pencil"),'run examples'),' or ', strong(icon("cube"),' load session data'),' and explore.',
                
                        hr(),
                        fluidRow(column(12,
                                        div(style="height: 60px;",
                                          fileInput("pathDbFile", label = tags$h5(icon("reply-all","fa-rotate-180"),strong("Enter your RNA secondary structures in "), HTML(' <a href="http://ultrastudio.org/en/Dot-Bracket_Notation" target="_blank">dot-bracket</a> ') ),
                                                    accept = c(
                                                      'text/plain',
                                                      '.db',
                                                      '.dot',
                                                      '.fasta'
                                                      )
                                                    )
                                        ),
                                        #br(),br(),
                                        #HTML('<code><textarea id="rnaDBinputTxt1" rows="13" cols="400" placeholder="Enter or paste RNA structures..."></textarea></code>'), 
                                        br(),
                                        div(id="validateStruct",""),
                                        br(),
                                        actionButton("go","Compute structural features",icon("play-circle")),
                                        br()
                                      )
                                ),
                        br(),
                        hr(),
                        fluidRow(column(6,
                                        strong(icon("pencil"),"Run examples"),": ",br(),
                                        conditionalPanel("input.bootstrap==0", 
                                        actionLink("ex_ss_linearRNA_pseudoknots", "229 struct. from 5S, tRNA and tmRNA (pseudoknots).",icon("arrow-circle-o-right")),br()
                                        ),
                                        actionLink("ex_ss_circularRNA", "21 struct. from viroids (circ. RNA).",icon("arrow-circle-o-right")),br(),
                                        conditionalPanel("input.bootstrap==0", 
                                        actionLink("ex_ss_linearRNA_g4", "88 struct. of 5S and struct. with G4.",icon("arrow-circle-o-right")),br()
                                        ),
                                        conditionalPanel("input.bootstrap==0", 
                                        actionLink("ex_ss_1000structures", "1186 struct. from 5S, HH, tRNA and 16s.",icon("arrow-circle-o-right")),br()
                                        ),
                                        
                                        actionLink("exInfo", "Source of structures",icon("info-circle")),br(),
                                        
                                        bsModal("modalExInfo", strong("Source of structures"), "exInfo", size = "medium",
                                                tags$ul(
                                                  tags$li( strong(HTML(' <a href="http://www.rnasoft.ca/strand/" target="_blank">RNASTRAND database</a> ')),
                                                    ' : 5S (5S ribosomal RNA), transfer RNA (tRNA), transfer messenger RNA (tmRNA), HH (Hammerhead Rizobyme) and 16S (16 Ribosomal RNA).'
                                                  ), 
                                                  tags$li(strong(HTML(' <a href="http://scottgroup.med.usherbrooke.ca/G4RNA/" target="_blank">G4RNA database</a> ')),
                                                    ' : Structures with g-quadruplexes (G4).'
                                                  ),
                                                  tags$li(strong("Giguère et al.(2014) Molecular Plant Pathology and Giguère et al. (2014) PLOS ONE"),
                                                    ' : Viroids (Circular RNA).'
                                                  ) 
                                                )
                                        )
                                      ),
                                 column(6,
                                        icon("cube"),strong("Session"),br(),
                                        "Save data of the current session or load data from previous session.",br(),  
                                        fluidRow(
                                                column(4,
                                                       downloadButton('saveSession', 'Save')
                                                ),
                                                column(8,
                                                         fileInput("loadSession", label = tags$h5(icon("upload"),"Load"),  accept = c('.R'))
                                                )

                                                )
                                        )
                        )
                  )
            ),
            column(6,
                  box(title = tagList(shiny::icon("share-alt","fa-rotate-180"), "Structural comparison parameters"), solidHeader = TRUE, status = "warning",collapsible = TRUE,
                      collapsed=TRUE,
                      fluidRow(
                        shinyjs::hidden(
                          column(6,
                                 radioButtons("distChoiceParam", label = "",
                                              choices = list("Super-n-motifs model" = 1, "RNAdistance" = 2,
                                                             "Base pair distance"=3), 
                                              selected = 1)
                          )
                        ),
                        column(12,
                               conditionalPanel(
                                 condition = "input.distChoiceParam == 1",
                                 fluidRow(column(12,h5(strong("The super-n-motifs model")))),
                                 #fluidRow(
                                   #column(6,
                                          selectInput("snm", label = tags$h6("Nb. of super-n-motifs"), 
                                                      choices = list("Auto." =0 ,"2" = 2, "3" = 3, "4" = 5, "6" = 6, "7" = 7, "8" = 8,"9"=9,
                                                                     "10"=10,"11"=11,"12"=12,"13"=13,"14"=14,"15"=15,"16"=16,"17"=17,"18"=18,"19"=19,"20"=20)
                                                      , selected = 0),
                                  #        ),
                                  # column(6,
                                          selectInput("max_n_motifs", label = tags$h6("Max. level of n-motifs"), 
                                                        choices = list("0-nmotifs" = 0, "1-nmotifs" = 1,
                                                                       "2-nmotifs" = 2), selected = 1)
                                  #        )
                                 #)
                                 ,
                                 actionLink("snmParam", "Super-n-motifs parameters",icon("info-circle")),br(),
                                 
                                 bsModal("modalsnmParam", strong("Super-n-motifs model parameters"), "snmParam", size = "medium",
                                         
                                           "The Super-n-motifs model performs efficient and accurate comparison of secondary structures. It is based on the ", 
                                           "super-n-motifs representation which considers secondary structure as combinations of motifs with their local context i.e. adjacent motifs.",
                                           "Motifs with their local context are called n-motifs. ex: S[HI] is a n-motif represented by a stem with hairpin and an internal loop." ,
                                           "S, S[HI] and S[HI]S[BI]S[BB] represent respectively a motif, a 1-motif and a 2-motif.",
                                           tags$ul(
                                             tags$li(strong("The number of super-n-motifs (Nb. of super-n-motifs)")," represents the maximum number of computed Super-n-motifs in the Super-n-motifs representation.",
                                                     "The automatic determination of the number of Super-n-motifs (", em('Auto.') ,") is based on the broken stick model."), 
                                             tags$li(strong("The maximum level of n-motifs (Max. level of n-motifs)")," represents the maximum level of n-motifs computed. For instance, if it set to 1-motifs, motifs and 1-motifs will be take into account.") 
                                           ),
                                          "Citation:..."
                                 )
                               )
                              
                        )
                      )
                  ),
                  box(title = tagList(shiny::icon("share-alt","fa-rotate-180"), "Clustering and boostrap parameters"), solidHeader = TRUE, status = "warning",collapsible = TRUE,
                      collapsed=TRUE,
                      fluidRow(
                        column(6,
                               div(id="step22",
                                   selectInput("HC", label = tags$h6("Hierarchical clustering"), 
                                               choices = list("UPGMA" = 0, "WARD" = 1,
                                                              "SLINK" = 2, "CLINK" = 3), selected = 0))),
                        column(6,
                               selectInput("maxClust", label = tags$h6("Max. nb. of clusters"), 
                                           choices = list("5" = 5, "10" = 10, "15" = 15, "20" = 25, "30" = 30)
                                           , selected = 10)
                        )
                      ),
                      selectInput("bootstrap", label = tags$h6("Nb. of bootstrap rep." ), 
                                  choices = list("0" = 0,"1000" = 1000), selected = 0)
                      ,
                      conditionalPanel(
                        condition = "input.bootstrap == 1000",
                        HTML('<p style="color:orange;font-weight:bold"><i class="fa fa-warning"></i> Bootstrap computation may take some times.</p>')
                        ),
                      
                      
                      actionLink("clusteringParam", "Clustering parameters",icon("info-circle")),
                      
                      bsModal("modalclusteringParam", strong("Clustering parameters"), "clusteringParam", size = "medium",
                              tags$ul(
                                tags$li(strong("Hierarchical clustering")," refers to the linkage function uses to compute the clustering. For further details on the silhouette coefficient, please see the function ", HTML(' <a href="https://stat.ethz.ch/R-manual/R-devel/library/stats/html/hclust.html" target="_blank">hclust</a>. ')
                                ), 
                                tags$li( strong("Maximum number of clusters (Max. nb. of clusters)")," defines the maximum number of clusters to be computed in order to find the best one."
                                ),
                                tags$li(strong("Number of bootstrap replications (Nb. bootstrap rep.)")," represent the number of times the bootstrap is performed. It is computed with the function ", HTML(' <a href="https://cran.r-project.org/web/packages/pvclust/index.html" target="_blank">pvclust</a>. '),
                                                'Bootstrap computation requires less than 50 structures as inputs and more than 5 computed super-n-motifs (Nb. of super-n-motifs parameter of the Super-n-motifs model). When the minimum number of computed super-n-motifs is not reach it is automatically set to 5.'
                                ) 
                              )
                          )
                      )
                  ,
                  box(width=12,
                  strong(icon("reply-all","fa-rotate-180"), HTML('StructureXpl<i class="fa fa-compass"></i>R')," supports structures of linear and circular RNA, pseudoknots and G-quadruplexes motifs"),
                  br(),
                  "RNA secondary structures must be in " ,HTML(' <a href="http://ultrastudio.org/en/Dot-Bracket_Notation" target="_blank">dot-bracket</a> ') ,' format',":",br(),
                  HTML('<code><textarea readonly id="rnaDBinputTxt2" rows="3" cols="200" placeholder="Enter RNA structures in dot-bracket...">>c_rnaID1
GAAAGGAAGGGGGAAAGGUUUGGAAAAGGGUUUGGGGUUGUUGGAAAAGGGGGGGGGGGGUUUUUUGG
.(((..{{.....AAA..)))..((((...aaa....}}.))..((((.++.++.++...))))))..</textarea></code>'),
                  h5(align="center",strong(textOutput("rna_prep1_Title"))),
                  htmlOutput("rna_prep1"),
                  br(),
                  icon("arrow-circle-o-right"),strong("Circular RNA")," require ",strong('c_')," at the beginning of the header : ",  strong('>c_rnaID1'),br(),
                  icon("arrow-circle-o-right"),strong("RNA Gquadruplexes (G4)")," are represented by the character", strong('+') , '. (See', HTML('<i class="fa fa-circle" style="color: #9E9E9E;"></i>'), ')', br(),
                  icon("arrow-circle-o-right"),strong("Pseudoknots")," are typically represented by", strong('<>'),', ',strong('{}'),', ',strong('[]'),
                  " or upper and lower case letter ", strong('Aa'),', ',strong('Bb'), ", etc.",'( See', HTML('<i class="fa fa-minus" style="color: red;"></i>'), ')'
                  )
                ,
                  box(width=12,
                    shiny::icon("warning"),strong(HTML('StructureXpl<i class="fa fa-compass"></i>R'),' require curated/validated structures'),
                    br(),
                    'Use MFE (Minimum Free Energy) structures or centroids of structure ensemble  with caution since they are error-prones and may lead to spurious structural patterns.'
                  )
            )
        
         ) 
    ),
    tabItem(tabName = "explore",
            fluidRow(
            conditionalPanel("input.go==0&&input.ex_ss_linearRNA_pseudoknots==0&&input.ex_ss_circularRNA==0&&input.ex_ss_linearRNA_g4==0&&input.ex_ss_1000structures==0&&output.sessionDataLoaded==0", 
                             tags$em(tags$h3("Please,", strong(" prepare your secondary structures "),"in the first place.",align = "center")),
                             fluidRow(column(12, align="center",
                             actionButton("buttonPrepare2","Prepare",icon("gears"),class="btn-primary",style='font-size:250%;color:white')))
            ),
            conditionalPanel("input.go>0||input.ex_ss_linearRNA_pseudoknots>0||input.ex_ss_circularRNA>0||input.ex_ss_linearRNA_g4>0||input.ex_ss_1000structures>0||output.sessionDataLoaded==1",
            tabBox( title = tagList(shiny::icon("compass"), "Explore"),
              id = "tabStructureXploR",width = 12,
              tabPanel(title = tagList(shiny::icon("check"), "Clustering quality"),value = "1",
                       fluidRow(column(12,
                                       
                                       actionLink("clusteringQualInfo", "Info. on clustering quality and silhouette coefficient",icon("info-circle")),br(),br(),
                                       bsModal("modalclusteringQualInfo", strong("Clustering quality and silhouette coefficient"), "clusteringQualInfo", size = "medium",
                                               'Global quality of clustering and cluster quality are associated with the computed silhouette coefficient( bounded [1,-1]) of the global clustering and each individual clusters.',
                                               'The quality is estimated as follows : ',
                                               tags$ul(
                                                 tags$li(strong("Very high quality")," when the computed coefficient is in the interval [1,0.7]. "
                                                 ), 
                                                 tags$li(strong("High quality")," when the computed coefficient is in  [0.7,0.5]. "
                                                 ),
                                                 tags$li(strong("Medium quality")," when the computed coefficient is in [0.5,0.3]. "
                                                 ),                                                  
                                                 tags$li(strong("Low quality")," when the computed coefficient is in [0.3,0]. "
                                                 ),                                                  
                                                 tags$li(strong("Very low quality")," when the computed coefficient is in [0,-1]. "
                                                 ),                                                  
                                                 tags$li(strong("Not applicable")," when the cluster size is too small, typically with cluster comprising less than 3 structures."
                                                 )                                                
                                               ),
                                               'For further details see the function ', HTML(' <a href="https://stat.ethz.ch/R-manual/R-devel/library/cluster/html/silhouette.html" target="_blank">silhouette</a>.'),br(),
                                               HTML('The ','<a href="https://ramnathv.github.io/rCharts/" target="_blank"> rCharts</a>',', ', '<a href="http://www.highcharts.com/" target="_blank"> Highcharts</a>' , ' and ' ,'<a href="https://rstudio.github.io/DT/" target="_blank"> DT</a>' )  ," libraries are used to visualize silhouette coefficients."
                                               
                                          ),
                                       div(id="step31",
                                          
                                          box( title = strong("Overview"), width = NULL, status = "primary", solidHeader = TRUE, 
                                           fluidRow(
                                             infoBoxOutput("nbClustersBox"),
                                             infoBoxOutput("qualClusteringBox"),
                                             infoBoxOutput("nbStructures")
                                            )
                                          )
                                       ),
                                       div(id="step32",
                                           
                                           fluidRow(
                                             column(3,
                                                    box(title = strong("Cluster quality"), width = NULL, solidHeader = TRUE, status = "primary",
                                                        uiOutput("cluterQualityBoxes")
                                                       )
                                                    )
                                                    ,
                                             column(9,
                                                    fluidRow(
                                                      column(6,
                                                             box(title = strong("Cluster silhouette coefficient"), width = NULL, solidHeader = TRUE, status = "warning",
                                                              uiOutput("hcontainer1") 
                                                             )
                                                      ),
                                                      column(6,
                                                             box(title = strong("Table of silhouette coefficients and sizes"), width = NULL, solidHeader = TRUE, status = "warning",
                                                               DT::dataTableOutput('clustInfo')
                                                             )
                                                      )
                                                    )
                                                    
                                            )
                                           )

                                         ) 
                                       )
                                )
                       ),
              tabPanel(title = tagList(icon("share-alt","fa-rotate-180"), "Cluster features"),value = "3", 
                       
                       actionLink("clusterFeat", "Info. on cluster Features",icon("info-circle")),br(),
                       bsModal("modalclusterFeat", strong("Clustering quality and silhouette coefficient"), "clusterFeat", size = "medium",
                               tags$ul(
                                 tags$li(
                                   strong("Structural variability of cluster (Struct. var.)")," yield information about how close to each other are the structures of a cluster . It is computed as the average distance between cluster members which is the intra-cluster distance. 
                                                          The rank represents the relative order of structural variabily of cluster from the one with the highest variability (1) to the lowest (n)."
                                 ), 
                                 tags$li(
                                   strong("Representative secondary structure of cluster (Representatives)")," is the structure that summarize structural features of a cluster. It is the closest structure to the cluster centroïde."
                                 ),
                                 tags$li(
                                   strong("Unusual secondary structures of cluster (Unusual structures)")," are structures with a rather different shape relative to the other members of the cluster. An unusual structure is considered as having a very low density relative to member of its cluster. Density around a structure is the mean distance between a structure and members of its cluster. "
                                 ),                                                  
                                 tags$li(
                                   strong("Top representatives regions (Top rep. regions) ")," are the n-motifs list that best describe structural features of cluster."
                                 )                                                
                               ),
                               HTML('The ','<a href="https://ramnathv.github.io/rCharts/" target="_blank"> rCharts</a>',', ', '<a href="http://www.highcharts.com/" target="_blank"> Highcharts</a>' , ' and ' ,'<a href="https://rstudio.github.io/DT/" target="_blank"> DT</a>' )  ," libraries are used to visualize cluster features."
                       ), br(),
                       
                       fluidRow(
                            column(4,
                                   box(title = strong("Cluster size"), width = NULL, solidHeader = TRUE, status = "primary",
                                     uiOutput("barChartClustSizeProp")
                                     #uiOutput("barChartClustSize")
                                     )
                                  ),
                            column(4,
                                   box(title = strong("Structural variability"), width = NULL, solidHeader = TRUE, status = "primary",
                                   uiOutput("barChartStructVar")
                                      )
                                   ),
                            column(4,
                                   box(title = strong("Length distribution"), width = NULL, solidHeader = TRUE, status = "primary",
                                   uiOutput("boxplotLengthDist")                                
                                   )
                                  )
                       ),
                       box(title = strong("Table of cluster features"), width = NULL, solidHeader = TRUE, status = "warning",

                           fluidRow(column(12,
                                           div(id="step51",

                                            DT::dataTableOutput('tbl')
                                           )
                                          )
                          )
                       )
              ),
              tabPanel(title = tagList(shiny::icon("line-chart"), "Features visualization"),value = "5", 
                    div(id="snmViz",
                      box(title = strong("Super-n-motifs representation of secondary structures"), width = NULL, solidHeader = TRUE, status = "primary",
                        fluidRow(
                         column(8,
                                
                                fluidRow(column(6,
                                                actionLink("infoSnmViz", "Info. on the super-n-motifs visualization",icon("info-circle")),br(),
                                                bsModal("modalinfoSnmViz", strong("Information on the super-n-motifs visualization"), "infoSnmViz", size = "medium",
                                                        HTML('The ','<a href="https://ramnathv.github.io/rCharts/" target="_blank"> rCharts</a> ' , ' and' ,'<a href="http://www.highcharts.com/" target="_blank"> Highcharts</a> ')  ," libraries are used to visualize the super-n-motifs representation."
                                                )
                                               ),
                                         column(6,align="right",actionButton("export-scatplotsnm-svg","SVG",icon("download"))
                                                )
                                         ),
                                 tags$em(tags$h5(strong(uiOutput("distHeader1Header2")),
                                                
                                                 shinyjs::hidden(actionLink("structDistInfo", "Info. on structural distance",icon("info-circle"))),br(),align="center"
                                                 )
                                         ),
                                bsModal("modalstructDistInfo", strong("Structural distance information"), "structDistInfo", size = "medium",
                                        strong("Structural distance (struct. a,struct. b) "),"is shown when two structures are selected.
                                                It represents the cosine dissimilarity, bounded [-1,1], between pairs of structures 
                                                in the space defined by the super-n-motifs retained."),

                                showOutput("scatplotsnm", "highcharts"),
                                
                                fluidRow( 
                                  
                                  column(4,""
                                        
                                      ),
                                  column(2,align="center",
                                         selectInput("snm_x", label = tags$strong(tags$h6("Super-n-motifs X")),
                                                       choices = list("1"=1,"2" = 2, "3" = 3, "4" = 5, "6" = 6, "7" = 7, "8" = 8,"9"=9,
                                                                      "10"=10,"11"=11,"12"=12,"13"=13,"14"=14,"15"=15)
                                                       , selected = 1)
                                         ),
                                column(2,align="center",
                                         selectInput("snm_y", label = tags$strong(tags$h6("Super-n-motifs Y")), 
                                                     choices = list("1"=1,"2" = 2, "3" = 3, "4" = 5, "6" = 6, "7" = 7, "8" = 8,"9"=9,
                                                                    "10"=10,"11"=11,"12"=12,"13"=13,"14"=14,"15"=15)
                                                     , selected = 2)
                                        ),
                                column(4,
                                       #strong("Control")," : ",
                                         tags$ul(
                                           tags$li(strong("Zoom"),": drag out a rectangle in the chart "
                                           ), 
                                           tags$li(strong("Show/hide elements of legend "),": click on legend element."
                                           )                                               
                                         )
                                       )
                               )
                               ),
                         column(4,
                                  fluidRow(
                                      column(12,
                                               
                                              fluidRow(column(6,h4(strong("Explained variability"))),
                                                       column(6,align="right",actionButton("export-varexplained-svg","SVG",icon("download")))
                                              ),
                                             showOutput("varExp", "highcharts"),
                                             tags$ul(
                                               tags$li(strong("Zoom"),": drag out a rectangle in the chart "
                                               ), 
                                               tags$li(strong("Show/hide elements of legend "),": click on legend element."
                                               )                                               
                                             )
                                             
                                      )
                                    )
                                  )
                               )
                          )
                        ),
                       fluidRow(
                         column(12,
                          box(title = strong("Visualize secondary structures"), width = NULL, solidHeader = TRUE, status = "warning",
                                    
                               actionLink("infoStructViz", "Info. on structure visualization",icon("info-circle")),br(),
                               bsModal("modalinfoStructViz", strong("Information on structure visualization"), "infoStructViz", size = "medium",
                                       HTML( '<a href="https://github.com/pkerpedjiev/fornac" target="_blank">Fornac</a> ')  ," is the tool used to visualize the secondary structures."
                               ),
                                fluidRow(column(3,checkboxInput("applyForce", label = "Move nucleotides (force layout)", value = FALSE)),
                                         column(3,checkboxInput("viznmotifs", label = "Visualize top rep. regions", value= FALSE)),
                                         column(3,checkboxInput("focusStructViz",label ="Focus on structure visualization",value=FALSE))
                                         ),
                             
                                fluidRow(
                                  column(6,
                                          div(id="headerss1_test",
                                          tags$em(tags$h5(strong(uiOutput("headerss1")),align="center"))),
                                          fluidRow(column(12,align="right",
                                                          downloadButton('exportdbStruct1', "Dot-bracket"),
                                                          actionButton("export-rna_ss1-svg","SVG",icon("download"))
                                                          )),
                                          fluidRow(column(12,htmlOutput("rna_ss1")))
                                          ),
                                  column(6,
                                          tags$em(tags$h5(strong(uiOutput("headerss2")),align="center")),
                                          fluidRow(column(12,align="right",
                                                          downloadButton('exportdbStruct2', "Dot-bracket"),
                                                          actionButton("export-rna_ss2-svg","SVG",icon("download"))
                                                          )),
                                          fluidRow(column(12,htmlOutput("rna_ss2")))
                                          )
                                  ),
                                br()
                                  ,
                               fluidRow(
                                 column(8, 
                                                  HTML('<i class="fa fa-circle" style="color: #36C3E0;"></i>'),
                                                  'Top representative regions', br(),
                                                  HTML('<i class="fa fa-circle" style="color: #9E9E9E;"></i>'),
                                                  'G-quadruplexes'
                                       ),
                                 column(4, 
                                        tags$ul(
                                          tags$li(strong("Zoom")," : scroll up/down or double click "
                                          ), 
                                          tags$li(strong("Reset/center")," : press 'c' key"
                                          )
                                        )
                                 )
                               )
                         )
                       )
                     )
                    ),
                    tabPanel(title = tagList(shiny::icon("sitemap"), "Cluster and structure hierarchy"),value = "2" ,
                        box(title = strong("Hierarchy"), width = NULL, solidHeader = TRUE, status = "primary",
                             fluidRow(column(12,
                                             actionLink("controlNodes", "Control over individual nodes",icon("info-circle")),br(),br(),
                                             bsModal("modalcontrolNodes", strong("Control over individual nodes"), "controlNodes", size = "medium",
                                                     tags$ul(
                                                       tags$li(
                                                         strong("Select all descendant, terminal, internal or incident branches of a node"),"by clicking on a node and selecting the corresponding option."
                                                       ),                                                  
                                                       tags$li(
                                                         strong("Reroot on a node"),"by clicking on a node and selecting 'Reroot on this node'."
                                                       )                                                
                                                     ),
                                                     HTML('The ','<a href="https://github.com/veg/phylotree.js" target="_blank"> phylotree.js</a>')," is the library used to visualize the hierarchy."
                                             )
                                             ,
                                             conditionalPanel(condition="input.bootstrap>0",
                                                              strong('Hierarchy with appromixately unbiased bootstrap values'),br(),br(),
                                                              
                                                              shinyjs::hidden(
                                                              selectInput("AU_bootstrap_values", label = tags$h6("Show bootstrap values"), 
                                                                          choices = list("Appromixately unbiased bootstrap values"=1)
                                                                          , selected = 1))   
                                                              
                                                              ),
                                             #treewidgetOutput("dendSS",width='100%',height='100%')
                                             HTML('
                                                  <!-- Brand and toggle get grouped for better mobile display -->
                                                  <div class="row">
                                                  
                                                  <div class="col-md-3">
                                                  
                                                  
                                                  <button type="button" class="btn btn-default btn-sm" data-direction = "vertical" data-amount = "1" title = "Expand vertical spacing">
                                                  <i class="fa fa-arrows-v" ></i>
                                                  </button>
                                                  <button type="button" class="btn btn-default btn-sm" data-direction = "vertical" data-amount = "-1" title = "Compress vertical spacing">
                                                  <i class="fa  fa-compress fa-rotate-135" ></i>
                                                  </button>
                                                  
                                                  <button type="button" class="btn btn-default btn-sm" id="sort_ascending" title="Sort deepest clades to the bototm">
                                                  <i class="fa fa-sort-amount-asc"></i>
                                                  </button>
                                                  <button type="button" class="btn btn-default btn-sm" id="sort_descending" title="Sort deepsest clades to the top">
                                                  <i class="fa fa-sort-amount-desc"></i>
                                                  </button>
                                                  <button type="button" class="btn btn-default btn-sm" id="sort_original" title="Restore original order">
                                                  <i class="fa fa-refresh"></i>
                                                  </button>
                                                  <div class="modal" id = \'newick_export_modal\'>
                                                  <div class="modal-dialog">
                                                  <div class="modal-content">
                                                  <div class="modal-body" id = \'newick_body\'>
                                                  <textarea id = \'nwk_export_spec\' autofocus = true placeholder = "" style = \'width: 100%; height: 100%\' rows = 20></textarea>
                                                  </div>
                                                  <div class="modal-footer">
                                                  <button type="button" class="btn btn-default" data-dismiss="modal">Close</button>
                                                  </div>
                                                  
                                                  </div><!-- /.modal-content -->
                                                  </div><!-- /.modal-dialog -->
                                                  </div><!-- /.modal -->
                                                  
                                                  
                                                  </div>
                                                  <div class="col-md-3">
                                                  <input type="text" id="branch_filter" class="form-control" placeholder="&#xf0b0; Filter branches on">
                                                  <label class="pull-right">Selected <span class="badge" id="selected_branch_counter">0</span> and filtered <span class="badge" id="selected_filtered_counter">0</span> branches</label>
                                                  
                                                  </div>
                                                  <div class="col-md-3 dropdown">
                                                  <button type="button" class="btn btn-default dropdown-toggle" data-toggle="dropdown">Selection 
                                                  <span class="caret"></span></button>
                                                  <ul class="dropdown-menu">
                                                  <li><a href="#" id="select_all">Select all</a></li>
                                                  <li><a href="#" id="select_all_internal">Select all internal nodes</a></li>
                                                  <li><a href="#" id="select_all_leaves">Select all leaf nodes</a></li>
                                                  <li><a href="#" id="clear_internal">Clear all internal nodes</a></li>
                                                  <li><a href="#" id="clear_leaves">Clear all leaves</a></li>
                                                  <li><a href="#" id="select_none">Clear selection</a></li>
                                                  </ul>
                                                  
                                                  </div>
                                                  <div class="col-md-3">
                                                  
                                                  <a href="#" data-toggle="modal" class="btn btn-default btn-sm" data-target="#newick_export_modal"> 
                                                  <i class="fa fa-download"> Newick</i>
                                                  </a>
                                                  
                                                  <button id="export-phylo-svg" type="button" class="btn btn-default btn-sm" data-amount="1" title="Expand vertical spacing">
                                                  <i class="fa fa-download"> SVG</i>
                                                  </button>
                                                  </div>
                                                  </div>
                                                  
                                                  ')
                                        ,
                                        fluidRow(column(12,align="center",htmlOutput("legendHierarchy")))
                                        ,
                                        div(id="step42",fluidRow(column(12,align="center",htmlOutput("dendSS2",class="tree-widget")))
                                        )
                                       )
                                      )
                                )
                             ),
                    tabPanel(title = tagList(shiny::icon("list"), "List of structures"),value = "4",
                         box(title = strong("Table of structure features"), width = NULL, solidHeader = TRUE, status = "primary",
                             fluidRow(column(12,
                                             fluidRow(column(6,align="left",
                                                             
                                                             actionLink("structFeat", "Info. on structure features",icon("info-circle")),br(),
                                                             bsModal("modalstructFeat", strong("Information on structure features"), "structFeat", size = "medium",
                                                                     tags$ul(
                                                                       tags$li(
                                                                         strong("Representative secondary structure of a cluster (Representatives)")," is the structure that summarize structural features of a cluster. It is structure which is the closest to the centroïde of the cluster."
                                                                       ),
                                                                       tags$li(
                                                                         strong("Length")," is the total number of nucleotides of a structure."
                                                                       ),
                                                                       tags$li(
                                                                         strong("Unusual secondary structures of a cluster (Unusual structures)")," are structures with a rather different shape relative to the other members of a cluster. An unusual structure is considered as having a very low density relative to member of its cluster. Density around a structure is the mean distance between a structure and members of its cluster."
                                                                       )                                                
                                                                     ),
                                                                     HTML('The ','<a href="https://rstudio.github.io/DT/" target="_blank">DT</a>')," library is used to generate the table of structure features."
                                                                          
                                                             )
                                             ),
                                             column(6,align="right",downloadButton('exportDB', "Filtered struct. in dot-bracket"))
                                             ),
                                             
                                             div(id="step61",
                                                 DT::dataTableOutput('tbOfStruct')
                                             )
                                        )
                             )
                          )
                    )
                  )
              )
            
          )
   # /tabItems
      )

    )
  )
  ),

#require to export phylotree in svg
HTML('                                               
      <img id="hyphy-chart-image"></img>
      <canvas id="hyphy-chart-canvas"></canvas>
     ')
) # /dashboardBody

dashboardPage(header, sidebar, body, skin = "blue")


