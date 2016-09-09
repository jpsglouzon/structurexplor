source("www/Functions/mainFunctions.R")
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


shinyServer(function(input, output,session) {
  
  
  options(shiny.maxRequestSize=30*1024^2)
  options(shiny.minified=TRUE)
  options(shiny.usecairo=FALSE)
  
  pathfile <- reactiveValues(data = NULL)
  data_res <- reactiveValues(data = NULL)
  optimalNbCluster <- reactiveValues(data = NULL)
  data_res_updatesetnmotifs <- reactiveValues(data = NULL)
  guidedTourActivated<- reactiveValues(data = FALSE)
  
  parsedRNA<-reactiveValues(data = NULL)
  parsingIsSuccessful<-reactiveValues(data = FALSE)
  
  session$onSessionEnded(function() {
    
    stopApp()
  })
  
  observe({
    #To activate/deactive the 'Compute structural patterns' button

    if(input$menu=="explore"&&!is.null(data_res$data)){
      shinyjs::show("clustConfig")
    }else 
      { shinyjs::hide("clustConfig") }
    
    if (!parsingIsSuccessful$data) {
      shinyjs::disable("go")
      a=2
    } else { shinyjs::enable("go")}
    
    if (is.null(input$click$header1) || nchar(input$click$header1)==0) {
      shinyjs::disable("exportdbStruct1")
      shinyjs::disable("export-rna_ss1-svg")
      
    } else { shinyjs::enable("exportdbStruct1")
             shinyjs::enable("export-rna_ss1-svg") }
    
     if (is.null(input$click$header2) || (nchar(input$click$header2)==0)) {
      shinyjs::disable("exportdbStruct2")
      shinyjs::disable("export-rna_ss2-svg")
     } else { shinyjs::enable("exportdbStruct2")
              shinyjs::enable("export-rna_ss2-svg")
     }
    
    if(input$focusStructViz==TRUE){
      shinyjs::hide("snmViz",anim=TRUE)
      
      #GA
      shinyjs::runjs(HTML("ga('send', 'event', 'Visualize structure', 'Focus on structure visualization','", input$viznmotifs ,"');")); 
      
    } else { shinyjs::show("snmViz",anim=TRUE)}
    
    if(is.null(data_res$data)){
      shinyjs::disable("saveSession")
    }else {shinyjs::enable("saveSession")}
    
    
  })
  
  observeEvent(input$bootstrap, {
    
    #nb struct maximum for bootstrap 

    if(!is.null(parsedRNA$data$headers))
    {
    if (length(parsedRNA$data$headers)>25 & input$bootstrap>0)
    {shinyjs::html("validateStruct", '<br><p style="color:orange;font-weight:bold"><i class="fa fa-warning"></i>Less than 25 structures is required. Please reduce the number of structures for bootstrap.</p>')
      parsingIsSuccessful$data=F}
    else
      { shinyjs::html("validateStruct", '<br><p style="color:green;font-weight:bold"><i class="fa fa-check"></i>Structures successfully parsed.</p>')
        parsingIsSuccessful$data=T}
    }
    
  })

  observeEvent(input$pathDbFile$datapath, {

    shinyjs::html("validateStruct", '<br><p style="font-weight:bold"><i class="fa fa-cog fa-spin"></i>Parsing structures. Please wait...</p>') 
    parsedRNA$data=parseDbFile(input$pathDbFile$datapath)
    
    #nb struct minimum
    if (length(parsedRNA$data$headers)<5)
    {shinyjs::html("validateStruct", '<br><p style="color:orange;font-weight:bold"><i class="fa fa-warning"></i>More than 5 structures is required. Please add more structures.</p>') 
      parsingIsSuccessful$data=F}
    #nb struct minimum
    if (length(parsedRNA$data$headers)>250)
    {shinyjs::html("validateStruct", '<br><p style="color:orange;font-weight:bold"><i class="fa fa-warning"></i>Less than 250 structures is required. Please reduce the number of structures.</p>') 
      parsingIsSuccessful$data=F}
    
    #structure with id identique
    else if(length(parsedRNA$data$headers)!=length(unique(parsedRNA$data$headers)))
    {shinyjs::html("validateStruct", '<br><p style="color:orange;font-weight:bold"><i class="fa fa-warning"></i>Headers of structures must be unique. Please check the headers of your structures.</p>') 
      parsingIsSuccessful$data=F}
    #check presence énergie
    else if(length(grep("[0123456789]",parsedRNA$data$structures))>0)
    {shinyjs::html("validateStruct", '<br><p style="color:orange;font-weight:bold"><i class="fa fa-warning"></i>Numerical characters are not supported. Please check your structures.</p>')  
      parsingIsSuccessful$data=F}
    else
    { shinyjs::html("validateStruct", '<br><p style="color:green;font-weight:bold"><i class="fa fa-check"></i>Structures successfully parsed.</p>')
      parsingIsSuccessful$data=T}
    
  })
  
  #Swith to tab 2.explore after clicking on run pipeline
  observeEvent(input$go, {
   
    pathfile$data=input$pathDbFile$datapath
    
    shinyjs::reset("hcExploreParam")
    
    setnbcluster_options <- list()
    setnbcluster_options[["Opt. nb. clusters"]] <- 0
    nbHeaders=length(parsedRNA$data$headers)
    if(nbHeaders<=30) #30 is the max of maxCluster and setnbCluster
    {
    setnbcluster_options<-c(setnbcluster_options,as.list(2:nbHeaders-1))
    updateSelectInput(session, "setnbcluster", choices = setnbcluster_options, selected = 0)
    }
    else{
      setnbcluster_options<-c(setnbcluster_options,as.list(2:30))
      updateSelectInput(session, "setnbcluster", choices = setnbcluster_options, selected = 0)
    }
  
    shinyjs::reset("setnmotifs")
    data_res_updatesetnmotifs$data=NULL  
    
    data_res$data=computeStructuralPatterns(pathfile$data,parsedRNA$data,input$distChoiceParam,input$snm,input$max_n_motifs,input$rnad,input$setnmotifs,input$maxClust,input$bootstrap,input$HC,input$setnbcluster)
    
    #si erreur repointer sur les données
    optimalNbCluster$data=data_res$data[[3]]$bestNbClusters
    
    updateTabItems(session, "menu","explore")
   
    shinyjs::show("clustConfig")


  })
  
  observeEvent(input$ex_ss_linearRNA_pseudoknots, {
    
    pathfile$data="www/Data/secStruc_linear_RNA_pseudoknots.db"
    shinyjs::reset("hcExploreParam")
    
    setnbcluster_options <- list()
    setnbcluster_options[["Opt. nb. clusters"]] <- 0
    setnbcluster_options<-c(setnbcluster_options,as.list(2:30))
    updateSelectInput(session, "setnbcluster", choices = setnbcluster_options, selected = 0)
    
    shinyjs::reset("setnmotifs")
    
    data_res_updatesetnmotifs$data=NULL  
    
    data_res$data=computeStructuralPatterns(pathfile$data,parseDbFile(pathfile$data),input$distChoiceParam,input$snm,input$max_n_motifs,input$rnad,input$setnmotifs,input$maxClust,input$bootstrap,input$HC,input$setnbcluster)
    optimalNbCluster$data=data_res$data[[3]]$bestNbClusters
      
    updateTabItems(session, "menu","explore")
    shinyjs::show("clustConfig")
    
    #update snm parameter in feature visualization
    setsnm_x_options <- list()
    ncol(data_res$data[[2]]$SuperMotif)
    setsnm_x_options<-c(setsnm_x_options,as.list(1:ncol(as.data.frame(data_res$data[[2]]$SuperMotif))))
    updateSelectInput(session, "snm_x", choices = setsnm_x_options, selected = 1)
    
    setsnm_y_options <- list()
    ncol(data_res$data[[2]]$SuperMotif)
    setsnm_y_options<-c(setsnm_y_options,as.list(1:ncol(as.data.frame(data_res$data[[2]]$SuperMotif))))
    

    
    updateSelectInput(session, "snm_y", choices = setsnm_y_options, selected = 2)
    
  }) 
  
  observeEvent(input$ex_ss_circularRNA, {
    
    
    pathfile$data="www/Data/secStruc_circular_RNA.db"
    shinyjs::reset("hcExploreParam")
    
    setnbcluster_options <- list()
    setnbcluster_options[["Opt. nb. clusters"]] <- 0
    setnbcluster_options<-c(setnbcluster_options,as.list(2:20))
    updateSelectInput(session, "setnbcluster", choices = setnbcluster_options, selected = 0)

    shinyjs::reset("setnmotifs")
    data_res_updatesetnmotifs$data=NULL  
    
  
    data_res$data=computeStructuralPatterns(pathfile$data,parseDbFile(pathfile$data),input$distChoiceParam,input$snm,input$max_n_motifs,input$rnad,input$setnmotifs,input$maxClust,input$bootstrap,input$HC,input$setnbcluster)
    optimalNbCluster$data=data_res$data[[3]]$bestNbClusters
    
    
    updateTabItems(session, "menu","explore")
    shinyjs::show("clustConfig")
    
    #update snm parameter in feature visualization
    setsnm_x_options <- list()
    ncol(data_res$data[[2]]$SuperMotif)
    setsnm_x_options<-c(setsnm_x_options,as.list(1:ncol(as.data.frame(data_res$data[[2]]$SuperMotif))))
    updateSelectInput(session, "snm_x", choices = setsnm_x_options, selected = 1)
    
    setsnm_y_options <- list()
    ncol(data_res$data[[2]]$SuperMotif)
    setsnm_y_options<-c(setsnm_y_options,as.list(1:ncol(as.data.frame(data_res$data[[2]]$SuperMotif))))
    
    
    updateSelectInput(session, "snm_y", choices = setsnm_y_options, selected = 2)
    

    
  
  }) 
  
  observeEvent(input$ex_ss_linearRNA_g4, {
    
    
    pathfile$data="www/Data/secStruc_linear_RNA_g4.db"
    shinyjs::reset("hcExploreParam")
    
    setnbcluster_options <- list()
    setnbcluster_options[["Opt. nb. clusters"]] <- 0
    setnbcluster_options<-c(setnbcluster_options,as.list(2:30))
    updateSelectInput(session, "setnbcluster", choices = setnbcluster_options, selected = 0)
    shinyjs::reset("setnmotifs")
    data_res_updatesetnmotifs$data=NULL  
    
    data_res$data=computeStructuralPatterns(pathfile$data,parseDbFile(pathfile$data),input$distChoiceParam,input$snm,input$max_n_motifs,input$rnad,input$setnmotifs,input$maxClust,input$bootstrap,input$HC,input$setnbcluster)
    optimalNbCluster$data=data_res$data[[3]]$bestNbClusters
    
    updateTabItems(session, "menu","explore")
    shinyjs::show("clustConfig")
    
    #update snm parameter in feature visualization
    setsnm_x_options <- list()
    ncol(data_res$data[[2]]$SuperMotif)
    setsnm_x_options<-c(setsnm_x_options,as.list(1:ncol(as.data.frame(data_res$data[[2]]$SuperMotif))))
    updateSelectInput(session, "snm_x", choices = setsnm_x_options, selected = 1)
    
    setsnm_y_options <- list()
    ncol(data_res$data[[2]]$SuperMotif)
    setsnm_y_options<-c(setsnm_y_options,as.list(1:ncol(as.data.frame(data_res$data[[2]]$SuperMotif))))
    updateSelectInput(session, "snm_y", choices = setsnm_y_options, selected = 2)
    
  }) 
  
  
  observeEvent(input$ex_ss_1000structures, {
    
  
    pathfile$data="www/Data/1000structures.db"
    shinyjs::reset("hcExploreParam")
    
    setnbcluster_options <- list()
    setnbcluster_options[["Opt. nb. clusters"]] <- 0
    setnbcluster_options<-c(setnbcluster_options,as.list(2:30))
    updateSelectInput(session, "setnbcluster", choices = setnbcluster_options, selected = 0)
    shinyjs::reset("setnmotifs")
    data_res_updatesetnmotifs$data=NULL  
    
    data_res$data=computeStructuralPatterns(pathfile$data,parseDbFile(pathfile$data),input$distChoiceParam,input$snm,input$max_n_motifs,input$rnad,input$setnmotifs,input$maxClust,input$bootstrap,input$HC,input$setnbcluster)
    optimalNbCluster$data=data_res$data[[3]]$bestNbClusters

   
        
    updateTabItems(session, "menu","explore")
    shinyjs::show("clustConfig")
    
    #update snm parameter in feature visualization
    setsnm_x_options <- list()
    ncol(data_res$data[[2]]$SuperMotif)
    setsnm_x_options<-c(setsnm_x_options,as.list(1:ncol(as.data.frame(data_res$data[[2]]$SuperMotif))))
    updateSelectInput(session, "snm_x", choices = setsnm_x_options, selected = 1)
    
    setsnm_y_options <- list()
    ncol(data_res$data[[2]]$SuperMotif)
    setsnm_y_options<-c(setsnm_y_options,as.list(1:ncol(as.data.frame(data_res$data[[2]]$SuperMotif))))
    updateSelectInput(session, "snm_y", choices = setsnm_y_options, selected = 2)
    
  
  })   
  
  
  #Session
  observeEvent(input$loadSession$datapath, {
    shinyjs::reset("hcExploreParam")
    shinyjs::reset("setnmotifs")
    data_res_updatesetnmotifs$data=NULL
    
    data_restemp=readRDS(input$loadSession$datapath)
    data_res$data=data_restemp$data
    
    setnbcluster_options <- list()
    setnbcluster_options[["Opt. nb. clusters"]] <- 0
    nbHeaders=length( data_res$data[[1]]$headers)
    if(nbHeaders<=30) #30 is the max of maxCluster and setnbCluster
    {
      setnbcluster_options<-c(setnbcluster_options,as.list(2:nbHeaders-1))
      updateSelectInput(session, "setnbcluster", choices = setnbcluster_options, selected = 0)
    }
    else{
      setnbcluster_options<-c(setnbcluster_options,as.list(2:30))
      updateSelectInput(session, "setnbcluster", choices = setnbcluster_options, selected = 0)
    }
  
   
    updateTabItems(session, "menu","explore")
    shinyjs::show("clustConfig")
    

  })
  
  output$sessionDataLoaded<- reactive({
    if (is.null(input$loadSession$datapath))
    {return(0) }
    else{return(1)}
  })
  outputOptions(output, 'sessionDataLoaded', suspendWhenHidden=FALSE)
  

  output$saveSession <- downloadHandler(
    
    
    filename = function() {
      paste('sessionData-', Sys.Date(), '.R', sep='')
    },
    content = function(fileConn) {
      saveRDS(data_res, file = fileConn)
    }
    
    
  )
  
  ##############
  
  
  #Modify nb. of cluster (Exploration parameters)
  observeEvent(input$hcExploreParam, {
    if(!is.null(data_res$data))
    {
      shinyjs::reset("setnmotifs")
      data_res_updatesetnmotifs$data=NULL  
      
      data_res$data=updateStructuralPatterns_with_setnbcluster(data_res$data[[1]],data_res$data[[2]],input$distChoiceParam,input$setnmotifs,input$maxClust,input$bootstrap,input$hcExploreParam,input$setnbcluster)
    
      #GA
      shinyjs::runjs(HTML("ga('send', 'event', 'Clust. Config', 'Hierarchical clustering','", input$hcExploreParam ,"');"));
      
    }
  })
  
  #Modify nb. of cluster (Exploration parameters)
  observeEvent(input$setnbcluster, {
      if(!is.null(data_res$data))
    {
      #shinyjs::reset("setnmotifs")
        data_res_updatesetnmotifs$data=NULL  
      data_res$data=updateStructuralPatterns_with_setnbcluster(data_res$data[[1]],data_res$data[[2]],input$distChoiceParam,input$setnmotifs,input$maxClust,input$bootstrap,input$hcExploreParam,input$setnbcluster)
      
      #GA
      shinyjs::runjs(HTML("ga('send', 'event', 'Clust. Config', 'Nb. of clusters','", input$setnbcluster ,"');"));
    }
    
  })
  
  #setnmotifs  (Exploration parameters)
  observeEvent(input$setnmotifs, {
      if(!is.null(data_res$data))
      {
        data_res_updatesetnmotifs$data=data_res$data  
        data_res_updatesetnmotifs$data=updateStructuralPatterns_with_setnmotifs(data_res$data,data_res$data[[1]],data_res$data[[2]],input$distChoiceParam,input$setnmotifs,input$maxClust,input$bootstrap,input$hcExploreParam,input$setnbcluster)
        data_res$data=data_res_updatesetnmotifs$data 
        
        #GA
        shinyjs::runjs(HTML("ga('send', 'event', 'Clust. Config', 'Top rep. regions','", input$setnmotifs ,"');")); 
        
        }
  })
  
  output$rna_prep1 <- renderUI({
    
    #print("") #allow to keep the graph visible
    rnaDBinputTxt2Parsed=strsplit(input$rnaDBinputTxt2, "\n")
    head=rnaDBinputTxt2Parsed[[1]][1]
    head=gsub('>','',head)
    seq=rnaDBinputTxt2Parsed[[1]][2]
    struct=rnaDBinputTxt2Parsed[[1]][3]
    
    #G4
    g4PoswithColors=""
    g4PosTemp=gregexpr("\\+", struct,ignore.case=FALSE)[[1]]
    if(g4PosTemp[1]!=-1){
      g4Pos=g4PosTemp[1:length(g4PosTemp)]
      g4PoswithColors=paste(g4Pos,":rgb(158,158,158)",sep="")
      g4PoswithColors=paste(g4PoswithColors,collapse=" ")
      struct=gsub("\\+", "\\.", struct)
    }
    
    #circRNA
    if(substr(head,0,2)=="c_"){struct=paste(struct,"*",sep="")} 

    atemp1=paste("
                  var container = new FornaContainer(\"#rna_prep1\",
                  {'applyForce': false, 'allowPanningAndZooming': true, 'initialSize':[400,220]});
                  var options = {
                  'structure': '",struct,"',
                  'sequence': '",seq,"'
                  };container.addRNA(options.structure, options);
                   var cs = new ColorScheme('",g4PoswithColors,"');
                   container.addCustomColors(cs.colorsJson);
                   container.changeColorScheme('custom'); 
                   ",
                sep="")

      ss=tags$script(HTML(atemp1))  
    
    return(ss)
  })
  
  #Title of the seconary structure in page prepare
  output$rna_prep1_Title<-renderText({
    
    rnaDBinputTxt2Parsed=strsplit(input$rnaDBinputTxt2, "\n")
    head=rnaDBinputTxt2Parsed[[1]][1]
    head=gsub('>','',head)
    
    return(paste("The secondary structure of ",head,sep=""))
  })
  
  
#Cluster quality tab
  
  
  output$nbStructures<- renderInfoBox({
    if (is.null(data_res$data)) return()
    
    dataPatterns=data_res$data
    infoBox(strong("Nb. of structures"),h3(strong(dim(dataPatterns[[1]])[1])),
       icon = icon("th"),  color = "teal",fill=TRUE
    )
  })
  
  output$nbClustersBox <- renderInfoBox({
    if (is.null(data_res$data)) return()
    dataPatterns=data_res$data
    
    #GA
    shinyjs::runjs(HTML("ga('send', 'event', 'Explore', 'Clustering Quality');"));  
    
    if (is.null(optimalNbCluster$data)){
      optimalNbCluster$data=dataPatterns[[3]]$bestNbClusters
    }
    
    if (optimalNbCluster$data==dataPatterns[[3]]$bestNbClusters)
    {
      infoBox(strong("Optimal nb. of clusters"),h3(strong(dataPatterns[[3]]$bestNbClusters)),
         icon = icon("th-large"),
        color = "aqua",fill=TRUE
      )
    }
    else
    {
      infoBox(strong("Nb. of clusters"),h3(strong(dataPatterns[[3]]$bestNbClusters)),
         icon = icon("th-large"),
        color = "aqua",fill=TRUE
      )     
    }

  })
  
  output$cluterQualityBoxes<- renderUI({
    dataPatterns=data_res$data
    if (is.null(data_res$data)) return()
    
    plot_output_list<-lapply(1:dataPatterns[[3]]$bestNbClusters, function(i) {
    output[[paste0('b', i)]] <- renderInfoBox({

      if(dataPatterns[[3]]$clustSize[i]<=2)
      {
        infoBox(strong(paste("Cluster",i)),'Not-app.',
              icon = icon("times"),
              color = 'black' ,fill=FALSE) 
      }
      else if (dataPatterns[[3]]$listSilhouetteCoefPerClusters[i]>0.7)
      {
        infoBox(strong("Cluster ",i),'Very high',
                icon = icon("arrow-up"),
                color = 'olive' ,fill=FALSE) 
        }
      else if (data_res$data[[3]]$listSilhouetteCoefPerClusters[i]>0.5)
      {
        infoBox(strong(paste("Cluster",i)),strong('High'),
                icon = icon("arrow-up","fa-rotate-45"),
                color = 'green' ,fill=FALSE) 
      }
      else if (data_res$data[[3]]$listSilhouetteCoefPerClusters[i]>0.3)
      {
      infoBox(strong(paste("Cluster",i)),strong('Medium'),
              icon = icon("arrow-right"),
              color = 'yellow' ,fill=FALSE) 
      }
      else if (data_res$data[[3]]$listSilhouetteCoefPerClusters[i]>0)
      {
      infoBox(strong(paste("Cluster",i)),strong('Low'),
              icon = icon("arrow-right","fa-rotate-45"),
              color = 'red' ,fill=FALSE) 
      }
      else if (data_res$data[[3]]$listSilhouetteCoefPerClusters[i]<=0)
      {
        infoBox(strong(paste("Cluster",i)),strong('Very low'),
                icon = icon("arrow-down"),
                color = 'black' ,fill=FALSE) 
      }
      
      })
    })
    plot_output_list
  })
  
  output$qualClusteringBox <- renderInfoBox({
    if (is.null(data_res$data)) return()
    dataPatterns=data_res$data
    
    globalClustsizeIsLessThan2=TRUE
     for (i in 1:length(dataPatterns[[3]]$clustSize)) 
     {
       if (dataPatterns[[3]]$clustSize[i]>2)
       {globalClustsizeIsLessThan2=FALSE}
     }
    
    if (globalClustsizeIsLessThan2==TRUE)
    {
      infoBox(strong("Global quality of clustering"),h3(strong("Not-applicable")),fill=FALSE,
              icon = icon("times"),
              color = "black"
      )
    }
    else if (dataPatterns[[3]]$avgSilhouetteCoef>0.7)
    {
      infoBox(strong("Global quality of clustering"),h3(strong("Very high")),fill=TRUE,
              icon = icon("arrow-up"),
              color = "olive"
      )
    }
    else if(dataPatterns[[3]]$avgSilhouetteCoef>0.5)
    {
    infoBox(strong("Global quality of clustering"),h3(strong("High")),fill=TRUE,
            icon = icon("arrow-up","fa-rotate-45"),
            color = "green"
    )
    
    }
    else if(dataPatterns[[3]]$avgSilhouetteCoef>0.3)
    {
    infoBox(strong("Global quality of clustering"),h3(strong("Medium")),fill=TRUE,
            icon = icon("arrow-right"),
            color = "yellow"
    )
    
    }
    else if(dataPatterns[[3]]$avgSilhouetteCoef>0)
    {
    infoBox(strong("Global quality of clustering"),h3(strong("Low")),fill=TRUE,
            icon = icon("arrow-right","fa-rotate-45"),
            color = "red"
    )
    
    }
    else if(dataPatterns[[3]]$avgSilhouetteCoef<=0)
    {
    infoBox(strong("Global quality of clustering"),h3(strong("Very low")),fill=TRUE,
            icon = icon("arrow-down"),
            color = "black"
    )
    }
  })
   
  output$hcontainer1 <- renderChart2({
    if (is.null(data_res$data)) return()
    dataPatterns=data_res$data
    
    datatemp1=dataPatterns[[3]]$avgSilhouetteCoef
    names(datatemp1)='Clustering'
    datatemp2=dataPatterns[[3]]$listSilhouetteCoefPerClusters
    names(datatemp2)=paste('Cluster', names(datatemp2))
    
    datatemp=signif(c(datatemp1,datatemp2),3)
        
    a <- rCharts:::Highcharts$new()
    a$chart(type = "column")
    a$xAxis( categories=names(datatemp), labels=list(rotation=-45))
    a$yAxis(title=list(text='Silhouette coef. (avg.)'),max=1)
    a$legend(enabled = FALSE)
    a$plotOptions(column=list(colorByPoint=TRUE))
    a$colors(c('#C0C0C0',dataPatterns[[3]]$colorClusters))
    
    a$data(data=as.vector(datatemp))

    
    a$tooltip(
      formatter = "#! function() {
              
              return '<table>'
              + '<center>'+this.x+'</center>'
              + '<br>'
              + '<center>'+this.y+'</center>'
              +'</table>';} !#"
    )
    
    a$set(width="100%",heigth='100%')
    
    return(a)
    
    
  })
  
  output$barChartClustSizeProp <- renderChart2({
    if (is.null(data_res$data)) return()
    dataPatterns=data_res$data
    
    #GA
    shinyjs::runjs(HTML("ga('send', 'event', 'Explore', 'Cluster features');"));  
    
    datatemp=signif((dataPatterns[[3]]$clustSize/sum(dataPatterns[[3]]$clustSize))*100,3)
    
    a <- rCharts:::Highcharts$new()
    a$chart(type = "column")
    a$xAxis( categories=paste('Cluster', names(datatemp)), labels=list(rotation=-45))
    a$yAxis(title=list(text='Number of structures (%)'), max=max(datatemp))
    a$legend(enabled = FALSE)
    a$plotOptions(column=list(colorByPoint=TRUE))
    a$data(data=datatemp)
    a$colors(dataPatterns[[3]]$colorClusters)
    a$tooltip(
      formatter = "#! function() {
              return '<table>'
              + '<center>'+this.x+'</center>'
              + '<br>'
              + '<center>'+this.y+'%</center>'
              +'</table>';} !#"
    )
    
    a$set(width="100%",heigth='100%')
    
    return(a)
  })
  
  output$barChartClustSize <- renderChart2({
    if (is.null(data_res$data)) return()
    dataPatterns=data_res$data
    
    datatemp=signif(dataPatterns[[3]]$clustSize,3)
    
    a <- rCharts:::Highcharts$new()
    a$chart(type = "column")
    a$xAxis( categories=paste('Cluster', names(datatemp)), labels=list(rotation=-45))
    a$yAxis(title=list(text='Number of structures'), max=max(datatemp))
    a$legend(enabled = FALSE)
    a$plotOptions(column=list(colorByPoint=TRUE))
    a$data(data=datatemp)
    a$colors(dataPatterns[[3]]$colorClusters)
    a$tooltip(
      formatter = "#! function() {
              return '<table>'
              + '<center>'+this.x+'</center>'
              + '<br>'
              + '<center>'+this.y+'</center>'
              +'</table>';} !#"
    )
    
    a$set(width="100%",heigth='100%')
    
    return(a)
  })
  
  output$barChartStructVar<- renderChart2({
    if (is.null(data_res$data)) return()
    dataPatterns=data_res$data
    
    datatemp=signif(dataPatterns[[3]]$confStabilit,3)
    
    a <- rCharts:::Highcharts$new()
    a$chart(type = "column")
    a$xAxis( categories=paste('Cluster', 1:length(datatemp)), labels=list(rotation=-45))
    a$yAxis(title=list(text='Intra-cluster distance'), max=max(datatemp))
    a$legend(enabled = FALSE)
    a$plotOptions(column=list(colorByPoint=TRUE))
    a$data(data=datatemp)
    a$colors(dataPatterns[[3]]$colorClusters)
    a$tooltip(
      formatter = "#! function() {
              return '<table>'
              + '<center>'+this.x+'</center>'
              + '<br>'
              + '<center>'+this.y+'</center>'
              +'</table>';} !#"
    )
    
    a$set(width="100%",heigth='100%')
    return(a)
  })
  
  output$boxplotLengthDist<- renderChart2({
    if (is.null(data_res$data)) return()
    dataPatterns=data_res$data
    
    bwtLength=c()
    
    for (i in 1:length(dataPatterns[[3]]$unikidxClustering)) 
    {
      clusCurrentInd=which(dataPatterns[[3]]$idxClustering %in% dataPatterns[[3]]$unikidxClustering[i])
      LengthDistribCurrentClust=dataPatterns[[1]]$lengthRNAs[clusCurrentInd]
      temp_boxplotLength=boxplot(LengthDistribCurrentClust, plot = FALSE)
      bwtLength=cbind(bwtLength,temp_boxplotLength$stat)
    }
    datatemp = setNames(as.data.frame(bwtLength),nm = NULL)
    
    a = Highcharts$new()
    a$chart(type = "boxplot")
    a$data(data = datatemp)
    
    a$xAxis( categories=paste('Cluster',1:dim(datatemp)[2]), labels=list(rotation=-45))
    a$yAxis(title=list(text='Number of nucleotides'))
    a$legend(enabled = FALSE)
    a$plotOptions(boxplot=list(colorByPoint=TRUE))
    a$colors(dataPatterns[[3]]$colorClusters)
    
    a$set(width="100%",heigth='100%')
    
    return(a)
  })
  
  
  output$clustInfo = DT::renderDataTable({
    if (is.null(data_res$data)) return()
    dataPatterns=data_res$data
    
    dataPatterns[[8]]
    })
  
  #Patterns tab  
  output$tbl = DT::renderDataTable({
    if (is.null(data_res$data)) return()
    dataPatterns=data_res$data
    dataPatterns[[9]]
  })
    
  #Tab of structures
    output$tbOfStruct = DT::renderDataTable(
      data_res$data[[4]][,1:5], server = FALSE,
      style = 'bootstrap',filter = 'top',
      rownames = FALSE, extensions = 'TableTools', options = list(searchHighlight = TRUE,search = list(regex = TRUE),
      dom = 'T<"clear">lfrtip',
      tableTools = list(sSwfPath = copySWF())
      )
    )
  
  #Dendrogram
#   output$dendSS <- renderTreewidget({
#     if (is.null(data_res$data)) return()
#     dataPatterns=data_res$data
#     
#     if (input$bootstrap>0 )
#     {
#       #dataPatterns[[3]]$resultClust %>% as.dendrogram %>%  color_branches(.,k=dataPatterns[[3]]$bestNbClusters,col=dataPatterns[[3]]$colorClusters) %>% plot(main = "Cluster dendrogram  (using pvclust and UPGMA) with AU/BP values (%)\n",horiz=T)
#       #return(dataPatterns[[3]]$resultClust %>% text)
#     }
#     else
#     {
#       
#       idxClusteringSortedAccordingToDendrogram=dataPatterns[[3]]$idxClustering[order.dendrogram(as.dendrogram(dataPatterns[[3]]$resultClust))]
#       headersStructure=paste(labels(dataPatterns[[3]]$resultClust)," cluster",idxClusteringSortedAccordingToDendrogram,sep="")
#       labels(dataPatterns[[3]]$resultClust)=headersStructure
#       return((treewidget(as.phylo(dataPatterns[[3]]$resultClust),browser = TRUE )))
#     }    
#     })
    
#Explore
  
  #Scatplot
  output$scatplotsnm <- renderChart2({
    if (is.null(data_res$data)) return()
    
    dataPatterns=data_res$data
    
    if (is.null(dataPatterns)) 
    {return(dataPatterns[[5]])}
    else{
      #GA
      shinyjs::runjs(HTML("ga('send', 'event', 'Explore', 'Feature visualization');"));   
      
      
      snm_x=as.numeric(input$snm_x)
      snm_y=as.numeric(input$snm_y)
      scatPlotSNM(dataPatterns[[2]]$SuperMotif[,snm_x],dataPatterns[[2]]$SuperMotif[,snm_y],dataPatterns[[3]]$bestNmotifsForClusters_coord[,snm_x],dataPatterns[[3]]$bestNmotifsForClusters_coord[,snm_y],snm_x,snm_y,dataPatterns[[3]]$bestNmotifsForClusters_pos,dataPatterns[[3]]$bestNmotifsForClusters_pos_namesNmotifs,dataPatterns[[2]]$singularValuesPercent,dataPatterns[[1]]$headers,dataPatterns[[3]]$idxClustering,dataPatterns[[3]]$colorClusters,dataPatterns[[3]]$repClustList,dataPatterns[[3]]$outliers,dataPatterns[[3]]$matnmPosBestnmotifsforForna,dataPatterns[[1]])
    
      
      
    }
  })
  
  output$distHeader1Header2<-renderUI({
    if(!is.null(input$click$header1)&&!(input$click$header1=="")&&!is.null(input$click$header2)&&!(input$click$header2==""))
    {    
      clusCurrentIndheader1=which(data_res$data[[1]]$headers %in% input$click$header1)#dissimilarity matrix based on cosine similarity -->test with output of snm (faster)
      clusCurrentIndheader2=which(data_res$data[[1]]$headers %in% input$click$header2)#dissimilarity matrix based on cosine similarity -->test with output of snm (faster)
      distTemp= signif(data_res$data[[3]]$matDissim[clusCurrentIndheader1,clusCurrentIndheader2],digits=3)
       shinyjs::show("structDistInfo")
      return(
        
        HTML(input$click$header1,' and ', input$click$header2,' selected.<br>',
        paste("Structural distance (",input$click$header1," , ",input$click$header2,") = ",distTemp,sep=""))
        )
      }
    else if (!is.null(input$click$header1)&&!(input$click$header1==""))
      {
      return(
        HTML(input$click$header1,' selected.<br>Ctrl key+ <i class="fa fa-hand-pointer-o"></i> to visualize a second secondary structure.' )
        
        )
      }
    else 
      {
      return(
        HTML('<i class="fa fa-hand-pointer-o"></i> to visualize a secondary structure')
      )
      }
  })
  
  
  #VarExplained
  output$varExp <- renderChart2({
    if (is.null(data_res$data)) return()
    dataPatterns=data_res$data
    
    dataPatterns[[6]]
    })
  
  #Secondary structures to explore
  output$headerss1<-renderUI({
    if(!is.null(input$click$header1)&&!(input$click$header1==""))
      {return(paste(input$click$header1,"|",input$click$length1,"nt.","|cluster", input$click$idxClustering1,sep=""))}
    else{return(HTML('<i class="fa fa-hand-pointer-o"></i> to visualize a secondary structure'))}
  })
  
  output$rna_ss1 <- renderUI({
    if (is.null(data_res$data)) return()
    if (is.null(input$click$header1)) return()
    
    dataPatterns=data_res$data
    if (input$click$header1=="") 
    {return()}
    else{
    

      struct1=input$click$structure1
      
        #Gquadruplexes color
        g4PoswithColors=""
        g4PosTemp=gregexpr("\\+", struct1,ignore.case=FALSE)[[1]]
        if(g4PosTemp[1]!=-1){
          g4Pos=g4PosTemp[1:length(g4PosTemp)]
          g4PoswithColors=paste(g4Pos,":rgb(158,158,158)",sep="")
          g4PoswithColors=paste(g4PoswithColors,collapse=" ")
          struct1=gsub("\\+", "\\.", struct1)
        }
        
        #circular RNA: add * at the end of the dotbracket 
        if(!is.null(input$click$header1)){
          if(substr(input$click$header1,0,2)=="c_"){struct1=paste(struct1,"*",sep="")} 
        }
        
        #Viz n-motif region
        if(input$viznmotifs==TRUE)
        {
          if (!is.null(data_res_updatesetnmotifs$data))
          {
            indTempHeader1=which(data_res_updatesetnmotifs$data[[1]]$headers %in% input$click$header1)#dissimilarity matrix based on cosine similarity -->test with output of snm (faster)
            nmotifsRegion1=data_res_updatesetnmotifs$data[[3]]$matnmPosBestnmotifsforForna[indTempHeader1]
          }
          else
          {
            indTempHeader1=which(data_res$data[[1]]$headers %in% input$click$header1)#dissimilarity matrix based on cosine similarity -->test with output of snm (faster)
            nmotifsRegion1=dataPatterns[[3]]$matnmPosBestnmotifsforForna[indTempHeader1]
          }
        }
        else{nmotifsRegion1=""}
        
        
        #Visualize structure
        atemp=paste("
                    var container = new FornaContainer(\"#rna_ss1\",
                    {'applyForce': ", tolower(input$applyForce), 
                    ",'allowPanningAndZooming': true, 'initialSize':[300,400]});
                    var options = {
                    'structure': '",
                    struct1,
                    "', 'sequence': '",
                    input$click$sequence1,
                    "'};container.addRNA(options.structure, options);
                     var cs = new ColorScheme('",
                    nmotifsRegion1," ",g4PoswithColors,
                     "');
                     container.addCustomColors(cs.colorsJson);
                     container.changeColorScheme('custom'); 
                     ",
                    sep="")
        
        a=tags$script(HTML(atemp))  
        return(a)
    }
    })
  
  output$headerss2<-renderUI({
    if(!is.null(input$click$header2)&&!(input$click$header2==""))
    {return(paste(input$click$header2,"|",input$click$length2,"nt.","|cluster", input$click$idxClustering2,sep=""))}
    else{return( HTML('Ctrl key+ <i class="fa fa-hand-pointer-o"></i> to visualize a second secondary structure.' ))}
  })
  
  output$rna_ss2 <- renderUI({
    if (is.null(data_res$data)) return()
    if (is.null(data_res$data)) return()
    if (is.null(input$click$header2)) return()
    dataPatterns=data_res$data
    
    if (is.null(dataPatterns)&&is.null(input$click$header2)&&(input$click$header2=="")) 
    {return(dataPatterns[[5]])}
    else{
      
      
      struct2=input$click$structure2
      #Gquadruplexes color
      g4PoswithColors=""
      g4PosTemp=gregexpr("\\+", struct2,ignore.case=FALSE)[[1]]
      if(g4PosTemp[1]!=-1){
        g4Pos=g4PosTemp[1:length(g4PosTemp)]
        g4PoswithColors=paste(g4Pos,":rgb(158,158,158)",sep="")
        g4PoswithColors=paste(g4PoswithColors,collapse=" ")
        struct2=gsub("\\+", "\\.", struct2)
      }
      
      #circular RNA: add * at the end of the dotbracket 
      if(!is.null(input$click$header2)){
        if(substr(input$click$header2,0,2)=="c_"){struct2=paste(input$click$structure2,"*",sep="")} 
      }
      
      #Viz n-motif region
      if(input$viznmotifs==TRUE)
      {
        if (!is.null(data_res_updatesetnmotifs$data))
        {
          indTempHeader2=which(data_res_updatesetnmotifs$data[[1]]$headers %in% input$click$header2)#dissimilarity matrix based on cosine similarity -->test with output of snm (faster)
          nmotifsRegion2=data_res_updatesetnmotifs$data[[3]]$matnmPosBestnmotifsforForna[indTempHeader2]
        }
        else
        {
          indTempHeader2=which(data_res$data[[1]]$headers %in% input$click$header2)#dissimilarity matrix based on cosine similarity -->test with output of snm (faster)
          nmotifsRegion2=dataPatterns[[3]]$matnmPosBestnmotifsforForna[indTempHeader2]
        }
        
      }
      else{nmotifsRegion2=""}
      
      
      
      #Visualize structure
      atemp=paste("
                    var container = new FornaContainer(\"#rna_ss2\",
                    {'applyForce': ", tolower(input$applyForce), 
                  ",'allowPanningAndZooming': true, 'initialSize':[300,400]});
                    var options = {
                    'structure': '",
                  struct2,
                  "', 'sequence': '",
                  input$click$sequence2,
                  "'};container.addRNA(options.structure, options);
                     var cs = new ColorScheme('",
                  nmotifsRegion2," ",g4PoswithColors,
                  "');
                     container.addCustomColors(cs.colorsJson);
                     container.changeColorScheme('custom'); 
                  ",
                  sep="")  
      
      a=tags$script(HTML(atemp))  
    }
    return(a)
  })
    

  
  observeEvent(input$buttonPrepare1, {
    shinydashboard::updateTabItems(session, "menu","prepare")
    #GA
    shinyjs::runjs(HTML("ga('send', 'event', 'button', 'Go to prepare page (button prepare 1)');"));
  })
  observeEvent(input$buttonExplore, {
    shinydashboard::updateTabItems(session, "menu","explore")
    #GA
    shinyjs::runjs(HTML("ga('send', 'event', 'button', 'Go to explore page');"));
  })
  observeEvent(input$buttonPrepare2, {
    shinydashboard::updateTabItems(session, "menu","prepare")
    #GA
    shinyjs::runjs(HTML("ga('send', 'event', 'button', 'Go to prepare page (button prepare 2)');"));
  })
  
  onclick("toggleAdvancedrnad", toggle(id = "advancedrnad", anim = TRUE))
  onclick("toggleAdvancedrnabpd", toggle(id = "advancedrnabpd", anim = TRUE))
  

  output$exportDB <- downloadHandler(
      filename = function() {
        paste('secStruct-', Sys.Date(), '.db', sep='')
      },
      content = function(fileConn) {
        
            #GA
            shinyjs::runjs(HTML("ga('send', 'event', 'button', 'Download Filtered struct. in dot-bracket');"));
        
            cat("", file=fileConn, append=FALSE, sep='')

            for (i in 1:length(input$tbOfStruct_rows_all)) 
            {
              temp=input$tbOfStruct_rows_all[i]
              #which(data_res$data[[4]] %in% data_res$data[[4]][i,1])
              
              cat(">", as.character(paste(data_res$data[[4]][temp,1],"|",data_res$data[[4]][temp,2],
                                          "|",data_res$data[[4]][temp,3],"|",data_res$data[[4]][temp,4],"|",data_res$data[[4]][temp,5],sep=""))
                  , "\n", file=fileConn, append=TRUE, sep='')
              cat(as.character(data_res$data[[4]][temp,6]), "\n", file=fileConn, append=TRUE, sep='')
              cat(as.character(data_res$data[[4]][temp,7]), "\n", file=fileConn, append=TRUE, sep='')
            
            }
      }
    )
  
  output$exportdbStruct1 <- downloadHandler(
    filename = function() {
      paste(input$click$header1,'-', Sys.Date(), '.db', sep='')
    },
    content = function(fileConn) {
      cat("", file=fileConn, append=FALSE, sep='')
      
        cat(">",as.character(input$click$header1), "\n", file=fileConn, append=TRUE, sep='')
        
        cat(as.character(input$click$sequence1), "\n", file=fileConn, append=TRUE, sep='')
        cat(as.character(input$click$structure1), "\n", file=fileConn, append=TRUE, sep='')
        
      
    }
  )
  
  output$exportdbStruct2 <- downloadHandler(
    filename = function() {
      paste(input$click$header2,'-', Sys.Date(),'.db', sep='')
    },
    content = function(fileConn) {
      cat("", file=fileConn, append=FALSE, sep='')
      
      cat(">",as.character(input$click$header2), "\n", file=fileConn, append=TRUE, sep='')
      
      cat(as.character(input$click$sequence2), "\n", file=fileConn, append=TRUE, sep='')
      cat(as.character(input$click$structure2), "\n", file=fileConn, append=TRUE, sep='')
      
      
    }
  )
  
  #Dendrogram###################################
  output$legendHierarchy <- renderUI({
  
    legend=''
    for (j in 1:data_res$data[[3]]$bestNbClusters) 
    {
      legend=paste(legend,'<i class="fa fa-minus" style="color:',data_res$data[[3]]$colorClusters[j],';"></i>','<b> Cluster ',j,'</b> ',sep="")
    }
    return(HTML(legend))
  })
  
  
  output$dendSS2 <- renderUI({
    
    data_res$data
    
    #GA
    shinyjs::runjs(HTML("ga('send', 'event', 'Explore', 'Cluster and structure hierarchy');"));  
    
    if(input$bootstrap>0)
    {
      if(input$AU_bootstrap_values==1)
      {values=data_res$data[[3]]$AU_values}
      else{values=data_res$data[[3]]$bootstrap_values}
      inputTree=write.tree(as.phylo.hclust.with.nodenames(data_res$data[[3]]$resultClust, nodenames=values)
                           ,tree.names=TRUE,digits=2)
    }else
    {
    inputTree = write.tree( as.phylo(data_res$data[[3]]$resultClust))
    }
    
    inputClust=rjson::toJSON( data_res$data[[3]]$idxClustering)

    a=rjson::toJSON(as.character(data_res$data[[3]]$unikidxClustering))
    b=rjson::toJSON(data_res$data[[3]]$colorClusters)
    

    
    atemp=paste('
// the global tree variable

                //var inputClust= {"PSTVd": 1, "TASVd": 1, "TCDVd": 1, "CSVd": 1, "CLVd": 1, "CCCVd": 1, "ASSVd": 1, "PBCVd": 1, "CVd-OS": 1, "CbVd-1": 1, "CbVd-2": 1, "CbVd-3": 1, "CCHMVd": 2, "ELVd": 2, "ASBVd": 2, "HSVd": 1, "CVdIV": 1, "CVDIII": 1, "CBLVd": 1, "CEVd": 1, "PLMVd": 2};

                //var inputTree="(((((((TASVd:0.06828169,CSVd:0.06828169)0.8170:0.05635267,CEVd:0.12463435)0.8290:0.06964112,(PSTVd:0.02617095,TCDVd:0.02617095)1.0000:0.16810452)0.6700:0.08751399,CLVd:0.28178947)0.8700:0.13005822,(HSVd:0.26704011,(CCCVd:0.20577326,CVdIV:0.20577326)0.6520:0.06126685)0.8790:0.14480757)0.8810:0.18306776,(((CbVd-2:0.08236826,CbVd-3:0.08236826)0.9300:0.10543001,CbVd-1:0.18779827)0.9970:0.29222981,(PBCVd:0.38569588,(CBLVd:0.30671661,(ASSVd:0.25961362,(CVd-OS:0.14517605,CVDIII:0.14517605)0.9870:0.11443757)0.5950:0.04710299)0.5660:0.07897927)0.7300:0.09433220)0.6860:0.11488736)0.9840:0.46199326,((CCHMVd:0.37098428,PLMVd:0.37098428)0.8850:0.20471288,(ELVd:0.37028264,ASBVd:0.37028264)0.8730:0.20541453)0.9840:0.48121154);";
                
                //var inputTree="(((((((TASVd:0.06828169,CSVd:0.06828169)0.8170:0.05635267,CEVd:0.12463435)0.8290:0.06964112,(PSTVd:0.02617095,TCDVd:0.02617095)1.0000:0.16810452)0.6700:0.08751399,CLVd:0.28178947)0.8700:0.13005822,(HSVd:0.26704011,(CCCVd:0.20577326,CVdIV:0.20577326)0.6520:0.06126685)0.8790:0.14480757)0.8810:0.18306776,(((CbVd-2:0.08236826,CbVd-3:0.08236826)0.9300:0.10543001,CbVd-1:0.18779827)0.9970:0.29222981,(PBCVd:0.38569588,(CBLVd:0.30671661,(ASSVd:0.25961362,(CVd-OS:0.14517605,CVDIII:0.14517605)0.9870:0.11443757)0.5950:0.04710299)0.5660:0.07897927)0.7300:0.09433220)0.6860:0.11488736)0.9840:0.46199326,((CCHMVd:0.37098428,PLMVd:0.37098428)0.8850:0.20471288,(ELVd:0.37028264,ASBVd:0.37028264)0.8730:0.20541453)0.9840:0.48121154);";

                var inputClust=',inputClust,'; 

                var inputTree="',inputTree,'";

                var tree;
                 
                // the dictionary mapping node names to clusters
                var clustering = {};
                
                // default scheme to color by date    
                //var coloring_scheme = d3.scale.category10();
var coloring_scheme = d3.scale.ordinal()
    .domain(',a,')
    .range(',b,');
                
                // this will be used to map bootstrap support values to edge thickness
                var bootstrap_scale = d3.scale.linear().domain ([0,0.5,0.7,0.9,0.95,1]).range ([1,2,3,4,5,6]).interpolate (d3.interpolateRound);
                


              //var control 
              var container_id = "#dendSS2";

              var width  = $(container_id).width()*0.90,
                  height = $(container_id).height()*0.90,
                  selection_set = [\'Foreground\'],
                  current_selection_name = $("#selection_name_box").val(),
                  current_selection_id = 0,
                  max_selections       = 10;

                  color_scheme = d3.scale.category10(),

                  selection_menu_element_action = "phylotree_menu_element_action";
    

                function edgeStyler (dom_element, edge_object) {
                //if ("bootstrap" in edge_object.target) {
                //dom_element.style ("stroke-width", bootstrap_scale (edge_object.target.bootstrap) + "pt");
                //}
                dom_element.style ("stroke", "cluster" in edge_object.target ? coloring_scheme (edge_object.target.cluster) : null);
                
                }
                
                function nodeStyler (dom_element, node_object) {
                if ("bootstrap" in node_object && node_object.bootstrap) { 
                var label = dom_element.selectAll (".bootstrap");
                if (label.empty()) {
                dom_element.append ("text").classed ("bootstrap", true).text (node_object.bootstrap).attr ("dx", ".3em").attr("text-anchor", "start").attr ("alignment-baseline", "middle");
                } else {
                if (tree.radial()) { // do not show internal node labels in radial mode
                label.remove();
                }
                }
                }
                }    
                
              
                  var svg=d3.select(container_id).append("svg")
                                 .attr("width",width)
                                 .attr("height",height)
                                 .attr("version","1.1")
                                 .attr("xmlns","http://www.w3.org/2000/svg")
                                 .attr("xmlns:xlink","http://www.w3.org/1999/xlink")

;

                function drawATree(newick) {

                tree = d3.layout.phylotree(container_id)
               //.svg(d3.select(container_id).append("svg"))
               .svg(svg)
                .separation (function (a,b) {return 0;})
                .count_handler (function (count) { 
                        $("#selected_branch_counter").text (function (d) {return count[current_selection_name];}); 
                        $("#selected_filtered_counter").text (count.tag);
                    }
                )
                .options({
                          \'left-right-spacing\'   : \'fit-to-size\', 
                          // fit to given size top-to-bottom
                         // \'top-bottom-spacing\'   : \'fit-to-size\',
               
                \'selectable\': true,
                \'collapsible\': true,
                \'transitions\' : true 

                })
                .size([height, width])
                .style_edges(edgeStyler)
                .style_nodes (nodeStyler)
                ;

tree.branch_length (function (n) {return undefined;});

                /* the next call creates the tree object, and tree nodes */
                tree(d3_phylotree_newick_parser(newick));

                
                // parse bootstrap support from internal node names
                _.each (tree.get_nodes(), function (node) {
                if (node.children) {
                node.bootstrap = parseFloat (node.name);
                }   
                });
                
                //tree.spacing_x (50).spacing_y (50);

                tree.spacing_x (23);
               // tree.spacing_y (50);

                if ($("#layout").prop("checked")) {
                tree.radial (true);
                }
                tree.placenodes().layout();
                
                
                // UI handlers
                $("#layout").on("click", function(e) {
                tree.radial($(this).prop("checked")).placenodes().update();
                });

                }



                function applyAnnotation (clustering) {
                //coloring_scheme.domain ([]); // reset the coloring scheme
                if (tree) {
                     // refresh cluster assignments for tips
                    /*_.each (tree.get_nodes(), function (node) {
                        if (node.name in clustering) {
                            node.cluster = clustering[node.name];
                        } else {
                            delete node.cluster;
                        }
                    });*/
    
                    
                    tree.traverse_and_compute (function (node) {
                        if (node.name in clustering) {
                                node.cluster = clustering[node.name];
                            } else {
                                delete node.cluster;
                                var children_clusters = _.keys(_.countBy (node.children, function (d) {
                                    return d.cluster;
                                }));
                                if (children_clusters.length == 1 && children_clusters[0]) {
                                    node.cluster =  children_clusters[0];
                                }
                            }
                        },
                        "post-order");
                    
                    tree.update();
                    
                    // update the legend
                    
                    d3.select ("#cluster-legend").selectAll (".row").remove();
                    var cluster_colors = d3.select ("#cluster-legend").selectAll (".row").data (coloring_scheme.domain().sort().map (function (d) {return [d];}));
                    cluster_colors.enter().append ("div").classed ("row", true);
                    cluster_colors.exit().remove();
                    cluster_colors = cluster_colors.selectAll ("span").data (function (d) {return d});
                    cluster_colors.enter().append ("span").classed ("cluster-text", true);
                    
                    cluster_colors.each (function (d) {
                    d3.select(this).style ("color", coloring_scheme (d), "important").classed ("cluster-text", true).text ("Cluster " + d);
                    });
                    
                 }
            }

//***************************************************************

                          //************************CONTROL**************************





$("#mp_label").on ("click", function (e) {
    tree.max_parsimony (true);
});
$ ("[data-direction]").on ("click", function (e) {
    var which_function = $(this).data ("direction") == \'vertical\' ? tree.spacing_x : tree.spacing_y;
    which_function (which_function () + (+ $(this).data ("amount"))).update();
}); 
$(".phylotree-layout-mode").on ("change", function (e) {
    if ($(this).is(\':checked\')) {
        if (tree.radial () != ($(this).data ("mode") == "radial")) {
            tree.radial (!tree.radial ()).placenodes().update ();
        }
    }
});         

       
function sort_nodes (asc) {
    tree.traverse_and_compute (function (n) {
            var d = 1;
            if (n.children && n.children.length) {
                d += d3.max (n.children, function (d) { return d["count_depth"];});
            }
            n["count_depth"] = d;
        }); 
        tree.resort_children (function (a,b) {
            return (a["count_depth"] - b["count_depth"]) * (asc ? 1 : -1);
        });
}


  $("#sort_original").on ("click", function (e) {
      
      drawATree(inputTree);
      applyAnnotation(inputClust);
    });


    $("#sort_ascending").on ("click", function (e) {
    sort_nodes (true);
    });
    $("#sort_descending").on ("click", function (e) {
    sort_nodes (false);
    });

    $("#save_tree").on ("click", function (e) {
      var tagged_tree=tree.get_newick (
            function (node) {
                var tags = [];
                selection_set.forEach (function (d) { if (node[d]) {tags.push(d)}; });
                if (tags.length) {
                    return "{" + tags.join (",") + "}";
                }
                return "";
            }
       );
      var blob = new Blob([tagged_tree], {type: "text/plain;charset=utf-8"});
      saveAs(blob, "phylowidget_tree.nwk");
    });

    $("#and_label").on ("click", function (e) {
    tree.internal_label (function (d) { return d.reduce (function (prev, curr) { return curr[current_selection_name] && prev; }, true)}, true);
    });
    $("#or_label").on ("click", function (e) {
    tree.internal_label (function (d) { return d.reduce (function (prev, curr) { return curr[current_selection_name] || prev; }, false)}, true);
    });
    $("#filter_add").on ("click", function (e) {
    tree.modify_selection (function (d) { return d.tag || d[current_selection_name];}, current_selection_name, false, true)
    .modify_selection (function (d) { return false; }, "tag", false, false);
    });
    $("#filter_remove").on ("click", function (e) {
    tree.modify_selection (function (d) { return !d.tag;});
    });
    $("#select_all").on ("click", function (e) {
    tree.modify_selection (function (d) { return true;});
    });
    $("#select_all_internal").on ("click", function (e) {
    tree.modify_selection (function (d) { return !d3_phylotree_is_leafnode (d.target);});
    });
    $("#select_all_leaves").on ("click", function (e) {
    tree.modify_selection (function (d) { return d3_phylotree_is_leafnode (d.target);});
    });
    $("#select_none").on ("click", function (e) {
    tree.modify_selection (function (d) { return false;});
    });
    $("#clear_internal").on ("click", function (e) {
    tree.modify_selection (function (d) { return d3_phylotree_is_leafnode (d.target) ? d.target[current_selection_name] : false;});
    });
    $("#clear_leaves").on ("click", function (e) {
    tree.modify_selection (function (d) { return !d3_phylotree_is_leafnode (d.target) ? d.target[current_selection_name] : false;});
    });
    $("#display_dengrogram").on ("click", function (e) {
    tree.options ({\'branches\' : \'step\'}, true);
    });


$("#branch_filter").on ("input propertychange", function (e) {  
   var filter_value = $(this).val();
   
   var rx = new RegExp (filter_value,"i");
   
  tree.modify_selection (function (n) { 
    return filter_value.length && (tree.branch_name () (n.target).search (rx)) != -1; 
   },"tag");
      
});


    var valid_id = new RegExp ("^[\\w]+$");
    $("#selection_name_box").on ("input propertychange", function (e) {  
    var name = $(this).val();         
    
    var accept_name = (selection_set.indexOf (name) < 0) &&
    valid_id.exec (name) ;
    
    d3.select ("#save_selection_button").classed ("disabled", accept_name ? null : true ); 
    });
    $("#selection_rename > a").on ("click", function (e) {
    d3.select ("#save_selection_button")
    .classed ("disabled",true)
    .on ("click", function (e) { // save selection handler
    var old_selection_name = current_selection_name;
    selection_set[current_selection_id] = current_selection_name = $("#selection_name_box").val();
    
    if (old_selection_name != current_selection_name) {
    tree.update_key_name (old_selection_name, current_selection_name);
    update_selection_names (current_selection_id);
    }
    send_click_event_to_menu_objects (new CustomEvent (selection_menu_element_action,
    {\'detail\' : [\'save\', this]}));
    });
    
    d3.select ("#cancel_selection_button")
    .classed ("disabled",false)
    .on ("click", function (e) { // save selection handler
    $("#selection_name_box").val(current_selection_name);
    send_click_event_to_menu_objects (new CustomEvent (selection_menu_element_action,
    {\'detail\' : [\'cancel\', this]}));
    });
    send_click_event_to_menu_objects (new CustomEvent (selection_menu_element_action,
    {\'detail\' : [\'rename\', this]}));
    e.preventDefault    (); 
    });
    $("#selection_delete > a").on ("click", function (e) {
    
    tree.update_key_name (selection_set[current_selection_id], null)
    selection_set.splice (current_selection_id, 1);
    
    if (current_selection_id > 0) {
    current_selection_id --;
    }
    current_selection_name = selection_set[current_selection_id];
    update_selection_names (current_selection_id)
    $("#selection_name_box").val(current_selection_name)
    
    
    send_click_event_to_menu_objects (new CustomEvent (selection_menu_element_action,
    {\'detail\' : [\'save\', this]}));
    e.preventDefault    ();
    
    });  
    
    $("#selection_new > a").on ("click", function (e) {
    
    d3.select ("#save_selection_button")
    .classed ("disabled",true)
    .on ("click", function (e) { // save selection handler
    current_selection_name = $("#selection_name_box").val();
    current_selection_id = selection_set.length;
    selection_set.push (current_selection_name);
    update_selection_names (current_selection_id);
    send_click_event_to_menu_objects (new CustomEvent (selection_menu_element_action,
    {\'detail\' : [\'save\', this]}));
    });
    
    d3.select ("#cancel_selection_button")
    .classed ("disabled",false)
    .on ("click", function (e) { // save selection handler
    $("#selection_name_box").val(current_selection_name);
    send_click_event_to_menu_objects (new CustomEvent (selection_menu_element_action,
    {\'detail\' : [\'cancel\', this]}));
    });
    
    send_click_event_to_menu_objects (new CustomEvent (selection_menu_element_action,
    {\'detail\' : [\'new\', this]}));
    e.preventDefault    ();
    
    });  
    function send_click_event_to_menu_objects (e) {
    $("#selection_new, #selection_delete, #selection_rename, #save_selection_name, #selection_name_box, #selection_name_dropdown").get().forEach (
    function (d) {
    d.dispatchEvent (e);
    }
    );
    }
    
    function update_selection_names (id, skip_rebuild) {
    skip_rebuild = skip_rebuild || false;
    id = id || 0;
    
    
    current_selection_name = selection_set[id];
    current_selection_id = id;
    
    if (!skip_rebuild) {
    d3.selectAll (".selection_set").remove();
    
    d3.select ("#selection_name_dropdown")
    .selectAll (".selection_set")
    .data (selection_set)
    .enter()
    .append ("li")
    .attr ("class", "selection_set")
    .append ("a")
    .attr ("href", "#")
    .text (function (d) { return d;})
    .style ("color", function (d,i) {return color_scheme(i);})
    .on ("click", function (d,i) {update_selection_names (i,true);});
    
    }
    
    
    d3.select ("#selection_name_box")
    .style ("color",  color_scheme(id))
    .property ("value", current_selection_name);
    
    tree.selection_label (selection_set[id]);
    }

function selection_handler_name_box (e) {
    var name_box = d3.select (this);
    switch (e.detail[0]) {
    case \'save\':
    case \'cancel\':
    name_box.property ("disabled", true)
    .style ("color",  color_scheme(current_selection_id));
    break;
    case \'new\':
    name_box.property ("disabled", false)
    .property ("value", "new_selection_name")
    .style ("color",  color_scheme(selection_set.length));
    break;
    case \'rename\':
    name_box.property ("disabled", false);
    break;
    }
    
  }
    function selection_handler_new (e) {
    var element = d3.select (this);
    $(this).data(\'tooltip\', false);         
    switch (e.detail[0]) {
    case \'save\':
    case \'cancel\':
    if (selection_set.length == max_selections) {
    element.classed ("disabled", true);
    $(this).tooltip ({\'title\' : \'Up to \' + max_selections + \' are allowed\', \'placement\' : \'left\'});
    } else {
    element.classed ("disabled", null);
    }
    break;
    default:
    element.classed ("disabled", true);
    break;
    
    }
    }
    function selection_handler_rename (e) {
    var element = d3.select (this);
    element.classed ("disabled", (e.detail[0] == "save" || e.detail[0] == "cancel") ? null : true);
    }
    function selection_handler_save_selection_name (e) {
    var element = d3.select (this);
    element.style ("display", (e.detail[0] == "save" || e.detail[0] == "cancel") ? "none" : null);
    }
    function selection_handler_name_dropdown (e) {
    var element = d3.select (this).selectAll (".selection_set");
    element.classed ("disabled", (e.detail[0] == "save" || e.detail[0] == "cancel") ? null : true);
    }
    function selection_handler_delete (e) {
    var element = d3.select (this);
    $(this).tooltip(\'destroy\');         
    switch (e.detail[0]) {
    case \'save\':
    case \'cancel\':
    if (selection_set.length == 1) {
    element.classed ("disabled", true);
    $(this).tooltip ({\'title\' : \'At least one named selection set <br> is required;<br>it can be empty, however\', \'placement\' : \'bottom\', \'html\': true});
    } else {
    element.classed ("disabled", null);
    }
    break;
    default:
    element.classed ("disabled", true);
    break;
    
    }}

    //************************DRAW TREE AND CLUSTER ANNOTATION**************************



      drawATree(inputTree);
      applyAnnotation(inputClust);


    $(\'#newick_export_modal\').on(\'show.bs.modal\', function (e) {
        $(\'textarea[id$="nwk_export_spec"]\').val(
          tree.get_newick (
            function (node) {
              var tags = [];
              selection_set.forEach (function (d) { if (node[d]) {tags.push(d)}; });
              if (tags.length) {
                return "{" + tags.join (",") + "}";
              }
              return "";
            }  
          )
        );
      })

                '
                ,sep="")
    
   # shinyjs::runjs(HTML(atemp));
    
   
    a=tags$script(HTML(atemp)) 
    
    
    
    return(a)
    
  })
  
  
  
})

