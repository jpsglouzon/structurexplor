#check R version
R_min_version = "3.1.9"
R_version = paste0(R.Version()$major, ".", R.Version()$minor)
if(compareVersion(R_version, R_min_version) < 0){
  stop("You do not have the latest required version of R installed (at least >=3.2).\n", 
       "Go to http://cran.r-project.org/ and update your version of R.")
}

#check for x86_64 version
CPUtype= "64"
CPUtype_current=substr(Sys.info()[['machine']], 5, 6)
if(CPUtype!=CPUtype_current){
  stop("StructurXploR is only design to run on 64 bit architecture\n")
}

#Create local library specific to structureXploR
#currentDir=dirname(sys.frame(1)$ofile)
#currentDirLib=paste(currentDir,"/localLib",sep="")
currentDirLib="localLib/"
dir.create(currentDirLib, showWarnings = F,mode="777")

#set structureXploR localLib as default library
bckpLibPaths=.libPaths()
.libPaths(currentDirLib) 

#Load/install packages
#require("tcltk")
#pb=tkProgressBar(title = "Load/install packages", label = "Load/install tcltk...",min = 0, max = 16, initial = 0, width = 500)

#setTkProgressBar(pb, 1, title = "Load/install packages", label = "Load/install shiny...")
#if (!require("shiny")){ install.packages("shiny",repo="http://cran.rstudio.com/");}

#setTkProgressBar(pb, 2, title = "Load/install packages", label = "Load/install shinydashboard...")
if (!require("shinydashboard")){ install.packages("shinydashboard",repo="http://cran.rstudio.com/");require("shinydashboard") }

#setTkProgressBar(pb, 3, title = "Load/install packages", label = "Load/install ape...")
if (!require("ape")){ install.packages("ape",repo="http://cran.rstudio.com/", lib=currentDirLib);require("ape") }

#setTkProgressBar(pb, 4, title = "Load/install packages", label = "Load/install rjson...")
if (!require("rjson")){ install.packages("rjson",repo="http://cran.rstudio.com/", lib=currentDirLib);require("rjson")}

#setTkProgressBar(pb, 4, title = "Load/install packages", label = "Load/install jsonlite...")
if (!require("jsonlite")){ install.packages("jsonlite",repo="http://cran.rstudio.com/", lib=currentDirLib);require("jsonlite") }

#setTkProgressBar(pb, 5, title = "Load/install packages", label = "Load/install pvclust...")
if (!require("pvclust")){ install.packages("pvclust",repo="http://cran.rstudio.com/", lib=currentDirLib);require("pvclust") }

#setTkProgressBar(pb, 6, title = "Load/install packages", label = "Load/install colorspace...")
if (!require("colorspace")){ install.packages("colorspace",repo="http://cran.rstudio.com/", lib=currentDirLib);}

#setTkProgressBar(pb, 7, title = "Load/install packages", label = "Load/install DT...")
if (!require("DT")){ install.packages("DT",repo="http://cran.rstudio.com/", lib=currentDirLib);require("DT")}

#setTkProgressBar(pb, 8, title = "Load/install packages", label = "Load/install cluster...")
if (!require("cluster" )){ install.packages("cluster",repo="http://cran.rstudio.com/", lib=currentDirLib); require("cluster" )}

#setTkProgressBar(pb, 9, title = "Load/install packages", label = "Load/install plyr...")
if (!require("plyr" )){ install.packages("plyr",repo="http://cran.rstudio.com/", lib=currentDirLib);require("plyr" ) }

#setTkProgressBar(pb, 10, title = "Load/install packages", label = "Load/install shinyjs...")
if (!require("shinyjs" )){ install.packages("shinyjs",repo="http://cran.rstudio.com/", lib=currentDirLib);require("shinyjs" )}

#setTkProgressBar(pb, 11, title = "Load/install packages", label = "Load/install devtools...")
if (!require("devtools")){ install.packages("devtools",repo="http://cran.rstudio.com/");}

#setTkProgressBar(pb, 12, title = "Load/install packages", label = "Load/install rCharts...")
if (!require("rCharts" )){withr::with_libpaths(new = currentDirLib, code=install_github('ramnathv/rCharts'));require("rCharts" )}

source("http://bioconductor.org/biocLite.R")
#setTkProgressBar(pb, 13, title = "Load/install packages", label = "Load/install Biostrings...")
if (!require(Biostrings ))
{ biocLite(pkgs=c("Biostrings"),suppressUpdates=T, lib=currentDirLib);require(Biostrings ) }

#setTkProgressBar(pb, 14, title = "Load/install packages", label = "Load/install BiocGenerics...")
if (!require(BiocGenerics ))
{ biocLite(pkgs=c("BiocGenerics"),suppressUpdates=T, lib=currentDirLib);require(BiocGenerics ) }

#setTkProgressBar(pb, 15, title = "Load/install packages", label = "Load/install shinyBS...")
if (!require("shinyBS" )){ install.packages("shinyBS",repo="http://cran.rstudio.com/", lib=currentDirLib);require("shinyBS" )}

#setTkProgressBar(pb, 15, title = "structureXploR", label = "Start structureXploR in your default browser...")
#Sys.sleep(1)
#close(pb)

print('Start structureXploR in your default browser...')
