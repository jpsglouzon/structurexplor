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

print("Load/install libraries. Please wait ...")

#install/start packrat
if (!require("packrat")){ install.packages("packrat",repo="http://cran.rstudio.com/");require("packrat") };
packrat::on();

#Restore all lib if shinydashboard does not exist
if (!require("shinydashboard")){ packrat::restore();};

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
library(BiocGenerics)
library(Biostrings)
library(shinyBS)

currentDir=dirname(sys.frame(1)$ofile)
print("Starting structureXploR in your defautl browser ...")
runApp(currentDir)
