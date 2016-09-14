#check for R version
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
  stop("StructurXploR is only design to run on 64 bit architecture.\n")
}

#detect plateform
platform=Sys.info()[['sysname']]
if(platform == "Linux") {
  typeForInstallPackages='source'
} else if (platform== "Windows"){
  typeForInstallPackages='binary'
} else if (platform== "Darwin"){
  typeForInstallPackages='binary'
} else
{ typeForInstallPackages='both'}

#Create local library specific to structureXploR
currentDir=dirname(sys.frame(1)$ofile)
setwd(currentDir)
currentDirLib=paste(currentDir,"/structurexplor_Lib",sep="")
dir.create(currentDirLib, showWarnings = F,mode="755")

#set structureXploR localLib as default library
bckpLibPaths=.libPaths()
.libPaths(currentDirLib) 

print("Loading libraries ...")

if (!require("devtools")){ install.packages("devtools",repo="http://cran.rstudio.com/"); require(devtools)}

if (!require("shiny")){ install.packages("shiny",repo="http://cran.rstudio.com/",type=typeForInstallPackages); require(shiny)}

if(packageVersion("shiny")<0.14){detach('package:shiny'); install.packages("shiny",repo="http://cran.rstudio.com/",type=typeForInstallPackages); require(shiny)} 

if (!require("shinydashboard")){ install.packages("shinydashboard",repo="http://cran.rstudio.com/"); require(shinydashboard)}
if (!require("ape")){ install.packages("ape",repo="http://cran.rstudio.com/"); require(ape)}
if (!require("rjson")){ install.packages("rjson",repo="http://cran.rstudio.com/"); require(rjson)}
if (!require("jsonlite")){ install.packages("jsonlite",repo="http://cran.rstudio.com/"); require(jsonlite)}
if (!require("pvclust")){ install.packages("pvclust",repo="http://cran.rstudio.com/"); require(pvclust)}
if (!require("colorspace")){ install.packages("colorspace",repo="http://cran.rstudio.com/"); require(colorspace )}
if (!require("DT")){ install.packages("DT",repo="http://cran.rstudio.com/"); require(DT )}
if(packageVersion("DT")<0.2){detach('package:DT'); install.packages("DT",repo="http://cran.rstudio.com/"); require(DT )} 

if (!require("cluster" )){ install.packages("cluster",repo="http://cran.rstudio.com/"); require(cluster )}
if (!require("plyr" )){ install.packages("plyr",repo="http://cran.rstudio.com/"); require(plyr )}
if (!require("shinyjs" )){ install.packages("shinyjs",repo="http://cran.rstudio.com/"); require(shinyjs )}
if (!require("rCharts" )){withr::with_libpaths(new = currentDirLib,code=install_github('ramnathv/rCharts')); require(rCharts )}
if (!require("shinyBS")){ install.packages("shinyBS",repo="http://cran.rstudio.com/"); require(shinyBS)}
source("http://bioconductor.org/biocLite.R")
if (!require(Biostrings ))
{biocLite(pkgs=c("Biostrings"),suppressUpdates=T); require(Biostrings)}
if (!require(BiocGenerics ))
{biocLite(pkgs=c("BiocGenerics"),suppressUpdates=T); require(BiocGenerics)}

source("www/Functions/mainFunctions.R")

#Run structureXploR app
print("Starting structureXploR in your default browser...")
runApp(currentDir)


