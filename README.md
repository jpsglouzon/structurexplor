# Structurexplor #

*A plateform for the exploration of structural features of RNA secondary structures.*

Discovering function-related structural features, such as the cloverleaf shape of transfer RNA secondary structures, is essential to understand RNA function. With this aim, we have developed a platform, named Structurexplor, to facilitate the exploration of structural features in populations of RNA secondary structures. It has been designed and developed to help biologists interactively search for, evaluate and select interesting structural features that can potentially explain RNA functions.

## Exploring your RNA secondary structures##

Go to [Structurexplor web application](http://structurexplor.dinf.usherbrooke.ca) and start exploring your data. 

*Structurexplor web app. has been successfully tested on modern browsers: Firefox 50, Chrome 55 and Safari 10.* 

### How to run Structurexplor on a local machine or deploy on a local server (Ubuntu 16.04) ###

1. Download and install the latest version of [R (>=3.2.3)](https://www.r-project.org)
for Ubuntu 16.04 also called Ubuntu Xenial Xerus.

2. Download and unzip the [Structurexplor repository](https://github.com/jpsglouzon/structurexplor/archive/master.zip)

3. Start R command line by opening the terminal and by typing 'R'.

4. At the R prompt, install and load the following packages using the following commands 
`packList=c('shiny','rmarkdown','devtools','shinydashboard','ape','rjson','jsonlite','pvclust','colorspace','DT','cluster','plyr','shinyjs','shinyBS','rCharts','BiocGenerics','Biostrings');`
`lapply(packList, install.packages, character.only = TRUE);`
`lapply(packList, require, character.only = TRUE)`

### Supplementary step to run Structurexplor on a local machine ###
* Start Structurexplor using
`shiny::runApp("path_to_Structurxplor") `

### Supplementary steps to run Structurexplor on a local server ### 

The middleware open source [ShinyProxy](http://www.shinyproxy.io/) has been used to deploy Structurexplor. 

* To install and configure ShinyProxy please follow the detailed procedure in the ['Getting started'](http://www.shinyproxy.io/getting-started/) menu.

* The instructions to deploy a local instance of Shinyproxy is provided in the ['Deploying Apps'](http://www.shinyproxy.io/deploying-apps/) menu.

## Post suggestions or issues ##
Please report suggestions or issues [here](https://github.com/jpsglouzon/structurexplor/issues).

## How to cite StructureXploR (soon) ##


## Licence ##

Structurexplor is released under the terms of the GNU GPL licence.
For further informations, please see the LICENCE file of the repository.



 