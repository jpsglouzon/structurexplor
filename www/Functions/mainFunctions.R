
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#####################################################
computeStructuralPatterns<-function(pathfile,resParseRNA,snm,max_n_motifs,setnmotifs,maxClust,HC,setnbcluster)
{
  patterns=NULL
  withProgress(message = 'Compute structural features', value = 0, { 
    allstep=9
    #incProgress(1/allstep, message = "Parse .db file...")   
    #resParseRNA <- parseDbFile(pathfile)
    incProgress(2/allstep, message = "Run the super-n-motifs model...")
    resComp <- runSNM (pathfile,snm,max_n_motifs) 
    resClus <- runClustering (resComp$SuperMotif,resComp$dissimSS,resComp$matDissim_SSbynm,resComp$SuperMotifnmotifs,resComp$matnmPos,setnmotifs,maxClust,HC,setnbcluster) 
    incProgress(4/allstep, message = "Generate datatable of structures...")
    #datass<-c()
    datass<-dataTableSS(resParseRNA$headers,resClus$idxClustering,resClus$repClustList,resParseRNA$sequences,resParseRNA$structures,resParseRNA$lengthRNAs,resClus$outliers)
    incProgress(5/allstep, message = "Generate datatable of Patterns...")
    dataPat=dataTablePattern(resClus$unikidxClustering,resClus$clustSize,resClus$confStabilit,resClus$repClustList,resClus$outliers,resClus$bestNmotifsForClusters_pos_namesNmotifs)
    incProgress(6/allstep, message = "Generate pattern visualization...")
    scatPlot<-scatPlotSNM(resComp$SuperMotif[,1],resComp$SuperMotif[,2],resClus$bestNmotifsForClusters_coord[,1],resClus$bestNmotifsForClusters_coord[,2],1,2,resClus$bestNmotifsForClusters_pos,resClus$bestNmotifsForClusters_pos_namesNmotifs,resComp$singularValuesPercent,resParseRNA$headers,resClus$idxClustering,resClus$colorClusters,resClus$repClustList,resClus$outliers,resClus$matnmPosBestnmotifsforForna,resParseRNA)###
    incProgress(7/allstep, message = "Generate explained variability plot...")
    variabilityExplained<-variability(resComp$singularValuesPercent,resComp$SuperMotif)
    incProgress(8/allstep, message = "Generate cluster hierarchy...")
    dendro1<-c()
    incProgress(9/allstep, message = "Generate info. on cluster quality...")
    dataclustInfo=dataTableQualClusters(resClus$unikidxClustering,resClus$clustSize,resClus$listSilhouetteCoefPerClusters)
    
    patterns=list(resParseRNA,resComp,resClus,datass,scatPlot,variabilityExplained,dendro1,dataclustInfo,dataPat)
    
    return(patterns)
  })
}

#####################################################
updateStructuralPatterns_with_setnbcluster<-function(resParseRNA,resComp,setnmotifs,maxClust,HC,setnbcluster)
{
  patterns=NULL
  withProgress(message = 'Update structural patterns', value = 0, { 
    allstep=7
    incProgress(1/allstep, message = "Run the clustering...")
    resClus <- runClustering (resComp$SuperMotif,resComp$dissimSS,resComp$matDissim_SSbynm,resComp$SuperMotifnmotifs,resComp$matnmPos,setnmotifs,maxClust,HC,setnbcluster) 
    incProgress(2/allstep, message = "Generate datatable of structures...")
    datass<-dataTableSS(resParseRNA$headers,resClus$idxClustering,resClus$repClustList,resParseRNA$sequences,resParseRNA$structures,resParseRNA$lengthRNAs,resClus$outliers)
    incProgress(3/allstep, message = "Generate datatable of Patterns...")
    dataPat=dataTablePattern(resClus$unikidxClustering,resClus$clustSize,resClus$confStabilit,resClus$repClustList,resClus$outliers,resClus$bestNmotifsForClusters_pos_namesNmotifs)
    incProgress(4/allstep, message = "Generate pattern visualization...")
    scatPlot<-scatPlotSNM(resComp$SuperMotif[,1],resComp$SuperMotif[,2],resClus$bestNmotifsForClusters_coord[,1],resClus$bestNmotifsForClusters_coord[,2],1,2,resClus$bestNmotifsForClusters_pos,resClus$bestNmotifsForClusters_pos_namesNmotifs,resComp$singularValuesPercent,resParseRNA$headers,resClus$idxClustering,resClus$colorClusters,resClus$repClustList,resClus$outliers,resClus$matnmPosBestnmotifsforForna,resParseRNA)###
    incProgress(5/allstep, message = "Generate explained variability plot...")
    variabilityExplained<-variability(resComp$singularValuesPercent,resComp$SuperMotif)
    incProgress(6/allstep, message = "Generate cluster hierarchy...")
    dendro1<-c()
    incProgress(7/allstep, message = "Generate info. on cluster quality...")
    dataclustInfo=dataTableQualClusters(resClus$unikidxClustering,resClus$clustSize,resClus$listSilhouetteCoefPerClusters)
    
    patterns=list(resParseRNA,resComp,resClus,datass,scatPlot,variabilityExplained,dendro1,dataclustInfo,dataPat)
    
    return(patterns)
  })
}

#####################################################
updateStructuralPatterns_with_setnmotifs<-function(patterns,resParseRNA,resComp,setnmotifs,maxClust,HC,setnbcluster)
{
  withProgress(message = 'Update structural patterns', value = 0, { 
    allstep=7

    incProgress(1/allstep, message = "Run the clustering...")
 
    resClus <- runClustering (resComp$SuperMotif,resComp$dissimSS,resComp$matDissim_SSbynm,resComp$SuperMotifnmotifs,resComp$matnmPos,setnmotifs,maxClust,HC,setnbcluster) 
    incProgress(2/allstep, message = "Generate datatable of structures...")
    datass<-dataTableSS(resParseRNA$headers,resClus$idxClustering,resClus$repClustList,resParseRNA$sequences,resParseRNA$structures,resParseRNA$lengthRNAs,resClus$outliers)
    incProgress(3/allstep, message = "Generate datatable of Patterns...")
    dataPat=dataTablePattern(resClus$unikidxClustering,resClus$clustSize,resClus$confStabilit,resClus$repClustList,resClus$outliers,resClus$bestNmotifsForClusters_pos_namesNmotifs)

    patterns=list(resParseRNA,resComp,resClus,datass,patterns[[5]],patterns[[6]],patterns[[7]],patterns[[8]],dataPat)
    
    return(patterns)
  })
}

#################### read and parse vienna########################
parseDbFile <- function (pathSSdnb) {
  
  fastaDataraw<-readBStringSet(pathSSdnb,"fasta")
  header=names(fastaDataraw)
  seqNstruct<-lapply(fastaDataraw, splitRNA2SeqNStruct)
  resultsParseRNA=c()
  resultsParseRNA$headers=header
  resultsParseRNA$sequences=unlist(lapply(seqNstruct, `[[`, 1))
  resultsParseRNA$structures=unlist(lapply(seqNstruct, `[[`, 2))
  resultsParseRNA$lengthRNAs=unlist(lapply(seqNstruct, `[[`, 3))
  
  resultsParseRNA=as.data.frame(resultsParseRNA)
  
  return(resultsParseRNA)
}

#####################################################
splitRNA2SeqNStruct<-function (x, ...){
  
  x=toString(x)
  seqTemp=subseq(x,start=1,end=width(x)/2)
  structTemp=subseq(x,start=(width(x)/2)+1,end=width(x))
  seqNstruct=list("sequence"=seqTemp[[1]],"structure"=structTemp[[1]],"length"=nchar(seqTemp[[1]]))
  
  return(seqNstruct)
} 

######################## 1.run supermotifs model #########################
runSNM <- function (pathSSdnb,nbSnm,maxNm) {

  platform=Sys.info()[['sysname']]
  if(platform == "Linux") {
    pathSNM="www/Functions/supernmotifs_ubuntu64_V1.3"
    a=system2("chmod", args=c("755",pathSNM))
  } else{
    print ("Unrecognized operating system. Structurexplor is supported by Ubuntu 64.")
    stopApp()
  }

  outputPath="www/Data/"
  
  a=system2(pathSNM, args=c("-i",pathSSdnb, "-o",outputPath,"-p","2","-s",nbSnm))
 
  #read super-n-motifs representation of SS
  pathMatSuperMotif<-paste(outputPath,"matSnmRep_SSbySnm.csv",sep="")
  SuperMotif <- read.table(pathMatSuperMotif,header = TRUE,row.names=1,sep = ",",comment.char = "")
  
  #dissim matrix
  pathdissimSS<-paste(outputPath,"matDissim_SSbySS.csv",sep="")
   dissimSS <- data.matrix(read.table(pathdissimSS, fill=T, sep = ",",  col.names=1:nrow(SuperMotif))) 
   dissimSS[is.na(dissimSS)]=0     
   dissimSS=rbind(rep(0,nrow(SuperMotif)),dissimSS)     
   dissimSS[upper.tri(dissimSS)] <- t(dissimSS)[upper.tri(dissimSS)]
   colnames(dissimSS)=row.names(SuperMotif)
   rownames(dissimSS)=row.names(SuperMotif)

  #read singular values
  pathSingularValues<-paste(outputPath,"singularValuesFull_supernmotifs.csv",sep="")
  singularValues <- read.csv(pathSingularValues,header = TRUE,sep = ",",comment.char = "")
  singularValues=singularValues[,-1]
  singularValuesPercent=(singularValues/sum(singularValues))*100  
  
  #read dissimilarity matrix SS*nm
  pathMatDissim_SSbynm<-paste(outputPath,"matDissim_SSbynm.csv",sep="")
  matDissim_SSbynmRaw <- read.csv(pathMatDissim_SSbynm,header = TRUE,sep = ",",comment.char = "")
  matDissim_SSbynm=matDissim_SSbynmRaw[,-1]
  row.names(matDissim_SSbynmRaw)=matDissim_SSbynmRaw[,1]
  matDissim_SSbynmRaw=matDissim_SSbynmRaw[,-1]

  #read super-n-motifs rep of all n-motifs
  pathMatSuperMotifnmotifs<-paste(outputPath,"matSnmRep_nmbySnm.csv",sep="")
  SuperMotifnmotifs <- read.table(pathMatSuperMotifnmotifs,header = TRUE,row.names=1,sep = ",",comment.char = "")  

  #read matnpos: all nmotifs
  pathMatnmPos<-paste(outputPath,"matnmPos.csv",sep="")
  matnmPos <- read.csv(pathMatnmPos,header = TRUE,sep = ",",check.names=FALSE,stringsAsFactors=FALSE,comment.char = "")
  colnames(matnmPos)=colnames(matnmPos)[-1]
  matnmPos=matnmPos[,-dim(matnmPos)[2]]

  results=c()
  results$SuperMotif=SuperMotif
  results$dissimSS=dissimSS
  results$singularValuesPercent=singularValuesPercent
  results$matDissim_SSbynm=matDissim_SSbynmRaw
  results$SuperMotifnmotifs=SuperMotifnmotifs
  results$matnmPos=matnmPos
  
  file.remove(pathMatSuperMotif,pathdissimSS,pathSingularValues,pathMatDissim_SSbynm,pathMatSuperMotifnmotifs,pathMatnmPos)
  
  return(results)
}

###################### 2. Clustering of secondary structures #########################

runClustering <- function (SuperMotifRaw,dissimSS,matDissim_SSbynm,SuperMotifnmotifs,matnmPos,setnmotifs,maxClust,methodHC,setnbcluster){
  
  if (methodHC==0) {methodHC="average"}
  if (methodHC==1) {methodHC="ward"}
  if (methodHC==2) {methodHC="single"}
  if (methodHC==3) {methodHC="complete"}
  
  maxClust=as.numeric(maxClust)
  
    SuperMotifRaw=t(SuperMotifRaw)
    SuperMotif=SuperMotifRaw
    
    if (maxClust>dim(SuperMotif)[2])
    {
      print("Maximum number of clusters (maxClust) is greater than the nb. of SS. maxClust is set so that maxClust = n.b of SS. -1")
      maxClust=dim(SuperMotif)[2]-1
    }
  
  ####### Determination of the number of clusters using silhouette coefficient #########

  result_1=hclust(as.dist(dissimSS), method=methodHC)
    
  silAvgWidth_allClustering=c()
  if(setnbcluster==0)
    {
      for (j in 2:maxClust) 
      {
        silAvgWidth_currentClustering=summary(silhouette(cutree(result_1,k=j),dissimSS))$avg.width
        silAvgWidth_allClustering<-c(silAvgWidth_allClustering,silAvgWidth_currentClustering)
      }
      #Silhouette of the best clustering
      avgSilhouetteCoef=max(silAvgWidth_allClustering) 
      #Best number of clusters : bestNbClusters
      bestNbClusters=which.max(silAvgWidth_allClustering)+1 
  }
  else
  {
    avgSilhouetteCoef=summary(silhouette(cutree(result_1,k=setnbcluster),dissimSS))$avg.width

    #Best number of clusters : bestNbClusters
    bestNbClusters=setnbcluster
  }
  
  #Optimal clustering based on the best number of clusters : idxClustering
  idxClustering=cutree(result_1,k=bestNbClusters) 
  
  #computed Silhouette
  computedSil=silhouette(idxClustering,dissimSS)
  
  #Silhouette coef of cluster of the best clustering
  bestSilouetteCoefClust=summary(silhouette(idxClustering,dissimSS))
  
  listSilhouetteCoefPerClusters=bestSilouetteCoefClust$clus.avg.widths
  clustSize=bestSilouetteCoefClust$clus.sizes
  
  #set colors of clusters
  colorClusters=gg_color_hue(as.numeric(bestNbClusters))
  #colorClusters=substr(colorClusters, 1, nchar(colorClusters)-2)
  #colorClusters=rainbow_hcl(as.numeric(bestNbClusters))

  #compute cluster representatives: repClustList
  unikidxClustering=unique(idxClustering);
  repClustList<-c()
  confStabilit<-c()
  outliers<-list()
  silDistribPerCluster<-c()
  lengthDistribPerCluster<-c()
  
  
    bestNmotifsForClusters_coord<-c()
    bestNmotifsForClusters_pos<-c()
    matnmPos_names=colnames(matnmPos)
    bestNmotifsForClusters_pos_namesNmotifs<-c()
    
    bestNmotifsForClusters_idx<-c()
    matnmPosBestnmotifsforForna<-c()
    
  length_unikidxClustering=length(unikidxClustering)
  for (i in 1:length_unikidxClustering) 
  {
    
    clusCurrentInd=which(idxClustering %in% unikidxClustering[i])#dissimilarity matrix based on cosine similarity -->test with output of snm (faster)
    names(clusCurrentInd)<-names(idxClustering[clusCurrentInd])

    if (length(clusCurrentInd)>1)
    {
      
      indRep=which(names(idxClustering) %in% names(which.min(colMeans(dissimSS[clusCurrentInd,clusCurrentInd]))));
      currentRep=rbind(idxClustering[indRep],indRep)
      row.names(currentRep)=c("idxClustering","indRep")
      #repClustList<-cbind(repClustList,currentRep)  
      
      silDistribPerCluster<-colMeans(dissimSS[clusCurrentInd,clusCurrentInd])

      temp_boxplotRes=boxplot(silDistribPerCluster, plot = FALSE)
      
      #Variability
      
      temp=dissimSS[clusCurrentInd,clusCurrentInd]
      temp[lower.tri(temp,T)] <- NA
      distancesOfCurrentClust <- temp[!is.na(temp)]
      confStabilit=c(confStabilit,mean(distancesOfCurrentClust))
      #confStabilit=c(confStabilit,temp_boxplotRes$stats[4,1]-temp_boxplotRes$stats[2,1])
      
      #Outliers detection: two methods
      
      outliersTemp=temp_boxplotRes$out[temp_boxplotRes$out >temp_boxplotRes$stat[5,1]]
      
      if (is.null(outliersTemp))
      {outliers[[i]]=NA}
      else
      {
        outliers[[i]]=outliersTemp
      }

    }else
    {
      indRep=clusCurrentInd
      currentRep=rbind(idxClustering[clusCurrentInd],indRep)
      row.names(currentRep)=c("idxClustering","indRep")
      #repClustList<-cbind(repClustList,currentRep)  
      confStabilit=c(confStabilit,NA)
      outliers[[i]]=NA
    }
    repClustList<-cbind(repClustList,currentRep)
    
    #Select n-motif representatives of clusters: n-motifs closest  to the representatives
    matDissim_SSbynm_clusterCurrent=matDissim_SSbynm[idxClustering==idxClustering[indRep],]
   
    
    matDissim_SSbynm_MeansnmotifsBycluster=colMeans(matDissim_SSbynm_clusterCurrent)
    bestNmotifsForClusters_ordered=order(matDissim_SSbynm_MeansnmotifsBycluster)
    
    #print(setnmotifs)
      for (j in 1:setnmotifs) 
        {  
        bestNmotifsForClusters_idxTemp=bestNmotifsForClusters_ordered[j]
        bestNmotifsForClusters_idx=rbind(bestNmotifsForClusters_idx,bestNmotifsForClusters_idxTemp)
        
        bestNmotifsForClusters_coordTemp=SuperMotifnmotifs[bestNmotifsForClusters_idxTemp,]
        bestNmotifsForClusters_coord=rbind(bestNmotifsForClusters_coord,bestNmotifsForClusters_coordTemp)
    
        bestNmotifsForClusters_posTemp=matnmPos[,bestNmotifsForClusters_idxTemp]
        bestNmotifsForClusters_pos=rbind(bestNmotifsForClusters_pos,bestNmotifsForClusters_posTemp)
    
        bestNmotifsForClusters_pos_namesNmotifsTemp=rbind(matnmPos_names[bestNmotifsForClusters_idxTemp],rep(idxClustering[indRep],length(matnmPos_names[bestNmotifsForClusters_idxTemp])))
        bestNmotifsForClusters_pos_namesNmotifs=cbind(bestNmotifsForClusters_pos_namesNmotifs,bestNmotifsForClusters_pos_namesNmotifsTemp)
        
        matnmPosBestnmotifsforFornaTemp=matnmPos[,bestNmotifsForClusters_idxTemp]
        matnmPosBestnmotifsforFornaTemp[idxClustering!=idxClustering[indRep]]="x"

        matnmPosBestnmotifsforFornaTemp=gsub('x',"",matnmPosBestnmotifsforFornaTemp)
        matnmPosBestnmotifsforFornaTemp=gsub('NA',"",matnmPosBestnmotifsforFornaTemp)
        matnmPosBestnmotifsforFornaTemp=gsub("\\|",":#36C3E0 ",matnmPosBestnmotifsforFornaTemp)
        
        matnmPosBestnmotifsforForna=paste(matnmPosBestnmotifsforFornaTemp,matnmPosBestnmotifsforForna,sep="")
      }
    
  }
  
  resultsClustering=c()
  resultsClustering$matDissim=dissimSS
  resultsClustering$repClustList=repClustList
  resultsClustering$idxClustering=idxClustering
  resultsClustering$unikidxClustering=unikidxClustering
  
  resultsClustering$bestNbClusters=bestNbClusters

  resultsClustering$resultClust=result_1
  resultsClustering$colorClusters=colorClusters
  #
  resultsClustering$avgSilhouetteCoef=avgSilhouetteCoef
  resultsClustering$listSilhouetteCoefPerClusters=listSilhouetteCoefPerClusters
  resultsClustering$clustSize=clustSize
  resultsClustering$colorClusters=colorClusters
  
  #
  resultsClustering$confStabilit=confStabilit
  resultsClustering$outliers=outliers
  
  #nmotifs
  resultsClustering$bestNmotifsForClusters_coord=bestNmotifsForClusters_coord
  resultsClustering$bestNmotifsForClusters_pos=bestNmotifsForClusters_pos
  resultsClustering$bestNmotifsForClusters_pos_namesNmotifs=bestNmotifsForClusters_pos_namesNmotifs
  resultsClustering$matnmPosBestnmotifsforForna=matnmPosBestnmotifsforForna
  
  return(resultsClustering)
}  

#################################### . dataTable of QualClusters. ######################################## 

dataTableQualClusters<-function(unikidxClustering,clustSize,listSilhouetteCoefPerClusters) {
 
  
  datass<-cbind(as.character(unikidxClustering),as.character(signif(listSilhouetteCoefPerClusters,3)),as.character(clustSize))
  datass=as.data.frame(datass)
  colnames(datass)<-c("Cluster id","Silhouette coef. (avg.)","Size")
  
  datassDatatable=datatable(datass,style = 'bootstrap',filter = 'top'
                            , rownames = FALSE,
                              extensions = 'Buttons', options = list(
                                searchHighlight = TRUE, search = list(regex = TRUE),
                                dom = 'Bfrtip',
                                buttons = list(list(extend='copy',filename='TableOfClusterSilhouetteCoefAndsize'), list(extend='csv',filename='TableOfClusterSilhouetteCoefAndsize'),list(extend='excel',filename='TableOfClusterSilhouetteCoefAndsize'))
                              )
  ) 
  return(datassDatatable) 
}

#################################### . dataTable of ss ######################################## 

dataTableSS<-function(headers,idxClustering,repClustList,sequences,structures,lengthRNAs,outliers) {
  
  lengthRNAs=lengthRNAs+0
  outliersList=rep("x",length(headers))
  
  for (i in 1:length(outliers)) 
  {
    indTemp<-c()
    indTemp=which(headers %in% names(outliers[[i]]))
    outliersList[indTemp]="Unusual structures"
  }
  
  repList=rep("x",length(headers))
  for (i in 1:dim(repClustList)[2]) 
  {repList[repClustList[2,i]]="Representative"}
  
  idxClustering=paste("cluster ",idxClustering,sep="")
  
  lengthRNAsString=toString(lengthRNAs)
  
  datass=cbind.data.frame(headers,idxClustering,repList,as.character(lengthRNAs),outliersList,sequences,structures)
  colnames(datass)<-c("Header","Cluster","Representatives","Length ","Unusual structures","Sequences","Structures")
  
  return(datass) 
}

#################################### . dataTable of patterns ######################################## 

dataTablePattern<-function(unikidxClustering,clustSize,confStabilit,repClustList,outliers,bestNmotifsForClusters_pos_namesNmotifs) {
  
  outlierStr<-c()
  for (i in 1:length(outliers)) 
  {outlierStr<-c(outlierStr,paste(names(outliers[[i]]), collapse = ', '))}
  
  RelConfVar=signif(confStabilit,digits=3)
  RelConfVarRank=rank(-RelConfVar)
  RelConfVar_RelConfVarRank=paste(RelConfVarRank," (",RelConfVar,")",sep="")
  
  bestNmotifsForClusters=c()
  for (i in unikidxClustering) 
  {
    bestNmotifsForCurrentClusters=bestNmotifsForClusters_pos_namesNmotifs[1,bestNmotifsForClusters_pos_namesNmotifs[2,]==i]
    bestNmotifsForCurrentClusters=paste(bestNmotifsForCurrentClusters, collapse = ', ')
    bestNmotifsForClusters=rbind(bestNmotifsForClusters,bestNmotifsForCurrentClusters)
  }

  datass<-cbind(unikidxClustering,clustSize,RelConfVar_RelConfVarRank, colnames(repClustList),outlierStr,bestNmotifsForClusters)
  datass=as.data.frame(datass)
  colnames(datass)<-c("Id cluster","Size","Struc. var. rank (value)","Representatives","Unusual structures","Top rep. regions")

    datassDatatable=datatable(datass,style = 'bootstrap',filter = 'top',
                            rownames = FALSE,
                            extensions = 'Buttons', options = list(
                              searchHighlight = TRUE, search = list(regex = TRUE),
                              dom = 'Bfrtip',
                              buttons = list(list(extend='copy',filename='TableOfClusterFeatures'), list(extend='csv',filename='TableOfClusterFeatures'),list(extend='excel',filename='TableOfClusterFeatures'))
                            )
                            
                      ) 

  return(datassDatatable) 
}


##################################### 3.Interactive Scatter Plot#########################################

scatPlotSNM<-function(x,y,xnmotifs,ynmotifs,indx,indy,bestNmotifsForClusters_pos,bestNmotifsForClusters_pos_namesNmotifs,singularValuesPercent,headers,idxClustering,colorClusters,repClustList,outliers,matnmPosBestnmotifsforForna,resultsParseRNA){
  
  coordXY=as.data.frame(cbind(x,y))
  
  Data_x_y_family_header=as.data.frame(cbind(coordXY,idxClustering,headers))
  Data_x_y_family_header$sequences=resultsParseRNA$sequences
  Data_x_y_family_header$structures=resultsParseRNA$structures
  Data_x_y_family_header$lengths=resultsParseRNA$lengthRNAs
  Data_x_y_family_header$topnmotifs=matnmPosBestnmotifsforForna
  
  coordCenters<-c()
  
  ##Scatter plot
  a <- rCharts:::Highcharts$new()
  a$chart(type = "scatter",zoomType="xy")
  unikidxClustering=unique(idxClustering)
  for (i in 1:length(unikidxClustering)) 
  {
    temp <- Data_x_y_family_header[Data_x_y_family_header$idxClustering == unikidxClustering[i],]
    a$series(data = toJSONArray2(temp, json = F), name = paste("cluster",unikidxClustering[i]),color=colorClusters[i],
             marker = list(symbol="circle",radius=5.5))
    
    coordCenters<-rbind(coordCenters,Data_x_y_family_header[repClustList[2,i],])
  }
  
  #centers
  a$series(data = toJSONArray2(coordCenters, json = F), name ="Representatives",
           marker = list(symbol="square",fillColor= '#000000',radius=4))  

  #outliers
  coordOutliers<-NULL

  for (i in 1:length(outliers)) 
  {
    indTemp<-c()
    indTemp=which(headers %in% names(outliers[[i]]))
    if (length(indTemp)!=0)
   { coordOutliers<-rbind(coordOutliers,cbind(Data_x_y_family_header[indTemp,],Data_x_y_family_header[indTemp,]))}
  }
  if (!is.null(coordOutliers))
  { a$series(data = toJSONArray2(coordOutliers, json = F), name ="Unusual structures",
           marker = list(symbol="diamond",fillColor= '#000000',radius=3)) }

  a$xAxis(title = list(style=list(color='#000000'),
    text = paste("Super-n-motif X:",format(singularValuesPercent[indx],digits=3, nsmall=3, decimal.mark="."),"% of explained variability.",sep=" ")))
  
  a$yAxis(title = list(style=list(color='#000000'),
    text = paste("Super-n-motif Y:",format(singularValuesPercent[indy],digits=3, nsmall=3, decimal.mark="."),"% of explained variability.",sep=" ")))
  
    a$tooltip(followPointer=FALSE,shadow=FALSE,
              useHTML = T,
              formatter = "#! function() {
              
              return '<table>'
              + '<center>'+this.point.headers+'</center>'
              + '<br>'
              + '<center>'+this.point.lengths+' nt. </center>'
              +'</table>';} !#"
    )
    
    
    a$plotOptions(
      series = list(
        animation=FALSE,
        allowPointSelect=TRUE,
        cursor = "pointer", 
        point = list(
          events = list(
            
            select = "#! function(event) {


                      var chart = this.series.chart;

                      var selectedPointsStr = '';
                      if (event.accumulate) {
                                selectedPoints.push(this);
                                if(selectedPoints.length>2){selectedPoints[0].select(false);selectedPoints.shift();}
                            } else {
                                selectedPoints = [this];
                                if(selectedPoints.length>2){selectedPoints[0].select(false);selectedPoints.shift();}
                            }

                      if(selectedPoints.length==1){
                            headerSelected1=selectedPoints[selectedPoints.length-1].headers
                            sequenceSelected1=selectedPoints[selectedPoints.length-1].sequences
                            structureSelected1=selectedPoints[selectedPoints.length-1].structures
                            lengthSelected1=selectedPoints[selectedPoints.length-1].lengths
                            idxClusteringSelected1=selectedPoints[selectedPoints.length-1].idxClustering
                            topnmotifsSelected1=selectedPoints[selectedPoints.length-1].topnmotifs
                            headerSelected2=''
                            sequenceSelected2=''
                            structureSelected2=''
                            lengthSelected2=''
                            idxClusteringSelected2=''
                            topnmotifsSelected2=''

                      }
                      else if (selectedPoints.length==2){
                            headerSelected1=selectedPoints[selectedPoints.length-2].headers
                            sequenceSelected1=selectedPoints[selectedPoints.length-2].sequences
                            structureSelected1=selectedPoints[selectedPoints.length-2].structures
                            lengthSelected1=selectedPoints[selectedPoints.length-2].lengths
                            idxClusteringSelected1=selectedPoints[selectedPoints.length-2].idxClustering
                            topnmotifsSelected1=selectedPoints[selectedPoints.length-2].topnmotifs

                            headerSelected2=selectedPoints[selectedPoints.length-1].headers
                            sequenceSelected2=selectedPoints[selectedPoints.length-1].sequences
                            structureSelected2=selectedPoints[selectedPoints.length-1].structures
                            lengthSelected2=selectedPoints[selectedPoints.length-1].lengths
                            idxClusteringSelected2=selectedPoints[selectedPoints.length-1].idxClustering
                            topnmotifsSelected2=selectedPoints[selectedPoints.length-1].topnmotifs
                      }

                  Shiny.onInputChange('click', {
                        header1: headerSelected1,
                        sequence1: sequenceSelected1,
                        structure1: structureSelected1,
                        length1: lengthSelected1,
                        idxClustering1: idxClusteringSelected1,
                        topnmotifs1: topnmotifsSelected1,

                        header2: headerSelected2,
                        sequence2: sequenceSelected2,
                        structure2: structureSelected2,
                        length2: lengthSelected2,
                        idxClustering2: idxClusteringSelected2,
                        topnmotifs2: topnmotifsSelected2,
                 })

                } !#"))
      )
    )
    
    a$addParams(dom = "chart")
    a$set(width="100%",heigth='100%')
  
  return(a) 
  }

##################################### 4. plot variability#######################################

variability<-function(singularValuesPercent,Supermotifs) {
  
  singularValuesPercent_cumsum=cumsum(as.numeric(singularValuesPercent))
  
  b <- rCharts:::Highcharts$new()
  b$chart(zoomType = "x",shadow=FALSE)
  
  b$series(data=t(signif(singularValuesPercent,digits=3)),type='column',name= 'Exp. var.', marker = list(symbol="square"),tooltip=list(valueSuffix='%'))
  b$series(data=t(signif(singularValuesPercent_cumsum,digits=3)),type='spline',name= 'Cum. exp. var.', marker = list(symbol="circle"),tooltip=list(valueSuffix='%'),yAxis = 1)
  
  
  b$legend(enabled = TRUE)

  b$xAxis(

          title = list(style=list(color='#000000'),text = "Super-n-motifs "),
          
          plotBands=list( 
            label=list(text="Super-n-motifs retained",rotation= 90,textAlign= 'left'),
            color= 'lightGrey', 
            from= dim(Supermotifs)[2], 
            to= 1 
          )
          
          )
  

  b$yAxis(list(
          
          list(title = list(style=list(color='#000000'),text = "Explained variability"), replace=F, labels=list(format='{value} %')),
          list(title = list(style=list(color='#000000'),text = "Cumulative explained variability"),opposite= TRUE,max=100, labels=list(format='{value} %'))
          
          )
  )
          
          #startOnTick= FALSE)
  
  b$tooltip(shared=TRUE,followPointer=FALSE,shadow=FALSE)
  
  b$plotOptions(
    
    spline = list(
      pointStart=1,animation=FALSE
    ),
    column = list(
      pointStart=1,animation=FALSE
    )
  )

  b$set(width='100%')

  return(b)
}

