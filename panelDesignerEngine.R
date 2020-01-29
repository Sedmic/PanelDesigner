
loadFiles <- function(){
  ssm <- read.csv(file = "SSM.csv",stringsAsFactors = FALSE);  rownames(ssm) <- ssm$SSM; ssm$SSM <- NULL
    ssm$BUV563 <- ssm$BUV563/10  #arbitrary way to get BUV563 back on scale
  fullAblist <- read.csv(file="fullAblist.csv",stringsAsFactors=FALSE) # list of marker and staining intensity
  mutEx <- read.csv(file="mutEx.csv",stringsAsFactors = FALSE)  #superficial approach to lineage-based mutual exclusivity, 1=mutex 
  rownames(mutEx) <- mutEx$X; mutEx$X <- NULL
  channels <- read.csv(file="Channels.csv",stringsAsFactors = FALSE)
  fluors <- read.csv(file="Fluors.csv",stringsAsFactors = FALSE); rownames(fluors) <- fluors$Fluorochrome
  rownames(channels) <- channels$Channel
  markerInfo <- list(ssm,fullAblist,mutEx,channels,fluors)
  names(markerInfo) <- c("ssm","fullAblist","mutEx","channels","fluors")
  print("loaded local files: ssm, fullAblist, mutEx, channels, fluors")
  return(markerInfo)
}

manualLoad <- function() {
  markerInfo <- loadFiles()
  selectedAb <- read.csv(file="selectedAblist.csv",stringsAsFactors = FALSE, strip.white = TRUE)  # matrix of desired antibodies
  # inventory <- inventoryLocal(markerInfo) 
  inventory <- inventoryCommercial(markerInfo)   # default is commercial inventory
  inventorytable <- convertInventory(inventory,markerInfo)
  numIter <- 10
}

panelDesigner <- function(numIter, selectedAb, inventorytable, markerInfo, updateProgress)
{
  library(foreach); library(doParallel)

  ssm <- markerInfo[[1]]
  mutEx <- markerInfo[[3]]
  channels <- markerInfo[[4]]; #print(channels, str(channels))
  fluors <- markerInfo[[5]]; 

  myPanel <- selectedAb
  myPanel$Priority <- sample(1:5,nrow(myPanel),replace = TRUE)
  chOptions <- unique(unlist(inventorytable)); chOptions <- chOptions[!is.na(chOptions)]
  if(length(chOptions) < nrow(selectedAb)) { print("not enough inventory options for this list"); stop() }
  
  inList <- subset(mutEx, rownames(mutEx) %in% selectedAb$Protein) #simplify mutEx to the selectedAb listing

  fluorToChannel <- assignFluorToChannel(fluors)   # assign fluor names to channels

  numOptions <- lapply(fluorToChannel, function(x) length(x))


  num_cores <- detectCores()-1  # how many cores are available on this CPU?
  cl<-makeCluster(num_cores)
  registerDoParallel(cl) # register the number of cores for foreach fnctions
  print(paste("using", num_cores,"CPU cores out of possible", num_cores+1))
  print("panelDesign started, calculation engine commencing")

  generatedPanels <- list(); penaltiesList <- list();
  z <- 1; loopCounter <- 1; loopIndex <- 0; assignIndex <- 0;
#  generatedPanels <- foreach(x=1:numIter, .packages=c("foreach"), .export="calculatePenalty_foreach") %dopar%  {
  while(length(generatedPanels) < numIter) {
    #loopIndex <- loopIndex + 1; print(c(loopIndex,assignIndex)) 
    #updateProgress(round(length(generatedPanels)/numIter,2),paste(100*round(length(generatedPanels)/numIter,2),"%"))
    
  ## **************** Assign antibodies to channels  ***************
  ## 1. check if there is an inventory entry for every antibody, if not then keep track for randomly pick
  ## 2. randomly pick one of the available colors for each antibody 
  ## 3. assign non-inventory antibodies to any color 
  ##           note that in dataframe, NA should only occur at the end of a line
    myPanel$Channel <- 0
    myPanel$Fluor <- 0
    badFlag <- 0; 
    myPanel <- myPanel[sample(nrow(myPanel)),]   # randomize antibodies 
    channelsRemaining <- rownames(channels) #names(fluorToChannel)
    if (!length(channelsRemaining)) {print("zero channels?!?!?"); badFlag <- 1;}
    #browser()
    for(y in 1:nrow(myPanel))
    {
      if(any(names(inventorytable) == myPanel$Protein[y]))   # is there an inventory entry for this antibody
      {
        #browser()
        lookUpAb <- which(names(inventorytable)==myPanel$Protein[y])
        inventoryOptions <- inventorytable[[lookUpAb]]
        chOptions <- array(character(),dim = length(inventoryOptions))
        for (m in 1:length(inventoryOptions))  # figure out options for distinct channels  
        {
          result <- lapply(fluorToChannel, function(x) grep(paste("^",inventoryOptions[m],"$",sep=""), x))
          if( ! (names(unlist(result)) %in% chOptions )) #if the channel isn't already part of the channel options then add it
            { chOptions[m] <- names(unlist(result)) } else { chOptions <- chOptions[-m] }
        }
        # now assign antibodies to channels      
        #print("assign antibodies to channels")
        if(length(chOptions)>0)
          {
            channelsPossible <- which(chOptions %in% channelsRemaining)
            if (length(channelsPossible)==0) { badFlag<-1; assignIndex <- assignIndex+1; break; } # cannot place Ab because no channels left for the available colors
            myPanel$Channel[y] <- sample(chOptions[channelsPossible],1)  # pick a channel
            fluorsForCandidate <- length(fluorToChannel[[myPanel$Channel[y]]])  # how many fluors are theoretically possible
            if ( fluorsForCandidate == 1 )    # if only one choice for fluor then great
            { 
              myPanel$Fluor[y] <- fluorToChannel[[myPanel$Channel[y]]]  
            } else if (fluorsForCandidate > 1)  #if possibility of multiple fluors for a channel
            {
              recallFluors <- which(fluorToChannel[[myPanel$Channel[y]]] %in% inventoryOptions) # num fluors are in inventory
              subsetInv <- fluorToChannel[[myPanel$Channel[y]]][recallFluors]
              if ( length(recallFluors) == 1) { myPanel$Fluor[y] <- fluorToChannel[[myPanel$Channel[y]]][recallFluors] }
              if ( length(recallFluors) > 1 ) { myPanel$Fluor[y] <- sample( subsetInv, 1) }
            }
            
            if(anyDuplicated(myPanel$Channel[1:y])) { badFlag<-2 }  # flag if double assigning a channel
          }  
      } else if ( !(myPanel$Protein[y] %in% names(inventorytable)))  # not listed in inventory?
      {  
        #browser()
        myPanel$Channel[y] <- sample(channelsRemaining, 1); print(c(myPanel$Protein[y], "...choosing randomly")) 
        myPanel$Fluor[y] <- sample( fluorToChannel[[myPanel$Channel[y]]],1 )
      }
      chHolder <- channelsRemaining ;  # browser()
      channelsRemaining <- channelsRemaining[-which(channelsRemaining %in% myPanel$Channel[y],arr.ind = TRUE)]
      if(length(channelsRemaining)==0) {  next }   # chRemaining <- chHolder
    }  
    print("channels all assigned") 
    #browser()
    
    #break out of this iteration because the antibody selection is flawed
    if (badFlag == 1)       {assignIndex <- assignIndex + 1 ; print("badFlag=1"); }
    if (badFlag == 2)       {print("double assigned Ab to channel"); }    
    
    ## *************** end assignment block ************************************  
    
    #print("assigned antibodies to channels based on inventory")
    #browser()
    if(badFlag == 0)  # if everything is good then continue to compute
    {
      myPanel$FluorName <- rownames(fluors[as.numeric(myPanel$Fluor),])
      myPanel$FluorIntensity <- fluors[as.numeric(myPanel$Fluor),2] #experim derived StainIndex per fluorochrome
      
      ssmSubset <- data.frame(); namesRow<- character()
      for (i in 1:nrow(myPanel))  # now make ssmSubset based on the fluors chosen for the panel
      {
        for (j in 1:nrow(myPanel))
        {
          lookupI <- which(rownames(ssm) %in% myPanel$FluorName[i])
          lookupJ <- which(rownames(ssm) %in% myPanel$FluorName[j])
          ssmSubset[i,j] <- ssm[lookupI,lookupJ]
        }
        namesRow[i] <- rownames(ssm[lookupI,])
      }; rownames(ssmSubset) <- namesRow; colnames(ssmSubset) <- namesRow    
      
      myPanel$PrimarySignal <- round(myPanel$FluorIntensity / myPanel$StainIntensity) #prefer dim markers on bright channels
      print("about to calculate penalty score")
      #system.time(Penalty <- calculatePenalty(myPanel,selectedAb,ssmSubset, inList))
      Penalty <- calculatePenalty_foreach(myPanel,selectedAb,ssmSubset, inList)
      print("i calculated ... something")
      myPanel$Penalty <- rowSums(Penalty) #calculate sum of penalties per Ab-channel
      myPanel$Score <- 50*myPanel$PrimarySignal - myPanel$Penalty  # take penalty away from primary signal
      myPanel$CumulScore <- sum(myPanel$Score)  
      if(is.nan(myPanel$CumulScore[1]) | is.infinite(myPanel$CumulScore[1]))  ## error! infinitiy or NaN in calculation
      {
        print(myPanel)
        print(Penalty)
        print(ssmSubset)
      }
      generatedPanels[[z]] <- myPanel; penaltiesList[[z]] <- Penalty  # sequentially add candidate panel to a list with scores
      z <- z+1
      print(paste0("the cumul score is:",  myPanel$CumulScore[1]))
      #if(length(generatedPanels) < numIter) {loopCounter=1} else if(length(generatedPanels) == numIter) {loopCounter=0}
      #return(myPanel)  # FOREACH VERSION
    }
    if (badFlag > 0)
    {
      myPanel$CumulScore <- 0
      #return(myPanel)   # FOREACH VERSION
    }
  }

  #browser()
  
  assignEfficiency <- 100*round(assignIndex / loopIndex, 2)
  print(paste("The Ab-to-channel assignment was unsuccessful", assignEfficiency, "% of the time."))
  library(beepr); beep()
  
  stopCluster(cl); registerDoSEQ()
  return(generatedPanels)
}  

scoreEvaluation <- function(generatedPanels) 
{
  scoreEval <- data.frame(Index=(seq(1:length(generatedPanels))), CumulScore=0)
  for(k in 1:(length(generatedPanels)))
  {
    scoreEval[k,2] <- round(generatedPanels[[k]]$CumulScore[1]) #pull CumulScore out of list to allow ordering in next step
  }
  return(scoreEval)
}

simpleResults <- function(generatedPanels) 
{
  scoreEval <- scoreEvaluation(generatedPanels)
  bestResults <- head(scoreEval[order(scoreEval$CumulScore,decreasing=TRUE),],n=10L)
  worstResults <- tail(scoreEval[order(scoreEval$CumulScore,decreasing=TRUE),],n=3L)
  df <- rbind(bestResults, worstResults)
  df2 <- as.data.frame(generatedPanels[[df$Index[1]]])
  for (i in 2:nrow(df))
  {
    df2 <- rbind(df2,generatedPanels[[df$Index[i]]])
    df2 <- rbind(df2,array(rep(0,11)))
  }
  return(as.matrix(df2))
}  


sanityCheck <- function(selectedAb, inventorytable, markerInfo) 
{
  ## verify that selectedAb is in inventory table
  unmatched <- which(! (selectedAb[,"Protein"] %in% names(inventorytable)))
  print(c("unmatched to inventory:", selectedAb[unmatched,]))
  #print("Every protein in selectedAb matches something in the Inventory!")
  return(unmatched)
}


calculatePenalty <- function(myPanel,selectedAb,ssmSubset, inList)
{
  Penalty <- data.frame(Penaltyi=double(),Penaltyj=double()) # temp dataframe to calculate penalty scores
  hiPriority <- 0; lowPriority <- 0;
  print("calculating penalties")
  for(i in 1:nrow(myPanel)) # where i is the row in the myPanel
  {
    hiPriority <- 0; lowPriority <- 0;
    if(myPanel[i,"Priority"] <= 2) {hiPriority <- 1} else {hiPriority<-0}
    if(myPanel[i,"Priority"] >= 4) {lowPriority <- 1} else {lowPriority<-0}
    for(j in 1:nrow(myPanel))  # where j represents all of the other antibodies in myPanel
    { 
      Penalty[i,j] <- 0.5*ssmSubset[i,j] * myPanel[j,2] * myPanel[j,7] - 
        myPanel[j,7]*inList[myPanel$Protein[i],myPanel$Protein[j]] - # bonus if mutually exclusive
        lowPriority*sum(ssmSubset[i,]) -   # bonus if lowPriority on messy ch  
        hiPriority*300/(0.1+sum(ssmSubset[i,])) +  # bonus if hiPriority on good ch
        lowPriority*( myPanel[j,2] * myPanel[j,7] ) # penalize low priority & bright on strong channel
      Penalty[i,j] <- round(Penalty[i,j])
      if(is.nan(Penalty[i,j]) | is.infinite(Penalty[i,j])) {brokeni<-i; brokenj<-j; stop()}
    }
  } 
  return(Penalty)
}

calculatePenalty_foreach <- function(myPanel,selectedAb,ssmSubset, inList)
{
  Penalty <- numeric() # temp dataframe to calculate penalty scores
  penaltyList <- list()
  #browser()
  #penaltyList <- foreach(x=1:nrow(myPanel), .combine=rbind ) %do%  {  # single-threaded penalty calculation
  penaltyList <- foreach(x=1:nrow(myPanel), .combine=rbind ) %dopar%  {  # multi-core penalty calculation
    hiPriority <- 0; lowPriority <- 0;
    if(myPanel[x,"Priority"] == "High") {hiPriority <- 1} else {hiPriority<-0}
    if(myPanel[x,"Priority"] == "Low") {lowPriority <- 1} else {lowPriority<-0}
    for(j in 1:nrow(myPanel))  # where j represents all of the other antibodies in myPanel
    { 
      Penalty[j] <- 0.5*ssmSubset[x,j] * myPanel[j,"StainIntensity"] * myPanel[j,"FluorIntensity"] -
        myPanel[j,"FluorIntensity"]*inList[myPanel$Protein[x],myPanel$Protein[j]]  - # bonus if mutually exclusive
        lowPriority*sum(ssmSubset[x,]) -   # bonus if lowPriority on messy ch  
        hiPriority*300/(0.1+sum(ssmSubset[x,])) +  # bonus if hiPriority on good ch
        lowPriority*( myPanel[j,"StainIntensity"] * myPanel[j,"FluorIntensity"] )
      Penalty[j] <- round(Penalty[j])
    }
    print(Penalty)
    return(Penalty)
  } 
  return(penaltyList)
}

