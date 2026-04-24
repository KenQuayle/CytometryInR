# Problem 1
  #Import data frame
  library(dplyr)
  fileLocation <- file.path("data")
  fileToImport <- list.files(path = fileLocation, pattern="Batter_Stats.csv", full.names=TRUE)
  CSV <- read.csv(file=fileToImport, check.names=FALSE)
  CSV <- CSV |> mutate("diff" = wOBA-xwOBA)

  # Define new function
  #' This function will return a text output indicating if a baseball player has overperformed 
  #'  or underperformed expectations based on xwOBA (expected weighted on-base average) 
  #' @param data A data frame object containing columns called "Name", "xwOBA", and "diff"
  #'  containing player name, xwOBA and wOBA-xwOBA respectively 
  #' @param player A character string with the name of a player that appears within data$Name
  #' 
  doTheMarinersSuck <- function(data, player){

    Name <- as.vector(data$Name)
    xwOBA <- as.vector(data$xwOBA)
    diff <- as.vector(data$diff)

    library(stringr)
    playerIndex <- which(str_detect(Name, player))
    firstMessage <- if(xwOBA[playerIndex] >= 0.340){paste0(player, " has been great")} 
                    else if(xwOBA[playerIndex] <= 0.300){paste0(player, " has sucked")} 
                    else{paste0(player, " has been average")}
    secondMessage <- if(diff[playerIndex] >= 0.01){paste0(" and has overperformed")} 
                     else if(diff[playerIndex] <= -0.01){paste0(" and has underperformed")}
                     else{paste0(" and has performed as expected")}
    
    return(paste0(firstMessage, secondMessage))
  }

  # Iterate output over the contents of CSV$Name
  library(purrr)
  Name <- as.vector(CSV$Name)
  map(.x=Name, .f=doTheMarinersSuck, data=CSV)

# Problem 2
  # Identify FCS files and create gating set
  fileLocation <- file.path("data")
  FCSFiles <- list.files(fileLocation, pattern=".fcs", full.names=TRUE)
  library(flowWorkspace)
  myCytoSet <- load_cytoset_from_fcs(files = FCSFiles, tranformation = FALSE, truncate_max_range = FALSE)
  myGatingSet <- GatingSet(myCytoSet)

  # Transform peak detectors for Pac Blue, PE, and APC
  peakDetectors <- c("V3-A", "B4-A", "R1-A")
  Biexponential <- flowjo_biexp_trans(channelRange=4096, maxValue=262144, pos=4.5, neg=0, widthBasis=-1000)
  myBiexTransform <- transformerList(peakDetectors, Biexponential)
  transform(myGatingSet, myBiexTransform)

  # Add gates
  library(openCyto)
  gateTemplate <- list.files(fileLocation, pattern="Gates.csv", full.names=TRUE)
  gatesToApply <- data.table::fread(gateTemplate)
  gt_gating(gatingTemplate(gatesToApply), myGatingSet)

  # Derive cell concentration
  #' This function will output cell concentration in cells per mL for a given gated population in an FCS file.
  #'
  #' @param x A GatingSet object
  #' @param subset The gate from which to retrieve cell counts
  #' @param dilutionFactor Amount by which to multiply the cell count to account for dilution of the sample
  #' 
  CellConcentration <- function(x, subset, dilutionFactor){
    Serial <- keyword(x)$'$CYTSN'
    FSCgain <- keyword(x)$'$P36V'
    TotalFileEvents <- keyword(x)$'$TOT'
    Volume <- keyword(x)$'$VOL'
    Volume <- as.numeric(Volume)
    targetSubset <- gs_pop_get_data(x, subset)
    Cells <- nrow(targetSubset)[[1]]
    Concentration <- (1000*Cells)/Volume
    Concentration <- round(Concentration*dilutionFactor, 0)
    Data <- data.frame(InstrumentSN=Serial, FSCgain=FSCgain, Volume=Volume, Cells=Cells, Concentration=Concentration)
    return(Data)
  }

  # Generate data frame with cell concentration and additional keywords
  Dataset <- map(.x=myGatingSet, .f=CellConcentration, subset="singlets", dilutionFactor=100) |> bind_rows()

# Problem 3
  # Identify FCS files and create gating set
  fileLocation <- file.path("data")
  FCSFiles <- list.files(fileLocation, pattern=".fcs", full.names=TRUE)
  library(flowWorkspace)
  myCytoSet <- load_cytoset_from_fcs(files = FCSFiles, tranformation = FALSE, truncate_max_range = FALSE)
  myGatingSet <- GatingSet(myCytoSet)

  # Derive flow rate in uL/min
  #' This function will output flow rate in uL/min for a FCS file or set of FCS files
  #' @param x A gating set object
  #' 
  flowRate <- function(x){
    library(lubridate)
    Volume <- as.numeric(keyword(x)$'$VOL')
    Start <- hms(keyword(x)$'$BTIM')
    End <- hms(keyword(x)$'$ETIM')
    Time <- as.numeric(as.duration(End-Start))
    Rate <- round(60*Volume/Time, 2)
    return(paste0(Rate, " uL/min"))
  }

  library(purrr)
  map(.x=myGatingSet, .f=flowRate)

  # The variation in flow rate appears to be relatively minor, all files are consistent with being recorded on 'high' flow rate on an Aurora.