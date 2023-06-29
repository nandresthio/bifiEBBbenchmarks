library(plyr)
library(stringr)

# Functions which take the experimental results and process them so that they can be plotted
# and analysed

# Basic functions which combines all of the experimental results given into a single
# dataframe. No processing is performed, only a storing of the results.
combineArrayResults <- function(runName, arrayStart, arrayEnd, jobsPerArray = 1, printInfo = TRUE){
  options(scipen=999)
  if(jobsPerArray == 1){combinedData <- read.table(paste0("data/clusterResults/", runName, "_arrayJob", arrayStart, ".txt"), header = TRUE, sep = " ", fill = TRUE)}
  else{combinedData <- read.table(paste0("data/clusterResults/", runName, "_arrayJob", (arrayStart - 1) * jobsPerArray + 1, "-", arrayStart * jobsPerArray, ".txt"), header = TRUE, sep = " ", fill = TRUE)}
  if(arrayEnd == arrayStart){return(combinedData)}
  for(i in (arrayStart + 1):arrayEnd){
    if(printInfo){cat(paste0("\rCombining files, done ", i, "/", arrayEnd))}
    if(jobsPerArray == 1){filename <- paste0("data/clusterResults/", runName, "_arrayJob", i, ".txt")}
    else{filename <- paste0("data/clusterResults/", runName, "_arrayJob", (i - 1) * jobsPerArray + 1, "-", i * jobsPerArray, ".txt")}
    
    if(!file.exists(filename)){
      print(paste0("Skipping file ", filename, " as it does not exist!"))
      next
    }
    newData <- read.table(filename, header = TRUE, sep = " ", fill = TRUE)
    if(nrow(newData) != 50){
      print(paste0("Weird, for this file have ", nrow(newData), " rows instead of 50!"))
      print(filename)
    }
    combinedData <- rbind.fill(combinedData, newData)
  }
  if(printInfo){cat(" - done.\n")}
  options(scipen=0)
  return(combinedData)
}


