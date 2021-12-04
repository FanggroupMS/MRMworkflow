#applying the total work-flow
packsneed <- c('xcms','magrittr','MSnbase','dplyr','tidyr',
               'ggplot2','tidyverse','ggpubr','miceadds',
               'rcdk','rcdklibs','randomForest','leaps',
               'caret','corrplot','mlr','Metrics','rio','hash')
lapply(packsneed, require, character.only = TRUE)

#outside function, file import
#define wrapper function for reading the mrm files
#funciton definition: read mrm file with default experiment setting and allow noise parameters imput  #done
#function param
filedat <- list.files(path = 'dmrm/', pattern = ".mzML$", full.names = TRUE)
pwid <- c(0.8,3)
snthresh <- 10
MSRT <- read.csv("MS1_RT.csv", col.names = c('MS1','RT'))
rtfile <- data.frame("row"=1:length(rt),"ms1"=MSRT$MS1, "predrt"=MSRT$RT)
delrt <- 1.5
group1 <- c(1:3)
group2 <- c(4:5)
peaksresult <- peaksdetect(filedat, rtfile,pwid, snthr, delrt, group1,group2)

#load new data for rt prediction  #done
newinput <- read.csv('OH_newdata_predictionRT.csv')
newsmilst <- newinput$smiles
a <- newsmilst[1:5]
newrt <- prediction(a,b='OH')

#load data for grouping based on transitions  #done
file <- import('grouping_example_data.csv')
num <- 50 
m <- 171.1 #for OH 171.1, for cooh 136.1
t <- 1000 #threshold in agilent experiment settings format 
groupresult <- prepare_data(file, num, m, t)
