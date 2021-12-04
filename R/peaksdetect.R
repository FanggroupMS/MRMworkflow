#loading packages
packsneed <- c('xcms','magrittr','MSnbase','dplyr','tidyr','ggplot2',
               'tidyverse','ggpubr')
lapply(packsneed, require, character.only = TRUE)

#outside function, file import
#define wrapper function for reading the mrm files
#funciton definition: read mrm file with default experiment setting and allow noise parameters imput
#function param
filedat <- list.files(path = 'dmrm/', pattern = ".mzML$", full.names = TRUE)
pwid <- c(0.8,3)
snthr<- 10
MSRT <- read.csv("MS1_RT.csv", col.names = c('MS1','RT'))
rtfile <- data.frame("row"=1:length(MSRT$RT),"ms1"=MSRT$MS1, "predrt"=MSRT$RT)
delrt <- 1.5
group1 <- c(1:3)
group2 <- c(4:5)
peaksresult <- peaksdetect(filedat, rtfile,pwid, snthr, delrt, group1,group2)
peaksresult

peaksdetect <- function(filedat, rtfile,pwid, snthr, delrt, group1,group2){
  dmrm <- readSRMData(filedat)
  #not considering rt correction
  cwp <- CentWaveParam(peakwidth = pwid, snthresh = snthr) #snthresh ratio adjust to 10, peakwidth to 0.5-3
  grouppeaks <- findChromPeaks(dmrm,cwp) #adjust the cwp for the best fit
  mspeaks <- chromPeaks(grouppeaks)
  rawpeaks <- as.data.frame(mspeaks)
  #add ms1,rt to the peakshandler table
  peaktable <- left_join(rawpeaks,rtfile, by='row')
  peaktable_nona <- peaktable %>% 
    na.omit() %>%
    mutate('deltart' = abs(rt-predrt)) %>%
    filter(deltart <= delrt)
  #method2: annova on sample 1-3, and 4-5
  #log2 foldchange
  anovtable <- peaktable %>%
    mutate('deltart' = abs(rt-predrt)) %>%
    filter(deltart <= 1.5) %>%
    select(column,row,into)%>% 
    pivot_wider(names_from = column, values_from =into, names_prefix='sample') %>%
    select(-c(row))
  rawpvalue <- apply(anovtable, 1, ttestpeaks, grp1 = group1, grp2 = group2)
  log2peaks <- log2(anovtable)
  ctrmean <- apply(log2peaks[, group1],1, mean)
  trmean <- apply(log2peaks[,group2],1,mean)
  foldchange <- ctrmean- trmean
  res1 <- cbind(foldchange,rawpvalue)
  res1 <- as.data.frame(res1)
  res1$MStransition <- rownames(res1)
  return(res1)
  }

ttestpeaks <- function(df, grp1, grp2) {
  x = df[grp1]
  y = df[grp2]
  x = as.numeric(x)
  y = as.numeric(y)  
  results = t.test(x, y)
  results$p.value
  return(results$p.value)
}
