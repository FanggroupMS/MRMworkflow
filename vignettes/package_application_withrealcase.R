#applying the total work-flow
packsneed <- c('xcms','magrittr','MSnbase','dplyr','tidyr',
               'ggplot2','tidyverse','ggpubr','miceadds',
               'rcdk','rcdklibs','randomForest','leaps',
               'caret','corrplot','mlr','Metrics','rio','hash','EnhancedVolcano')
lapply(packsneed, require, character.only = TRUE)

#outside function, file import
#define wrapper function for reading the mrm files
#funciton definition: read mrm file with default experiment setting and allow noise parameters imput  #done
#function param
filedat <- list.files(path = 'datasets/', pattern = ".mzML$", full.names = TRUE)
pwid <- c(0.8,3)
snthr<- 10
MSRT <- read.csv("MS1_RT.csv", col.names = c('MS1','RT'))
rtfile <- data.frame("row"=1:length(MSRT$RT),"ms1"=MSRT$MS1, "predrt"=MSRT$RT)
delrt <- 1.5
group1 <- c(1:3)
group2 <- c(4:5)
peaksresult <- peaksdetect(filedat, rtfile,pwid, snthr, delrt, group1,group2)
peaksresult
volcanoplot(peaksresult)
dev.copy(png,filename="volcano.png")
dev.off ()

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



########################################

#load smiles of molecules from OH compounds
rtdata <- read.csv("OH_RT_zj.csv")
smilst <- rtdata$smiles
rt <- rtdata$RT
#model start training
datadescs <- getdesc(smilst)
norm_data <- getnorm(datadescs,rt)
OH_model <- getmodel(norm_data, picname = 'OH_data')

#load smiles of molecules from COOH compounds
rtdata2 <- read.csv("COOH_RT_zj.csv")
smilst2 <- rtdata2$smiles
rt2 <- rtdata2$RT
#model start training
datadescs2 <- getdesc(smilst2)
norm_data2 <- getnorm(datadescs2,rt2)
COOH_model <- getmodel(norm_data2, picname = 'COOH_data')

#save model settings
OH_descs <- colnames(datadescs) #descripotrs before norm
OH_descnorsettings <- getnormsettings(datadescs)
COOH_descs <- colnames(datadescs2) #descriptors before norm
COOH_descnorsettings <- getnormsettings(datadescs2)

save(OH_model,file = 'OH_model.RData')
save(OH_descs,file = 'OH_descs.RData')
save(OH_descnorsettings,file = 'OH_normsettings.RData')
save(COOH_model,file = 'COOH_model.RData')
save(COOH_descs,file = 'COOH_descs.RData')
save(COOH_descnorsettings,file = 'COOH_normsettings.RData')
