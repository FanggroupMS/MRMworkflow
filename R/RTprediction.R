#get the rt prediction and the ie prediciton
#get the function for risk socres weight
#get the function for ms confidence scores weight
#------------------------rt prediction module------import rt model--------------------
#load all packages at once
packsneed <- c('rcdk','rcdklibs','randomForest','leaps','caret','corrplot',
               'mlr','dplyr','Metrics','ggpubr','ggplot2','miceadds')
lapply(packsneed, require, character.only = TRUE)

#get prediction from new data
prediction <- function(input,b){
  descs <- getdesc_nopp(input) #get mds without any data imputation
  #find a way to load colnames data, settings,and model
  load.Rdata(paste0(b,'_descs.RData'),'descname')
  load.Rdata(paste0(b,'_normsettings.RData'),'settings')
  load.Rdata(paste0(b,'_model.RData'),'predmodel')
  descsdata <- descs[,descname] #select mds to prepare same dimensions for normalzation settings
  descsdata <- predict(settings,descsdata) #get normalized mds for modeling
  predictionrt <- predict(predmodel, newdata = descsdata) #model select the mds as it trained
  output <- cbind('smiles'=input,'predictionRT'= predictionrt)
  return(output)
}

# load new data
# newinput <- read.csv('OH_newdata_predictionRT.csv')
# newsmilst <- newinput$smiles
# a <- newsmilst[1:5]
# newrt <- prediction(a,b='OH')

#define functions for molecular descriptors without preprocess
getdesc_nopp <- function(input){
  mols <- parse.smiles(input)
  descNames <- unique(unlist(sapply(get.desc.categories(),get.desc.names)))
  descs_tot <- data.frame()
  for (mol in mols) {
    descs_temp <- eval.desc(mol,descNames,verbose = FALSE)
    if (dim(descs_tot)[1] == 0) {
      descs_tot <- descs_temp
    }
    descs_tot <- rbind(descs_tot,descs_temp)
  }
  descs_start <- descs_tot[-c(1),]
  return(descs_start)
}

#define functions for molecular descriptors with preprocess
getdesc <- function(smilist) {
  mols <- parse.smiles(smilist)
  descNames <- unique(unlist(sapply(get.desc.categories(), get.desc.names))) #get 51 types of mds
  descs_tot <- data.frame()
  for (mol in mols) {
    descs_temp <- eval.desc(mol,descNames,verbose = FALSE)
    if (dim(descs_tot)[1] == 0) {
      descs_tot <- descs_temp
    }
    descs_tot <- rbind(descs_tot,descs_temp)
  }
  descs_start <- descs_tot[-c(1),]
  descs_start2<- descs_start[,!apply(descs_start, 2,function(x) any(is.na(x)))]  #220 left
  descs_clean2 <- descs_start2[,!apply(descs_start2,2,function(x) length(unique(x))==1)] #166 left
  r2 <- which(cor(descs_clean2)^2 >.9, arr.ind = TRUE)
  r2 <- r2[r2[,1]>r2[,2],]
  descs_output <- descs_clean2[,-unique(r2[,2])]
  return(descs_output)
}

#normalization by scaling with mean for each column
getnorm <- function(a,b){
  input_norm <- predict(preProcess(a), a)
  input_norm$RT <- b
  return(input_norm)
}

#function for getting model
getmodel <- function(input,picname){
  #data splitting
  set.seed(123)
  smp_size <- round(nrow(input)*0.75)
  index <- sample(seq_len(nrow(input)), size = smp_size)
  train <- input[index, ]
  test <- input[-index, ]
  #rfe feature selection
  subsets <- c(10, 15, 20, 25,30,50)
  ctrl <- rfeControl(functions = rfFuncs, 
                     method = 'repeatedcv', 
                     repeats = 10, 
                     verbose = FALSE)
  y<-train$RT
  x<-subset(train,select = -c(RT))
  lmprofile <- rfe(x,y,
                   sizes = subsets,
                   rfeControl = ctrl)
  #random forest modeling
  selmds <- predictors(lmprofile)
  x <- subset(x,select = c(selmds))
  rf.model = randomForest(y~., data=x, 
                          mtry=9,
                          ntree =500,
                          importance = TRUE)
  #output model performance with prediction error
  test_y<-test$RT
  test_x<-subset(test,select = -c(RT))
  prd_rt = predict(rf.model, newdata = test_x)
  #data prep for plot
  prediction_table <- as.data.frame(cbind(prd_rt, test_y),row.names = NULL)
  #plot 
  plot(lmprofile, type=c('g','o'), main='RMSE of 5cv', lwd=1.0, cex = 1.0)
  plot(rf.model)
  p <- ggscatter(prediction_table, x = "prd_rt", y = "test_y", 
                 add = "reg.line", conf.int = TRUE, 
                 cor.coef = TRUE, cor.method = "pearson",
                 xlab = "RT prediction (mins)", ylab = "RT experiment (mins)")
  p+ geom_text(x=9,y=13, label=paste('mae =',round(mae(prd_rt,test_y),2)),size = 5) #need to change x,y positions
  ggsave(paste (picname,"prediction vs experiment.png"))
  return(rf.model)
}

#get normalization settings
getnormsettings <- function(input){
  normsettings <- preProcess(input)
  return(normsettings)
}

#function for features selection after normalized data
getseldesc <- function(input){
  set.seed(123)
  smp_size <- round(nrow(input)*0.75)
  index <- sample(seq_len(nrow(input)), size = smp_size)
  train <- input[index, ]
  test <- input[-index, ]
  #rfe feature selection
  subsets <- c(5,10, 15, 20, 25,30,50)
  ctrl <- rfeControl(functions = rfFuncs, 
                     method = 'repeatedcv', 
                     repeats = 10, 
                     verbose = FALSE)
  y<-train$RT
  x<-subset(train,select = -c(RT))
  lmprofile <- rfe(x,y,
                   sizes = subsets,
                   rfeControl = ctrl)
  seldesc <- predictors(lmprofile)
  return(seldesc)
}

# #load smiles of molecules from OH compounds
# rtdata <- read.csv("OH_RT_zj.csv")
# smilst <- rtdata$smiles
# rt <- rtdata$RT
# #model start training
# datadescs <- getdesc(smilst)
# norm_data <- getnorm(datadescs,rt)
# OH_model <- getmodel(norm_data, picname = 'OH_data')

# #load smiles of molecules from COOH compounds
# rtdata2 <- read.csv("COOH_RT_zj.csv")
# smilst2 <- rtdata2$smiles
# rt2 <- rtdata2$RT
# #model start training
# datadescs2 <- getdesc(smilst2)
# norm_data2 <- getnorm(datadescs2,rt2)
# COOH_model <- getmodel(norm_data2, picname = 'COOH_data')

# #save model settings
# OH_descs <- colnames(datadescs) #descripotrs before norm
# OH_descnorsettings <- getnormsettings(datadescs)
# COOH_descs <- colnames(datadescs2) #descriptors before norm
# COOH_descnorsettings <- getnormsettings(datadescs2)

# save(OH_model,file = 'OH_model.RData')
# save(OH_descs,file = 'OH_descs.RData')
# save(OH_descnorsettings,file = 'OH_normsettings.RData')
# save(COOH_model,file = 'COOH_model.RData')
# save(COOH_descs,file = 'COOH_descs.RData')
# save(COOH_descnorsettings,file = 'COOH_normsettings.RData')
