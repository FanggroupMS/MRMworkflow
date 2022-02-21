# CILMRM-workflow

To Install from R console:

````
install.packages("devtools", dependencies=TRUE)
library(devtools) 

install_github("YANGJJ93MS/MRMworkflow")
library(CILMRM)
````


To Install R packages required:

````
install.packages(c('xcms','magrittr','MSnbase','dplyr','tidyr',
                   'ggplot2','tidyverse','ggpubr','miceadds',
                   'rcdk','rcdklibs','randomForest','leaps',
                   'caret','corrplot','mlr','Metrics','rio','hash'))
                   
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
BiocManager::install('EnhancedVolcano')
````


This package requires R version 4.1.2.
