packsneed <- c('xcms','magrittr','MSnbase','dplyr','tidyr','ggplot2','tidyverse','ggpubr',"ggrepel","rio")
lapply(packsneed, require, character.only = TRUE)

#function for screening peaks with both light and heavy signals
data_process = function (raw_data, wrong_data, output) {
  # read through data
  category_list = unique(raw_data$ID2)
  for (category in category_list) {
    # examine the same row
    data = filter(raw_data, ID2==category)
    # see if legal 12,or 13
    data_12_line_num = c()
    data_13_line_num = c()
    for (i in 1:nrow(data)) {
      row_data = data[i, ]
      if (grepl(12, row_data$ID)) {
        data_12_line_num = c(data_12_line_num, i)
      }
      if (grepl(13, row_data$ID)) {
        data_13_line_num = c(data_13_line_num, i)
      }
    }
    count_12 = length(data_12_line_num)
    count_13 = length(data_13_line_num)
    if (count_12 == 0 | count_13 == 0) {
      # wrong data if 0 is inside
      wrong_data = rbind(wrong_data, data)
    } else if (count_12 > count_13) { # take smaller value from rows
      ret = data_filter(data, data_12_line_num, data_13_line_num, 
                               wrong_data, output)
      wrong_data = ret$w
      output = ret$o
    } else {
      ret = data_filter(data, data_13_line_num, data_12_line_num, 
                               wrong_data, output)
      wrong_data = ret$w
      output = ret$o
    }
  }
  return(list(w=wrong_data, o=output))
}

#function for tag the light and heavy isotope amount the data, add the 4th ID for further matching
data_filter = function (data, num_1, num_2, wrong, output) {
  # set num_2 as starts
  num = 1
  for(i in num_2) { # split
    chosen_row = data[i, ]
    quantity = 0
    for (j in num_1) {
      compare_row = data[j, ]
      # condition
      deltart2 =abs(chosen_row$rt-compare_row$rt)
      deltamaxo = chosen_row$maxo/compare_row$maxo
      if (deltart2 <= 0.1 & deltamaxo >= 0.25 & deltamaxo <= 4) {
        # take row satisfy these three rules
        # not considering 1 vs 3
        quantity = quantity + 1
        #add a new 4th name for later cross matching
        compare_row$ID3 = paste(compare_row$ID2, "-", num, sep="")
        chosen_row$ID3 = paste(compare_row$ID2, "-", num, sep="")
        output = rbind(output, compare_row)
        output = rbind(output, chosen_row)
      }
    }
    num = num + 1
    if (quantity != 1) {
      wrong = rbind(wrong, data)
      return(list(w=wrong, o=output))
    }
  }
  return(list(w=wrong, o=output))
} 

get_category = function (data) {
  #define wrong data
  wrong_data = data.frame()
  #define right data
  output = data.frame()
  for(num in 1:6) {
    w_data = data.frame()
    o_data = data.frame()
    datalist = filter(data, column == num)
    ret = data_process(datalist, w_data, o_data) #easy to go wrong
    wrong_data = rbind(wrong_data, ret$w)
    output = rbind(output, ret$o)
  }
  return(list(w=wrong_data, o=output))
}

#set function for screening the extract matches between the 12 and 13 peaks 
kill = function (data) {
  for (i in 1:6) {
    waiting_list = c()
    list = filter(data, column == i)
    for (category in unique(list$ID3)) {
      if (nrow(count(list, ID3 == category)) == 1) {
        if (count(list, ID3 == category)$n[1] != 2) {
          waiting_list = c(waiting_list, category)
        }
      } else {
        if (count(list, ID3 == category)$n[2] != 2) {
          waiting_list = c(waiting_list, category)
        }
      }
    }
    print(waiting_list)
    if (!is.null(waiting_list)) {
      for (id3 in waiting_list) {
        rt13 = filter(list, ID3 == id3 & grepl("13", ID))$rt[1]
        rt12 = filter(list, ID3 == id3 & grepl("12", ID))$rt
        final_12 = 0
        x = 0
        for (rt in rt12) {
          if (abs(rt - rt13) > x) {
            final_12 = rt
            x = abs(rt - rt13)
          }
        }
        final_data = rbind(filter(list, ID3 == id3 & grepl("13", ID))[1,],
                           filter(list, ID3 == id3 & grepl("12", ID) & rt == final_12))
        # find row id in data
        line_num = c()
        for (id3 in waiting_list) {
          line_num = c(line_num, which(data$ID3 == id3 & data$column == i))
        }
        for (num in line_num) {
          data = data[-line_num[1],]
        }
        data = rbind(data, final_data)
      }
    }
  }
  return(data)
}

#remove peakid with invalide samples
#take the finaldata format from the export function
peaksep = function (data) {
  final = data.frame()
  peak = unique(data$ID3)
  peak2 = c()
  for (id3 in peak) {
    nacount1 = 0
    nacount2 = 0
    #find out the bad rows in each treatment group
    for (i in 1:3) {
      if (nrow(filter(data, column == i & ID3 == id3)) == 0) {
        nacount1 = nacount1 + 1 
        next
      }
    }
    for (i in 4:6) {
      if (nrow(filter(data, column == i & ID3 == id3)) == 0) {
        nacount2 = nacount2 + 1 
        next
      }
    }
    if (nacount1 <= 1 & nacount2 <=1){
      peak2 = c(peak2,id3)
    } else{
      next
    }
  }
  return(peak2)
}

#transfer and combined all the data in the row and samples format
export_func = function (data,group1,group2) {
  finaldat = data.frame()
  # peak = peaklist #unique peak id
  peak = unique(data$ID3)
  head = c("ID3",group1,group2)
  #collect peaks according to peakid
  for (id3 in peak) {
    temprow = c(id3)
    #find out the bad rows in each treatment group
    for (i in 1:6) {
      if (nrow(filter(data, column == i & ID3 == id3)) == 0) {
        temprow[i+1] = NA
        next
      }
      temprow[i+1] = filter(data, column == i & ID3 == id3)$maxo
    }
    finaldat = rbind(finaldat, temprow)
  }
  names(finaldat) = head
  # rownames(final) <- final$ID3
  # final <- subset(final,select=-c(ID3))
  # export(finaldat, 'final_intensitydata_concentration.xlsx')
  return(finaldat)
}

#ttest analysis
ttestpeaks <- function(df, group1, group2) {
  x = df[group1]
  y = df[group2]
  x = as.numeric(x)
  y = as.numeric(y)  
  results = t.test(x, y,
                   alternative=c("two.sided"),
                   paired=FALSE,var.equal=FALSE)
  results$p.value
  return(results$p.value)
}

#define function removing nas with more than 2 samples in one group
filterna <- function(dataset){
  dataset2 = data.frame()
  for (i in 1:nrow(dataset)){
    nacounts = sum(is.na(dataset[i,]))
    if (nacounts >=1 ){
      napos = which(is.na(dataset[i,]))
      if(napos <= 3){
        temprow = dataset[i,1:3]
        temprow = temprow[-napos]
        temprowmean = mean(as.numeric(temprow))
        dataset[i,napos] = temprowmean
      }else if(napos >= 4){
        temprow = dataset[i,4:6]
        temprow = temprow[-napos]
        temprowmean = mean(as.numeric(temprow))
        dataset[i,napos] = temprowmean
      }
    }else{
      next
    }
  }
  for(j in 1:ncol(dataset)){
    dataset[,j] = as.numeric(dataset[,j])
  }
  dataset2 = dataset
  return(dataset2)
}

#difine foldchange function
foldchange <- function(dataset,group1,group2){
  log2peaks <- log2(dataset)
  mean1 <- apply(log2peaks[, group1],1, mean)
  mean2 <- apply(log2peaks[,group2],1,mean)
  Log2Foldchange <- mean1-mean2
  return(Log2Foldchange)
}

#define the peak detection function and peak table generator
peaksdetect <- function(input,rtdata,pwid, snthr, delrt, intenmin,peakres){
  dmrm <- readSRMData(input)
  cwp <- CentWaveParam(peakwidth = pwid, snthresh = snthr)#snthresh ratio adjust to 10, peakwidth to 0.01-30
  grouppeaks <- findChromPeaks(dmrm,cwp)
  mspeaks <- chromPeaks(grouppeaks)
  rawpeaks <- as.data.frame(mspeaks)
  rawpeaks['Mz1'] <- precursorMz(grouppeaks)[mspeaks[, "row"], ][,1]
  rawpeaks['Mz2'] <- productMz(grouppeaks)[mspeaks[, "row"], ][,1]
  peaktable <- left_join(rawpeaks,rtdata, by='row')
  #peak screening
  #using the Predicted RT for the filtering  #rt <=1.5
  #using the intensity threshold filtering   maxo >= intenmin
  peaktable <- peaktable %>%
    mutate('deltart'= abs(rt-predrt)) %>%
    filter(deltart <= delrt & maxo >= intenmin) %>%
    select(rt,maxo,row,column,Mz1,Mz2,ID, deltart,ID2)
  # save(peaktable, file="concentration_method1.RData")
  export(peaktable, file = paste(peakres,".csv",sep = ""))
  return(peaktable)
}

#define the peak table editor
peakrefine <- function(input,group1,group2){
  res = get_category(input)
  output = res$o
  rownames(output) = seq(1, nrow(output), 1)
  output = filter(output, grepl("12", ID)) # only consider isotope-12 peaks
  # export(ret$w, 'wrong.xlsx') # export wrong data # wrong_data = res$w
  # export(output, 'finaldata.xlsx') #for double check
  #retain peaks for further analysis
  # print(output)
  retain_peak = peaksep(output)
  # print(retain_peak)
  final_output = filter(output,ID3%in%retain_peak)
  final_table = export_func(final_output,group1,group2) #group1 and group2 as parameters
  rownames(final_table) = final_table$ID3
  testtable = subset(final_table,select = -c(ID3))
  # export(testtable, 'finaltestsampledata.xlsx') #for mannual check
  return(testtable)
}

#define the function for differential express and graphs output
diffexp = function(input,group1, group2,resoutput){
  pvalue <- apply(input, 1, ttestpeaks, grp1 = group1, grp2 = group2)
  #replace nas in each group with group means
  input <- filterna(input)
  log2FoldChange <- foldchange(dataset=input,group1=group1,group2=group2)
  res1 <- cbind(log2FoldChange,pvalue)
  res1 <- as.data.frame(res1)
  res1$diffexpressed <- "NO" #add a column of NAs #add label for down and up regulation
  # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
  res1$diffexpressed[res1$log2FoldChange > 0.6 & res1$pvalue < 0.05] <- "UP"
  # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  res1$diffexpressed[res1$log2FoldChange < -0.6 & res1$pvalue < 0.05] <- "DOWN"
  res1$delabel <- NA
  res1$symbol <- rownames(res1)
  res1$delabel[res1$diffexpressed != "NO"] <- res1$symbol[res1$diffexpressed != "NO"]
  #volcanoplot
  p <- ggplot(data=res1, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) +
    geom_point(size=3) + 
    theme(text = element_text(size=15, family = 'TT Arial',face='bold'),
          axis.text.x = element_text(size=15, family = 'TT Arial',face='bold'),
          axis.text.y = element_text(size=15, family = 'TT Arial',face='bold')) +
    theme_minimal() +
    geom_text_repel() +
    scale_color_manual(values=c("blue", "black", "red")) +
    geom_vline(xintercept=c(-0.6, 0.6), col="red") +
    geom_hline(yintercept=-log10(0.05), col="red")
  print(p) #return the plot
  ggsave(filename = paste(resoutput,".tiff",sep = ""), units="in", width = 8, height=6, dpi=400)
  return(res1)
}
