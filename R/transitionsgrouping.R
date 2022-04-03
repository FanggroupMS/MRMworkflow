library(rio)
library(dplyr)
library(hash)

#if compound is OH, deltams1 light=1, heavy =2
#if compound is COOH, deltams1 light =1, heavy= 3
#parameters: file= original transitions, num =, m =,
transitiongroup = function(file, num, m) {
  # read data
  source_data = file
  # calculate the precursor M+H+ mass
  source_data = mutate(source_data, MS1 = MS1 + 1.0073) 

  # round number upwords
  source_data = mutate(source_data, RTgroup=ceiling(RT), 
                       MSgroup=sprintf("%0.2f", ceiling(MS1 * 10) / 10))
  #count the number of compounds with same MS1 mass
  table = count(source_data, MSgroup)
  #grouping compounds into n files (n = count of compounds with same MS1 mass)
  export_num = slice(table, which.max(table$n))
  #find the folder for storing output files
  time_now = round(as.numeric(as.POSIXct(Sys.time())))
  dir_name = paste("./temp_", time_now, sep = "")
  dir.create(dir_name)
  #the folder for output result
  result_dir_name = paste("./result_", time_now, sep = "")
  dir.create(result_dir_name)
  
  export(table, paste(result_dir_name, "/statistics.xlsx", sep = ""))
  build_intermediate_xlsx(source_data, table, export_num$n, num, dir_name, m, 
                          result_dir_name)
}

build_intermediate_xlsx = function(source_data, table, export_num, num, 
                                   dir_name, m, result_dir_name) {
  #loop over every output files
  for (i in 1:as.numeric(export_num)) {
    #output data format
    result = data_frame()
    names(result) = c("casrn", "smiles", "RT", "MW", "MS1", "RTgroup", "MSgroup")
    #form each result for each ms group
    for (line in filter(table, n >= i)$MSgroup) {
      result = rbind(result, slice(filter(source_data, MSgroup == line), i))
    }
    #output the final results
    result = arrange(result, RTgroup)
    export(result, paste(dir_name, "/result", i, ".xlsx", sep = ""), which = "sheet1")
  }
  get_random_result(dir_name, num, m, result_dir_name)
}

format_export_data = function (path, i, data, m, sheet_name) {
  data_length = dim(data)[1]
                    
  output = data.frame(
    "Compound Group" = rep("neg", data_length),
    "Compound Name" = data$casrn,
    "ISTD?" = rep("FALSE", data_length),
    "Precursor Ion" = data$MS1,
    "MS1 Res" = rep("Unit", data_length),
    "Product Ion" = rep(m, data_length),
    "MS2 Res" = rep("Unit", data_length),
    "Primary" = rep("TRUE", data_length),
    "Trigger" = rep("FALSE", data_length) ,
    "Ret Time (min)" = data$RT,
    "Delta Ret Time" = rep("2", data_length),
    "Fragmentor" = rep("166", data_length) ,
    "Collision Energy" = rep("12", data_length),
    "Cell Accelerator Voltage" = rep("4", data_length) ,
    "Polarity" = rep("Positive", data_length),
    "Trigger Entrance Delay (cycles)" = rep("0", data_length),
    "Trigger Delay (cycles)" = rep("0", data_length),
    "Trigger Window" = rep("0", data_length),
    "IsLogicEnabled" = rep("FALSE", data_length),
    "Trigger Logic Flag" = rep("AND", data_length),
    "Trigger Ratio" = rep("1", data_length),
    "Trigger Ratio Window" = rep("1", data_length),
    "Ignore MRM" = rep("FALSE", data_length), check.names = F
  )
  #output final files with experiment settings
  export(output, paste(path, "/result", i, ".xlsx", sep = ""), 
         which = sheet_name)
}

get_random_result = function(dir_name, num, m, result_dir_name) {
  i = 1
  for (file in list.files(dir_name)) {
    print(paste("---processing", i, "files---"))
    if (file == "statistics.xlsx") {
      next
    }
    data = import(paste(dir_name, "/", file, sep = ""))
    # RTgroup counting
    temp_table = arrange(data, RTgroup) %>% count(RTgroup) %>% 
      filter(n > num) %>% mutate(left_count = n - num)
    if (dim(temp_table)[1] == 0) {
      # if no further match for the MS group, then continue the searching
      format_export_data(result_dir_name, i, data, m, "sheet1")
      i = i + 1
      next
    }
    # store the row from the data out of the RT groups into a new datafile
    left_data = data_frame()
    names(left_data) = c("casrn", "smiles", "RT", "MW", "MS1", "RTgroup", "MSgroup")
    for (line in temp_table$RTgroup) {
      # take tables for splitting, eliminate those not wanted
      left_num = filter(temp_table, RTgroup==line)$left_count
      RT_data = filter(data, RTgroup == line)
      left_data = rbind(left_data, RT_data[sample(nrow(RT_data), left_num, 
                                                  replace = F),])
    }
    data = setdiff(data, left_data)
    # order data
    data = arrange(data, RTgroup)
    left_data = arrange(left_data, RTgroup)
    # export data table
    format_export_data(result_dir_name, i, data, m, "sheet1")
    format_export_data(result_dir_name, i, left_data, m, "sheet2")
    i = i + 1
  }
  # delet temp files
  unlink(dir_name, recursive = TRUE)
}

execute = function () {
  print("choose files input")
  file = import(file.choose())
  print("input num")  #number of transitions allowed in a RT group
  num = as.numeric(readline())
  print("input m")   #product ion value
  m = readline()
  transitiongroup(file,num,m)
}
