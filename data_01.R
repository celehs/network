rm(list=ls())
library("visNetwork")
library("RColorBrewer")
library(htmlwidgets)
library(base)
library(data.table)
library(dplyr)

setwd("/Users/clara-lea/Documents/GitHub/network.vis")
source("myfunction_new.R")

clean.fun=function(feature_desc){
  list.remove="LOINC:|Other lab:|Shortname:|ShortName:| preservative free|,,|, ,"
  name.new=gsub(list.remove, "", feature_desc)
  name.new=tolower(name.new)
  #text.rm=gsub("\\(([^()]*)\\)|.", "\\1", name.new)
  #name.new=gsub(text.rm, "", name.new)
  name.new=gsub("[()]", "", name.new)
  if(nchar(name.new) > 50){
    temp = unlist(strsplit(name.new, split = " "))
    nm1=temp[1]
    nm2=temp[2]
    nm3=temp[3]
    nm4=temp[length(temp)]
    nm = paste0(nm1, " ", nm2, " ", nm3, "...", nm4)
  }else(nm=name.new)
  nm
}

# upload dictionnary results
load("dict_combine_uniform_molei.Rdata")
dict.combine=dict.combine[,c("feature_id", "feature_desc")]
dict.combine$feature_desc=unlist(lapply(1:dim(dict.combine)[1], function(jj) clean.fun(dict.combine$feature_desc[jj])))
colnames(dict.combine)=c("Variable","Description")

# upload network data
interest.list=list.files("network_data")[-1]
interest.list=gsub("tab_|\\[|\\]|.csv", "",interest.list)

data.res <- vector("list", length(interest.list))
names(data.res) <- interest.list
for (ii in 1:length(interest.list)){
  interest=interest.list[ii]
  interest.label=dict.combine[dict.combine$Variable==interest,"Description"]
  res=read.csv(paste0("network_data/tab_[", interest,"]", ".csv"))
  res[is.na(res)]=0
  res[,1]=as.character(res[,1])
  res[,2]=as.character(res[,2])
  
  res=res[,c("Variable", "Description", "Integrative_Estimator_RPDR", "Integrative_Estimator_VA")]
  colnames(res)[-c(1:2)]=c("int_RPDR", "int_VA")
  #
  res['int01'] <-rep(0,length(res$int_RPDR))
  res$int01[which(res$int_RPDR!=0 | res$int_VA!=0)] <- 1
  res <- subset(res, select=c("Variable", "Description","int01"))
  
  data.res[[interest]] <- res
}

saveRDS(data.res, 'data_res01.rds')


