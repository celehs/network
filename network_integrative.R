rm(list=ls())
library("visNetwork")
library("RColorBrewer")
library(htmlwidgets)
library(base)
library(data.table)
library(dplyr)

setwd("/n/data1/hsph/biostat/celehs/cb334/feature_selection/result/network_result")
source("../../code/myfunction_new.R")

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

load("/n/data1/hsph/biostat/celehs/ch263/feature_selection/dictionary/result/dict_combine_uniform_molei.Rdata")
dict.combine=dict.combine[,c("feature_id", "feature_desc")]
dict.combine$feature_desc=unlist(lapply(1:dim(dict.combine)[1], function(jj) clean.fun(dict.combine$feature_desc[jj])))
colnames(dict.combine)=c("Variable","Description")

interest.list=list.files("../network_data")
interest.list=gsub("tab_|\\[|\\]|.csv", "",interest.list)

main.kern=function(ii){
  interest=interest.list[ii]
  interest.label=dict.combine[dict.combine$Variable==interest,"Description"]
  res=read.csv(paste0("../network_data/tab_[", interest,"]", ".csv"))
  res[is.na(res)]=0
  res[,1]=as.character(res[,1])
  res[,2]=as.character(res[,2])

  res=res[,c("Variable", "Description", "Integrative_Estimator_RPDR", "Integrative_Estimator_VA")]
  colnames(res)[-c(1:2)]=c("int_RPDR", "int_VA")
  res1=res[res[,"int_RPDR"]!=0&res[,"int_VA"]!=0,c(1,2,3,4)]
  res1$B=rowMeans(res1[,c(3,4)])
  res1=res1[,c("Variable", "Description", "B")]
  res2=res[,c("Variable","Description","int_RPDR")]
  colnames(res2)[3]="B"
  res3=res[,c("Variable","Description","int_VA")]
  colnames(res3)[3]="B"
  res2=res2[which(res2$Variable%in%res1$Variable!=1),]
  res2=res2[which(res2$B!=0),]
  res3=res3[which(res3$Variable%in%res1$Variable!=1),]
  res3=res3[which(res3$B!=0),]
  colnames(res1)=colnames(res2)=colnames(res3)=c("feature_id", "feature_desc", "B")
  net=net.fun(interest, interest.label, res1,res2,res3)
  saveWidget(net, file = paste0("network_[",gsub(":",".",interest),"].html"))
}


for (ii in 1:length(interest.list)){
  print(ii)
  tryCatch(main.kern(ii),error=function(e) NA)
}

