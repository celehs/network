rm(list=ls())
library(wordcloud)
library(RColorBrewer)
library(dplyr)
library(data.table)
setwd("/n/data1/hsph/biostat/celehs/ch263/feature_selection/revision/wordcloud/")

interest = c("U07.1")
label <- "COVID19"
  res.rpdr=fread("/n/data1/hsph/biostat/celehs/ch263/feature_selection/revision/covid/cos_covidmgb_final.csv")[,c(2,3,4)]
  res.va=fread("/n/data1/hsph/biostat/celehs/ch263/feature_selection/revision/covid/cos_covidva_final.csv")
  colnames(res.rpdr)=colnames(res.va)=c("feature_id",  "coef","feature_desc")
  res.rpdr=res.rpdr[which(is.na(res.rpdr$coef)!=1 & res.rpdr$coef!=0 & res.rpdr$feature_id!="ICD10:U07.1"),]
  res.va=res.va[which(is.na(res.va$coef)!=1 & res.va$coef!=0 & res.va$feature_id!="U07.1"),]
  res.va=res.va[res.va$coef>quantile(res.va$coef, 0.98),]
  res.rpdr=res.rpdr[res.rpdr$coef>quantile(res.rpdr$coef, 0.98),]
  
  set.seed(124)
  pdf(file=paste0("result/wordcloud-RPDR-", label,".pdf"),width=12,height=12)
  wordcloud(words = res.rpdr$feature_desc, freq = c(as.numeric(res.rpdr$coef))*10000, random.order=F,
            colors=brewer.pal(8, "Dark2"),scale = c(3.5, 0.2),rot.per = 0)
  dev.off()
  
  set.seed(124)
  pdf(file=paste0("result/wordcloud-VA-", label,".pdf"),width=12,height=12)
  wordcloud(words = res.va$feature_desc, freq = c(as.numeric(res.va$coef))*10000, random.order=F,
            colors=brewer.pal(8, "Dark2"),scale = c(3.5, 0.2),rot.per = 0)
  dev.off()

 

