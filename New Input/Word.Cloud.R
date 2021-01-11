rm(list=ls())
library(wordcloud)
library(RColorBrewer)
library(dplyr)
library(data.table)
setwd("/n/data1/hsph/biostat/celehs/ch263/feature_selection/revision/wordcloud/")

interest.list = c("PheCode:714.1","PheCode:411.4","PheCode:296.2", "PheCode:250.2", 
                  "PheCode:250.1","PheCode:555.2","PheCode:555.1", "PheCode:335")

code.label <- data.frame(rbind(
  c(interest.list[1],'Rheumatoid arthritis'),
  c(interest.list[2],'Coronary Artery Disease'),
  c(interest.list[3],'Depression'),
  c(interest.list[4],'Type 2 diabetes'),
  c(interest.list[5],'Type 1 diabetes'),
  c(interest.list[6],'Ulcerative colitis'),
  c(interest.list[7],'Regional enteritis'),
  c(interest.list[8],'Multiple sclerosis')))

for(ii in c(1:dim(code.label)[1])){
  interest=as.character(code.label[ii,1])
  label=as.character(code.label[ii,2])
  tab.all=fread(paste0("/n/data1/hsph/biostat/celehs/ch263/feature_selection/revision/network/result/tab_[", interest,"]", ".csv"),data.table=F)
  res.rpdr=tab.all[,c("Variable", "Description","Local_Estimator_RPDR")]
  res.va=tab.all[,c("Variable", "Description", "Local_Estimator_VA")]
  colnames(res.rpdr)=colnames(res.va)=c("feature_id", "feature_desc", "coef")
  res.rpdr=res.rpdr[which(is.na(res.rpdr$coef)!=1 & res.rpdr$coef!=0),]
  res.va=res.va[which(is.na(res.va$coef)!=1 & res.va$coef!=0),]
  
  set.seed(1234)
  pdf(file=paste0("result/wordcloud-RPDR-", label,".pdf"),width=12,height=12)
  wordcloud(words = res.rpdr$feature_desc, freq = c(as.numeric(res.rpdr$coef))*10000, random.order=F,
            colors=brewer.pal(8, "Dark2"),scale = c(3.5, 0.2),rot.per = 0)
  dev.off()
  
  set.seed(1234)
  pdf(file=paste0("result/wordcloud-VA-", label,".pdf"),width=12,height=12)
  wordcloud(words = res.va$feature_desc, freq = c(as.numeric(res.va$coef))*10000, random.order=F,
            colors=brewer.pal(8, "Dark2"),scale = c(3.5, 0.2),rot.per = 0)
  dev.off()
}
 

