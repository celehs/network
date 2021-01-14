rm(list=ls())
library(visNetwork)
library(RColorBrewer)
library(htmlwidgets)
library(base)
library(data.table)
library(dplyr)
library(network)
library(igraph)
library(igraphdata)
library(Matrix)
library(stringr)

mainpath = "D:/xxiong/program/network/network.vis-main"
setwd(mainpath)

clean.fun=function(feature_desc){
  list.remove="LOINC:|Other lab:|Shortname:|ShortName:| preservative free|,,|, ,"
  name.new=gsub(list.remove, "", feature_desc)
  name.new=tolower(name.new)
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
interest.list=list.files("network_data")
interest.list=gsub("tab_|\\[|\\]|.csv", "",interest.list)
interest.list=gsub("_", ":",interest.list)

for (i in 1:length(interest.list)) {
  interest = interest.list[[i]]
  interest1=gsub(":", "_",interest)
  res = read.csv(paste0("network_data/tab_[",interest1,"].csv"))
  dict.combine = rbind(dict.combine, data.frame("Variable"=res$Variable[!res$Variable %in% dict.combine$Variable],
                                                "Description"=res$Description[!res$Variable %in% dict.combine$Variable]))
}

dict.combine$label = sub(pattern = "(^\\D*)", "", dict.combine$Variable,fixed = FALSE)
Cap = substr(dict.combine$Variable,1,1)
dict.combine$Cap = Cap
dict.combine$label = paste0(Cap,dict.combine$label)
dict.combine$label[Cap=="S"] = sub(pattern = "(^\\D*[:])", "", dict.combine$Variable[Cap=="S"],fixed = FALSE)
n.node = nrow(dict.combine)

# update description
dict.combine$Description_s = dict.combine$Description
dict.combine$Description_s[Cap=="L"] = sub(pattern = "(^.*group:)", "",
                                           dict.combine$Description[Cap=="L"],fixed = FALSE)


name = dict.combine$Variable


add.edge <- function(interest, n.node, name){
  vec.all = vec.RPDR.local = vec.RPDR.inte =
    vec.VA.local = vec.VA.inte = matrix(0,1,n.node)
  colnames(vec.all) = colnames(vec.RPDR.local) = colnames(vec.VA.local) =
    colnames(vec.RPDR.inte) = colnames(vec.VA.inte) = name
  interest1=gsub(":", "_",interest)
  res = read.csv(paste0("network_data/tab_[",interest1,"].csv"))
  is.all = as.character(res$Variable[which(abs(res$Integrative_Estimator_RPDR)>=1e-5 | abs(res$Integrative_Estimator_VA)>=1e-5)])
  all = apply(res[abs(res$Integrative_Estimator_RPDR)>=1e-5 | abs(res$Integrative_Estimator_VA)>=1e-5,
                  c("Integrative_Estimator_RPDR","Integrative_Estimator_VA")],1,function(x){max(abs(x))})

  is.RPDR.local = as.character(res$Variable[which(abs(res$Local_Estimator_RPDR)>=1e-5)])
  RPDR.local = res[abs(res$Local_Estimator_RPDR)>=1e-5, "Local_Estimator_RPDR"]

  is.RPDR.inte = as.character(res$Variable[which(abs(res$Integrative_Estimator_RPDR)>=1e-5)])
  RPDR.inte = res[abs(res$Integrative_Estimator_RPDR)>=1e-5, "Integrative_Estimator_RPDR"]

  is.VA.local = as.character(res$Variable[which(abs(res$Local_Estimator_VA)>=1e-5)])
  VA.local = res[abs(res$Local_Estimator_VA)>=1e-5, "Local_Estimator_VA"]

  is.VA.inte = as.character(res$Variable[which(abs(res$Integrative_Estimator_VA)>=1e-5)])
  VA.inte = res[abs(res$Integrative_Estimator_VA)>=1e-5, "Integrative_Estimator_VA"]

  vec.all[,is.all] = all
  vec.RPDR.local[,is.RPDR.local] = RPDR.local
  vec.RPDR.inte[,is.RPDR.inte] = RPDR.inte
  vec.VA.local[,is.VA.local] = VA.local
  vec.VA.inte[,is.VA.inte] = VA.inte
  a = list(list(1,2,3,4,5))
  a[[1]] = list(list(vec.all),list(vec.RPDR.local),list(vec.VA.local),
                list(vec.RPDR.inte),list(vec.VA.inte))
  return(a)
}

edge.matrix.list = sapply(interest.list, add.edge, n.node=n.node, name=name)

edge.matrix.all = sapply(edge.matrix.list, function(x){unlist(x[[1]])},
                         simplify = TRUE)
edge.matrix.RPDR.local = sapply(edge.matrix.list, function(x){unlist(x[[2]])},
                         simplify = TRUE)
edge.matrix.VA.local = sapply(edge.matrix.list, function(x){unlist(x[[3]])},
                          simplify = TRUE)
edge.matrix.RPDR.inte = sapply(edge.matrix.list, function(x){unlist(x[[4]])},
                                simplify = TRUE)
edge.matrix.VA.inte = sapply(edge.matrix.list, function(x){unlist(x[[5]])},
                              simplify = TRUE)
rownames(edge.matrix.all) = rownames(edge.matrix.RPDR.local) = rownames(edge.matrix.VA.local) =
  rownames(edge.matrix.RPDR.inte) = rownames(edge.matrix.VA.inte) = dict.combine$Variable


sum.nonzero = function(x){
  return(sum(x!=0))
}


edge.matrix.all = edge.matrix.all[apply(edge.matrix.all, 1, sum.nonzero)!=0,]
edge.matrix.RPDR.local = edge.matrix.RPDR.local[apply(edge.matrix.RPDR.local, 1, sum.nonzero)!=0,]
edge.matrix.VA.local = edge.matrix.VA.local[apply(edge.matrix.VA.local, 1, sum.nonzero)!=0,]
edge.matrix.RPDR.inte = edge.matrix.RPDR.inte[apply(edge.matrix.RPDR.inte, 1, sum.nonzero)!=0,]
edge.matrix.VA.inte = edge.matrix.VA.inte[apply(edge.matrix.VA.inte, 1, sum.nonzero)!=0,]

edge.matrix.all = as(as.matrix(edge.matrix.all), "dgCMatrix")
edge.matrix.RPDR.local = as(as.matrix(edge.matrix.RPDR.local), "dgCMatrix")
edge.matrix.VA.local = as(as.matrix(edge.matrix.VA.local), "dgCMatrix")
edge.matrix.RPDR.inte = as(as.matrix(edge.matrix.RPDR.inte), "dgCMatrix")
edge.matrix.VA.inte = as(as.matrix(edge.matrix.VA.inte), "dgCMatrix")
edge.list = edge.full.list = list(edge.matrix.all, edge.matrix.RPDR.local, edge.matrix.VA.local,
                                  edge.matrix.RPDR.inte, edge.matrix.VA.inte)
Cap = Cap[!duplicated(dict.combine$Description)]
dup.name = dict.combine$Variable[duplicated(dict.combine$Description)]
dict.combine=dict.combine[!duplicated(dict.combine$Description),]

for (i in 1:5) {
  edge.matrix = edge.list[[i]]
  row.id = na.omit(match(dup.name,rownames(edge.matrix)))
  col.id = na.omit(match(dup.name,colnames(edge.matrix)))
  edge.matrix = edge.matrix[-row.id,]
  edge.matrix = edge.matrix[,-col.id]
  nodes = unique(c(colnames(edge.matrix),rownames(edge.matrix)))
  n.nodes = length(nodes)
  edge.full.matrix = Matrix(0, n.nodes, n.nodes)
  colnames(edge.full.matrix) = rownames(edge.full.matrix) = nodes
  edge.full.matrix[match(colnames(edge.matrix), nodes),
                    match(rownames(edge.matrix), nodes)] = t(edge.matrix)
  edge.full.list[[i]] = edge.full.matrix
  edge.list[[i]] = edge.matrix
}


n.node = nrow(dict.combine)
dict.combine$group = " "

## PheCode group assignment
load("label data/pheinfo.rda")
P.dict = dict.combine[dict.combine$Cap=="P",]
p.code = sub(pattern="P","",P.dict$label)
P.dict$group = pheinfo$group[match(p.code,pheinfo$phecode)]
table(P.dict$group)
sum(is.na(P.dict$group))
aa = paste0(P.dict$Variable[is.na(P.dict$group)],": ",
       P.dict$Description[is.na(P.dict$group)])
sum(is.na(P.dict$group))
a = table(P.dict$group)
print(length(a))




## CCS group assignment
ccs_map = read.csv("label data/prmlabel-09.csv")
ccs_group = read.csv("label data/ccs_multi_pr_tool_2015.csv")
ccs_group = ccs_group[,c(2,3)]
colnames(ccs_group) = c("cap","group")
ccs_group$cap = ccs_group$cap %>%
  str_replace_all(c("\'" = ""))
ccs_group$group = ccs_group$group %>%
  str_replace_all(c("\'" = ""))
ccs_group = ccs_group[!duplicated(ccs_group$cap),]

C.dict = dict.combine[dict.combine$Cap=="C",]
colnames(ccs_map) = c("ccs.id","ccs.des")
ccs_map$ccs.id = ccs_map$ccs.id %>%
  str_replace_all(c("\'" = "", " " = ""))
ccs_map$ccs.des = ccs_map$ccs.des %>%
  str_replace_all(c("\'" = "", ";" = ",","\\("="","\\)"="","\\[.*]"="")) %>%
  str_to_lower() %>%
  str_trim(side = "both")
c.code = str_match_all(C.dict$Description, "(^[a-zA-Z0-9, -]*)(...)?([a-zA-Z0-9, ]*)")
c.group = sapply(c.code, function(code){
  if(!is.na(code[[3]])){
    pattern = paste0("^",code[[2]],".*",code[[4]],"$")
  }else{
    pattern = paste0("^",code[[2]])
  }
  ccs_map$ccs.id[which(is.na(str_match(ccs_map$ccs.des, pattern))!=TRUE)]
})

c.null = c.mul = c()
jj = ii = c()
c.group.na = rep(1, length(c.code))
for (i in 1:length(c.code)) {
  if(length(c.group[[i]])==0){
    c.null = c(c.null,c.code[[i]][[1]])
    ii = c(ii,i)
    c.group.na[i] = NA
  }else{
    if(length(c.group[[i]])>1){
      c.mul = c(c.mul,c.code[[i]][[1]])
      jj = c(jj,i)
      c.group.na[i] = ccs_group$group[match(str_extract(c.group[[i]][length(c.group[[i]])],"^[0-9]*"),ccs_group$cap)]
    }else{
      c.group.na[i] = ccs_group$group[match(str_extract(c.group[[i]],"^[0-9]*"),ccs_group$cap)]
    }
  }
}
C.dict$group = c.group.na

aa = paste0(C.dict$Variable[ii], ": ", c.null)
paste0(C.dict$Variable[ii], ": ", c.null)
sum(is.na(C.dict$group))
a = table(C.dict$group)
print(length(a))


## LOINC group assignment
loinc_map = read.csv("label data/MultiAxialHierarchy.csv")
L.dict = dict.combine[dict.combine$Cap=="L",]
l.code = str_replace(L.dict$label,"L","")
loinc_map$CODE = str_replace_all(loinc_map$CODE,"LP","")
loinc_map$PATH_TO_ROOT = str_replace_all(loinc_map$PATH_TO_ROOT,"LP","")
l.group = as.character(loinc_map$PATH_TO_ROOT[match(l.code, loinc_map$CODE)])
l.group.list = str_split(l.group, "\\.")
n.group.max = max(sapply(l.group.list, function(x){length(x)}))
for (j in 1:length(l.group.list)) {
  k = 0
  for (i in 1:n.group.max) {
    if(length(l.group.list[[j]])>=i){
      L.dict[j,paste0("H",i)] = loinc_map$CODE_TEXT[match(l.group.list[[j]][i],loinc_map$CODE)]
      k = k + 1
    }else{
      L.dict[j,paste0("H",i)] = loinc_map$CODE_TEXT[match(l.group.list[[j]][k],loinc_map$CODE)]
    }
  }
}
L.dict$group = L.dict$H4
L.dict$group = str_replace_all(L.dict$group,"\\(.*\\)","")
aa = paste0(L.dict$Variable[is.na(L.dict$group)],": ",
            L.dict$Description[is.na(L.dict$group)])
#write.csv(aa,"label data/L.csv")
#write.csv(L.dict,"label data/L data.csv")
paste0(L.dict$Variable[is.na(L.dict$group)],": ",
       L.dict$Description[is.na(L.dict$group)])
sum(is.na(L.dict$group))


a=unique(L.dict$H4)
length(a)
## RXNORM group assignment
rxnorm_map = read.csv("label data/rxnorm_class_all.csv")
rxnorm_map$VA = str_to_lower(rxnorm_map$VA)
r.group = str_split(rxnorm_map$VA,"\\|")
r.group.ex = sapply(r.group, function(x){
  sapply(x, str_match, pattern = "[a-zA-Z-\\. /]*")
})

R.dict = dict.combine[dict.combine$Cap=="R",]
#R.dict.write = R.dict[,c(1,2)]
#R.dict.write$RXCUI = str_remove_all(R.dict$label,"R")
#write.csv(R.dict.write,"label data/R data.csv")
R.dict$group = sapply(r.group.ex, function(x){
  if(length(x)==1 & x[[1]]!=""){
    x
  }else{
    if(length(x)>1){
      x[which.min(str_length(x))]
    }else{
      NA
    }
  }
})

R.dict$Variable[is.na(R.dict$group)]
sum(is.na(R.dict$group))
a = table(R.dict$group)
print(length(a))
#write.csv(R.dict,"label data/R data.csv")
aa = paste0(R.dict$Variable[is.na(R.dict$group)],": ",
            R.dict$Description[is.na(R.dict$group)])
#write.csv(aa,"label data/R.csv")


Cap = substr(dict.combine$Variable,1,1)
dict.combine$group[Cap=="P"] = P.dict$group
dict.combine$group[Cap=="C"] = C.dict$group
dict.combine$group[Cap=="R"] = R.dict$group
dict.combine$group[Cap=="L"] = as.character(L.dict$group)
#dict.combine$group[Cap=="O"] = "VA Lab Code"
#dict.combine$group[Cap=="S"] = "VA Lab Group"

dict.combine$Description_s = dict.combine$Description_s_nosym = dict.combine$Description
dict.combine$Description_s[Cap=="L"] =
  dict.combine$Description_s_nosym[Cap=="L"] = sub(pattern = "(^.*group:)", "",
                                           dict.combine$Description[Cap=="L"],fixed = FALSE)



label = str_replace_all(dict.combine$Description_s[str_length(dict.combine$Description_s)>20], "\\.\\.\\."," ...")
label = str_replace_all(label,"/","/ ")
label_split = str_split(label," ")
label_com = sapply(label_split, function(x){
  y = ""
  i = 1
  k = 0
  while(i <= length(x)){
    y = paste(y, x[i])
    k = k + str_length(x[i])
    if(k>=15 & i!=length(x)){
      y = paste0(y,"\n")
      k = 0
    }
    i = i + 1
  }
  return(str_trim(y,side = "both"))
})


dict.combine$Description_s[str_length(dict.combine$Description_s)>20] = label_com

dict.combine$group_label = dict.combine$group
dict.combine$group_label[is.na(dict.combine$group_label)|dict.combine$group_label==" "] =
  dict.combine$Description_s[is.na(dict.combine$group_label)|dict.combine$group_label==" "]
dict.combine$group[is.na(dict.combine$group)|dict.combine$group==" "] =
  dict.combine$Description_s_nosym[is.na(dict.combine$group)|dict.combine$group==" "]

label = dict.combine$group[str_length(dict.combine$group_label)>20]
label = str_replace_all(label,"/","/ ")
label_split = str_split(label," ")
label_com = sapply(label_split, function(x){
  y = ""
  i = 1
  k = 0
  while(i <= length(x)){
    y = paste(y, x[i])
    k = k + str_length(x[i])
    if(k>=15 & i!=length(x)){
      y = paste0(y,"\n")
      k = 0
    }
    i = i + 1
  }
  return(str_trim(y,side = "both"))
})


dict.combine$group_label[str_length(dict.combine$group_label)>20] = label_com

group.total = dict.combine[,c("group","group_label","Cap")]
group.total = group.total[!duplicated(group.total$group),]


edge.full.list.2 = list(list(1,1),list(1,1),list(1,1),list(1,1),list(1,1))
for (i in 1:5) {
  edge.full.list.2[[i]][[1]] = edge.full.list[[i]]
  edge.full.list.2[[i]][[2]] = dict.combine$group[match(colnames(edge.full.list[[i]]),
                                                                   dict.combine$Variable)]
}
edge.full.list = edge.full.list.2

n.col = ncol(dict.combine)

dict.combine$degree.VA.inte = dict.combine$degree.RPDR.inte =
  dict.combine$degree.VA.local =  dict.combine$degree.RPDR.local =
  dict.combine$degree.both = 0

for (method in 1:5) {
  x = dict.combine$Variable
  x.neighbor = sapply(x, function(xx){
    if(xx %in% colnames(edge.full.list[[method]][[1]])){
      sum(edge.full.list[[method]][[1]][xx,]!=0) +
        sum(edge.full.list[[method]][[1]][,xx]!=0) -
        sum(edge.full.list[[method]][[1]][,xx] != 0 &
              edge.full.list[[method]][[1]][xx,] != 0)
    }else{
      0
    }
  })
  dict.combine[,method + n.col] = x.neighbor
}

group.total = dict.combine[,c("group","group_label","Cap")]
group.total = group.total[!duplicated(group.total$group),]


a<-ls()


rm(list=a[which(a!='edge.full.list' & a !='dict.combine' & a != 'edge.list' & a != 'group.total' & a != 'n.col')])
rm(a)


save.image("shiny/data/edge_matrix.RData")





###replicate
a = table(dict.combine$Description)
tail(sort(a),100)
b = dict.combine[dict.combine$Description=="protein",]
View(b)

save.image("network_shiny/data/edge_matrix.RData")
