rm(list=ls())
library("visNetwork")
library("Matrix")

mainpath = "D:/xxiong/program/network/network.vis-main"
setwd(mainpath)
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

dict.combine$label = sub(pattern = "(^\\D*)", "", dict.combine$Variable,fixed = FALSE)
Cap = substr(dict.combine$Variable,1,1)
dict.combine$label = paste0(Cap,dict.combine$label)
dict.combine$label[Cap=="S"] = sub(pattern = "(^\\D*[:])", "", dict.combine$Variable[Cap=="S"],fixed = FALSE)
n.node = nrow(dict.combine)

# update description
dict.combine$Description_s = dict.combine$Description
dict.combine$Description_s[Cap=="L"] = sub(pattern = "(^.*group:)", "", 
                                           dict.combine$Description[Cap=="L"],fixed = FALSE)

# upload network data
interest.list=list.files("network_data")
interest.list=gsub("tab_|\\[|\\]|.csv", "",interest.list)
interest.list=gsub("_", ":",interest.list)
name = dict.combine$Variable

add.edge <- function(interest, n.node, name){
  vec.all = vec.RPDR.local = vec.RPDR.inte = 
    vec.VA.local = vec.VA.inte = matrix(0,1,n.node)
  colnames(vec.all) = colnames(vec.RPDR.local) = colnames(vec.VA.local) =
    colnames(vec.RPDR.inte) = colnames(vec.VA.inte) = name
  interest1=gsub(":", "_",interest)
  res = read.csv(paste0("network_data/tab_[",interest1,"].csv"))
  is.all = as.character(res$Variable[which(res$Integrative_Estimator_RPDR!=0 | res$Integrative_Estimator_VA!=0)])
  is.RPDR.local = as.character(res$Variable[which(res$Local_Estimator_RPDR!=0)])
  is.RPDR.inte = as.character(res$Variable[which(res$Integrative_Estimator_RPDR!=0)])
  is.VA.local = as.character(res$Variable[which(res$Local_Estimator_VA!=0)])
  is.VA.inte = as.character(res$Variable[which(res$Integrative_Estimator_VA!=0)])
  vec.all[,is.all] = vec.RPDR.local[,is.RPDR.local] = vec.RPDR.inte[,is.RPDR.inte] = 
    vec.VA.local[,is.VA.local] = vec.VA.inte[,is.VA.inte] = 1
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

edge.matrix.all = edge.matrix.all[apply(edge.matrix.all, 1, sum)!=0,]
edge.matrix.RPDR.local = edge.matrix.RPDR.local[apply(edge.matrix.RPDR.local, 1, sum)!=0,]
edge.matrix.VA.local = edge.matrix.VA.local[apply(edge.matrix.VA.local, 1, sum)!=0,]
edge.matrix.RPDR.inte = edge.matrix.RPDR.inte[apply(edge.matrix.RPDR.inte, 1, sum)!=0,]
edge.matrix.VA.inte = edge.matrix.VA.inte[apply(edge.matrix.VA.inte, 1, sum)!=0,]

edge.matrix.all = as(as.matrix(edge.matrix.all), "dgCMatrix")
edge.matrix.RPDR.local = as(as.matrix(edge.matrix.RPDR.local), "dgCMatrix")
edge.matrix.VA.local = as(as.matrix(edge.matrix.VA.local), "dgCMatrix")
edge.matrix.RPDR.inte = as(as.matrix(edge.matrix.RPDR.inte), "dgCMatrix")
edge.matrix.VA.inte = as(as.matrix(edge.matrix.VA.inte), "dgCMatrix")
edge.list = edge.full.list = list(edge.matrix.all, edge.matrix.RPDR.local, edge.matrix.VA.local,
                                  edge.matrix.RPDR.inte, edge.matrix.VA.inte)
for (i in 1:5) {
  edge.matrix = edge.list[[i]]
  nodes = unique(c(colnames(edge.matrix),rownames(edge.matrix)))
  n.nodes = length(nodes)
  edge.full.matrix = Matrix(0, n.nodes, n.nodes)
  colnames(edge.full.matrix) = rownames(edge.full.matrix) = nodes
  edge.full.matrix[match(colnames(edge.matrix), nodes),
                    match(rownames(edge.matrix), nodes)] = t(edge.matrix)
  edge.full.list[[i]] = edge.full.matrix
  
}


n.node = nrow(dict.combine)
set.seed(10101)
dict.combine$group = sample(1:30,n.node,replace = TRUE)
edge.full.list.2 = list(list(1,1),list(1,1),list(1,1),list(1,1),list(1,1))
for (i in 1:5) {
  edge.full.list.2[[i]][[1]] = edge.full.list[[i]]
  edge.full.list.2[[i]][[2]] = paste0("G",dict.combine$group[match(colnames(edge.full.list[[i]]),
                                                                   dict.combine$Variable)])
}
edge.full.list = edge.full.list.2

dict.combine$degree.VA.inte = dict.combine$degree.RPDR.inte = 
  dict.combine$degree.VA.local =  dict.combine$degree.RPDR.local = 
  dict.combine$degree.both = 0

for (method in 1:5) {
  x = dict.combine$Variable
  x.neighbor = sapply(x, function(xx){
    if(xx %in% colnames(edge.full.list[[method]][[1]])){
      sum(edge.full.list[[method]][[1]][xx,]==1) + 
        sum(edge.full.list[[method]][[1]][,xx]==1) - 
        sum(edge.full.list[[method]][[1]][,xx] == 1 & 
              edge.full.list[[method]][[1]][xx,] == 1)
    }else{
      0
    }
  })
  dict.combine[,method + 5] = x.neighbor
}




a<-ls()
rm(list=a[which(a!='edge.full.list' & a !='dict.combine' & a != 'edge.list')])
rm(a)
mainpath = "D:/xxiong/Github/network.vis"
setwd(mainpath)
save.image("shiny/data/edge_matrix.RData")





###replicate
a = table(dict.combine$Description)
tail(sort(a),100)
b = dict.combine[dict.combine$Description=="protein",]
View(b)
