rm(list=ls())
library("visNetwork")
library("RColorBrewer")
library(htmlwidgets)
library(base)
library(data.table)
library(dplyr)
library(network)
library(igraph)
library(igraphdata)
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


# create the edge adjacency matrix
dict.combine$label = sub(pattern = "(^\\D*)", "", dict.combine$Variable,fixed = FALSE)
Cap = substr(dict.combine$Variable,1,1)
dict.combine$label = paste0(Cap,dict.combine$label)
dict.combine$label[Cap=="S"] = sub(pattern = "(^\\D*[:])", "", dict.combine$Variable[Cap=="S"],fixed = FALSE)
n.node = nrow(dict.combine)

# upload network data
interest.list=list.files("network_data")
interest.list=gsub("tab_|\\[|\\]|.csv", "",interest.list)
interest.list=gsub("_", ":",interest.list)

# node = dict.combine$Variable
# a = node%in%interest.list
# sum(a)
name = dict.combine$Variable

add.edge <- function(interest, n.node, name){
  vec.all = vec.RPDR = vec.VA = matrix(0,1,n.node)
  colnames(vec.all) = colnames(vec.RPDR) = colnames(vec.VA) = name
  interest1=gsub(":", "_",interest)
  res = read.csv(paste0("network_data/tab_[",interest1,"].csv"))
  is.all = as.character(res$Variable[which(res$Integrative_Estimator_RPDR!=0 & res$Integrative_Estimator_VA!=0)])
  is.RPDR = as.character(res$Variable[which(res$Integrative_Estimator_RPDR!=0 & res$Integrative_Estimator_VA==0)])
  is.VA = as.character(res$Variable[which(res$Integrative_Estimator_RPDR==0 & res$Integrative_Estimator_VA!=0)])
  vec.all[,is.all] = vec.RPDR[,is.RPDR] = vec.VA[,is.VA] = 1
  a = list(list(1,2,3))
  a[[1]] = list(list(vec.all),list(vec.RPDR),list(vec.VA))
  return(a)
}


edge.matrix.list = sapply(interest.list, add.edge, n.node=n.node, name=name)

edge.matrix.all = sapply(edge.matrix.list, function(x){unlist(x[[1]])},
                         simplify = TRUE)
edge.matrix.RPDR = sapply(edge.matrix.list, function(x){unlist(x[[2]])},
                         simplify = TRUE)
edge.matrix.VA = sapply(edge.matrix.list, function(x){unlist(x[[3]])},
                          simplify = TRUE)
rownames(edge.matrix.all) = rownames(edge.matrix.RPDR) = rownames(edge.matrix.VA) = 
  dict.combine$Variable

edge.matrix.all = edge.matrix.all[apply(edge.matrix.all, 1, sum)!=0,]
edge.matrix.RPDR = edge.matrix.RPDR[apply(edge.matrix.RPDR, 1, sum)!=0,]
edge.matrix.VA = edge.matrix.VA[apply(edge.matrix.VA, 1, sum)!=0,]

edge.matrix.all = as(as.matrix(edge.matrix.all), "dgCMatrix")
edge.matrix.RPDR = as(as.matrix(edge.matrix.RPDR), "dgCMatrix")
edge.matrix.VA = as(as.matrix(edge.matrix.VA), "dgCMatrix")

edge.list = edge.full.list = list(edge.matrix.all, edge.matrix.RPDR, edge.matrix.VA)
for (i in 1:3) {
  edge.matrix = edge.list[[i]]
  nodes = unique(c(colnames(edge.matrix),rownames(edge.matrix)))
  n.nodes = length(nodes)
  edge.full.matrix = Matrix(0, n.nodes, n.nodes)
  colnames(edge.full.matrix) = rownames(edge.full.matrix) = nodes
  edge.full.matrix[match(colnames(edge.matrix), nodes),
                    match(rownames(edge.matrix), nodes)] = t(edge.matrix)
  
  # names = names.0 = colnames(edge.full.matrix)
  # names.1 = sub(pattern = "(^\\D*)", "", names,fixed = FALSE)
  # Cap = substr(names,1,1)
  # names.0 = paste0(Cap,names.1)
  # names.0[Cap=="S"] = sub(pattern = "(^\\D*[:])", "", names[Cap=="S"], fixed = FALSE)
  # colnames(edge.full.matrix) = names.0
  # 
  # names = names.0 = rownames(edge.full.matrix)
  # names.1 = sub(pattern = "(^\\D*)", "", names,fixed = FALSE)
  # Cap = substr(names,1,1)
  # names.0 = paste0(Cap,names.1)
  # names.0[Cap=="S"] = sub(pattern = "(^\\D*[:])", "", names[Cap=="S"], fixed = FALSE)
  # rownames(edge.full.matrix) = names.0
  # 
  edge.full.list[[i]] = edge.full.matrix
  
}


n.node = nrow(dict.combine)
set.seed(10101)
dict.combine$group = sample(1:10,n.node,replace = TRUE)
edge.full.list.2 = list(list(1,1),list(1,1),list(1,1))
for (i in 1:3) {
  edge.full.list.2[[i]][[1]] = edge.full.list[[i]]
  edge.full.list.2[[i]][[2]] = paste0("G",dict.combine$group[match(colnames(edge.full.list[[i]]),
                                                                   dict.combine$Variable)])
}
edge.full.list = edge.full.list.2


a<-ls()
rm(list=a[which(a!='edge.full.list' & a !='dict.combine' & a != 'edge.list')])
rm(a)


method=3
root.node=c(10,20,30,40,50)

vis.InOutNodes <- function(root.node, method){
  n.root = length(root.node)
  edge.matrix = t(edge.list[[method]])[root.node,]
  edge.matrix.full = edge.full.list[[method]]
  edges = data.frame()
  for (i in 1:nrow(edge.matrix)) {
    to=from = rownames(edge.matrix)[i]
    newto = colnames(edge.matrix)[which(edge.matrix[i,]==1)]
    i.loc = root.node[i]
    newfrom = rownames(edge.matrix.full)[which(edge.matrix.full[,i.loc]==1)]
    newtofrom = intersect(newto,newfrom)
    newto = setdiff(newto, newtofrom)
    newfrom = setdiff(newfrom,newtofrom)
    if(length(newto)!=0){
      edges = rbind(edges, data.frame(from=as.character(rep(from,length(newto))),
                                      to=as.character(newto),
                                      color.color="#66DD00",
                                      color.highlight="#00AA00",
                                      color.hover="#00AA00",
                                      #arrows="to",
                                      length=250,
                                      hover=TRUE,
                                      width = 1,
                                      title = paste0(as.character(rep(from,length(newto))),"<b> &rarr; </b>",
                                                     as.character(newto))))
      
    }
    if(length(newfrom)!=0){
      edges = rbind(edges, data.frame(from=as.character(newfrom),
                                      to=as.character(rep(to,length(newfrom))),
                                      color.color="#33CCFF",
                                      color.highlight="#0066FF",
                                      color.hover="#0066FF",
                                      #arrows="to",
                                      length=350,
                                      hover=TRUE,
                                      width=1,
                                      title = paste0(as.character(newfrom),"<b> &rarr; </b>",
                                                     as.character(rep(to,length(newfrom))))))
      
    }
    if(length(newtofrom)!=0){
      edges = rbind(edges, data.frame(from=as.character(newtofrom),
                                      to=as.character(rep(to,length(newtofrom))),
                                      color.color="#FFFF00",
                                      color.highlight="#FFBB00",
                                      color.hover="#FFBB00",
                                      #arrows="middle",
                                      length=150,
                                      hover=TRUE,
                                      width=2,
                                      title= paste0(as.character(newtofrom),"<b> &harr; </b>",
                                                    as.character(rep(to,length(newtofrom))))))
      
    }
  }
  
  nodes = data.frame(id = unique(c(rownames(edge.matrix),
                                   as.character(edges$from),
                                   as.character(edges$to))))
  n.nonroot = nrow(nodes) - n.root
  nodes$id = as.character(nodes$id)
  DegreeofNodes = cbind(rep(0,nrow(nodes)),rep(0,nrow(nodes)))  #1:from(out degree); 2:to(in degree)
  DegreeofNodes[,1] = apply(edge.matrix.full,1,sum)[match(nodes$id,rownames(edge.matrix.full))]
  DegreeofNodes[,2] = apply(edge.matrix.full,2,sum)[match(nodes$id,colnames(edge.matrix.full))]
  
  
  title = dict.combine$Description[match(nodes$id,dict.combine$Variable)]
  
  nodes = cbind(nodes, data.frame(title = paste0("<b>",title,"</b>",
                                                 "<br>",nodes$id,
                                                 "<br>Out degree:",DegreeofNodes[,1],
                                                 "<br>In degree:",DegreeofNodes[,2]),
                                  label = dict.combine$label[match(nodes$id,dict.combine$Variable)],
                                  group = sub(pattern = "([:.0-9-]*$)", "", nodes$id,fixed = FALSE),
                                  shape = c(rep("star",n.root),rep("box",n.nonroot)),
                                  font.size = c(rep(40,n.root),rep(15,n.nonroot)),
                                  font.background = c(rep("#FFFFF0",n.root),rep("",n.nonroot)),
                                  borderWidth = c(rep(2,n.root),rep(1,n.nonroot)),
                                  borderWidthSelected = c(rep(3,n.root),rep(2,n.nonroot))))
  nodes$group = as.character(nodes$group)
  group.colors = c("#87CEFA","#98F898","#FFB6C1","#FFD700","#EE82EE","#D2B48C") #blue green pink yellow purple tan
  group.names = sort(as.character(unique(nodes$group)))
  nodes.colors = nodes$id
  group.legend = list()
  
  for (i in 1:length(group.names)) {
    group.legend[[i]] = list(label = group.names[i], shape = "ellipse", 
                                        color = group.colors[i])
  }

  ledges <- data.frame(color = c("#66DD00", "#33CCFF", "#FFFF00"),
                       label = c("out", "in", "both"), arrows =c("to", "from",""),
                       width = 2)  
  for (i in 1:length(group.names)) {
    nodes.colors[which(nodes$group==group.names[i])] = group.colors[i] 
  }
  nodes$color = nodes.colors
  nodes$color[1:n.root] = "#FF0000"
  
  nodes = rbind(nodes[((n.root+1):nrow(nodes)),],nodes[(1:n.root),])
  
  visNetwork(nodes, edges, width = "100%") %>%
    visNodes(color = list(background = "lightblue",
                          border = "darkblue",
                          highlight = "yellow"),
             shadow = list(enabled = TRUE, size = 10)) %>%
    visEdges(
      #physics = FALSE, 
             smooth = FALSE,
             hoverWidth = 2.5) %>%
    visOptions(highlightNearest = 
                 list(enabled = T, degree = 1, hover = T),
               selectedBy = "group",
               collapse = TRUE) %>%
    visLegend(width = 0.1, position = "right", main = "Group")  %>%
    visLegend(addNodes = group.legend, addEdges = ledges,
      useGroups = FALSE) %>%
    visInteraction(hover = TRUE) %>%
    visLayout(randomSeed = 10) # to have always the same network
}


vis.InOutNodes(root.node=c(10,20,50,15,19,300), method=1)


ma = matrix(c(1:6),2,3,byrow = FALSE)
layout(ma)
par(pin = c(1.7,1.7),mai=c(.5,.5,.3,.3),
    omi = c(0, 0, 1, 0),mgp=c(2,1,0))
for (i in 1:3) {
  a = apply(edge.full.list[[i]],1,sum)
  hist(a)
  a = apply(edge.full.list[[i]],2,sum)
  hist(a)
}

library("Matrix")
edge.matrix.full <- as(as.matrix(edge.full.list[[1]]), "dgCMatrix")

 # nodes.neighbor <- function(g,nodes){
#   neighbors = c()
#   for (i in 1:length(nodes)) {
#     neighbors = c(neighbors, get.neighborhood(g,i,"combined"))
#   }
#   return(unique(neighbors))
# }
# 
# 
# nodes.diglevel <- function(root.node, level, g){
#   neibors = c()
#   sub.nodes.all = c(root.node)
#   for (node.id in root.node) {
#     neibors = get.neighborhood(g,node.id,"combined")
#     sub.nodes = c(node.id)
#     sub.nodes.new = unique(c(sub.nodes, neibors))
#     if(level > 1){
#       for (j in 2:level) {
#         new.nodes = setdiff(sub.nodes.new, sub.nodes)
#         sub.nodes = sub.nodes.new
#         neibors = nodes.neighbor(g, new.nodes)
#         sub.nodes.new = unique(c(sub.nodes, neibors))
#       }
#     }
#     sub.nodes.all = unique(c(sub.nodes.all,sub.nodes.new))
#   }
#   return(sort(sub.nodes.all))
# }
# 
# a=nodes.diglevel(1:2,3,g)
