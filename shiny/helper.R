vis.InOutNodes <- function(root.node, method, edge.list, edge.full.list, dict.combine){
  n.root = length(root.node)
  if(n.root==1){
    edge.matrix = t(as.matrix(t(edge.list[[method]])[root.node,]))
    rownames(edge.matrix) = rownames(edge.full.list[[method]][[1]])[root.node]
  }else{
    edge.matrix = t(edge.list[[method]])[root.node,]
  }
  edge.matrix.full = edge.full.list[[method]][[1]]
  node.group = edge.full.list[[method]][[2]]
  node.used = c()
  edges.dataframe = data.frame()
  for (i in 1:nrow(edge.matrix)) {
    edges = data.frame()
    to = from = rownames(edge.matrix)[i]
    newto = colnames(edge.matrix)[which(edge.matrix[i,]!=0)]
    i.loc = root.node[i]
    newfrom = rownames(edge.matrix.full)[which(edge.matrix.full[,i.loc]!=0)]
    newto = setdiff(newto, node.used)
    newfrom = setdiff(newfrom, node.used)
    newtofrom = intersect(newto,newfrom)
    newtofrom.value = abs(edge.matrix.full[newtofrom,i.loc])
    newtofrom[rank(newtofrom.value)] = newtofrom
    newtofrom.value = sort(newtofrom.value)
    
    newto = setdiff(newto, newtofrom)
    newto.value = abs(edge.matrix.full[i.loc,newto])
    newto[rank(newto.value)] = newto
    newto.value = sort(newto.value)
    
    newfrom = setdiff(newfrom,newtofrom)
    newfrom.value = abs(edge.matrix.full[newfrom,i.loc])
    newfrom[rank(newfrom.value)] = newfrom
    newfrom.value = sort(newfrom.value)
    length.total = c()
    if(length(newto)!=0){
      #color.opa = round(ifelse((newto.value-min(newto.value))*20+0.1>1,1,
      #                   (newto.value-min(newto.value))*20+0.1),2)
      #length = abs(newto.value)^(-1.15)*1
      length = abs(newto.value)^(-1.1)*0.8
      length = length - min(length) + 5
      a.min = min(500,quantile(length,.85))
      a.max = ifelse(a.min==500,600,min(600,quantile(length,.95)))
      length[length>a.min] = sort(runif(sum(length>a.min), min = a.min, max = a.max),decreasing = TRUE)
      a.min = 5
      a.max = min(400,quantile(length,.6))
      length[length<a.max] = sort(runif(sum(length<a.max), min = a.min, max = a.max),decreasing = TRUE)
      length.total = c(length.total,length)
      edges = rbind(edges, data.frame(from=as.character(rep(from,length(newto))),
                                      to=as.character(newto),
                                      length=length,
                                      title = paste0(as.character(rep(from,length(newto))),"<b> &rarr; </b>",
                                                     as.character(newto))))
      
    }
    if(length(newfrom)!=0){
      #color.opa = round(ifelse((newfrom.value-min(newfrom.value))*20+0.1>1,1,
      #                         (newfrom.value-min(newfrom.value))*20+0.1),2)
      #length=abs(newfrom.value)^(-1.15)*1
      length=abs(newfrom.value)^(-1.1)*0.8
      length = length - min(length) + 5
      a.min = min(500,quantile(length,.85))
      a.max = ifelse(a.min==500,600,min(600,quantile(length,.95)))
      length[length>a.min] = sort(runif(sum(length>a.min), min = a.min, max = a.max),decreasing = TRUE)
      a.min = 5
      a.max = min(400,quantile(length,.6))
      length[length<a.max] = sort(runif(sum(length<a.max), min = a.min, max = a.max),decreasing = TRUE)
      length.total = c(length.total,length)
      edges = rbind(edges, data.frame(from=as.character(newfrom),
                                      to=as.character(rep(to,length(newfrom))),
                                      length=length,
                                      title = paste0(as.character(newfrom),"<b> &rarr; </b>",
                                                     as.character(rep(to,length(newfrom))))))
      
    }
    if(length(newtofrom)!=0){
      #color.opa = round(ifelse((newtofrom.value-min(newtofrom.value))*10+0.1>1,1,
      #                         (newtofrom.value-min(newtofrom.value))*10+0.1),2)
      #length=abs(newtofrom.value)^(-1.15)*1
      length=abs(newtofrom.value)^(-1.1)*0.8
      length = length - min(length) + 5
      a.min = min(500,quantile(length,.85))
      a.max = ifelse(a.min==500,600,min(600,quantile(length,.95)))
      length[length>a.min] = sort(runif(sum(length>a.min), min = a.min, max = a.max),decreasing = TRUE)
      a.min = 5
      a.max = min(400,quantile(length,.6))
      length[length<a.max] = sort(runif(sum(length<a.max), min = a.min, max = a.max),decreasing = TRUE)
      length.total = c(length.total,length)
      edges = rbind(edges, data.frame(from=as.character(newtofrom),
                                      to=as.character(rep(to,length(newtofrom))),
                                      length=length,
                                      title= paste0(as.character(newtofrom),"<b> &harr; </b>",
                                                    as.character(rep(to,length(newtofrom))))))
      
    }
    node.used = c(node.used, rownames(edge.matrix)[i])
    color.opa = round(ifelse((1-(length.total-min(length.total))/
                         (max(length.total)-min(length.total)))^5+0.2>1,1,
                       (1-(length.total-min(length.total))/
                         (max(length.total)-min(length.total)))^5+0.2),2)
    plot(length.total,color.opa)
    edges = cbind(edges, data.frame(color.color=paste0("rgba(128,128,128,", color.opa,")"),
                                    color.opa = color.opa,
                                    color.highlight="rgba(205,104,57)",
                                    color.hover="rgba(205,104,57)",
                                    smooth = TRUE,
                                    hover=TRUE,
                                    width=1.5,
                                    hoverWidth=2.5,
                                    selectionWidth=2.5))
    
    edges.dataframe = rbind(edges.dataframe,edges)
  }
  edges = edges.dataframe
  a = edges$length[edges$from %in% rownames(edge.matrix) &
                     edges$to %in% rownames(edge.matrix)]
  edges$length[edges$from %in% rownames(edge.matrix) &
               edges$to %in% rownames(edge.matrix)] = sapply(a, function(x){max(x,100*min(10,nrow(edge.matrix)))})
  edges$width[edges$from %in% rownames(edge.matrix) &
                 edges$to %in% rownames(edge.matrix)] = 4
  edges$hoverWidth[edges$from %in% rownames(edge.matrix) &
                edges$to %in% rownames(edge.matrix)] = 
    edges$selectionWidth[edges$from %in% rownames(edge.matrix) &
                       edges$to %in% rownames(edge.matrix)] = 4
  edges$color.color = as.character(edges$color.color)
  edges$color.highlight = as.character(edges$color.highlight)
  edges$color.hover = as.character(edges$color.hover)
  edges$color.color[edges$from %in% rownames(edge.matrix) &
                     edges$to %in% rownames(edge.matrix)] = "#0000CD"
  edges$color.hover[edges$from %in% rownames(edge.matrix) &
                      edges$to %in% rownames(edge.matrix)] = "#0000CD"
  edges$color.highlight[edges$from %in% rownames(edge.matrix) &
                      edges$to %in% rownames(edge.matrix)] = "#0000CD"
  edges$smooth[edges$from %in% rownames(edge.matrix) &
                      edges$to %in% rownames(edge.matrix)] = FALSE
  nodes = data.frame(id = unique(c(rownames(edge.matrix),
                                   as.character(edges$from),
                                   as.character(edges$to))))
  n.nonroot = nrow(nodes) - n.root
  nodes$id = as.character(nodes$id)
  
  title = dict.combine$Description[match(nodes$id,dict.combine$Variable)]
  


  nodes = cbind(nodes, data.frame(label = dict.combine$Description_s[match(nodes$id,dict.combine$Variable)],
                                  group = dict.combine$group_label[match(nodes$id,dict.combine$Variable)],
                                  Cap = dict.combine$Cap[match(nodes$id,dict.combine$Variable)],
                                  shape = c(rep("star",n.root),rep("box",n.nonroot)),
                                  mass = c(rep(7,n.root),rep(3,n.nonroot)),
                                  font.size = c(rep(30,n.root),rep(22,n.nonroot)),
                                  size =  c(rep(30,n.root),rep(20,n.nonroot)),
                                  font.background = c(rep("#FFFFF0",n.root),rep("",n.nonroot)),
                                  borderWidth = c(rep(2,n.root),rep(1,n.nonroot)),
                                  borderWidthSelected = c(rep(3,n.root),rep(2,n.nonroot))))
  
  nodes$group = as.character(nodes$group)
  nodes$group[1:n.root] = as.character(nodes$label[1:n.root])
  group.colors = c("#87CEFA","#98F898","#FFB6C1","#FFD700","#EE82EE","#D2B48C") #blue green pink yellow purple tan
  group.colors.rgb = c("rgba(135,206,250,","rgba(152,248,152,","rgba(255,182,193,",
                       "rgba(255,215,0,","rgba(238,130,238,","rgba(210,180,140,")
  group.names = c("C","L","O","P","R","S")
  group.names.long = c("CCS","LOINC"," VA Lab Code","PheCode","RXNORM","VA Lab Group")
  nodes.colors = nodes.colors.select = nodes$id
  group.legend = list(1,1,1,1,1,1,
                      1,1,1,1,1,1)
  group.legend[[1]] = list(label = "Node:",shape="box",color="rgba(0,0,0,0)",
                           physics = FALSE,size=10)
  k=1
  for (i in c(1,2,4,5,3,6)) {
    group.legend[[k+1]] = list(label = group.names.long[i], shape = "box", 
                             color = group.colors[i], physics = FALSE)
    k = k + 1
  }
  group.legend[[8]] = list(label = "Group:",shape="box",color="rgba(0,0,0,0)",
                           physics = FALSE,size=10)
  k=1
  for (i in c(1,2,4,5,3,6)) {
    group.legend[[k+8]] = list(label = group.names.long[i], shape = "ellipse",
                               color = group.colors[i], physics = FALSE)
    k = k + 1
  }
  ledges <- data.frame(color = c("#0000CD","rgba(128,128,128)","rgba(205,104,57)"),
                      label = c("target-target", "target-other","target-other\n(selected)"), 
                      width = 3, arrows = c("","",""),
                      physics = FALSE,font.align = c("top","top","horizontal"),font.size=15)  
  for (i in 1:length(group.names)) {
    nodes.colors[which(nodes$Cap==group.names[i])] = group.colors.rgb[i]
    nodes.colors.select[which(nodes$Cap==group.names[i])] = group.colors[i]
  }

  nodes$color.opa = sapply(1:length(nodes$id), function(i){max(edges$color.opa[nodes$id[i]==edges$from],
                                                               edges$color.opa[nodes$id[i]==edges$to])})
  nodes$color.background = paste0(nodes.colors,nodes$color.opa,")")
  nodes$color.highlight.background = paste0(nodes.colors,"1)")
  nodes$color.hover.background = paste0(nodes.colors,"1)")
  
  nodes$color.border = nodes$color.background
  nodes$color.highlight.border = nodes$color.highlight.background
  nodes$color.hover.border = nodes$color.hover.background
  
  #nodes$color.background[1:n.root] = nodes$color.highlight.background[1:n.root] =
  #  nodes$color.hover.background[1:n.root] =  nodes$color.border[1:n.root] = "#FF0000"
  nodes$color.highlight.border[1:n.root] = nodes$color.hover.border[1:n.root] = "#FF0000"
 
  nodes$font.size = sapply(nodes$color.opa*20,function(x){max(min(x,18),5)})
  nodes$font.size[str_length(nodes$label)>40] = 
    nodes$font.size[str_length(nodes$label)>40]*40/
    str_length(nodes$label[str_length(nodes$label)>40])
  nodes$font.size[str_length(nodes$label)>40] = sapply(nodes$font.size[str_length(nodes$label)>40],
                                                       function(x){max(x,5)})
  nodes$font.size[1:n.root] = 20
  
  group.target = unique(nodes$group[1:n.root])
  group.other.C = setdiff(unique(nodes$group[nodes$Cap=="C"]),group.target)
  group.other.L = setdiff(unique(nodes$group[nodes$Cap=="L"]),group.target)
  group.other.O = setdiff(unique(nodes$group[nodes$Cap=="O"]),group.target)
  group.other.P = setdiff(unique(nodes$group[nodes$Cap=="P"]),group.target)
  group.other.R = setdiff(unique(nodes$group[nodes$Cap=="R"]),group.target)
  group.other.S = setdiff(unique(nodes$group[nodes$Cap=="S"]),group.target)
  group.other = list(group.other.C,group.other.L,group.other.O,
                     group.other.P,group.other.R,group.other.S)
  text.group.color = list(NA,NA,NA,NA,NA,NA)
  for (i in 1:6) {
    if(length(group.other[[i]])!=0){
      text.group.color[[i]] = paste0(group.colors.rgb[i], 
                                     sapply(group.other[[i]], function(group){max(nodes$color.opa[nodes$group==group])}),")")
      
    }
   }
  group.other = unlist(group.other)
  text.group.color = na.omit(unlist(text.group.color))
  nodes$title = paste0("<b>",title,"</b>",
                 "<br>ID:",nodes$id,
                 "<br>Group:", nodes$group,
                 "<br>Degree:", dict.combine[match(nodes$id,dict.combine$Variable), method+8])
  #nodes = rbind(nodes[((n.root+1):nrow(nodes)),],nodes[(1:n.root),])
  
  return(list(nodes, edges, list(group.legend,ledges), list(group.target, group.other, text.group.color)))
}

