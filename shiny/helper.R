
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
                                      color.color="CCCCCC",
                                      color.highlight="FFFF00",
                                      color.hover="FFFF00",
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
                                      color.color="CCCCCC",
                                      color.highlight="FFFF00",
                                      color.hover="FFFF00",
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
                                      color.color="CCCCCC",
                                      color.highlight="FFFF00",
                                      color.hover="FFFF00",
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
  
  title = dict.combine$Description_s[match(nodes$id,dict.combine$Variable)]
  

 
  nodes = cbind(nodes, data.frame(label = dict.combine$Description_s[match(nodes$id,dict.combine$Variable)],
                                  group = stringr::str_match(string = nodes$id, pattern = "(?<ID>.*):.*")[,2],
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
  group.legend[[1]] = list(label = "group", shape = "circle", 
                           color = "#D1BBFF", physics = FALSE,
                           font.size = 20,borderWidth=2)
  group.legend[[2]] = list(label = "*group", shape = "circle", 
                           color = "#9955FF", physics = FALSE,
                           font.size = 20,borderWidth=2)
  
  for (i in 1:length(group.names)) {
    group.legend[[i+2]] = list(label = group.names[i], shape = "ellipse", 
                             color = group.colors[i], physics = FALSE)
  }
  
  ledges <- data.frame(color = c("#66DD00", "#33CCFF", "#FFFF00"),
                      label = c("out", "in", "both"), 
                      arrows =c("to", "from",""), width = 3,
                      physics = FALSE,font.align = "top",font.size=20)  
  for (i in 1:length(group.names)) {
    nodes.colors[which(nodes$group==group.names[i])] = group.colors[i] 
  }
  nodes$color = nodes.colors
  nodes$color[1:n.root] = "#FF0000"
  nodes$group = node.group[match(nodes$id,colnames(edge.matrix.full))]
  group.target = unique(nodes$group[1:n.root])
  group.other = setdiff(unique(nodes$group),group.target)
  
  
  
  nodes$title = paste0("<b>",title,"</b>",
                 "<br>ID:",nodes$id,
                 "<br>Group:", nodes$group,
                 "<br>Degree:", dict.combine[match(nodes$id,dict.combine$Variable), method+5])
  nodes = rbind(nodes[((n.root+1):nrow(nodes)),],nodes[(1:n.root),])
  
  return(list(nodes,edges,list(group.legend,ledges),list(group.target,group.other)))
}

