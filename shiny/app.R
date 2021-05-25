library(visNetwork)
library(shiny)
library(shinyBS)
library(shinyWidgets)
library(Matrix)
library(shinydashboard)
library(shinycssloaders)
library(DT)
library(igraph)
library(stringr)
load("data/edge_matrix.RData")
source("helper.R")


ui <- dashboardPage(
  dashboardHeader(title = "visNetwork",
                  disable = FALSE
  ),
  ## Sidebar content
  dashboardSidebar(
    width = "200pt",
    disable = FALSE,
    sidebarMenu(
      menuItem("Network", tabName = "Network", icon = icon("dashboard")),
      selectInput("select", label = "Select data from:", 
                  choices = list("Both" = 1, "RPDR local" = 2,"VA local" = 3,
                                 "RPDR integrative" = 4,"VA integrative" = 5), 
                  selected = 1,width = '100%'),
      
      checkboxGroupInput("inCheckboxGroup2", "Candidates",
                         c("")),
    fluidRow(column(1),
      column(3,actionButton("goButton", "Show", 
                                   icon = tags$i(class = "far fa-play-circle",
                                                 style="font-size: 10px"), 
                                   class = "btn-success")),
             column(3,actionButton('refresh', 'Unselect', 
                                   icon = tags$i(class = "fa fa-refresh", 
                                                 style="font-size: 10px")))), 
    checkboxInput("cluster", "Cluster by groups*", value = FALSE),
    conditionalPanel(condition =  "input.current_node_id!=null",
                     hr(),
                     h5("Clicked node info:"),
                     tableOutput("shiny_return_subnode_table"),
                     actionButton("addButton", "Add to candidates",class="btn-primary active"),
                     br(),
                     h5("Connected node(s) info:"),
                     tableOutput("shiny_returntable")),
    conditionalPanel(condition =  "input.current_node_id==null",
                     hr(),
                     h4(" Try to click on a non-center node."))

    )
  ),
  
  dashboardBody(
    bsModal(
      id = "instruction", title = "Instruction", trigger = FALSE,
      size = "small",
      p('1. Please refer to "Possible input" box to see appropriate values.'),
      p('2. If you select "cluster by groups", double click any group node (ellipse) to unfold or
                                  any expanded individual node (square) to fold the group it belongs to. 
                                  Click "Reinitialize clustering" botton to fold all groups.'),
      p('3. Click any individual node (ellipse) to see its neighbors.'),
      p('4. The maximum number of input nodes has been set to be 50. However, less than 10 nodes
                    are recommended considering the possessing time and network density.')
    ),
    tabItems(
      tabItem(tabName = "Network",
              fluidRow(
                column(1,
                       dropdown(
                         h4("Possible inputs"),
                         DT::dataTableOutput("table"),
                         style = "unite", icon = icon("table"),
                         status = "primary", width = "350px",
                         animate = animateOptions(
                           enter = animations$fading_entrances$fadeInDown,
                           exit = animations$fading_exits$fadeOutUp,
                           duration = 0
                         )
                       )
                ),
              column(1,       
                     dropdown(
                       selectInput("Focus", label = "Focus on node:",
                                   choices = NA,width = "100%"),
                       sliderInput("scale_id","Focus scale:", min=1,max=10,value=5,width = "100%"),
                       sliderInput("slider_h", "Graph height:",
                                   min = 100, max = 1500, value = 750,width = "100%"),
                       style = "unite", icon = icon("eye"),right = FALSE,
                       status = "primary", width = "200px",
                       animate = animateOptions(
                         enter = animations$fading_entrances$fadeInDown,
                         exit = animations$fading_exits$fadeOutUp,
                         duration = 0
                       ))
                     ),
              column(8),
              column(1,downloadBttn(outputId = "downloadData", 
                                      label = "Download nodes as .csv",
                                      style = "material-circle",
                                      color = "default")),
              column(1,actionBttn("instruct", "Help", style = "material-circle",
                                    icon = icon("question")))
              ), 
              fluidRow(
                column(width = 12,
                       uiOutput("network")
                )
                
              )

      )
    )
  )
)





server <- function(input, output, session) {
  method = reactive({input$select})
  cluster = reactive({input$cluster})
  
  observeEvent(input$instruct, {
    toggleModal(session, "instruction", toggle = "open")
  })
  
  output$network <- renderUI({

    shinycssloaders::withSpinner(
      visNetworkOutput("network_proxy_nodes", 
                       height =  paste0(input$slider_h,"px")), type = 6
    )
  })

  output$table <- DT::renderDataTable(DT::datatable({
    #method = as.numeric(method())
    data.frame("Node id"=colnames(edge.list[[1]]),
               "Description"=dict.combine$Description[match(colnames(edge.list[[1]]),dict.combine$Variable)])
    #"Degree"=dict.combine[match(colnames(edge.list[[1]]),dict.combine$Variable),method+n.col],
    #"Group"=dict.combine$group[match(colnames(edge.list[[1]]),dict.combine$Variable)])
    
  }, rownames = FALSE,options = list(
    pageLength = 5
  )))
  
  proxy = dataTableProxy('table')
  
  observe({
    input$refresh
    isolate({
      reloadData(
        proxy,
        resetPaging = TRUE,
        clearSelection = c("all"))
    })
  })
  
  
  observe({
    s = input$table_rows_selected
    if(length(s)!=0){
      method = as.numeric(method())
      x = colnames(edge.list[[method]])[s[seq(1,min(50,length(s)),by=1)]]
      x.neighbor = sapply(x, function(xx){
        sum(edge.full.list[[method]][[1]][xx,]!=0) + 
          sum(edge.full.list[[method]][[1]][,xx]!=0) - 
          sum(edge.full.list[[method]][[1]][,xx] != 0 & 
                edge.full.list[[method]][[1]][xx,] != 0)
      })
      x.name = dict.combine$Description[match(x,dict.combine$Variable)]
      x.neighbor = paste0(x.name," (" ,x.neighbor," degrees)")
    }else{
      x = x.name = x.neighbor = character(0)  # Can use character(0) to remove all choices
    }
    
    # Can also set the label and select items
    updateCheckboxGroupInput(session, "inCheckboxGroup2",
                             label = paste(length(x), " candidate nodes:"),
                             choiceValues = x,
                             choiceNames = x.neighbor,
                             selected = x
    )
  })
  
  
  id <- NULL
  
  observeEvent(input$goButton, {
    # If there's currently a notification, don't add another
    if (!is.null(id))
      return()
    if(length(input$inCheckboxGroup2)>=10){
      id <<- showNotification(paste("You choose ", length(input$inCheckboxGroup2)," nodes. It will take a while to finish plotting..."), 
                              duration = max(5,3*length(input$inCheckboxGroup2)), type = "message")
      
    }
    # Save the ID for removal later
  })
  
  output$network_proxy_nodes <- renderVisNetwork({
    input$goButton
    s = isolate(input$inCheckboxGroup2)
    if(length(s)!=0){
      method = as.numeric(method())
      input.correct = s[seq(1,min(50,length(s)),by=1)]
      root.node = match(input.correct,rownames(edge.full.list[[method]][[1]]))
      # root.node = 1:5
      draw.data = vis.InOutNodes(root.node=root.node, method=method,
                                 edge.list, edge.full.list, dict.combine)
      nodes = draw.data[[1]]
      edges = draw.data[[2]]
      group.legend = draw.data[[3]][[1]]
      ledges = draw.data[[3]][[2]]
      group.target = draw.data[[4]][[1]]
      group.other = draw.data[[4]][[2]]
      text.group.color = draw.data[[4]][[3]]
      if(cluster()==TRUE){
        edges$length[edges$length<300] =  as.numeric(edges$length[edges$length<300])*2
        if(sum(edges$length>400)>0){
          edges$length[edges$length>400] =  runif(sum(edges$length>400),400,450)
        }
        nodes$mass[1:length(root.node)]=40
        a = edges$length[edges$from %in% nodes$id[1:length(root.node)] &
                         edges$to %in% nodes$id[1:length(root.node)]]
        edges$length[edges$from %in% nodes$id[1:length(root.node)] &
                       edges$to %in% nodes$id[1:length(root.node)]] = sapply(a, function(x){max(x,300*min(10,length(root.node)))})
        
        visNetwork(nodes, edges, width = "100%",height = "100%") %>%
          visNodes(color = list(background = "lightblue",
                                border = "darkblue",
                                highlight = "yellow"),
                   shadow = list(enabled = TRUE, size = 10)) %>%
          visEdges(
            physics = TRUE,
            smooth = FALSE,
            hoverWidth = 2.5) %>%
          visOptions(highlightNearest =
                       list(enabled = T, degree = 1, hover = T,
                            hideColor = "rgba(200,200,200,0.2)"),
                     selectedBy = list(`variable` = "Cap_label",
                                       `multiple` = TRUE,
                                       `main` = "Select by group"),
                     #clickToUse = TRUE,
                     collapse = FALSE) %>%
          visLegend(width = 0.09, position = "right",
                    addNodes = group.legend,
                    addEdges = ledges,
                    useGroups = FALSE, zoom = FALSE,
                    stepX = 150, stepY = 75,ncol=1) %>%
          visInteraction(hover = TRUE,
                         navigationButtons = TRUE) %>%
          visPhysics(barnesHut = list("avoidOverlap"=0.1)) %>%
          visIgraphLayout(layout = "layout_nicely",physics = TRUE,
                          smooth = TRUE) %>%
          
          visEvents(selectNode = "function(nodes) {
              Shiny.onInputChange('current_node_id', nodes);
              ;}") %>%

          visClusteringByGroup(groups = c(group.target,group.other),
                               label = "Group:\n",
                               scale_size = TRUE,
                               shape = c(rep("ellipse",length(group.target)),
                                         rep("ellipse",length(group.other))),
                               color = c(rep("#9955FF",length(group.target)),
                                         text.group.color),
                               force = TRUE)%>%
          
          visLayout(randomSeed = 10) # to have always the same network
      }else{
        group.legend = list(group.legend[[1]],group.legend[[2]],group.legend[[3]],
                            group.legend[[4]],group.legend[[5]],group.legend[[6]],group.legend[[7]])
        visNetwork(nodes, edges, width = "100%",height = "100%") %>%
          visNodes(color = list(background = "lightblue",
                                border = "darkblue",
                                highlight = "yellow"),
                   shadow = list(enabled = TRUE, size = 10)) %>%
          visEdges(
            #physics = FALSE,
            smooth = FALSE,
            hoverWidth = 2.5) %>%
          visOptions(highlightNearest =
                       list(enabled = T, degree = 1, hover = T,
                            hideColor = "rgba(200,200,200,0.2)"),
                     selectedBy = list(`variable` = "Cap_label",
                                       `multiple` = TRUE,
                                       `main` = "Select by group"),
                     #clickToUse = TRUE,
                     collapse = FALSE) %>%
          visLegend(width = 0.09, position = "right",
                    addNodes = group.legend,
                    addEdges = ledges,
                    useGroups = FALSE, zoom = FALSE,
                    stepX = 150, stepY = 75,ncol=1) %>%
          visInteraction(hover = TRUE,
                         navigationButtons = TRUE) %>%
          visPhysics(barnesHut = list("avoidOverlap"=0.7)) %>%
          visIgraphLayout(layout = "layout_nicely",physics = TRUE,
                          smooth = TRUE) %>%
          visEvents(selectNode = "function(nodes) {
              Shiny.onInputChange('current_node_id', nodes);
              ;}") %>%

          visLayout(randomSeed = 10) # to have always the same network
      }
    }else{
      visNetwork(data.frame(), data.frame(), width = "100%",
                 main = paste("Try to click some rows in",tagList(icon("table")),"to specify your nodes"))
    }
    
    
  })
  

  output$shiny_return_subnode_table <- renderTable({
    if(length(input$current_node_id$nodes[[1]])!=0 & length(input$inCheckboxGroup2)!=0){
      node_now = input$current_node_id$nodes[[1]]
      edge.ma.now = edge.full.list[[as.numeric(method())]][[1]]
      loc.node_now = match(node_now, colnames(edge.ma.now))
      if(is.na(loc.node_now)==FALSE){
        re = data.frame("info" = c("id","label","degree","sub_group"),
                        "value"= c(node_now,
                                   dict.combine$Description[match(node_now,dict.combine$Variable)],
                                   dict.combine[match(node_now,dict.combine$Variable), 1+n.col],
                                   dict.combine$group[match(node_now,dict.combine$Variable)]))
      }
    }else{
      "You haven't clicked on a node."
    }
  })
  
  output$shiny_returntable <- renderTable({
    input$goButton
    isolate({
      s = input$inCheckboxGroup2
      })
    if(length(input$current_node_id$nodes[[1]])!=0 & length(s)!=0){
      node_now = input$current_node_id$nodes[[1]]
      if(node_now %in% s){
        data.frame("message"="Try to click on a non-center node!")
      }else{
        edge.ma.now = edge.full.list[[as.numeric(method())]][[1]]
        loc.node_now = match(node_now, colnames(edge.ma.now))
        if(is.na(loc.node_now)==FALSE){
          node.id = edge.ma.now[loc.node_now,]!=0 | edge.ma.now[,loc.node_now]!=0
          connected.node_now = colnames(edge.ma.now)[node.id]
          connected.node_now = connected.node_now[connected.node_now %in% s]
          nn = length(connected.node_now)
          re = data.frame("info"=rep(c("id","label","degree","cos_Sim","---"),nn),
                          "value"=rep(NA,nn*5))
          for(i in 1:nn){
            re[i*5-4,2] = connected.node_now[i]
            re[i*5-3,2] = dict.combine$Description[match(connected.node_now,dict.combine$Variable)][i]
            re[i*5-2,2] = dict.combine[match(connected.node_now,dict.combine$Variable), as.numeric(method())+n.col][i]
            re[i*5-1,2] = max(edge.ma.now[loc.node_now,connected.node_now][i],
                              edge.ma.now[connected.node_now,loc.node_now][i])
            re[i*5,2] = "---"
          }
          re
        }
      }
    }
  })
  
  observe({
    input$addButton
    isolate({
      method = as.numeric(method())
      s = input$table_rows_selected
      if(length(s)!=0 & length(input$current_node_id$nodes[[1]])!=0){
        node_now = input$current_node_id$nodes[[1]]
        edge.ma.now = edge.full.list[[method]][[1]]
        loc.node_now = match(node_now, colnames(edge.ma.now))
        if(is.na(loc.node_now)==FALSE & !(loc.node_now %in% s)){
          node_now_name = colnames(edge.ma.now)[loc.node_now]
          x = input$inCheckboxGroup2
          if(!(node_now_name %in% input$inCheckboxGroup2)){
            x = c(node_now_name, x)
            x.neighbor = sapply(x, function(xx){
              sum(edge.full.list[[method]][[1]][xx,]!=0) + 
                sum(edge.full.list[[method]][[1]][,xx]!=0) - 
                sum(edge.full.list[[method]][[1]][,xx] != 0 & 
                      edge.full.list[[method]][[1]][xx,] != 0)
            })
            x.name = dict.combine$Description[match(x,dict.combine$Variable)]
            x.neighbor = paste0(x.name," (" ,x.neighbor," degrees)")
            updateCheckboxGroupInput(session, "inCheckboxGroup2",
                                     label = paste(length(x), " candidate nodes:"),
                                     choiceValues = x,
                                     choiceNames = x.neighbor,
                                     selected = x
            ) 
          }
        }
      }
    })
  })
  observe({
    s = input$table_rows_selected
    if(length(s)!=0){
      method = as.numeric(method())
      x = dict.combine$Description_s[match(colnames(edge.list[[method]])[s[seq(1,min(50,length(s)),by=1)]],
                                           dict.combine$Variable)]
      x = c(NA,x)
      updateSelectInput(session, "Focus","Choose one node to focus on:",
                        choices = x, selected = NA)
      
    }
    
  })
  observe({
    if(!is.na(input$Focus)){
      id = dict.combine$Variable[match(input$Focus,dict.combine$Description_s)]
      visNetworkProxy("network_proxy_nodes") %>%
        visFocus(id = id, scale = input$scale_id/10)
    }
  })
  

  output$downloadData <- downloadHandler(
    filename = function() {
      paste("node.csv")
    },
    content = function(path) {
      s = input$inCheckboxGroup2
      height = "900px"
      if(length(s)!=0){
        method = as.numeric(method())
        input.correct = s[seq(1,min(50,length(s)),by=1)]
        root.node = match(input.correct,rownames(edge.full.list[[method]][[1]]))
        draw.data = vis.InOutNodes(root.node=root.node, method=method,
                                   edge.list, edge.full.list, dict.combine)
        nodes = draw.data[[1]]
        edges = draw.data[[2]]
        # target = data.frame("m1" = match(edges$from, s),
        #                     "m2" = match(edges$to, s))
        # target[!is.na(target$m1),1] = paste0("node",target[!is.na(target$m1),1])
        # target[!is.na(target$m2),2] = paste0("node",target[!is.na(target$m2),2])
        # target[is.na(target$m1),1]="other"
        # target[is.na(target$m2),2]="other"
        # edges$direction = paste0(target$m1," -> ",target$m2)
        file = edges[,c(1,2,3,4)]
        file$from.term = nodes$label[match(file$from,nodes$id)]
        file$to.term = nodes$label[match(file$to,nodes$id)]
        file = file[,c(3,1,5,2,6,4)]
        file = file[order(file$corvalue,decreasing = TRUE),]
      }else{
        file = data.frame("Warning"="Try to click some rows in the 'Possible inputs' box to specify your nodes!")
      }
      write.csv(file,path,row.names = FALSE)
    }
  )

  

  
  
  
  
}



shinyApp(ui = ui, server = server)

