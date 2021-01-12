library(visNetwork)
library(shiny)
library(Matrix)
library(shinydashboard)
library(shinycssloaders)
library(DT)
source("helper.R")
load("data/edge_matrix.RData")



ui <- dashboardPage(
  dashboardHeader(title = "visNetwork"),
  ## Sidebar content
  dashboardSidebar(
    disable = TRUE,
    sidebarMenu(
      menuItem("Network", tabName = "Network", icon = icon("dashboard")),
      menuItem("Instructions", icon = icon("th"), tabName = "possibleinput")
    )
  ),
  
  dashboardBody(
    fluidRow(
      column(
        width = 12,
        box(collapsible=TRUE,status = "info",
            width = NULL, title = "Selectable inputs",
            DT::dataTableOutput("table"))
      )
    ),
    fluidRow(
      column(width = 3,
             box(width = NULL,title = "Arguments",
                 checkboxGroupInput("inCheckboxGroup2", "Input checkbox 2",
                                    c("")),
                 actionButton("goButton", "Show", icon = tags$i(class = "fas fa-project-diagram", style="font-size: 10px"), class = "btn-success"),
                 actionButton('refresh', 'Unselect', icon = tags$i(class = "fa fa-refresh", style="font-size: 10px")),
                 hr(),
                 selectInput("select", label = "Select data from station(s)", 
                             choices = list("Both" = 1, "RPDR local" = 2,"VA local" = 3,
                                            "RPDR integrative" = 4,"VA integrative" = 5), 
                             selected = 1,width = '75%'),
                 checkboxInput("cluster", "Cluster by groups*", value = TRUE)
             )
      ),
      column(width = 9,
             box(title = "Network visualization",
                 width = NULL, shinycssloaders::withSpinner(
                   visNetworkOutput("network_proxy_nodes", height = "500px"), type = 6
                 ))
             
      )
      
    ),
    fluidRow(
      column(width = 3),
      column(
        width = 9,
        box(collapsible=TRUE,collapsed = FALSE,
            width = NULL, title = "Neighbors of",
            textOutput("shiny_returntext",container = h4),
            actionButton("addButton", "Add to candidates", icon = tags$i(class = "fas fa-project-diagram", style="font-size: 10px"), class = "btn-success"),
            
            hr(),
            br(),
            DT::dataTableOutput("shiny_return"))
      )
    ),
    fluidRow(
      column(
        width = 12,
        hr(),
        strong("Instruction:"),
        p('1. Please refer to "Possible input" box to see appropriate values.
                                  Input NA if you do not need more inputs.'),
        p('2. If you select "cluster", double click any group node (purple circle) to unfold or
                                  any individual node (ellipse) to fold the group it belongs to. 
                                  Click "Reinitialize clustering" botton to fold all groups.'),
        p('3. Click any individual node (ellipse) to see its neighbors.'),
        p('4. The maximum number of input nodes has been set to be 50. However, less than 10 nodes
                    are recommended considering the possessing time.')
        
      )
    )
  )
)





server <- function(input, output, session) {
  method = reactive({input$select})
  cluster = reactive({input$cluster})
  
  output$table <- DT::renderDataTable(DT::datatable({
    method = as.numeric(method())
    data.frame("Node id"=colnames(edge.list[[1]]),
               "Description"=dict.combine$Description[match(colnames(edge.list[[1]]),dict.combine$Variable)],
               "Degree"=dict.combine[match(colnames(edge.list[[1]]),dict.combine$Variable),method+5],
               "Group"=dict.combine$group[match(colnames(edge.list[[1]]),dict.combine$Variable)])
    
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
        sum(edge.full.list[[method]][[1]][xx,]==1) + 
          sum(edge.full.list[[method]][[1]][,xx]==1) - 
          sum(edge.full.list[[method]][[1]][,xx] == 1 & 
                edge.full.list[[method]][[1]][xx,] == 1)
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

      if(cluster()==TRUE){
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
                     collapse = FALSE) %>%
          visLegend(width = 0.1, position = "right",
                    addNodes = group.legend,
                    useGroups = FALSE, zoom = FALSE,
                    stepX = 150, stepY = 75,ncol=1) %>%
          visInteraction(hover = TRUE) %>%
          visEvents(selectNode = "function(nodes) {
              Shiny.onInputChange('current_node_id', nodes);
              ;}") %>%
          visClusteringByGroup(groups = c(group.target,group.other),
                               label="",
                               scale_size = TRUE,
                               shape = c(rep("circle",length(group.target)),
                                         rep("circle",length(group.other))),
                               color = c(rep("#9955FF",length(group.target)),
                                         rep("#D1BBFF",length(group.other))),force = TRUE)%>%
          visLayout(randomSeed = 10) # to have always the same network
      }else{
        visNetwork(nodes, edges, width = "100%",
                   main = "Network Visualization") %>%
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
                     collapse = FALSE) %>%
          visLegend(width = 0.1, position = "right",
                    addNodes = group.legend,
                    useGroups = FALSE, zoom = FALSE,
                    stepX = 150, stepY = 75,ncol=1) %>%
         
          visInteraction(hover = TRUE) %>%
          visEvents(selectNode = "function(nodes) {
              Shiny.onInputChange('current_node_id', nodes);
              ;}") %>%
          visLayout(randomSeed = 10) # to have always the same network
      }
    }else{
     visNetwork(data.frame(), data.frame(), width = "100%",
                main = "Try to click some rows in the above table to specify your nodes")
    }
    
    
  })
  

  output$shiny_return <-
    DT::renderDataTable(DT::datatable({
      if(length(input$current_node_id$nodes[[1]])!=0){
        node_now = input$current_node_id$nodes[[1]]
        edge.ma.now = edge.full.list[[as.numeric(method())]][[1]]
        loc.node_now = match(node_now, colnames(edge.ma.now))
        if(is.na(loc.node_now)==FALSE){
          node.id = edge.ma.now[loc.node_now,]==1 | edge.ma.now[,loc.node_now]==1
          connected.node_now = colnames(edge.ma.now)[node.id]
          data.frame("Neighbor id"=connected.node_now,
                     "Description"=dict.combine$Description[match(connected.node_now,dict.combine$Variable)],
                     "Degree"=dict.combine[match(connected.node_now,dict.combine$Variable), as.numeric(method())+5],
                     "Group"=dict.combine$group[match(connected.node_now,dict.combine$Variable)]
          )
        }else{
          data.frame()
        }
      }else{
        data.frame()
      }
    }, rownames = FALSE))
  
  output$shiny_returntext <- renderText({
    if(length(input$current_node_id$nodes[[1]])!=0){
      node_now = input$current_node_id$nodes[[1]]
      edge.ma.now = edge.full.list[[as.numeric(method())]][[1]]
      loc.node_now = match(node_now, colnames(edge.ma.now))
      if(is.na(loc.node_now)==FALSE){
        node_now_name = colnames(edge.ma.now)[loc.node_now]
        paste0(dict.combine$Description[match(node_now_name,dict.combine$Variable)],":")
      }else{
        "Try to click on an individual node instead of a group circle."
      }
    }else{
      "You haven't clicked on a node."
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
        if(is.na(loc.node_now)==FALSE){
          node_now_name = colnames(edge.ma.now)[loc.node_now]
          x = input$inCheckboxGroup2
          x = c(node_now_name, x)
          x.neighbor = sapply(x, function(xx){
            sum(edge.full.list[[method]][[1]][xx,]==1) + 
              sum(edge.full.list[[method]][[1]][,xx]==1) - 
              sum(edge.full.list[[method]][[1]][,xx] == 1 & 
                    edge.full.list[[method]][[1]][xx,] == 1)
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
    })
  })
  
}



shinyApp(ui = ui, server = server)

