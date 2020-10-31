rm(list=ls())

# Library For shiny App
options(rgl.useNULL=TRUE)
if (!require("shiny")){
 install.packages("shiny")
  library("shiny")
}

# Library For SNF Algo
options(rgl.useNULL=TRUE)
if (!require("SNFtool")){
 install.packages("SNFtool")
 library("SNFtool")
 }

# Library For Force Directed Graph
options(rgl.useNULL=TRUE)
if (!require("networkD3")){
 install.packages("networkD3")
 library("networkD3")
}

# Library For Correlation Matrix
options(rgl.useNULL=TRUE)
if (!require("corrplot")){
  install.packages("corrplot")
  library("corrplot")
} 

# Library For Clustering Techniques
options(rgl.useNULL=TRUE)
if (!require("cluster")){
  install.packages("cluster")
  library("cluster")
}
  
ui <- pageWithSidebar( 

  headerPanel("Similarity Network Fusion"),
  
  sidebarPanel(
   
      h3("Select Parameters"),
      tags$div( class ="well", style = "color: green; border-radius: 16px;
                border-width: 12px: border-color: blue; border-style: solid",
      br(),
      selectInput("SelectCancer", "Cancer Type:",
                        c("Breast" = "Breast",
                          "Colon" = "Colon",
                          "GBM" = "GLIO",
                          "Kidney" = "Kidney",
                          "Lung" = "Lung")
                  ),
                  
      numericInput('neighbors', 'Number Of Neighbors (K)', 20,
                   min = 1, max = 30),
      numericInput('hyper', 'Hyper Parameter (alpha), usually (0.3 to 0.8)', 0.5,
                   min = 0.1, max = 1),
      numericInput('iterations', 'Number Of Iterations (T),  usually (10 to 20)', 20,
                   min = 1, max = 50)
      ),
      tags$div( class ="well",
      actionButton(inputId= "SNF", label ="RUN SNF", 
                   style = "color: white;
                   font-size: 20px;
                   background-color: #6495ED; 
                   position: relative; 
                   left: 3%;
                   height: 45px;
                   width: 200px;
                   text-align:center;
                   text-indent: -2px;
                   border-radius: 6px;
                   border-width: 2px"), 
      p("Click the button and find the estimated clusters from main panel.")
      ),
      tags$div( class ="well",
      h3("Choose Clusters and % for Network visualization"),
      numericInput('clusters', 'Number Of Clusters', 0,
                   min = 0, max = 9),
      hr(),

      selectInput("ClusterMethodChoice", "Clustering Technique:",
                  c("Spectral Clustering" = "spectralClustering",
                    "Hierarchical Clustering" = "hclust",
                    "Agglomerative Clustering" = "agnes",
                    "Partitioning Around Medoids" = "pam",
                    "Clustering Large Applications" = "clara")
                  ),
      
      actionButton(inputId= "ClusterMethod", label ="Cluster Data", 
                   style = "color: white;
                   font-size: 20px;
                   background-color: #6495ED; 
                   position: relative; 
                   left: 3%;
                   height: 45px;
                   width: 200px;
                   text-align:center;
                   text-indent: -2px;
                   border-radius: 6px;
                   border-width: 2px"), 
      p("Click the button for visualizing node-link interactions under 'Detailed View' tab."),
      p("Select top 'X' % of interactions to visualize, using the slider below."),
      sliderInput("toppercentage",
                  "Top % of Interactions",
                  min = 1,
                  max = 100,
                  value = 5)                   
      
),
tags$div( class ="well",
      h3("Downloading Data"),
      br(),
      downloadButton("downloadData", "Download Plotted Node Links"),
      br(),
      br(),
      downloadButton("downloadClusters", "Download Clusters"),
      br(),
      br(),
      downloadButton("downloadCompleteData", "Download Complete SNF Matrix")
)
      
    ),
    mainPanel(
    
      tabsetPanel(type = "tabs",
                        tabPanel("StepsToRunSNF",
                                 tags$style(type='text/css', '#clusters {background-color: rgba(255,255,0,0.40); color: green;font-size: 18px;}'),
                                 tags$style(type='text/css', '#summary {background-color: rgba(0,255,0,0.40); color: blue;font-size: 18px;}'),
                                 tags$style(type='text/css', '#text,text2 { color: black;font-size: 20px;}'),  
                                 #helpText('Steps to run SNF: '),
                                 br(),
                                 htmlOutput("text1"),
                                 br(),
                                 verbatimTextOutput("summary"),
                                 br(),
                                 htmlOutput("text"),
                                 br(),
                                 htmlOutput("text2"),
                                 br(),
                                 verbatimTextOutput("clusters", placeholder = TRUE)),
                                 tabPanel("Overview", plotOutput("corrPlot", width = "100%", height = "700px"),  plotOutput("hclust", width = "100%", height = "700px")), 
                                 tabPanel("Detailed View", forceNetworkOutput(outputId = "net", height="500px"))
                                 

    )
))

server <- function(input, output, session) {

DATA <- reactive({ input$SelectCancer })  
toReadPath <- reactive({  paste("D:/Research/SNF-NewCode/SNF-Data/InputData/",DATA(),"/", sep="")  })
toWritePath <- reactive({  paste("D:/Research/SNF-NewCode/SNF-Data/InputData/",DATA(),"/", sep="")  })
  
Methylation <- reactive({  paste(DATA(),"_Methy_Expression",sep="")  })
mRNA <- reactive({ paste(DATA(),"_Gene_Expression",sep="") })
miRNA <- reactive({ paste(DATA(),"_Mirna_Expression",sep="") })  

Data1 <- reactive({
  
if ((DATA() == "Colon") || (DATA() == "Kidney") ) {
     dat<-  read.table(paste(toReadPath(),Methylation(),".csv", sep =""),sep =",",header = TRUE)
     dat <- dat[-1] 
 }
  
else {
     dat <- read.table(paste(toReadPath(),Methylation(),".txt", sep =""),header = TRUE)
 }  
})
observe({ print("Complete Data Size") })
observe({ print(dim(Data1())) })
#Data1 <-  reactive({ read.table(paste(toReadPath(),Methylation(),".txt", sep =""),header = TRUE) }) 
Data2 <-  reactive({ read.table(paste(toReadPath(),mRNA(),".txt", sep =""),header = TRUE) }) 
Data3 <-  reactive({ read.table(paste(toReadPath(),miRNA(),".txt", sep =""),header = TRUE) }) 
  
## Calculate the pair-wise distance;  
Dist1 <- reactive({ dist2(as.matrix(t(Data1())),as.matrix(t(Data1()))) }) 
Dist2 <- reactive({ dist2(as.matrix(t(Data2())),as.matrix(t(Data2()))) }) 
Dist3 <- reactive({ dist2(as.matrix(t(Data3())),as.matrix(t(Data3()))) }) 

## Hyperparameters;
K <- reactive({ (input$neighbors) })  
alpha <- reactive({ (input$hyper) })  
iter <- reactive({ (input$iterations) })

## col and row details
Features <- reactive({ ncol(Data2()) })
names <- reactive({ colnames(Data2()) })

W1 <- reactive({  affinityMatrix(Dist1(), K(), alpha())  })
W2 <- reactive({  affinityMatrix(Dist2(), K(), alpha())  })
W3 <- reactive({  affinityMatrix(Dist3(), K(), alpha()) })

output$text1 <- renderUI({ HTML(paste("<h3><b> I) Steps to run SNF:</b></h3>")) })

output$summary <- renderText({ HTML(paste('The Number of Features (N) of the dataset is,', Features())) })

output$text <- renderUI({  
  str1 <- paste("<u><b>Select the following parameters from the side panel</b></u>", "1. <i>Number of neighbors (K)</i> -- <b>K = N/C</b>, where C = Number of clusters.</br> &thinsp; &thinsp; &thinsp; &thinsp; &thinsp; &thinsp; &thinsp; &thinsp; &thinsp; &thinsp; &thinsp; &thinsp; &thinsp; &thinsp; &thinsp; &thinsp; &thinsp; &thinsp; &thinsp; &thinsp; &thinsp; &thinsp; &thinsp; &thinsp; &thinsp; &thinsp; &thinsp; &thinsp;<b>K = N/10</b> If C is not known ", "2. <i>Hyper parameter (alpha)</i> -- Normally in the range of 0.3 to 0.8", "3. <i>Number of iterations for the network to be fused (T)</i> -- Normally in the range of 10 to 20.","&thinsp; &thinsp; &thinsp;Empirically found, iterations of 20 is enough for converge.", sep="<br/>")
  str2 <- paste("<b>Click &quotRUN SNF&quot button</b>")
  HTML(paste(str1, str2, sep = '<br/><br/>'))
})


output$text2 <- renderUI({
  str1 <- (paste("<h3><b>II) Identifying &quotNumber of Clusters (K)&quot, to find the interactions accordingly: </b></h3>"))
  str2 <- paste("<h3>Computed Clusters using &quotEigen gap&quot (C1, C12) and &quotRotate cost&quot (C2, C22) methods:</h3")
  HTML(paste(str1, str2, sep = '<br/>'))
})

snfoutput <- eventReactive(input$SNF, {
  W <- reactive({ SNF <- SNF(list(W1(),W2(),W3()), K(), iter()) * 2 })
  estimatedClusters <- reactive({ estimateNumberOfClustersGivenGraph(W(), 2:15) })
  DATAOutput <- list(W = W, estimatedClusters = estimatedClusters)
})


output$clusters <- renderText({
  DATAOutput <- snfoutput()
  estimatedClusters <- DATAOutput$estimatedClusters
  observe({ print("Computed Clusters using Eigen gap (C1, C12) and Rotate cost (C2, C22) methods ") })
  observe({ print(estimatedClusters()) })
  str1 <- paste(" C1 -", estimatedClusters()[[1]])
  str2 <- paste(" C12 -", estimatedClusters()[[2]])
  str3 <- paste(" C2 -", estimatedClusters()[[3]])
  str4 <- paste(" C22 -", estimatedClusters()[[4]])
  HTML(paste(str1,str2,str3,str4, sep = ','))
})



labels <- eventReactive(input$ClusterMethod, {
  DATAOutput <- snfoutput()
  W <- DATAOutput$W
  if (input$ClusterMethodChoice=="spectralClustering") {
    cluster <- spectralClustering(W(), input$clusters) 
  }
  else if (input$ClusterMethodChoice =="hclust") {
    toCluster <- as.matrix(W())
    d <- dist(toCluster, method = "euclidean")
    fit <- hclust(d, method="average") 
    cl <- cutree(fit, input$clusters) # cut tree for clusters
  }
  else if (input$ClusterMethodChoice =="agnes") {
    toCluster <- as.matrix(W())
    d <- dist(toCluster, method = "euclidean")
    fit <- agnes(d, method="ward") 
    cl <- cutree(fit, input$clusters) # cut tree for clusters
  }
  else if (input$ClusterMethodChoice =="pam") {
    toCluster <- as.matrix(W())
    cl <- pam(toCluster, input$clusters, diss = FALSE) 
    cl$cluster
  }
  else if (input$ClusterMethodChoice =="clara") {
    toCluster <- as.matrix(W())
    cl <- clara(toCluster, input$clusters, metric = "euclidean", stand = FALSE)
    cl$cluster
  }
  # Commented Fuzzy Clustering, as this type of technique is not appropriate for patient Similarity network
  #else if (input$ClusterMethodChoice =="fanny") {
    #toCluster <- as.matrix(W())
    #cl <- fanny(toCluster, input$clusters, diss = FALSE)
    #cl$cluster
  #}
})

Clusters <- reactive ({ 
      nodeIDs <- names()
      group <- labels()
      d <- data.frame(nodeIDs, group)    
  })

#Plotting corrplot for all cluster techniques
output$corrPlot <- renderPlot({
  DATAOutput <- snfoutput()
  M <- (as.matrix(DATAOutput$W()))
  orderedRows <- Clusters()[with(Clusters(), order(group)), ]
  rownames <- as.numeric(rownames(orderedRows))
  W1 <- M[,c(rownames)]
  W2 <- W1[c(rownames),]
  colnames(W2) <- c(rownames)
  rownames(W2) <- c(rownames)
  n <- length(unique(Clusters()$group))
  col <- colorRampPalette(c("white", "gray", "black"))
  #corrplot(MAT, tl.col="black", tl.cex = 0.5, cl.cex = 1, method = "square", is.corr = TRUE, cl.lim = c(-1, 1))
 
  newW2 <- (W2 - min(W2))/(max(W2) - min(W2))*10
 
  diag(newW2) <- 1
  M <- range(newW2)
  corrplot(newW2, method = "shade", tl.col = "black", tl.cex = 0.6, cl.cex = 1, is.corr = TRUE, cl.lim = c(M[1], M[2]), title ="SNF Integrated Matrix")
})

#Plotting dendogram for h-clust and agglom
observeEvent(input$ClusterMethod,{ 
  
  if (input$ClusterMethodChoice =="hclust") {
  DATAOutput <- snfoutput()
  M <- (as.matrix(DATAOutput$W()))
  d <- dist(M, method = "euclidean")
  fit <- hclust(d, method="ward.D2")
  output$hclust <-  renderPlot({
    plot(fit) # display dendogram
    rect.hclust(fit, input$clusters, border="red") # draw dendogram with red borders around the clusters
    })
  }
  else if (input$ClusterMethodChoice =="agnes") {
    DATAOutput <- snfoutput()
    M <- (as.matrix(DATAOutput$W()))
    fit <- agnes(M, method="ward", stand = FALSE)
    output$hclust <-  renderPlot({
      plot(fit) # display dendogram
      rect.hclust(fit, input$clusters, border="red") # draw dendogram with red borders around the clusters
    })
  }
  else {
    output$hclust <-  renderPlot({
      plot <- ""
    })
  }
})

#Commented Heatmap Code
#output$hclust <-  renderPlot({

#  if (input$ClusterMethodChoice =="hclust") {
    #plot(onClick()) # display dendogram
    #rect.hclust(onClick(), input$clusters, border="red") # draw dendogram with red borders around the clusters 
 #   }

  #  my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 20)
  #normalize <- function(X) X / rowSums(X)
  
  #  heatmap.2(1-M,
  #            cellnote = 1-M,  # same data set for cell labels
  #            main = "Similarity between gene pairs after network fusion", # heat map title
  #            notecol="black",      # change font color of cell labels to black
  #            density.info="none",  # turns off density plot inside color legend
  #            trace="none",         # turns off trace lines inside the heat map
  #            margins =c(12,9),     # widens margins around plot
  #            col=my_palette,       # use on color palette defined earlier
  #            #breaks=col_breaks,    # enable color transition at specified limits
  #            dendrogram="row",     # only draw a row dendrogram
  #            Colv="NA")            # turn off column clustering
  #  
#})


subsetdata <- reactive({ 
  DATAOutput <- snfoutput()
  W <- DATAOutput$W
  linksData <- NULL
  for (i in 1:Features())
  {
    value = source = target = NULL
    for (j in 1:Features())
    {
      if (i != j & i > j)
      {
        source <- rbind(source, i-1)
        target <- rbind(target, j-1)
        value <- rbind (value, W()[i,j])
      }
      
    }
    newsetValues <- cbind(source, target, value)
    linksData <- rbind(linksData, newsetValues)
  }
  source <- (linksData[,1])
  target <- (linksData[,2])
  value <- linksData[,3]
  a <- data.frame(source, target, value)
})  
observe({ print(dim(subsetdata())) })
newlinks <- reactive ({ dat <-  subset(subsetdata(), subsetdata()$value > quantile(subsetdata()$value, prob = 1 - input$toppercentage/100)) })
vals <- reactive({ quantile(subsetdata()$value, prob = 1 - input$toppercentage/100) })
observe({ print(dim(newlinks())) })

#Plotting force directed network using D3
output$net <- renderForceNetwork({
  p <- forceNetwork(
    Links  = newlinks(), Nodes   = Clusters(),
    Source = "source", Target  = "target",
    Value  = "value",  NodeID  = "nodeIDs",
    Group  = "group",  opacity = 1, zoom = TRUE, legend = TRUE, bounded = TRUE, linkColour = "black",
    fontSize = 20)  
})

output$downloadData <- downloadHandler(
  filename = function() {
    paste("SelectedFinalData-", Sys.Date(), ".csv", sep="")
  },
  content = function(file) {
    write.csv(newlinks(), file, row.names=FALSE, col.names=names)
  }
)

output$downloadClusters <- downloadHandler(
  filename = function() {
    paste("ClusterData-", Sys.Date(), ".csv", sep="")
  },
  content = function(file) {
    write.csv(Clusters(), file, row.names = FALSE, col.names=c("GeneSymbol", "Clusters"))
  }
)

output$downloadCompleteData <- downloadHandler(
  filename = function() {
    paste("SNF_ComplteMatrix-", Sys.Date(), ".csv", sep="")
  },
  content = function(file) {
    DATAOutput <- snfoutput()
    W <- DATAOutput$W
    write.csv(W(), file, row.names = FALSE, col.names=c("source", "target", "value"))
  }
)
}

shinyApp(ui = ui, server = server)


