library(shiny); library(DT)
library(foreach); library(doParallel)
source("panelDesignerEngine.R")
source('Inventory.R', encoding = 'UTF-8')


server <- function(input, output) {
  markerInfo <- loadFiles()
  output$mytable = DT::renderDataTable(markerInfo[[2]],server=TRUE)
  #output$x4 = renderPrint({
    #s = input$mytable_rows_selected
    #if (length(s)) {   cat('These', length(s) ,'antibodies were selected:\n\n');    cat(markerInfo[[2]][s,1], sep = ', ') }
  #  })
  
  output$antibodies <- renderUI({
    s = input$mytable_rows_selected
    if (length(s)) {
      antibodyList <- lapply(1:length(input$mytable_rows_selected), function(i) {
        list(tags$p(tags$u(h4(renderText(markerInfo[[2]][s[i],1])))),
             column(3,
           sliderInput(paste0("intensity",i), label = "Stain Intensity:  Lo 1 to 5 Hi", value=3, min=1, max=5, ticks=FALSE),
           radioButtons(paste0("priority", i), label = "Priority?", choices=c("High","Normal","Low"), selected="Normal"))  
          )
        })
      }
    })

  observeEvent(input$button, { isolate({
    
    inputLabelsIntensity <- c(); inputLabelsPriority <- c()
    userModSelectedAb <- data.frame(Protein = markerInfo[[2]][input$mytable_rows_selected,1]) #protein names
    for (i in 1:length(input$mytable_rows_selected))
    {
      inputLabelsIntensity[i] <- paste0("intensity",i)
      inputLabelsPriority[i] <- paste0("priority",i)
      userModSelectedAb$StainIntensity[i] <- input[[inputLabelsIntensity[i]]]
      userModSelectedAb$Priority[i] <- input[[inputLabelsPriority[i]]] 
    }
    
    #showModal(modalDialog(renderTable(userModSelectedAb), easyClose = TRUE))
    selectedAb <- userModSelectedAb
    
    progress <- shiny::Progress$new(min=0, max=1);  progress$set(message = "Calculating", detail="progress", value=0); 
      on.exit(progress$close())
    updateProgress <- function(x,detail) {        progress$set(value=x,detail=detail)    }
    
    if (input$inventory == "local")
    { inventory <- inventoryLocal(markerInfo); print("success Local inventory") } else if (input$inventory == "commercial") 
    { inventory <- inventoryCommercial(markerInfo); print("success Commercial inventory")}
    inventory <- convertInventory(inventory,markerInfo)
    
    #selectedAb <- markerInfo[[2]][input$mytable_rows_selected,]
    #browser()
    unmatched <- sanityCheck(selectedAb,inventory,markerInfo)
    #if(length(unmatched) > 0)   {   showModal(modalDialog(title = "Not in inventory", "The following antibodies were not found in the catalog so they'll be randomly assigned: ",
    #    print(selectedAb[unmatched]),        easyClose = TRUE      ))    }  
    #browser()
    if (input$cytometer == "LSR E")
    { 
      holdingMatrix <- data.frame(channels=markerInfo[[4]][-c(5,7:13,15,21),],stringsAsFactors = FALSE)
      rownames(holdingMatrix) <- holdingMatrix$channels
      markerInfo[[4]] <- holdingMatrix  ; print("using LSR E panel")
    } 
    
    generatedPanels <- panelDesigner(input$numIter, selectedAb, inventory, markerInfo, updateProgress)
    scoreEval <- scoreEvaluation(generatedPanels)
    output$resultSummary = renderDataTable({    
      datatable(scoreEval, selection=list(mode="single", selected=1), options = list(searching=FALSE)) })
    output$selectedResult = DT::renderDataTable({
      datatable(generatedPanels[[input$resultSummary_rows_selected]][c(1,5:9)],selection="none",
                options = list(paging = FALSE,searching=FALSE))
      })
    generatedPanels.df <- simpleResults(generatedPanels) 
    output$downloadData <- downloadHandler(
      filename = function() { paste("generatedPanels",Sys.Date(),".csv",sep="")}, 
      content = function(file) { write.csv(generatedPanels.df, file)}
    )
  }) })
  
}

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      radioButtons("inventory", "Choose your fluorochrome inventory:", c("local","commercial")),
      radioButtons("cytometer", "Choose which cytometer you are using:", c("X50", "LSR E")),
      numericInput("numIter","How many iterations do you want to run:",value =1,min=1,max=1000,step=1),
      verbatimTextOutput('x4'), 
      uiOutput("antibodies"),
      actionButton("button", "Calculate")
      ),
    mainPanel(
      tabsetPanel(
        tabPanel("Setup", 
                 column(3, DT::dataTableOutput('mytable'))
                 ),
        tabPanel("Results", 
                 column(3,DT::dataTableOutput("resultSummary"), downloadButton('downloadData','Download')),
                 column(4,DT::dataTableOutput("selectedResult"))
                 )
      )
    )
  )
)

shinyApp(ui = ui, server = server)