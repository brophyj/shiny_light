library(shiny)
library(DT)

# Sample data
initial_data <- data.frame(
  PMID = c("12345678"),
  Title = c("Effect of Drug X on Cardiovascular Outcomes"),
  Journal = c("Cardiology Today"),
  Year = c("2022"),
  Authors = c("Smith J, Doe A"),
  AbstractConclusion = c("Drug X reduced cardiovascular events significantly."),
  ArmSizes = c("100; 120"),
  EventCounts = c("10; 15"),
  PrimaryOutcome = c("Major adverse cardiovascular events"),
  EffectSize = c("HR 0.85"),
  ConfidenceInterval = c("0.75–0.95"),
  stringsAsFactors = FALSE
)

ui <- fluidPage(
  titlePanel("Cardiovascular RCT Data Entry"),
  sidebarLayout(
    sidebarPanel(
      textInput("pmid", "PMID"),
      textInput("title", "Title"),
      textInput("journal", "Journal"),
      textInput("year", "Year"),
      textInput("authors", "Authors"),
      textAreaInput("abstract", "Abstract Conclusion", "", rows = 3),
      textInput("arms", "Arm Sizes (e.g., 100; 120)"),
      textInput("events", "Event Counts (e.g., 10; 15)"),
      textInput("outcome", "Primary Outcome"),
      textInput("effect", "Effect Size (e.g., HR 0.85)"),
      textInput("ci", "Confidence Interval (e.g., 0.75–0.95)"),
      actionButton("add", "Add Entry"),
      actionButton("delete", "Delete Selected Row"),
      downloadButton("save", "Save to CSV")
    ),
    mainPanel(
      DTOutput("table")
    )
  )
)

server <- function(input, output, session) {
  data <- reactiveVal(initial_data)
  selected <- reactiveVal(NULL)
  
  observeEvent(input$add, {
    # Validate numeric fields
    arm_valid <- all(grepl("^\\d+$", unlist(strsplit(input$arms, ";\\s*"))))
    event_valid <- all(grepl("^\\d+$", unlist(strsplit(input$events, ";\\s*"))))
    
    if (!arm_valid || !event_valid) {
      showModal(modalDialog(
        title = "Validation Error",
        "Arm Sizes and Event Counts must be numeric values separated by semicolons.",
        easyClose = TRUE
      ))
      return()
    }
    
    new_row <- data.frame(
      PMID = input$pmid,
      Title = input$title,
      Journal = input$journal,
      Year = input$year,
      Authors = input$authors,
      AbstractConclusion = input$abstract,
      ArmSizes = input$arms,
      EventCounts = input$events,
      PrimaryOutcome = input$outcome,
      EffectSize = input$effect,
      ConfidenceInterval = input$ci,
      stringsAsFactors = FALSE
    )
    data(rbind(data(), new_row))
  })
  
  output$table <- renderDT({
    datatable(data(), selection = "single", editable = TRUE)
  })
  
  observeEvent(input$table_cell_edit, {
    info <- input$table_cell_edit
    df <- data()
    df[info$row, info$col] <- info$value
    data(df)
  })
  
  observeEvent(input$delete, {
    sel <- input$table_rows_selected
    if (!is.null(sel)) {
      df <- data()
      df <- df[-sel, ]
      data(df)
    }
  })
  
  output$save <- downloadHandler(
    filename = function() {
      paste0("cardiovascular_rct_results_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(data(), file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
