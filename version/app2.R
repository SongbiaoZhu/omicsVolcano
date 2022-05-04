library(shiny)
library(shinydashboard)
library(DT)
library(ggpubr)
library(ggplot2)

vchoices <- 1:4
names(vchoices) <- c("foldchange",
                     "pvalue")

header <- dashboardHeader(title = "Volcano visualization")

side <- dashboardSidebar(sidebarMenu(
  menuItem(
    "Home",
    tabName = "Home",
    icon = icon("home"),
    badgeColor = "green"
  ),
  menuItem(
    "Volcano",
    tabName = "Volcano",
    icon = icon("bar-chart-o")
  )
))

body <- dashboardBody(tabItems(
  tabItem(
    tabName = "Home",
    h2("Welcome to shinyVolcano!"),
    br(),
    h4("You can upload your data and visualize it with volcano online."),
    h4("You can also download the volcano."),
    
    br(),
    h3("Required data format:"),
    img(
      src = "example_data_format_volcano.png",
      width = 512,
      height = 369
    ),
    br(),
    br(),
    h3("Useful links: "),
    h4(
      a("Protein Chemistry & Proteomics Facility",
        href = "http://phoenix.tsinghua.edu.cn/index.php?c=show&id=93")
    ),
    h4(
      a("蛋白质化学与组学平台实验方法",
        href = "http://phoenix.tsinghua.edu.cn/index.php?c=show&id=584")
    ),
    
    br(),
    h3("Contact information:"),
    h4("songbiao_zhu@163.com")
  ),
  tabItem(tabName = "Volcano",
          fluidRow(
            column(
              width = 3,
              fluidRow(
                fileInput(
                  inputId = 'file1',
                  label = 'Choose tab-delimited File',
                  accept = c('.txt'),
                  buttonLabel = "Browse...",
                  placeholder = "No file selected"
                )
              ),
              fluidRow(
                actionButton(inputId ="goButton", 
                             label = "Click to upload data !")
              ),
              fluidRow(
                selectInput(
                  inputId = "variable_X",
                  label = "X axis:",
                  choices = vchoices
                )
              ),
              fluidRow(
                selectInput(
                  inputId = "variable_Y",
                  label = "Y axis:",
                  choices = vchoices
                )
              )
            ),
            column(width = 9,
                   fluidRow(plotOutput("p_volcano")),
                   fluidRow(dataTableOutput("dat_raw")))
          ))
))







shinyApp(
  ui = dashboardPage(header,
                     side,
                     body),
  server = function(input, output, session) {
    dataset <- reactive({
      inFile <- input$file1
      if (!is.null(inFile))
        read.delim(inFile$datapath,
                   header = TRUE,
                   stringsAsFactors = FALSE)
    })
    
    info <- eventReactive(input$goButton, {
      inFile <- input$file1
      if (is.null(inFile))
        return(NULL)
      isolate(f <- read.delim(
        inFile$datapath,
        header = TRUE,
        stringsAsFactors = FALSE
      ))
      f
    })
    
    observe({
      f <- info()
      if (!is.null(f)) {
        updateSelectInput(session, "variable_X", choices = names(f))
        updateSelectInput(session, "variable_Y", choices = names(f))
      }
      
    })
    
    output$p_volcano <- renderPlot({
      ggplot(data = info(), aes(
        x = (!!!input$variable_X),
        y = (!!!input$variable_Y),
      )) +
        geom_point() +
        ggpubr::theme_pubr()
    })
    
    output$dat_raw <- renderDataTable({
      datatable(dataset())
    })
    
  }
)
