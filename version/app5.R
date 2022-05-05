library(shiny)
library(shinydashboard)
library(DT)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(magrittr)

manual_colors <- c("#0000FF", "grey15", "#FF0000")

header <- dashboardHeader(title = "Omics Volcano")

side <- dashboardSidebar(sidebarMenu(
  menuItem(
    "Home",
    tabName = "Home",
    icon = icon("home",
                verify_fa = FALSE),
    badgeColor = "green"
  ),
  menuItem(
    "Volcano",
    tabName = "Volcano",
    icon = icon("bar-chart-o",
                verify_fa = FALSE)
  )
))

body <- dashboardBody(tabItems(
  tabItem(
    tabName = "Home",
    h2("Welcome to omicsVolcano!"),
    br(),
    h4(
      "A volcano plot is a kind of scatter plot which represents differential expression of genes or proteins."
    ),
    h4("We typically find the fold change on the x-axis and the p-value on the y-axis."),
    h4("You can upload your data, visualize with volcano plot online and download the volcano plot."),
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
              align = "center",
              fluidRow(
                radioButtons(
                  "data_file_type",
                  h3("Try example data or upload your data"),
                  c("Example Data" = "examplecounts",
                    "Upload Data" = "upload"),
                  selected = "examplecounts"
                ),
                conditionalPanel(
                  condition = "input.data_file_type=='upload'",
                  fileInput(
                    inputId = 'file1',
                    label = 'Choose tab-delimited File',
                    accept = c('.txt'),
                    buttonLabel = "Browse...",
                    placeholder = "No file selected"
                  )
                )
              ),
              fluidRow(
                selectInput(
                  inputId = "variable_X",
                  label = "X axis:",
                  choices = "foldchange"
                )
              ),
              fluidRow(
                selectInput(
                  inputId = "variable_Y",
                  label = "Y axis:",
                  choices = "pvalue"
                )
              ),
              fluidRow(
                textInput(
                  inputId = "axisTitle_X",
                  label = "X axis title:",
                  value = "Log2(Fold change)"
                )
              ),
              fluidRow(
                textInput(
                  inputId = "axisTitle_Y",
                  label = "Y axis title:",
                  value = "-log10(P)"
                )
              ),
              fluidRow(
                numericInput(
                  inputId = "fc_cut",
                  label = "fold change cutoff",
                  value = 1.5
                )
              ),
              fluidRow(
                numericInput(
                  inputId = "p_cut",
                  label = "p value cutoff",
                  value = 0.05
                )
              ),
              fluidRow(
                downloadButton("downloadData1",
                               "Download table as csv"),
                downloadButton("downloadPlot1",
                               "Download plot as pdf")
              ),
            ),
            column(width = 9,
                   fluidRow(plotOutput("p_volcano")),
                   fluidRow(dataTableOutput("table_analysis")))
          ))
))


shinyApp(
  ui = dashboardPage(header,
                     side,
                     body),
  server = function(input, output) {
    dataset <- reactive({
      isolate({
        validate(need((input$data_file_type == "examplecounts") |
                        (!is.null(input$filedata)),
                      message = "Please upload data file"
        ))
        if (input$data_file_type == "examplecounts") {
          inFileData <- read.delim("www/example_data.txt",
                                   header = TRUE,
                                   stringsAsFactors = FALSE)
        } else {
          inFile <- input$file1
          if (!is.null(inFile))
            inFileData <- read.delim(inFile$datapath,
                                     header = TRUE,
                                     stringsAsFactors = FALSE)
        }
        return(inFileData)
      })
    })
    
    data <- reactive(
      dataset() %>%
        dplyr::mutate(log2FC = log2(fold_change)) %>%
        dplyr::mutate(neglogP = -log10(p_value)) %>%
        dplyr::mutate(is_dep = p_value < input$p_cut &
                        abs(log2FC) >= log2(input$fc_cut)) %>%
        dplyr::mutate(regulated = ifelse(
          is_dep == TRUE,
          ifelse(log2FC > log2(input$fc_cut),
                 "Up",
                 "Down"),
          "None"
        ))
    )
    
    volcano <- reactive({
      ggplot(data = data(), aes(x = log2FC,
                                y = neglogP,
                                color = regulated)) +
        geom_point(size = 1.0, alpha = 0.5) +
        scale_x_continuous(name = input$axisTitle_X) +
        scale_y_continuous(name = input$axisTitle_Y) +
        geom_hline(
          yintercept = -log10(input$p_cut),
          linetype = "dashed",
          color = "grey",
          size = 0.5
        ) +
        geom_vline(
          xintercept = log2(input$fc_cut),
          linetype = "dashed",
          color = "grey",
          size = 0.5
        ) +
        geom_vline(
          xintercept = log2(1 / input$fc_cut),
          linetype = "dashed",
          color = "grey",
          size = 0.5
        ) +
        scale_color_manual(values = manual_colors,
                           name = NULL) +
        ggpubr::theme_pubr()
    })
    
    output$p_volcano <- renderPlot({
      volcano()
    })
    
    output$table_analysis <- renderDataTable({
      datatable(data())
    })
    
    output$downloadData1 <- downloadHandler(
      filename = function() {
        paste('Analysis_result_', Sys.Date(), '.csv', sep = '')
      },
      content = function(con) {
        write.csv(data(), con,
                  row.names = FALSE)
      }
    )
    
    output$downloadPlot1 <- downloadHandler(
      filename = function() {
        paste0("volcano_plot_", Sys.Date(), ".pdf")
      },
      
      content = function(file) {
        ggsave(
          file,
          volcano(),
          device = "pdf",
          width = 8,
          height = 8,
          units = "in"
        )
      }
    )
    
    
  }
)
