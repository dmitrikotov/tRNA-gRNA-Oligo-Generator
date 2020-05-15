library(shiny)
library(Biostrings)

# Define UI for data upload app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("tRNA-gRNA Oligo Generator"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Select a file ----
      fileInput("file1", "Choose CSV File",
                multiple = TRUE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      
      # Horizontal line ----
      tags$hr(),
      
      # Input: Checkbox if file has header ----
      checkboxInput("header", "Header", TRUE),
      
      # Input: Select separator ----
      radioButtons("sep", "Separator",
                   choices = c(Comma = ",",
                               Semicolon = ";",
                               Tab = "\t"),
                   selected = ","),
      
      # Input: Select quotes ----
      radioButtons("quote", "Quote",
                   choices = c(None = "",
                               "Double Quote" = '"',
                               "Single Quote" = "'"),
                   selected = '"'),
      
      # Horizontal line ----
      tags$hr(),
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Data file ----
      tableOutput("contents")
      
    )
    
  )
)

Aar1_tRNA = "TATCACCTGCCCCC"
Aar1_rep  = "TATCACCTGCCCCA"
tRNA = "TGCACCAGCCGGGAATCG"
rep = "GTTTTAGAGCTAGAAATAGC"

oligos <- function(x){
  x <- DNAString(x)
  rep_guide = substr(x,9,20)
  if (substr(as.character(x),1,4) == "GGTG") stop("Guide generates an AarI cut site")
  tRNA_guide = reverseComplement(substr(x,1,12))
  tRNA_oligo = paste(Aar1_tRNA,tRNA_guide,tRNA, sep = "")
  rep_oligo = paste(Aar1_rep,rep_guide,rep, sep = "")
  df <- data.frame(tRNA_oligo, rep_oligo)
  return(df)
}


# Define server logic to read selected file ----
server <- function(input, output, Aar1_rep, Aar1_tRNA, tRNA, rep) {
  
  output$contents <- renderTable({
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.
    
    req(input$file1)
    
    df <- read.csv(input$file1$datapath,
                   header = input$header,
                   sep = input$sep,
                   quote = input$quote)

   results <- as.data.frame(matrix(0, nrow(df), 3))
   colnames(results) <- c("gene","tRNA_oligo","rep_oligo")
   
   for (i in 1:nrow(df)) {
     results[i,2] <- as.character(oligos(df[i,2])[1,1])
     results[i,3] <- as.character(oligos(df[i,2])[1,2])
     results[i,1] <- as.character(df[i,1])
   }
   
   write.csv(as.data.frame(results),file="tRNA-gRNA oligo.csv",quote=F)
  })
  
}

# Create Shiny app ----
shinyApp(ui, server)