## Module for file input

# Gooitzen Zwanenburg, g.zwanenburg@uva.nl, June 2018
# Version: 1.0
#
# - Reads data- and design files
# - Checks if the input from the files is valid
# - Cleans the data: removes NA's
# - Summarizes the data

# Contents
#   readFileUI: dashboard layout
#   readFiles:  module function
#     getData:    fuction to read file
#     ReadInput:  reactive function
#     CleanData:  reactive function

#
# UI
#

readFilesUI <- function(id) {
  ns <- NS(id)
  tagList(
    box(width = 4, status = "primary",
      div(style="height: 80px",
        fileInput(ns("data.file"),                        # input data file,
                  "Choose data file",           
                  accept = c("csv",
                    "comma-separated-values",
                    ".csv"
                  )
        )
      ),
      checkboxInput(ns("Xheader"), "Header", TRUE),       # header yes/no      
      tags$hr(style="border-color: black;"),              # horizontal rule
      div(style="height: 80px", 
        fileInput(ns("design.file"),                      # input design file
                  "Choose design file",        
                  accept = c("csv",
                    "comma-separated-values",
                    ".csv"
                  )
        )
      ),
      checkboxInput(ns("Fheader"), "Header", TRUE),       # header yes/no
      tags$hr(style="border-color: black;")
    ),
    box(width = 5, status = "primary",
         title = "Data summary",
         tableOutput(ns("datasummary"))
    ),
    box(width = 3, status = "primary",
       title = "Factor levels",
       tableOutput(ns("factor.levels"))
    )
  )
}
#
# Server
#

readFiles <- function(input, output, session, DAT.session) {
  # - Reads data- and design file
  # - Checks consistency of data- and design files
  # - Removes rows with NA's
  # - Summarizes the data
  # Returns:
  #  ReadInput: (reactive) List with datafile, designfile, summary and warnings
  #
  
  
  getData <- function(file, header, type) {
    # Reads the data file, errors caught with tryCatch
    # Input: 
    #  - file: character, full file name of data file
    #  - header: boolean, data file with/withoun header
    #  - type: character, selects for data design file to read
    # Output:
    #  - on success: df.data: data frame with data
    #  - on error/warning: 1
    
    if(type == "data") {
      colClasses <- "numeric"       # Data are numeric
    } else if (type == "design") {
      colClasses <- "factor"        # Design columns are factors
    }
    res <- tryCatch(
      {
        read.csv(file,
                 header = header,
                 colClasses = colClasses, na.strings = "")
      },
      error = function(e) {
        return(1)
      },
      warning = function(w) {
        return(1)
      }
    )
    return(res)
  }  
  
  # Make empty summary dataframe
  data.summary           <- data.frame(c("Number of observations", 
                                         "Number of variables",
                                         "Number of NaN's",
                                         "Number of removed rows",
                                         "Number of factors",
                                         "Max number of observations / cell",
                                         "Min number of observations / cell",
                                         "Number of empty cells",
                                         "Is data balanced?")
  )
  
  ReadInput <- reactive({
    req(input$data.file)
    data                       <- input$data.file
    dataPath                   <- data$datapath
    raw.data                   <- getData(dataPath, input$Xheader, "data")
    validate(need(class(raw.data) == "data.frame", 
                         paste("Data file could not be read or",
                               "did not have the proper format.",
                               "Check if Header-setting is correct")))
    
    req(input$design.file)
    design                     <- input$design.file
    designPath                 <- design$datapath
    raw.design                 <- getData(designPath, input$Fheader, "design")
    validate(
      need(raw.design != 1, paste("Design file could not be read or", 
                                  "did not have the proper format")),
      need(sum(is.na(raw.design)) == 0, paste("Design file has empty fields"))
    )
    
    raw.dim                    <- dim(raw.data)
    # Remove rows with NA and empty fields from data file
    number.removed.rows        <- 0
    number.nan                 <- 0
    data.na                    <- as.data.frame(sapply(raw.data, function(x)
                                          suppressWarnings(as.numeric(
                                          as.character(x)))))
    number.nan                 <- sum(is.na(data.na))
    df.data                    <- raw.data[complete.cases(data.na), ]  
    number.removed.rows        <- sum(!complete.cases(data.na))
    validate(
      need(dim(df.data)[1] > 2, "Too few rows left after NaN's were removed"),
      need(dim(raw.data)[1] == dim(raw.design)[1], 
           paste("Design and data files have\n", 
                 "different numbers of rows. \n",
                 "Incorrect Header checkbox value?"))
    )
  
    df.design             <- raw.design[complete.cases(data.na), ] # design
  # Save headers ...
    if(input$Xheader) {
      DAT.session$data.headers   <- colnames(df.data)
    } else if (!input$Xheader) {
        colnames(df.data)        <- paste("Var", 
                                        as.character(seq(ncol(df.data))), 
                                        sep = "_")
        DAT.session$data.headers <- colnames(df.data) 
    } 

    if(input$Fheader) {
      DAT.session$design.headers <- colnames(df.design)
    } else if (!input$Fheaders) {
        colnames(df.design)        <- paste("Factor", 
                                          as.character(seq(ncol(df.design))), 
                                          sep = "_") 
        DAT.session$design.headers <- colnames(df.design)
    }
  
    DAT.session$n.factors  <- dim(df.design)[2]
    design                 <- lapply(seq(ncol(df.design)), function(x) 
                                as.numeric(df.design[, x]))
    
    # number of levels for each factor
    num.levels              <- unlist(lapply(design, max))
    level.numbers           <- cbind(DAT.session$design.headers, num.levels)
    colnames(level.numbers) <- c("Factor", "Levels")
    
    # Check for empty cells in the design
    #
    # This does not work check last column of table for 0.This is done a few lines down
    # Empty cells may arise because of 
    design.table <- tryCatch(  # Not sure why this is in a tryCatch
      {
        data.frame(table(df.design))
      },
      error = function(e) {
        return(1)
      }
    )
    # This is to catch design files that contain numeric data.
    validate(need(class(design.table) == "data.frame", 
                  paste("Number of factors is too large.\n", 
                        "Please check design file."))
    )
    empty.cells <- sum(design.table[, "Freq"] == 0) 
    cell.min  <- min(design.table[design.table[,"Freq"] > 0, "Freq"])
    cell.max  <- max(design.table[design.table[,"Freq"] > 0, "Freq"])
    if(cell.min == cell.max) {
      balanced <- "Yes"
    } else {
      balanced <- "No"
    }
    
    data.summary[1,2]      <- dim(raw.data)[1]
    data.summary[2,2]      <- dim(raw.data)[2]
    data.summary[3,2]      <- number.nan
    data.summary[4,2]      <- number.removed.rows
    data.summary[5,2]      <- DAT.session$n.factors
    data.summary[6,2]      <- cell.max
    data.summary[7,2]      <- cell.min
    data.summary[8,2]      <- empty.cells
    data.summary[9,2]      <- balanced
    colnames(data.summary) <- c("", "")  

  return(list(df.data, df.design, data.summary, level.numbers, raw.data, raw.design))
})
  
  output$datasummary   <- renderTable(ReadInput()[[3]])
  output$factor.levels <- renderTable(ReadInput()[[4]])
  
  return(ReadInput)
  
}
