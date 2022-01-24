## Module for data preporcessing

# Gooitzen Zwanenburg, g.zwanenburg@uva.nl, June 2018
# Version: 1.0
#
# - Balances, scales and centers the data

# Contents
#   preProcessUI: dashboard layout
#   preProcess: module function
#     balancingAct: balances data
#     scaling:  scale data
#     PreProcessData: reactive function
#
# UI
#

#' Title
#'
#' @param id 
#'
#' @return
#' @export
#'
#' @examples
preProcessUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      box(width = 3, status = "primary",
          radioButtons(ns("balance"), "Balancing",
                       c("No balancing" = "nobalancing",
                         "Replace by means" = "means",
                         "Trim" = "trim", 
                         "Fill" = "fill",
                         "Fill positive" = "fill.pos"),
                       selected = "nobalancing"
          ),      
          tags$hr(style="border-color: black;"),
          uiOutput(ns("removeempties"))
      ),
      box(width = 3, status = "primary",
          radioButtons(ns("transform"), "Transformation",
                       c("No transfromation" = "notransformation",
                         "Square root" = "squareroot",
                         "Log" = "log"),
                       selected = "notransformation"
          ),
          tags$hr(style="border-color: black;"),
          checkboxInput(ns("center"), "Center", TRUE)
      ),
      box(width = 3, status = "primary",
          radioButtons(ns("scale"), "Scaling",
                       c("No scaling" = "noscaling",
                         "Standardize" = "standardize",
                         "Pareto" = "pareto"),
                       selected = "noscaling"
          )
      ),
      column(width = 3,
             fluidRow(
        box(width = 12, status = "primary",
            div(style="display: inline-block;vertical-align:top; width: 100px;",
                radioButtons(ns("plottype"), 
                             "Plot",
                             c("Data"   = "data",
                               "Scores" = "scores",
                               "Biplot" = "biplot",
                               "Loadings" = "loadings"),
                             selected = "data"
                )
            )
        )
             ),
        fluidRow(
          conditionalPanel(paste0("input['", ns("plottype"), "'] != 'loadings' "),
            box(width = 12, status = "primary",
                div(style="display: inline-block;vertical-align:top; width: 100px;",
                    uiOutput(ns("colorlevels"))
                )
            )
          )
        )
  #      fluidRow(
  #        conditionalPanel(paste0("input['", ns("plottype"), "'] == 'biplot' "),
  #          box(width = 12, status = "primary",
  #              div(style="display: inline-block;vertical-align:top; width: 100px;",
  #                  uiOutput(ns("biplotscaling"))
  #              )
  #          )
  #        )
  #      )
      )
    ),
    fluidRow(
      box(width = 12, status = "primary",
          plotOutput(ns("pca.plot"))
      )
    )
  )
}
#
# Server
#

preProcess <- function(input, output, session, data, DAT.session) {
  # - Balances, scales and centers the data
  # - Plots score plot of whole data set
  # Returns
  # 
  #

  balancingAct <- function(design.df, data.df) {
    # Balance cleaned data 1. as is, i.e. balanced or unbalanced
    #                      2. trim: remove rows to minimum cell size
    #                      3. means: replace cells by their means, each cell has a single value
    #                      4. fill: cells are filled from normal distribution
    #                      5. fill positive: as fill but with only positive numbers
    # Args:
    #   design.df: design dataframe
    #   data.df: data dataframe
    #
    # Returns:
    #   design.df: balanced design dataframe
    #   balanced.df: balanced data dataframe
    # Note: balanced.df may not actually be balanced if 'nobalancing' is 
    #       selected or if the 'fill' option does not succeed.
    
    # merge design and data
    df                    <- cbind(design.df, data.df)
    # Set up empty dataframe for balanced data
    df.balanced           <- data.frame(matrix(nrow=0, 
                                               ncol=length(colnames(df))))
    colnames(df.balanced) <- colnames(df)
    # Find multiplicity of the rows
    tmp                   <- ddply(design.df, colnames(design.df), nrow)
    cell_numbers          <- tmp[, ncol(tmp)]     
    cell_max              <- max(cell_numbers)    
    cell_min              <- min(cell_numbers)   
    unique.design.df      <- unique.data.frame(design.df)   
    unique.rows           <- rownames(unique.design.df)
    
    # select balancing method
    if(input$balance == "nobalancing") {         
      df.balanced      <- df
    }
    # all cells minimum number of measurements
    else if (input$balance == "trim") {
      df.balanced <- suppressMessages(do.call(rbind, lapply(unique.rows, function(x)
        df[sample(rownames(match_df(design.df, design.df[x, ])),
                  cell_min), ])))
    }
    # Replace with cell means
    else if (input$balance == "means") {
      df.balanced <- suppressMessages(do.call(rbind, lapply(unique.rows, function(x)
        cbind(design.df[x, ], as.data.frame.list(colMeans(
          data.df[rownames(match_df(design.df, design.df[x, ])),  ]))))))
    }
    
    
    # fill all cells to maximum number of obserservations per cell (cell_max)
    else if (input$balance == "fill" || input$balance == "fill.pos") {
  
      validate(need(cell_min >= 2, paste("To fill, each cell must have",
                                         "at least two measurments")))
      for (r in unique.rows) {
        suppressMessages(m <- rownames(match_df(design.df, design.df[r, ])))
        df.balanced        <- rbind(df.balanced, df[m, ])
        cell.mean          <- colMeans(data.df[m, ])
        cell.sd            <- as.data.frame.list(sapply(data.df[m, ], sd))
        num.to.add         <- cell_max - length(m)
        if (num.to.add > 0) {
          for (i in seq(num.to.add)) {
            to.add         <- as.data.frame.list(mapply(function(x, y)
                                {rnorm(1, mean=x, sd=y)}, x=cell.mean, y=cell.sd))
            if(input$balance == "fill.pos") {
              to.add     <- abs(to.add)
            }
            df.balanced  <- rbind(df.balanced, cbind(design.df[r, ], to.add))
          }
        }
      }
    }
    design.df <- df.balanced[, (1 : ncol(design.df))]
    data.df   <- df.balanced[, -(1 : ncol(design.df))]
    
    return(list(design.df, data.df))
  }
  
  transforming <- function(data.df) {
    # Transforms the data matrix
    # Args:
    #  data.df: (Balanced) dataframe
    # Returns:
    #  data.scaled: scaled dataframe
    
    # No Transform
    if(input$transform == "notransformation") {
      data.transformed <- data.df
    # Square root scaling
    } else if (input$transform == "squareroot") {
      data.transformed <- as.data.frame(sapply(
        data.df, function(x) sign(x)*sqrt(abs(x)))
      )
    # 10 log scaling
    } else if (input$transform == "log") {
    # Replace 0's by 10 times machine precision
      tol <- .Machine$double.eps * 10     
      data.log    <- data.df
      data.log[abs(data.log) < tol] <- tol
      data.transformed <- as.data.frame(sapply(
        data.log, function(x) log10(abs(x + 1)))
      )
    }
    return(data.transformed)
  }
  
  scaling <- function(data.df) {
    # Scales the data matrix
    # Args:
    #  data.df: (Balanced) dataframe
    # Returns:
    #  data.scaled: scaled dataframe
    
    # No scaling
    if(input$scale == "noscaling") {
      data.scaled <- data.df
      # Standardize
    } else if (input$scale == "standardize") {
      data.scaled <- as.data.frame(scale(
                       data.df, center = TRUE, scale = TRUE)
                     )
      # Pareto scaling
    } else if (input$scale == "pareto") {
      # first center the data
      data.centered <- sweep(data.df, 2, colMeans(data.df), "-")
      data.scaled   <- as.data.frame(t(
                         t(data.centered) / sqrt(apply(data.centered, 2, sd)))
                     )
    }
    return(data.scaled)
  }
  
 PreProcessData <- reactive({
    # Balances, scales and centers data.
    # Removes selected rows
    # Returns:
    #  List with designfile, balanced datafile and scaled datafile

  df.data          <- data()[[1]]
  df.design        <- data()[[2]]
  raw.data         <- data()[[5]]

  # Balance data
  df.balanced      <- balancingAct(df.design, df.data)
  design.balanced  <- df.balanced[[1]]
  data.balanced    <- df.balanced[[2]]
 
  data.transformed <- transforming(data.balanced)
  
  # Center data
  if(input$center) {
    data.transformed   <- sweep(data.transformed, 2, 
                                colMeans(data.transformed), "-")
  }
  
  # Scale balanced data
  data.scaled      <- scaling(data.transformed)
  
  return(list(design.balanced, data.balanced, data.scaled))
  }) 
 
  observe({
     
    plotScores <- function(scores, levels, perc.explnd, selected.factor) {
      g <- ggplot(data = scores, aes(x = PC_1, y = PC_2,
                                     col = levels)
      )
      g <- g + geom_point(size = 3)
      g <- g + labs(title = "Scores",
                    x = paste("PC", "1", 
                              "(", sprintf("%.2f", perc.explnd[1]),"%", ")",   
                              sep = " "),
                    y = paste("PC", "2",  
                              "(", sprintf("%.2f", perc.explnd[2]),"%", ")",  
                              sep = " " )
      )
     # g <- g + stat_ellipse()
      g <- g + scale_color_discrete(name = selected.factor)
      return(g)
    }
    
    plotBiplot <- function(u, d, v, levels, selected.factor, 
                           pc1 = 1, pc2 = 2, type = "1") {
      # Plot US vs V
      if(type == "1") {
        scores             <- as.data.frame(u %*% diag(d))[, c(pc1, pc2)]
        loadings           <- as.data.frame(v[, c(pc1, pc2)])
      # Plot U*sqrt(nrow -1) vs VS/sqrt(nrow - 1) Gabriel: PCA biplot
      } else if (type == "2") {
        scores           <- as.data.frame(u*sqrt(nrow(u) - 1))[, c(pc1, pc2)]
        loadings         <- as.data.frame( (v %*% diag(d))/
                                             sqrt(nrow(u) - 1))[, c(pc1, pc2)]
      }
      
      perc.explnd        <- 100*d[c(pc1, pc2)]^2 / sum((d)^2 )
      colnames(scores)   <- paste("PC", c(pc1, pc2), sep = "")
      colnames(loadings) <- paste("L", c(pc1, pc2), sep = "")
      max.scores         <- max(sapply(seq(nrow(scores)), function(x) sum(scores[x, ]^2)))
      max.loadings       <- max(sapply(seq(nrow(loadings)), function(x) sum(loadings[x, ]^2)))
      lambda             <- sqrt(max.loadings)/sqrt(max.scores)
      scaled.scores      <- scores
      scaled.loadings    <- loadings/lambda
      
      g <- ggplot(data = scaled.scores, aes_string(x = "PC1", 
                                           y = "PC2", 
                                           col = levels))
      g <- g + geom_point(size = 3)
      g <- g + labs(title = "Biplot",
                     x = paste("PC", pc1, 
                               "(", sprintf("%.2f", perc.explnd[1]),"%", ")",   
                               sep = " "),
                     y = paste("PC", pc2,  
                               "(", sprintf("%.2f", perc.explnd[2]),"%", ")",  
                               sep = " " )
                   )
      
      g <- g + geom_segment(data = scaled.loadings,
                              aes(x = 0, y = 0, xend = L1, yend = L2),
                              arrow = arrow() , color = "red", size = 1
                           )
      
      g <- g + geom_text_repel(data = scaled.loadings, 
                         aes(x = L1, y = L2, 
                             label = colnames(DAT.session$data)), 
                          col = "black",
                         show.legend = FALSE)
      g <- g + scale_color_discrete(name = selected.factor)
     
      return(g)
    }
    
    plotPCA <- function(pca.data, design, levels, selected.factor, plottype) {
      svd.data             <- svd(pca.data)
      u                    <- svd.data$u
      d                    <- svd.data$d
      v                    <- svd.data$v
      scores               <- as.data.frame(u %*% diag(d))
      loadings             <- as.data.frame(v[, c(1,2)])
      colnames(scores)     <- paste("PC",
                                    as.character(seq(ncol(scores))), sep = "_")
      colnames(loadings)   <- paste("PC",
                                    as.character(seq(ncol(loadings))), sep = "_")
      perc.explnd          <- 100*d[c(1, 2)]^2 / sum((d)^2 )
      
      if(plottype == "data") {
        data.df            <- cbind(pca.data, design)
   
        data.long          <- gather_(data.df, "variable", "values", 
                                      colnames(pca.data), factor_key = TRUE)
        unique.levels    <- unique(data.long[, selected.factor])
        color.levels     <- factor(data.long[, selected.factor], levels = unique.levels)
        p
        
        g <- ggplot(data = data.long, aes(x = variable, y = values, 
                                          col = color.levels)) 
        g <- g + geom_jitter(width = 0.2, height = 0)
        g <- g + scale_color_discrete(name = selected.factor)
      }

      if(plottype == "scores") {
        g <- plotScores(scores, levels, perc.explnd, selected.factor ) 
        
      } else if(plottype == "biplot") {
        g <- plotBiplot(u, d, v, levels, selected.factor, type = 1) 
        
      } else if (plottype == "loadings") {
        plot.loadings           <- cbind.data.frame(colnames(DAT.session$data), loadings)
        colnames(plot.loadings) <- c("Variable", "L1", "L2")
        plot.loadings           <- gather(plot.loadings, Loading, value, "L1", "L2")
        
        plot.title              <- "Loadings"
         
        
        g <- ggplot(data = plot.loadings, 
                    aes(x = factor(Variable), y = value, fill = Loading)) +
          geom_bar(stat = "identity", 
                   position = position_dodge(preserve = "total")) +
          labs(title = plot.title, x = "Variable", y = "Loadings")
    
      }
      g <- g + theme(legend.text = element_text(size = 14), 
                     legend.title = element_text(size = 16),
                     axis.title = element_text(size = 14),
                     axis.text.x = element_text(size = 14),
                     axis.text.y = element_text(size = 14),
                     plot.title = element_text(size = 16, 
                                               face="bold", hjust = 0.5)
                    )
      return(g)
    }
 

    data.design        <- PreProcessData()[[1]]
    data.scaled        <- PreProcessData()[[3]]    

    color.by           <- colnames(data.design)
    output$colorlevels <- renderUI({ 
                            selectInput(session$ns("colorlevels"), 
                            "Color by ",
                            choices = color.by,
                            selected = input$colorlevels)
    }) 
    output$biplotscaling <- renderUI({
      radioButtons(session$ns("biplotscaling"), 
                   "Biplot scaling",
                   c("Plot U vs VS" = "1",
                     "plot US vs V" = "2"),
                   selected = "1")
    })
  
    req(input$colorlevels)
    index            <- which(color.by == input$colorlevels)
    unique.levels    <- unique(data.design[, index])
    color.levels     <- factor(data.design[, index], levels = unique.levels)
    output$pca.plot  <- renderPlot(plotPCA(data.scaled, data.design, color.levels, 
                                          input$colorlevels, input$plottype))
    DAT.session$design <- data.design
    DAT.session$data   <- data.scaled
  })
}



