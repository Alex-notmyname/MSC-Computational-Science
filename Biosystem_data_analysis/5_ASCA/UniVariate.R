## Module for uni-variate data exploration

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
uniVariateUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      box(status = "primary", width = 3,
          selectInput(ns("xvar"), 
                      "X-axis", 
                      choices = "", selected = "")
      ),
      box(status = "primary", width = 3,
          selectInput(ns("yvar"), 
                      "Variables", 
                      choices = "", 
                      selected = "", 
                      multiple = TRUE)
      ),
      box(status = "primary", width = 3,
          selectInput(ns("colorvar"), 
                      "Color by",
                      choices = "",
                      selected = "")
      )
    ),
    fluidRow(
      box(width = 12, status = "primary",
          plotOutput(ns("univariate.plot"))
      )
    )
  )
}
#
# Server
#

uniVariate <- function(input, output, session, DAT.session) {
  
  MakeUnivariateMenu <- observe({
    design.df <- DAT.session$design 
    data.df   <- DAT.session$data 
    factors   <- colnames(design.df)
    vars      <- colnames(data.df)
    updateSelectInput(session, "xvar", 
                      choices = factors,
                      selected = factors[1])
    updateSelectInput(session, "yvar",
                      choices = vars,
                      selected = vars[1])
    updateSelectInput(session, "colorvar",
                      choices = factors,
                      selected = factors[1])
  }) 
  
  PlotUnivariate       <- reactive({
    design.df          <- DAT.session$design 
    data.df            <- DAT.session$data 
    combined           <- cbind.data.frame(data.df, design.df)
   
    req(input$xvar, input$yvar, input$colorvar)
    data.long          <- gather(combined, variable, values, input$yvar, 
                                 factor_key = TRUE)
    
    # Order of the x-axis is in which variables appear in design file
    ordered.xvar       <- unique(design.df[, input$xvar])
    ordered.colorvar   <- unique(design.df[, input$colorvar])
    data.long$xvar     <- factor(data.long[ ,input$xvar], 
                                 levels = ordered.xvar)
    data.long$colorvar <- factor(data.long[ ,input$colorvar], 
                                 levels = ordered.colorvar)
    
    p <- ggplot(data = data.long, aes(x = xvar, y = values))
    p <- p + geom_boxplot(aes(fill = colorvar), alpha = 0.5, width = 0.4)
    p <- p + geom_point(aes(color = colorvar), 
                        size = 2, position = position_dodge(width = 0.4))
    p <- p + stat_summary(fun.y = "mean", geom = "point", size = 5, 
                            aes(color = colorvar),
                            position = position_dodge(width = 0.4))
    if(input$xvar != input$colorvar) {
    p <- p + stat_summary(fun.y = "mean", geom = "line", 
                          aes(x = xvar, y = values,
                          color = colorvar, group = colorvar), 
                          size = 1, position = position_dodge(width = 0.4)
                      )
    }
    p <- p + scale_fill_discrete(name = input$colorvar, 
                                 breaks = as.character(ordered.colorvar))
    p <- p + scale_color_discrete(name = input$colorvar, 
                                  breaks = as.character(ordered.colorvar))
    p <- p + theme(legend.text = element_text(size = 14), 
                   legend.title = element_text(size = 16),
                   axis.title = element_text(size = 14),
                   axis.text.x = element_text(size = 14),
                   axis.text.y = element_text(size = 14),
                   plot.title = element_text(size = 16, 
                                             face="bold", hjust = 0.5)
    )
    p <- p + labs(x = input$xvar)
    p <- p + facet_wrap( ~ variable)
    return(p)
  })
  
  output$univariate.plot <- renderPlot(PlotUnivariate())
}



