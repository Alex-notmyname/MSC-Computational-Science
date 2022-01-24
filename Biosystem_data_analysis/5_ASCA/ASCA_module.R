# ASCA Model

buildModelUI <- function(id) {
  ns <- NS(id)
  tagList(
    tabItem("buildmodel",
            fluidRow(
              box(width = 3, status = "primary",
                  textInput(ns("include.factors"), "Include factors", "")
              ),
              box(width = 3, status = "primary",
                  textInput(ns("include.interactions"), "Include interactions", "")
              ),
              box(width = 3, status = "primary",
                  textInput(ns("include.combinations"), "Combine terms", "")
              )
            ),
            fluidRow(
              box(width = 9, status = "primary", title = "Model",
                  tableOutput(ns("model"))
              )
            ),
            fluidRow(
              box(width = 9, status = "primary", title = "Variances",
                  tableOutput(ns("variances"))
              )
            )
    )
  )
}

makeFactorsUI <- function(id) {
  ns <- NS(id)
  fluidPage(
    fluidRow(
      box(status = "primary", width = 3,
          selectInput(ns("choose.factor"), 
                      "Factor", 
                      choices = "", selected = ""),
          conditionalPanel(condition = paste0("input['", ns("fpt"), "'] == 'biplot' |
                                               input['", ns("fpt"), "'] == 'scores' "),
            tags$hr(style="border-color: black;"),
            checkboxInput(ns("fplotprojctns"), "Projections", TRUE)
          )
      ),
      box(status = "primary", width = 3,
          radioButtons(ns('fpt'), "Plot", 
                       c("Scores" = "scores", 
                         "Levels"   = "levels",
                         "Biplot"   = "biplot",
                         "Loadings" = "loadings"),
                       selected   = "scores")
      ),
      box(width = 3, status = "primary",
          div(style="display: inline-block;vertical-align:top; width: 100px;",
              selectInput(ns("factor.pc1"), 
                          "First PC",
                          choices  = "",
                          selected = "")
          ),
          div(style=
                "display: inline-block;vertical-align:top; width: 20px; ",""),
          div(style="display: inline-block;vertical-align:top; width: 100px;",
              selectInput(ns("factor.pc2"), 
                          "Second PC",
                          choices  = "",
                          selected = "")
          )
      )
  #    conditionalPanel(condition = paste0("input['", ns("fpt"), "'] == 'biplot' "),
  #      box(width = 3, status = "primary", 
  #         div(style="display: inline-block;vertical-align:top; width: 100px;",
  #             uiOutput(ns("fbiplot"))
  #         )
  #      )
  #    )
    ),
    fluidRow(
      box(width = 12, status = "primary",
          plotOutput(ns("factors"))
      )
    )
  )
}

makeInteractionsUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(   
      box(status = "primary", width = 3,
          selectInput(ns("choose.interaction"), 
                      "Interaction", 
                      choices = "", selected = ""),
          conditionalPanel(condition = paste0("input['", ns("ipt"), "'] == 'biplot' |
                                               input['", ns("ipt"), "'] == 'scores' "),
                           tags$hr(style="border-color: black;"),
                           checkboxInput(ns("iplotprojctns"), "Projections", TRUE)
          )
      ),
      box(status = "primary", width = 3,
          radioButtons(ns('ipt'), 
                       "Plot", 
                       c("Scores" = "scores", 
                         "Levels"   = "levels",
                         "Biplot"  = "biplot",
                         "Loadings" = "loadings"),
                       selected   = "scores")
      ),
      box(width = 3, status = "primary",
          div(style="display: inline-block;vertical-align:top; width: 100px;",
              selectInput(ns("interaction.pc1"), 
                          "First PC",
                          choices  = "",
                          selected = "")
          ),
          div(style="display: inline-block;vertical-align:top; width: 20px;",""),
          div(style="display: inline-block;vertical-align:top; width: 100px;",
              selectInput(ns("interaction.pc2"), 
                          "Second PC",
                          choices  = "",
                          selected = "")
          )
      ),
      conditionalPanel(condition = paste0("input['", ns('ipt'), "'] == 'levels' "),
                       box(width = 3, 
                           selectInput(ns("ixvalues"),
                                       "Factor on x-axis",
                                       choices  = "",
                                       selected = ""),
                           selectInput(ns("igroupfactor"),
                                       "Color by factor",
                                       choices  = "",
                                       selected = "")
                       )
      ),
      conditionalPanel(condition = paste0("input['", ns('ipt'), "'] == 'scores' |
                                           input['", ns('ipt'), "'] == 'biplot'"),
                       box(width = 3, 
                           selectInput(ns("icolor"),
                                       "Color",
                                       choices  = "",
                                       selected = ""),
                           selectInput(ns("ishape"),
                                       "Shape",
                                       choices  = "",
                                       selected = "")
                       )
      )
#      conditionalPanel(condition = paste0("input['", ns("ipt"), "'] == 'biplot' "),
#                       box(width = 3, status = "primary", 
#                           div(style="display: inline-block;vertical-align:top; width: 100px;",
#                               uiOutput(ns("ibiplot"))
#                           )
#                       )
#      )
    ),
    fluidRow(
      box(width = 12, status = "primary",
          plotOutput(ns("interactions"))
      )
    )
  )
}

makeCombinationsUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(          
      box(status = "primary", width = 3,
          selectInput(ns("choose.combination"), 
                      "Combination", 
                      choices = "", selected = ""),
          conditionalPanel(condition = paste0("input['", ns("cpt"), "'] == 'biplot' |
                                               input['", ns("cpt"), "'] == 'scores' "),
                           tags$hr(style="border-color: black;"),
                           checkboxInput(ns("cplotprojctns"), "Projections", TRUE)
          )
      ),
      box(status = "primary", width = 3,
          radioButtons(ns('cpt'), "Plot", 
                       c("Scores" = "scores", 
                         "Levels"   = "levels",
                         "Biplot"   = "biplot",
                         "Loadings" = "loadings"),
                       selected   = "scores")
      ),
      box(width = 3, status = "primary",
          div(style="display: inline-block;vertical-align:top; width: 100px;",
              selectInput(ns("combination.pc1"), 
                          "First PC",
                          choices = "",
                          selected = "")
          ),
          div(style=
                "display: inline-block;vertical-align:top; width: 20px;", ""),
          div(style=
                "display: inline-block;vertical-align:top; width: 100px;",
              selectInput(ns("combination.pc2"), 
                          "Second PC",
                          choices = "",
                          selected = "")
          )
      ),
      conditionalPanel(paste0("input['", ns("cpt"), "'] == 'levels' "),
                       box(width = 3, 
                           selectInput(ns("cxvalues"),
                                       "Factor on x-axis",
                                       choices = "",
                                       selected = ""),
                           selectInput(ns("cgroupfactor"),
                                       "Color by factor",
                                       choices = "",
                                       selected = "")
                       )
      ),
      conditionalPanel(condition = paste0("input['", ns('cpt'), "'] == 'scores' | 
                                           input['", ns('cpt'), "'] == 'biplot'"),
                       box(width = 3, 
                           selectInput(ns("ccolor"),
                                       "Color",
                                       choices  = "",
                                       selected = ""),
                           selectInput(ns("cshape"),
                                       "Shape",
                                       choices  = "",
                                       selected = "")
                       )
      )
#     conditionalPanel(condition = paste0("input['", ns("cpt"), "'] == 'biplot' "),
#                       box(width = 3, status = "primary", 
#                           div(style="display: inline-block;vertical-align:top; width: 100px;",
#                               uiOutput(ns("cbiplot"))
#                           )
#                       )
#      )
    ),
    fluidRow(
      box(width = 12, status = "primary",
          plotOutput(ns("combination"))
      )
    )
  )
}

makeResidualsUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      box(status = "primary", width = 3,
          radioButtons(ns("rpt"), "Plot", 
                       c("Scores"      = "scores", 
                         "Loadings"    = "loadings"),
                       selected        = "scores")
      ),
      box(width = 3, status = "primary",
          div(style="display: inline-block;vertical-align:top; width: 100px;",
              selectInput(ns("residual.pc1"), 
                          "First PC",
                          choices = "",
                          selected = "")
          ),
          div(style=
                "display: inline-block;vertical-align:top; width: 20px;", ""),
          div(style="display: inline-block;vertical-align:top; width: 100px;",
              selectInput(ns("residual.pc2"), 
                          "Second PC",
                          choices = "",
                          selected = "")
          )
      ),
      conditionalPanel(paste0("input['", ns("rpt"), "'] == 'scores' "),
                       box(width = 3, status = "primary",
                           selectInput(ns("colorlevels"), 
                                       label = "Color by factor",
                                       choices = "",
                                       selected = "")
                       )
      )
    ),
    fluidRow(
      box(width = 12, status = "primary",
          plotOutput(ns("residuals"))
      )
    )
  )
}




#
# Server
#

ASCAmodel <- function(input, output, session, DAT.session) {
  # List of all terms in the model.
  # TermList has three elements: TermList$factors, TermList$interactions,
  # and TermList$combinations
  
  MakeTermsList <- function(input.string, fic) {
    # Function make a list of terms to include in the model
    # The function recognizes factor, interaction and combination input
    # Input is checked on validity and for duplicate terms.
    # interactions of combinations
    # Args:
    # Input: 
    #  input.string: user input
    #  fic: character: f factor, i interaction, c combination
    # Returns: checked.terms: list with terms to include
    
    ### Auxillary functions ###
    CheckInteraction <- function(st) {
      # Gets factors from an interaction
      interacting.factors <- as.numeric(unlist(str_extract_all(st, "[0-9]+")))
      # check for duplicate (eg. 2:2) and out of range factors
      if(any(duplicated(interacting.factors)) | 
         any(interacting.factors > DAT.session$n.factors) |
         any(interacting.factors <= 0)) {
        return(NULL)
      } else {
        return(paste(sort(interacting.factors), collapse = ":"))
      }
    }
    
    CheckFactors <- function(st) {
      # Makes a comma separated list of factors from a range n-m
      tmp.list    <- str_extract_all(st, "[0-9]+")
      from        <- as.numeric(tmp.list[[1]][1])
      to          <- as.numeric(tmp.list[[1]][2])
      if(length(tmp.list[[1]]) == 1) {
        to <- from
      }
      fac.seq <- seq(from = from, to = to)
      return(fac.seq[(fac.seq <= DAT.session$n.factors) & (fac.seq > 0)])
    }
    
    CheckCombinations <- function(st) {
      # Check individual combinations
      tmp.terms <- unique(unlist(sapply(st, function(x) 
        str_split(x, "\\+"), USE.NAMES = FALSE)))
      
      factor.or.ineraction <- function(s) {
        # check interaction or factor
        if(grepl(":", s, fixed = TRUE)) {
          res <- CheckInteraction(s) 
        } else {
          res <- CheckFactors(s)
        }
        return(res)
      }
      
      combination <- unlist(lapply(tmp.terms, factor.or.ineraction))
      return(paste(combination, collapse = "+"))    
    }
    
    ### Main code of MakeTermsList() ###
    
    # Replace space as separater by ","
    # Do this twice as the replacement works on pairs 
    # Allow for letter typo's, they are removed later on
    tmp <- gsub("([a-zA-Z0-9]) +([a-zA-Z0-9])", "\\1,\\2", input.string)
    tmp <- gsub("([a-zA-Z0-9]) +([a-zA-Z0-9])", "\\1,\\2", tmp)
    
    # Remove remaining spaces
    tmp <- gsub(" +", "", tmp)
    # combinations
    if(fic == "c") {
      # get individal combinations
      combinations   <- unique(strsplit(tmp, ",", fixed = TRUE)[[1]])
      tmp.terms      <- unique(unlist(lapply(combinations, function(x) 
        str_extract_all(x,
                        "^([0-9]+)(:[0-9]+)*(\\+[0-9]+(:[0-9]+)*)+$")[[1]])
      )
      )
      
      checked.terms  <- lapply(tmp.terms, CheckCombinations)
      
      # interactions
    } else if(fic == "i") { 
      # extract interactions separated by ","
      interactions   <- unique(strsplit(tmp, ",", fixed = TRUE)[[1]])
      # in case some non-numerics made it through
      tmp.terms      <- unique(unlist(lapply(interactions, function(x) 
        str_extract_all(x,
                        "^[0-9]+(:[0-9]+)+$")[[1]])
      )
      )
      checked.terms  <- as.list(unique(unlist(lapply(tmp.terms, 
                                                     CheckInteraction))))
      
      # factors
    } else if(fic == "f") {
      factors        <- unique(strsplit(tmp, ",", fixed = TRUE)[[1]])
      tmp.terms      <- unique(unlist(lapply(factors, function(x)
        str_extract_all(x, 
                        "^[0-9]+([-][0-9])+$|^[0-9]+$")[[1]])
      )
      )
      checked.terms  <- as.list(
        as.character(sort(unique(
          unlist(
            lapply(
              tmp.terms, CheckFactors), 
            use.names = FALSE)
        )
        )
        )
      )
    }
    else {    # No numbers
      checked.terms  <- list()
    }
    
    return(checked.terms)
  }
  
  MakeModelTerms <- function(f.list, i.list, c.list) {
    # Makes the model from factors, interactions and combinations.
    # In the model, the combination list is leading:
    # if terms appear in both the combination list (c.list) and 
    # factor (f.list) or interaction list (i.list) they are remvoed 
    # from the latter two.
    # Args:
    #  f.list: list with factors
    #  i.list: list with interactions
    #  c.list: list with combinations
    # Returns:
    #  list(f.list, i.list, c.list): list with updated factors, interactions
    #  and combinations
    
    # check for overlapping interactions
    if(length(c.list) > 0 & length(i.list) > 0) { 
      # get interactions out
      ints     <- lapply(c.list, function(x) 
        unlist(str_extract_all(x, "[0-9]+(:[0-9]+)+"))
      )
      tmp.list <- lapply(i.list, function(x)
        sapply(ints, function(y) any((x %in% y)))
      )
      mask     <- lapply(tmp.list, function(x) Reduce("|", x))
      i.list   <- i.list[!unlist(mask)]
    }
    
    # check for overlapping factors
    if(length(c.list) > 0 & length(f.list) > 0) { 
      # remove interactions
      facs     <- lapply(c.list, function(x) gsub("[0-9]+(:[0-9]+)+", "", x))
      tmp.list <- lapply(facs, 
                         function(x) unlist(str_extract_all(x, "[0-9]+"))
      )
      tmp.list <- lapply(f.list, function(x) 
        sapply(tmp.list, function(y) any((x %in% y)))
      )
      mask     <- lapply(tmp.list, function(x) Reduce("|", x))
      f.list   <- f.list[!unlist(mask)]
    }
    return(list(f.list, i.list, c.list))
  }
  
  MakeCells <- function(col.nrs) {
    # Select rows that belong to the same cell
    # Args:
    #  col.nrs: vector with columns numbers from the design matrix. 
    #           The columns are combined to make cells
    # Returns:
    #  list: filled.cells: row assignments for each level
    #        levels: vector with for each row the levels
    
    
    # create empty dataframe with same size as the data matrix
    # Make selected columns of the design matrix numeric
    design        <- lapply(col.nrs, function(x) 
      as.numeric(DAT.session$design[, x])
    )
    design.matrix <- sapply(design, cbind)
    
    # number of levels for each column
    num.levels    <- sapply(design, function(x) length(unique(x)))
    
    # all level combinations
    i.grid        <- as.matrix(expand.grid(lapply(num.levels, seq)))
    
    # Get matching rows from the design matrix
    matching.rows <- lapply(seq(nrow(i.grid)), function(x) 
      which(sapply(seq(nrow(design.matrix)), function(y) 
        identical(design.matrix[y, ], as.numeric(i.grid[x, ])))
      )
    )
    
    filled.cells  <- matching.rows[sapply(matching.rows, function(x) 
      length(x) > 0)]
    row.levels    <- vector(length = nrow(design.matrix))
    
    for(i in seq(length(filled.cells))) {
      for(j in filled.cells[i]) {
        row.levels[j] <- i
      }
    }
    return(list(filled.cells, as.factor(row.levels)))
  }
  
  MakeCellAverages <- function(col.nrs) {
    # Makes cell averages for all columns over selected rows
    # Args:
    #  col.nrs: vector with columns numbers from the design matrix. 
    #           The columns are combined to make cells and averages over the
    #           cells are calculated.
    # Calls: MakeCells
    # Returns:
    #  cell.averages: data frame with cell averages
    
    # create empty dataframe with same size as the data matrix
    cell.averages <- as.data.frame(matrix(integer(0),
                                          nrow = dim(DAT.session$data)[1],
                                          ncol = dim(DAT.session$data)[2]) 
    )
    filled.cells  <- MakeCells(col.nrs)[[1]]
    
    # collect the cell average in selected rows
    for(x in filled.cells) {
      # Note: use rbind to get a matrix out
      cell.averages[x, ] <- as.data.frame(rbind(colMeans(DAT.session$data[x, ])))
    }
    return(cell.averages)
  }
  
  MakeFactorMatrix <- function(factor) {
    matrix               <- MakeCellAverages(as.numeric(factor))
    overall.means.matrix <- matrix(rep(1, nrow(DAT.session$data)))%*%
      colMeans(DAT.session$data)
    factor.matrix        <- matrix - overall.means.matrix
    return(factor.matrix)
  }
  
  MakeInteractionMatrix <- function(interaction) {
    # Makes interaction matrix from factors
    # Args:
    #  interaction: interaction term (character string)
    # Returns:
    #  interaction.matrix: list: matrix of the interaction
    
    # Make all possible combinations with 'factors':
    # A, B, c, AB, BC, AC, ABC
    
    factors              <- as.numeric(unlist(str_extract_all(
      interaction, "[0-9]+"))
    )
    combinations         <- lapply(seq(length(factors)), function(x) 
      combn(factors, x)
    )
    
    # Make matrix with level averages for all combinations
    matrices             <- lapply(combinations, function(x) 
                              lapply(seq(ncol(x)), function(y) 
                                MakeCellAverages(x[, y])
                              )
                             )
    
    # Make interaction matrix by adding matrices 
    # with level averages with the right sign
    # e.g. ABC - AB - AC - BC + A + B + C
    interaction.matrix   <- Reduce("+", 
                                   lapply(seq(
                                     from = length(combinations), to = 1, by = -1), 
                                     function(x) ((-1)^(length(combinations) - x))*
                                       Reduce("+", matrices[[x]])
                                   )
    )
    overall.means.matrix <- matrix(rep(1, nrow(DAT.session$data)))%*%
      colMeans(DAT.session$data)
    
    # Add or subtract overall means from the interaction matrix
    interaction.matrix   <- interaction.matrix + 
      ((-1)^(length(factors)))*overall.means.matrix
    return(interaction.matrix)
  }
  
  MakeCombinationMatrix <- function(combination) {
    
    MakeMatrix <- function(str) {
      if(grepl(":", str)) {
        matrix <- MakeInteractionMatrix(str)
      } else {
        matrix <- MakeFactorMatrix(str)
      }
      return(matrix)
    }
    
    # split on '+'
    tmp.split <- unlist(strsplit(combination, "\\+" ))
    matrix <- Reduce("+", lapply(tmp.split, function(x) MakeMatrix(x)))
    return(matrix)
  }
  
  MakeResiduals <- function(model) {
    # Make matrix with residuals
    # Args:
    #  factors: factor matrices
    #  interactions: interaction matrices
    #  combinations: combination matrices
    #  data.df: scaled data matrix
    # Returns:
    #  residuals: matrix with residuals
    
    residuals   <- DAT.session$data
    residuals   <- Reduce("-", init = residuals, model$means$matrices)
    residuals   <- Reduce("-", init = residuals, model$factors$matrices)
    residuals   <- Reduce("-", init = residuals, model$interactions$matrices)
    residuals   <- Reduce("-", init = residuals, model$combinations$matrices)
    return(list(residuals))
  }
  
  
  
  ##################################
  
  MakeModel <- reactive({
    req(DAT.session$data)
    ASCA                       <- list()
    
    # Terms
    f.list                     <- MakeTermsList(input$include.factors, "f")
    i.list                     <- MakeTermsList(input$include.interactions, "i")
    c.list                     <- MakeTermsList(input$include.combinations, "c")
    model.terms                <- MakeModelTerms(f.list, i.list, c.list)
    factors                    <- model.terms[[1]]
    interactions               <- model.terms[[2]]
    combinations               <- model.terms[[3]] 
  
    ## Model ##
    
    ASCA$data$ssq              <- lapply(list(DAT.session$data), function(x) 
      sum(x^2)
    )
    # Overall means
    ASCA$means$matrices        <- list(matrix(rep(1, nrow(DAT.session$data)))%*%
                                         colMeans(DAT.session$data))
    ASCA$means$ssq             <- lapply(ASCA$means$matrices, function(x) 
      sum(x^2)
    )
    
    # Factors
    ASCA$factors$terms         <- factors
    ASCA$factors$matrices      <- lapply(factors, function(x) 
                                    MakeFactorMatrix(x)
                                  )
    ASCA$factors$levels        <- lapply(factors, function(x) 
                                    MakeCells(as.numeric(x))[[2]]
                                  )
    ASCA$factors$ssq           <- lapply(ASCA$factors$matrices, function(x) 
                                    sum(x^2)
                                  )
    # Interactions
    ASCA$interactions$terms    <- interactions
    ASCA$interactions$matrices <- lapply(interactions, function(x) 
                                    MakeInteractionMatrix(x)
                                  )
    ASCA$interactions$levels   <- lapply(interactions, function(x)
                                    MakeCells(as.numeric(unlist(
                                    str_extract_all(x, "[0-9]+")))
                                  )[[2]]
    )
    ASCA$interactions$ssq      <- lapply(ASCA$interactions$matrices, 
                                    function(x) sum(x^2)
                                  )
    # Combinations
    ASCA$combinations$terms    <- combinations
    ASCA$combinations$matrices <- lapply(combinations, function(x) 
                                    MakeCombinationMatrix(x)
                                  )
    ASCA$combinations$levels   <- lapply(combinations, function(x)
                                    MakeCells(as.numeric(unlist(
                                    str_extract_all(x, "[0-9]+")))
                                  )[[2]]
    )
    ASCA$combinations$ssq      <- lapply(ASCA$combinations$matrices, 
                                    function(x) sum(x^2)
                                  )
    # Residuals
    ASCA$residuals$terms       <- lapply(seq(ncol(DAT.session$design)), 
                                    function(x) as.character(x))
    ASCA$residuals$matrices    <- MakeResiduals(ASCA)
    ASCA$residuals$ssq         <- lapply(ASCA$residuals$matrices, function(x) 
                                    sum(x^2)
                                  )
    
    ## Model table ##
    max.length                 <- max(length(factors), length(interactions), 
                                      length(combinations))
    col.factors                <- c(unlist(factors), rep("", 
                                                         max.length - length(factors))
    )
    col.interactions           <- c(unlist(interactions), rep("", 
                                                              max.length - length(interactions))
    )
    col.combinations           <- c(unlist(combinations), rep("", 
                                                              max.length - length(combinations))
    )
    model.text                 <- cbind.data.frame(col.factors, col.interactions, 
                                                   col.combinations)
    colnames(model.text)       <- c("Factors", "Interactions", "Combinations")
    
    
    ## Variance table
    ssq.names                  <- c("Data", "Overall Means", 
                                    unlist(factors), 
                                    unlist(interactions), 
                                    unlist(combinations), 
                                    "Residuals")
    ssq.values                 <- c(unlist(ASCA$data$ssq), 
                                    unlist(ASCA$means$ssq), 
                                    unlist(ASCA$factors$ssq), 
                                    unlist(ASCA$interactions$ssq), 
                                    unlist(ASCA$combinations$ssq),
                                    unlist(ASCA$residuals$ssq))
    ssq.percentages            <- (ssq.values / unlist(ASCA$data$ssq))*100
    ssq.table                  <- cbind.data.frame(ssq.names, 
                                                   ssq.values, 
                                                   ssq.percentages)
    colnames(ssq.table)        <- c("Source", "Sum of Squares", 
                                    "Percentage of variation")
    DAT.session$models$ASCA    <- ASCA
    return(list(model.text, ssq.table))
    
  })

  SetPCs <- function(n.levels) {
    # Set the values for the principal components, PC1 an PC2 that 
    # are plotted against each other
    # Args:
    #  n.levels (numeric): maximum number of levels
    # Returns: 
    #   list:
    #     choices.pc1 (character): choices for first principal component
    #     selected.pc1 (character): selected first principal component
    #     choices.pc2 (character): choices for second principal component
    #     selected.pc1 (character): selected second principal component
    
    if(n.levels == 2) {                           # two levels, single PC
      
      choices.pc1  <- "1"
      selected.pc1 <- "1"
      choices.pc2  <- "None"
      selected.pc2 <- "None"
    } else if (n.levels > 2) {                    # more than 2 levels
      choices.pc1  <- sapply(seq(n.levels - 1), as.character)
      selected.pc1 <- choices.pc1[1]
      choices.pc2  <- sapply(seq(n.levels - 1), as.character)
      selected.pc2 <- choices.pc2[2]
    }
    return(list(choices.pc1, selected.pc1, choices.pc2, selected.pc2))
  }
  
  
  SelectPlotType <- function(pc1, pc2, wtp, plot.factor, plot.type, 
                             plot.projections = "yes", x.values = 1,
                             group.by = 1, biplot.type = 1) {
    # Makes scores, loadings and projections matrices and
    # percentage explained for pc1 and pc2 
    # Args:
    #   pc1 (numeric): first principal component to plot
    #   pc2 (numeric): second principal component to plot (can be 0)
    #   wtp (character): what to plot, factor, interaction, combination or
    #                    residuals
    #   plot.factor: (character): which factor/interaction/combination 
    #                   to plot
    # Returns:
    #   All matrices returned are one (level plots) or two columns
    #   list:
    #     scores (dataframe):  scores
    #     loadings (dataframe): loadings
    #     projections (dataframe): projections
    #     perc.explained (numeric): percentage explained variance by pc's for factor
    
    # scores from svd-function are normalized
    index             <- which(DAT.session$models$ASCA[[wtp]]$terms == 
                                 plot.factor)
    svd.list          <- svd(DAT.session$models$ASCA[[wtp]]$matrices[[index]])
    u                 <- svd.list$u
    d                 <- svd.list$d
    v                 <- svd.list$v
    scores.matrix     <- u %*% diag(d)
    residuals         <- as.matrix(
                          DAT.session$models$ASCA$residuals$matrices[[1]])
    if(pc2 != 0) {
      scores      <- scores.matrix[ , c(pc1, pc2)]
      loadings    <- v[ , c(pc1, pc2)]
      projections <- residuals %*% loadings + scores
    } else if(pc2 == 0) {
      # add column of zeros to plot y = 0
      scores      <- cbind(scores.matrix[, pc1], 
                           rep(0, nrow(scores.matrix)))
      loadings    <- cbind(v[, pc1], rep(0, nrow(v)))
      projections <- residuals %*% loadings + scores
    }
    
    perc.explained <- 100*d^2 / sum((d)^2 )
    pca.list          <- list(as.data.frame(scores),
                              as.data.frame(loadings),
                              as.data.frame(diag(d)),
                              as.data.frame(projections), 
                              perc.explained)
    
    if(plot.type == "scores") {
      p <- PlotScores(pc1, pc2, wtp, plot.factor, plot.projections, pca.list)
    } else if (plot.type == "loadings") {
        p <- PlotLoadings(pc1, pc2, wtp, plot.factor, pca.list)
    } else if (plot.type == "levels") {
        p <- PlotLevels(pc1, wtp, pca.list, x.values, group.by)
    } else if (plot.type == "biplot") {
        p <- PlotBiplot(pca.list, wtp, plot.factor, plot.projections, 
                        pc1 = pc1, pc2 = pc2)
    }
    
    return(p)
  }
  
  Order.levels <- function(design.factor) {
    # Order the design levels as the appear in the design file.
    # Args:
    #  design.factor: ("character"): factor in design matrix
    # Returns:
    #  ordered.levels: (factor): levels ordered as in design file
    
    unique.levels  <- unique(DAT.session$design[, as.numeric(design.factor)])
    ordered.levels <- factor(DAT.session$design[, as.numeric(design.factor)], 
                             levels = unique.levels)
    return(ordered.levels)
  }
  
  
  PlotScores <- function(pc1, pc2, wtp, plot.factor, plot.projections, 
                         pca.list) {
    # Makes the score plots
    # Args: 
    #   pc1 (numeric): first principal component to plot
    #   pc2 (numeric): second principal component to plot (can be 0)
    #   wtp (character): what to plot, factor, interaction, combination or
    #                    residuals
    #   plot.factor (character):  plot factor: which factor, interaction or
    #                             combination to plot
    #   plot.projections (character): if 'yes' projections are plotted
    #   pca.list: (list): scores, loadings, projections and precentage explained
    # Returns:
    #   g (ggplot object): scores plot
    
    scores            <- pca.list[[1]]
    loadings          <- pca.list[[2]]
    eigen.values      <- pca.list[[3]]
    projections       <- pca.list[[4]]
    perc.explnd       <- pca.list[[5]]
  
    if(wtp == "factors") {
      color.to.plot   <- Order.levels(plot.factor)
      color.header    <- colnames(DAT.session$design)[as.numeric(plot.factor)]
      shape.to.plot   <- as.factor(1)
      shape.header    <- ""
    } else if(wtp == "interactions") {
        # icolor is factor with color levels; ishape factor with shape levels
        req(input$icolor, input$ishape)
        color.to.plot <- Order.levels(input$icolor)
        shape.to.plot <- Order.levels(input$ishape)
        color.header  <- colnames(DAT.session$design)[as.numeric(input$icolor)]
        shape.header  <- colnames(DAT.session$design)[as.numeric(input$ishape)]
    } else if(wtp == "combinations") {
        # ccolor is factor with color levels; cshape factor with shape levels
        req(input$ccolor, input$cshape)
        color.to.plot <- Order.levels(input$ccolor)
        shape.to.plot <- Order.levels(input$cshape)
        color.header  <- colnames(DAT.session$design)[as.numeric(input$ccolor)]
        shape.header  <- colnames(DAT.session$design)[as.numeric(input$cshape)]
    } else if(wtp == "residuals") {
        index         <- which(colnames(DAT.session$design) == input$colorlevels)
        color.to.plot <- Order.levels(index)
        shape.to.plot <- as.factor(1)
        color.header  <- input$colorlevels
        shape.header  <- ""
    }
    # Plot averages
    g <- ggplot(data = scores, aes(x = scores[, 1], y = scores[, 2],
                                   col = color.to.plot,
                                   shape = shape.to.plot ))
    g <- g + geom_point(size = 6)
    
    # Plot projections. Note, the projections dataframe has only two columns
    if(plot.projections == "yes") {
      g <- g + geom_point(data = projections, aes(x = projections[, 1], 
                            y  = projections[, 2],
                            col = color.to.plot, 
                            shape = shape.to.plot), size = 3)
    }

    x.label <- paste("PC", as.character(pc1),
                       "(", sprintf("%.2f", perc.explnd[pc1]),"%", ")", 
                       sep = " ")
    if(pc2 == 0) {
      y.label <- ""
    } else if (pc2 != 0) {
      y.label <- paste("PC", as.character(pc2),  
                       "(", sprintf("%.2f", perc.explnd[pc2]),"%", ")",  
                       sep = " " )
    }
    

    g <- g + labs(title = "Scores", 
                  x = x.label,
                  y = y.label)
    if(wtp == "factors" | wtp == "residuals") {
      g <- g + guides(shape = FALSE)
    }
     
    g <- g + scale_color_discrete(name = color.header)
    g <- g + scale_shape_discrete(name = shape.header)
    g <- g + theme(legend.text = element_text(size = 14), 
                   legend.title = element_text(size = 16),
                   axis.title = element_text(size = 14),
                   plot.title = element_text(size = 16, 
                   face="bold", hjust = 0.5)
                  )
    return(g)
  }
  
  PlotLoadings <- function(pc1, pc2, wtp, plot.factor, pca.list) {
    # Makes the loadings plot
    
    loadings                  <- pca.list[[2]]
    # Make dataframe with first column the variables,

    if(pc2 == 0) {
      plot.loadings             <- cbind.data.frame(colnames(DAT.session$data), 
                                                    loadings)
      colnames(plot.loadings) <- c("Variable", "PC")
      plot.title              <- paste("Loadings PC", as.character(pc1))
      p <- ggplot(data = plot.loadings, 
                  aes(x = factor(Variable), y = PC) ) +
        geom_bar(stat = "identity", fill = "#F8766D")   +
        labs(title = plot.title, x = "Variable", y = "Loadings")
    } else {
      req(pc1 != pc2)
      plot.loadings             <- cbind.data.frame(colnames(DAT.session$data), 
                                                    loadings)
      colnames(plot.loadings) <- c("Variable", as.character(pc1), as.character(pc2))
      plot.loadings           <- gather(plot.loadings, Loading, value, 
                                        c(as.character(pc1), as.character(pc2)))
      plot.title              <- paste("Loadings PC", as.character(pc1), "and PC",
                                       as.character(pc2))
      p <- ggplot(data = plot.loadings, 
                  aes(x = factor(Variable), y = value, fill = Loading)) 
      p <- p + geom_bar(stat = "identity", position = position_dodge(preserve = "total")) 
      p <- p + labs(title = plot.title, x = "Variable", y = "Loadings")
    }
    
    p <- p + theme(legend.text = element_text(size = 14), 
                   legend.title = element_text(size = 16),
                   axis.title = element_text(size = 14),
                   axis.text.x = element_text(size = 14, angle = 90),
                   axis.text.y = element_text(size = 14),
                   plot.title = element_text(size = 16, 
                                             face="bold", hjust = 0.5)
    )
    return(p)
  }
  
  PlotLevels <- function(pc, wtp, pca.list, x.values, group.factor) {
    # Plots PC-value against plot.levels and group by group.levels
    scores       <- pca.list[[1]]
    perc.explnd  <- pca.list[[5]]
    x.values     <- Order.levels(x.values)
    group.by     <- Order.levels(group.factor)
    group.header <- colnames(DAT.session$design)[as.numeric(group.factor)]
    if(wtp == "factors") {

      g          <- ggplot(data = scores, aes(x = x.values, y = scores[, 1],
                                     col = "#F8766D" , group = 1))
      g          <- g + geom_line(size = 1)
      g          <- g + geom_point(size = 4)
      g          <- g + theme(legend.position = "none")
    }
    else {

      g          <- ggplot(data = scores, aes(x = x.values, y = scores[, 1],
                                     col = group.by, group = group.by))
      g          <- g + geom_line(size = 1)
      g          <- g + geom_point(size = 4)
    }
    g            <- g + scale_color_discrete(name = group.header)
    g            <- g + labs(title = "PC vs levels", 
                             x = "Levels", 
                             y = paste("PC", as.character(pc),  
                                       "(", sprintf("%.2f", perc.explnd[pc]),"%", ")",  
                                       sep = " " )
                             )
    
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
  
  
  PlotBiplot <- function(pca.list, wtp, plot.factor, plot.projections, 
                         pc1 = 1, pc2 = 2) {
  
    scores             <- pca.list[[1]]
    loadings           <- pca.list[[2]]
    perc.explnd        <- pca.list[[5]]
    colnames(loadings) <- c("L1", "L2")
    max.scores         <- max(sapply(seq(nrow(scores)), function(x) sum(scores[x, ]^2)))
    max.loadings       <- max(sapply(seq(nrow(loadings)), function(x) sum(loadings[x, ]^2)))
    # Scale length of loadings
    lambda             <- (sqrt(max.loadings)/sqrt(max.scores))
    scaled.loadings    <- loadings/lambda
    residuals          <- as.matrix(
                            DAT.session$models$ASCA$residuals$matrices[[1]])

    g <- PlotScores(pc1, pc2, wtp, plot.factor, plot.projections, pca.list)
   
    g <- g + labs(title = "Biplot",
                  x = paste("PC", pc1, 
                            "(", sprintf("%.2f", perc.explnd[pc1]),"%", ")",   
                            sep = " "),
                  y = paste("PC", pc2,  
                            "(", sprintf("%.2f", perc.explnd[pc2]),"%", ")",  
                            sep = " " )
    )
    
    g <- g + geom_segment(data = scaled.loadings,
                          aes(x = 0, y = 0, xend = L1, yend = L2),
                          arrow = arrow() , color = "red", size = 1,
                          inherit.aes = FALSE
    )
    
    g <- g + geom_text_repel(data = scaled.loadings, 
                            aes(x = L1, y = L2, 
                                 label = colnames(DAT.session$data)), 
                             col = "black",
                            inherit.aes = FALSE,
                             show.legend = FALSE)
    
    return(g)
  }
  
  
  ###########################################
  
  
  MakePlotTypeMenu <- observe({
    # Makes the menu items for the objects (factor, interaction, combination, residual) 
    # to be plotted
    
    req(DAT.session$data, DAT.session$models$ASCA)
    factors          <- unlist(DAT.session$models$ASCA$factors$terms)
    interactions     <- unlist(DAT.session$models$ASCA$interactions$terms)
    combinations     <- unlist(DAT.session$models$ASCA$combinations$terms)
    if(length(factors) > 0) {
      updateSelectInput(session, "choose.factor", 
                        choices = factors, selected = factors[1])
    } else if(length(factors) == 0) {
        updateSelectInput(session, "choose.factor", 
                          choices = "", selected = "")
      
    }
    if(length(interactions) > 0) {
      updateSelectInput(session, "choose.interaction", 
                        choices = interactions, selected = interactions[1])
    } else if(length(interactions) == 0) {
        updateSelectInput(session, "choose.interaction", 
                          choices = "", selected = "")
    }
    if(length(combinations) > 0) {
      updateSelectInput(session, "choose.combination", 
                        choices = combinations, selected = combinations[1])
    } else if(length(combinations) == 0) {
        updateSelectInput(session, "choose.combination", 
                          choices = "", selected = "")
    }
  })
  
  UpdateFactorsMenu <- observe ({
    # Update factor menu with pc's to choose from, depending on the 
    # factor selected in the plot menu in MakePlotMenu
    
    req(DAT.session$data, 
        DAT.session$models$ASCA,
        input$choose.factor,
        length(DAT.session$models$ASCA$factors$terms) > 0)

    plot.levels  <- DAT.session$design[, as.numeric(input$choose.factor)]
    # The number of levels depends on the rank of the data matrix and the number 
    # of experimental levels
    n.levels     <- min(max(as.numeric(plot.levels)), 
                          ncol(DAT.session$data), 
                          nrow(DAT.session$data)
                        ) 
    pcs          <- SetPCs(n.levels)
    updateSelectInput(session, "factor.pc1", 
                      choices = pcs[[1]], 
                      selected = pcs[[2]])
    
#    if(input$fpt == "biplot") {
#      output$fbiplot <- renderUI({
#        radioButtons(session$ns("fbiplot"), 
#                     "Biplot scaling",
#                     c("Plot U vs VS" = "1",
#                       "plot US vs V" = "2"),
#                     selected = "1")
#      })
#    }
    # for "levels", plot the scores of a single pc against the factor levels
    if(input$fpt == "levels") {   
      updateSelectInput(session, "factor.pc2", 
                        choices = "None", 
                        selected = "None")
      updateSelectInput(session, "flevelfactor", 
                        choices = input$choose.factor,
                        selected = input$choose.factor)
    } else { 
      updateSelectInput(session, "factor.pc2", 
                        choices = pcs[[3]], 
                        selected = pcs[[4]])
    }
  })
  
  PlotFactor <- reactive ({
    # Determine which PC's can be plotted for factor and call 
    # "SelectPlotType"
    
    req(input$choose.factor,
        DAT.session$data, 
        DAT.session$models$ASCA,
        length(DAT.session$models$ASCA$factors$terms) > 0,
        input$factor.pc1, 
        input$factor.pc2)
    
    if(input$fplotprojctns) {
      plotprojections <- "yes"
    } else {
      plotprojections <- "no"
    }
    pc1   <- as.numeric(input$factor.pc1)
    pc2   <- input$factor.pc2
    if(pc2 == "None") {
      pc2 <- 0
    } else {
      pc2 <- as.numeric(input$factor.pc2)
    }
    if(input$fpt == "levels") {
      req(input$choose.factor)
      p   <- SelectPlotType(pc1, pc2, "factors", 
                            input$choose.factor, 
                            input$fpt,
                            x.values = input$choose.factor,
                            group.by = input$choose.factor)
   
    } else if(input$fpt == "biplot") {
      p   <- SelectPlotType(pc1, pc2, "factors", input$choose.factor, 
                            input$fpt, plot.projections = plotprojections)
    }  else {
      p   <- SelectPlotType(pc1, pc2, "factors", input$choose.factor, 
                            input$fpt, plot.projections = plotprojections)
    } 
    return(p)
  })
  
  UpdateInteractionsMenu <- observe ({
    # Update interacion menu with pc's to choose from, depending on the 
    # interaction selected in the plot menu in MakePlotMenu
    req(DAT.session$data, 
        DAT.session$models$ASCA,
        input$choose.interaction,
        length(DAT.session$models$ASCA$interactions$terms) > 0 )
    index          <- which(
                        DAT.session$models$ASCA$interactions$terms == 
                        input$choose.interaction
                      )
    req(length(index) > 0)
    plot.levels    <- DAT.session$models$ASCA$interactions$levels[[index]]
    n.levels       <- min(max(as.numeric(plot.levels)), 
                          ncol(DAT.session$data), 
                          nrow(DAT.session$data)) 
    pcs            <- SetPCs(n.levels)
    updateSelectInput(session, "interaction.pc1", choices = pcs[[1]], 
                      selected = pcs[[2]])
    
#    if(input$ipt == "biplot") {
#      output$ibiplot <- renderUI({
#        radioButtons(session$ns("ibiplot"), 
#                     "Biplot scaling",
#                     c("Plot U vs VS" = "1",
#                       "plot US vs V" = "2"),
#                     selected = "1")
#      })
#    }
    if(input$ipt == "levels") {   
      updateSelectInput(session, "interaction.pc2", choices = "None", 
                        selected = "None")
      # Get factors involved
      choice       <- unlist(str_extract_all(
                        DAT.session$models$ASCA$interactions$terms[[index]],
                        "[0-9]+")
                      )
      # Factor on x-axis
      updateSelectInput(session, "ixvalues", choices = choice, 
                        selected = choice[1])
      # Factor to group
      updateSelectInput(session, "igroupfactor", choices = choice, 
                        selected = choice[2])
    } else if(input$ipt =="scores") {
      # Get factors involved
        choice       <- unlist(str_extract_all(
                          DAT.session$models$ASCA$interactions$terms[[index]],  
                          "[0-9]+")
                        )
        updateSelectInput(session, "icolor", choices = choice, 
                        selected = choice[1])
        updateSelectInput(session, "ishape", choices = choice, 
                        selected = choice[2])
        updateSelectInput(session, "interaction.pc2", choices = pcs[[3]], 
                          selected = pcs[[4]])
    } else {
      updateSelectInput(session, "interaction.pc2", choices = pcs[[3]], 
                        selected = pcs[[4]])
    }
  })
  
  PlotInteractions <- reactive ({
    req(DAT.session$data, 
        DAT.session$models$ASCA,
        input$choose.interaction,
        length(DAT.session$models$ASCA$interactions$terms) > 0,
        length(DAT.session$models$ASCA$interactions$levels) > 0)
    index          <- which(
                        DAT.session$models$ASCA$interactions$terms == 
                        input$choose.interaction
                      )
    plot.levels    <- DAT.session$models$ASCA$interactions$levels[[index]]
    req(input$interaction.pc1, 
        input$interaction.pc2)
    
    if(input$iplotprojctns) {
      plotprojections <- "yes"
    } else {
        plotprojections <- "no"
    }
    pc1            <- as.numeric(input$interaction.pc1)
    pc2            <- input$interaction.pc2
    if(pc2 == "None") {
      pc2 <- 0
    } else {
      pc2 <- as.numeric(input$interaction.pc2)
    }
    if(input$ipt == "levels") {
      req(input$ixvalues, input$igroupfactor)
      p            <- SelectPlotType(pc1, pc2, "interactions", 
                                     input$choose.interaction, 
                                     input$ipt, 
                                     x.values = input$ixvalues, 
                                     group.by = input$igroupfactor)
    } else if (input$ipt == "biplot") {
      p            <- SelectPlotType(pc1, pc2, "interactions", 
                                     input$choose.interaction, 
                                     input$ipt, 
                                     plot.projections = plotprojections)
    } else {
      p            <- SelectPlotType(pc1, pc2, "interactions", 
                                     input$choose.interaction, 
                                     input$ipt, 
                                     plot.projections = plotprojections)
    }
    return(p)
  })
  
  UpdateCombinationsMenu <- observe ({
    req(DAT.session$data, 
        DAT.session$models$ASCA,
        input$choose.combination,
        length(DAT.session$models$ASCA$combinations$terms) > 0)
    index          <- which(
                            DAT.session$models$ASCA$combinations$terms == 
                            input$choose.combination
                           )
    req(length(index) > 0)
    plot.levels    <- DAT.session$models$ASCA$combinations$levels[[index]]
    n.levels       <- min(max(as.numeric(plot.levels)), 
                          ncol(DAT.session$data), 
                          nrow(DAT.session$data)) 
    pcs            <- SetPCs(n.levels)
    updateSelectInput(session, "combination.pc1", choices = pcs[[1]], 
                      selected = pcs[[2]])
    
 #   if(input$cpt == "biplot") {
 #     output$cbiplot <- renderUI({
 #       radioButtons(session$ns("cbiplot"), 
#                     "Biplot scaling",
#                     c("Plot U vs VS" = "1",
#                       "plot US vs V" = "2"),
#                     selected = "1")
#      })
#    }
    
    if(input$cpt == "levels") {  # for "levels", plot a single pc against the levels
      updateSelectInput(session, "combination.pc2", choices = "None", 
                        selected = "None")
      choice       <- unique(unlist(str_extract_all(
                        DAT.session$models$ASCA$combinations$terms[[index]],
                        "[0-9]+"))
                      )
      updateSelectInput(session, "cxvalues", choices = choice, 
                        selected = choice[1])
      updateSelectInput(session, "cgroupfactor", choices = choice, 
                        selected = choice[2])
    } 
    else if(input$cpt == "scores") {
      # Get factors involved
      choice       <- unique(unlist(str_extract_all(
                        DAT.session$models$ASCA$combinations$terms[[index]],  
                        "[0-9]+"))
                      )
      updateSelectInput(session, "ccolor", choices = choice, 
                        selected = choice[1])
      updateSelectInput(session, "cshape", choices = choice, 
                        selected = choice[2])
      updateSelectInput(session, "combination.pc2", choices = pcs[[3]], 
                        selected = pcs[[4]])
    } else {
      updateSelectInput(session, "combination.pc2", choices = pcs[[3]], 
                        selected = pcs[[4]])
    }
  })
  
  PlotCombination <- reactive ({
    req(DAT.session$data, 
        DAT.session$models$ASCA,
        input$choose.combination,
        length(DAT.session$models$ASCA$combinations$terms) > 0)
    index          <- which(
      DAT.session$models$ASCA$combinations$terms == 
        input$choose.combination
    )
    req(length(index) > 0)
    plot.levels    <- DAT.session$models$ASCA$combinations$levels[[index]]
    req(input$combination.pc1, input$combination.pc2)
    
    if(input$cplotprojctns) {
      plotprojections <- "yes"
    } else {
      plotprojections <- "no"
    }
    
    pc1            <- as.numeric(input$combination.pc1)
    pc2            <- input$combination.pc2
    if(pc2 == "None") {
      pc2 <- 0
    } else {
      pc2 <- as.numeric(input$combination.pc2)
    }
    if(input$cpt == "levels") {
      req(input$cxvalues, input$cgroupfactor)
      
      # Order levels according to appearence in the design file
      design.levels  <- DAT.session$design[, as.numeric(input$cxvalues)]
      unique.levels  <- unique(design.levels)
      plot.levels    <- factor(design.levels, levels = unique.levels)
      
      design.levels  <- DAT.session$design[, as.numeric(input$cgroupfactor)]
      unique.levels  <- unique(design.levels)
      group.by       <- factor(design.levels, levels = unique.levels)
      p           <- SelectPlotType(pc1, pc2, "combinations", 
                                    input$choose.combination, 
                                    input$cpt,
                                    x.values = input$cxvalues,
                                    group.by = input$cgroupfactor)
    } else if(input$cpt == "biplot") {
      p           <- SelectPlotType(pc1, pc2, "combinations", 
                                    input$choose.combination,
                                    input$cpt,
                                    plot.projections = plotprojections)
    } else {
      p           <- SelectPlotType(pc1, pc2, "combinations", 
                                    input$choose.combination,
                                    input$cpt,
                                    plot.projections = plotprojections)
    }
    return(p)
  })
  
  UpdateResiduals <- observe ({
    force(MakeModel())   # update ASCA
    req(DAT.session$data, 
        DAT.session$models$ASCA, 
        input$rpt,
        DAT.session$n.factors)
    n.levels       <- min(ncol(DAT.session$data), nrow(DAT.session$data)) 
    pcs            <- SetPCs(n.levels)
    n.factors      <- DAT.session$n.factors
  
    updateSelectInput(session, "colorlevels", 
                      "Color by:",
                      choices  = colnames(DAT.session$design[, seq(n.factors)]),
                      selected = colnames(DAT.session$design[, seq(n.factors)])[1])
    updateSelectInput(session, 
                      "residual.pc1", 
                      choices = pcs[[1]], 
                      selected = pcs[[2]])
    updateSelectInput(session, "residual.pc2", choices = pcs[[3]], 
                        selected = pcs[[4]])
  })
  
  PlotResiduals <- reactive ({
    req(DAT.session$data, 
        DAT.session$design,
        DAT.session$models$ASCA,
        input$colorlevels,
        length(DAT.session$models$ASCA$residuals$matrices) > 0)
    n.levels       <- min(ncol(DAT.session$data), nrow(DAT.session$data))
    pcs            <- SetPCs(n.levels)
    req(input$residual.pc1, input$residual.pc2)
    pc1            <- as.numeric(input$residual.pc1)
    pc2            <- input$residual.pc2
    if(pc2 == "None") {
      pc2 <- 0
    } else {
      pc2 <- as.numeric(input$residual.pc2)
    }
      p            <- SelectPlotType(pc1, pc2, "residuals", "1", 
                                     input$rpt, 
                                     plot.projections = "no")
    return(p)
  })
  
  output$factors       <- renderPlot(print(PlotFactor()))
  output$interactions  <- renderPlot(print(PlotInteractions()))
  output$combination   <- renderPlot(print(PlotCombination()))
  output$residuals     <- renderPlot(print(PlotResiduals()))
  output$model         <- renderTable(MakeModel()[[1]])
  output$variances     <- renderTable(MakeModel()[[2]])
  
}