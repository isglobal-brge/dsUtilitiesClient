#' Create Analysis of Variance Table
#'
#' Fit an analysis of variance model by calling \code{ds.glm} on data from single or multiple sources on the server-side.
#'
#' @usage
#' ds.aov(
#'   formula = NULL,
#'   data = NULL,
#'   type = "split",
#'   ANOVAtype = "I",
#'   datasources = NULL
#' )
#'
#' @param formula A valid formula specifying the model, this needs to be compatible with \code{ds.glm}.
#' @param data A character string specifying the name of the data frame that contains all of the variables in the formula.
#' @param type A character string that represents the type of analysis to carry out. This can be set as 'combine', 'combined', 'combines', 'split', 'splits', 's', 'both' or 'b'.
#' @param ANOVAtype A character string specifying the ANOVA type (I, II & III) if two-way.
#' @param datasources A list of DSConnection-class objects obtained after login. If the datasources argument is not specified the default set of connections will be used: see datashield.connections_default.
#'
#' @return ANOVA table for each connection if type is split, or combined ANOVA table if type is combined.
#' @author MV
#'
#' @examples TODO
ds.aov <- function(formula=NULL, data=NULL, type = "split", ANOVAtype = "I", datasources = NULL){
  # Check input
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }

  # Ensure datasources is a list of DSConnection-class
  if(!(is.list(datasources) && all(unlist(lapply(datasources, function(d) {methods::is(d,"DSConnection")}))))){
    stop("The 'datasources' were expected to be a list of DSConnection-class objects", call.=FALSE)
  }

  # Verify that 'formula' is set
  if(is.null(formula)){
    stop(" Please provide a valid formula!", call.=FALSE)
  }
  formula <- stats::as.formula(formula)

  # Ff the argument 'data' is set, check that the data frame is defined (i.e. exists) on the server site
  if(is.null(data)){
    stop("data=NULL. Please provide the name of a matrix or dataframe!", call.=FALSE)
  }else{
    defined <- dsBaseClient:::isDefined(datasources, data)  # What  do you do with the result of this?
  }

  # To hold results
  anova.tables <- list()

  if (type == "both" || type == "split" || type == "splits" || type == "s"){
    # Make a call for each server
    for (i in 1:length(datasources)){
      anova.table.server <- ds.makeAovTable(formula = formula, data = data, ANOVAtype = ANOVAtype, datasources = datasources[i])
      anova.tables <- c(list(anova.table.server), anova.tables)
    }
  } else if (type == "both" || type == "combine" || type == "combined" || type == "c") {
    # Combined call
    anova.table.server <- ds.makeAovTable(formula = formula, data = data, ANOVAtype = ANOVAtype, datasources = datasources)
    anova.tables <- c(list(anova.table.server), anova.tables)
  }

  return (anova.tables)
}

# This function should be private, do I still document with roxygen?
# It does the ANOVA computations, and creates the table, the connection object is passed onto ds.RSS
# For computing two-way types I, II (and later III) I used the formulas in
# https://mcfromnz.wordpress.com/2011/03/02/anova-type-iiiiii-ss-explained/
ds.makeAovTable <- function(formula=NULL, data=NULL, ANOVAtype = "I", datasources=NULL){
  # Split the formula into variables
  pred <- formula[[2]]
  vars <- labels(terms(formula, keep.order = TRUE))

  # Total degrees of freedom and sum sq. are the same in all cases
  DF.T <- ds.length(x = paste0(data, "$", pred), datasources = datasources)
  DF.T <- DF.T[[length(DF.T)]] - 1
  SS.T <- ds.RSS(formula = paste0(pred, " ~ 1"), data = data, datasources = datasources)

  # One-way ANOVA
  if(length(vars) == 1){
    # Residual sum squares from model
    R.A <- ds.RSS(formula = formula, data = data, datasources = datasources)

    # Sum squares
    SS.A <- SS.T - R.A
    SS.E <- SS.T - SS.A

    # Degrees of freedom
    DF.A <- length(ds.levels(x = paste0(data, "$", vars[1]), datasources = datasources)[[1]]$Levels) - 1
    DF.E <- DF.T - DF.A

    # Mean sum squares
    MS.A <- SS.A / DF.A
    MS.E <- SS.E / DF.E

    # F values
    F.A <- MS.A / MS.E

    # P values
    P.A <- pf(F.A, DF.A, DF.E, lower.tail = FALSE)

    # Making the ANOVA table
    anova.table.server <- data.frame(Df = c(DF.A, DF.E, DF.T))
    anova.table.server$'Sum Sq' <- c(SS.A, SS.E, SS.T)
    anova.table.server$'Mean Sq' <- c(MS.A, MS.E, NA)
    anova.table.server$'F value' <- c(F.A, NA, NA)
    anova.table.server$'Pr(>F)' <- c(P.A, NA, NA)
    rownames(anova.table.server) <- c(vars[1], "Residuals", "Total")

  } else {
    # Two-way ANOVAs

    if(ANOVAtype == "I" || ANOVAtype == "1" || ANOVAtype == 1){
      # Residual sum squares from model
      R.A <- ds.RSS(formula = paste0(pred, " ~ ", vars[1]), data = data, datasources = datasources)
      R.AB <- ds.RSS(formula = paste0(pred, " ~ ", vars[1], " + ", vars[2]), data = data, datasources = datasources)
      R.AxB <- ds.RSS(formula = paste0(pred, " ~ ", vars[1], " * ", vars[2]), data = data, datasources = datasources)

      # Sum squares
      SS.A <- SS.T - R.A
      SS.B <- R.A - R.AB
      SS.AB <- R.AB - R.AxB
      SS.E <- SS.T - SS.A - SS.B - SS.AB
    }
    else if (ANOVAtype == "II" || ANOVAtype == "2" || ANOVAtype == 2){
      # Residual sum squares from model
      R.A <- ds.RSS(formula = paste0(pred, " ~ ", vars[1]), data = data, datasources = datasources)
      R.B <- ds.RSS(formula = paste0(pred, " ~ ", vars[2]), data = data, datasources = datasources)
      R.AB <- ds.RSS(formula = paste0(pred, " ~ ", vars[1], " + ", vars[2]), data = data, datasources = datasources)
      R.AxB <- ds.RSS(formula = paste0(pred, " ~ ", vars[1], " * ", vars[2]), data = data, datasources = datasources)

      # Sum squares
      SS.A <- R.B - R.AB
      SS.B <- R.A - R.AB
      SS.AB <- R.AB - R.AxB
      SS.E <- R.AB - SS.AB
    } else if (ANOVAtype == "III" || ANOVAtype == "3" || ANOVAtype == 3) {
      # not implemented
    } else {
      # Do something
    }
    # The rest is the same for all ANOVA types

    # Degrees of freedom
    DF.A <- length(ds.levels(x = paste0(data, "$", vars[1]), datasources = datasources)[[1]]$Levels) - 1
    DF.B <- length(ds.levels(x = paste0(data, "$", vars[2]), datasources = datasources)[[1]]$Levels) - 1
    DF.AB <- DF.A * DF.B
    DF.E <- (DF.T + 1)  - (DF.A + 1)*(DF.B + 1)

    # Mean sum squares
    MS.A <- SS.A / DF.A
    MS.B <- SS.B / DF.B
    MS.AB <- SS.AB / DF.AB
    MS.E <- SS.E / DF.E

    # F values
    F.A <- MS.A / MS.E
    F.B <- MS.B / MS.E
    F.AB <- MS.AB / MS.E

    # P values
    P.A <- pf(F.A, DF.A, DF.E, lower.tail = FALSE)
    P.B <- pf(F.B, DF.B, DF.E, lower.tail = FALSE)
    P.AB <- pf(F.AB, DF.AB, DF.E, lower.tail = FALSE)

    # Making the ANOVA table
    anova.table.server <- data.frame(Df = c(DF.A, DF.B, DF.AB, DF.E, DF.T))
    anova.table.server$'Sum Sq' <- c(SS.A, SS.B, SS.AB, SS.E, SS.T)
    anova.table.server$'Mean Sq' <- c(MS.A, MS.B, MS.AB, MS.E, NA)
    anova.table.server$'F value' <- c(F.A, F.B, F.AB, NA, NA)
    anova.table.server$'Pr(>F)' <- c(P.A, P.B, P.AB, NA, NA)
    rownames(anova.table.server) <- c(vars[1], vars[2], paste0(vars[1], ":", vars[2]), "Residuals", "Total")
  }

  return (anova.table.server)
}

# This function should be private, do I still document with roxygen?
# It fits a ds.glm object, and uses the coefficients to predict and compute the residual sum squares.
# I reconstruct the predictions from the coefficients, could be better.
# I think this should be a DS function in general, but this version is not good for that.
ds.RSS <- function(formula=NULL, data=NULL, datasources = NULL){ # private?
  # Split the formula into variables
  formula <- stats::as.formula(formula)
  pred <- formula[[2]]
  vars <- labels(terms(formula, keep.order = TRUE))

  # Fit the glm object, and save its intercept on the server
  fit <- ds.glm(formula, data = data, family = "gaussian", datasources = datasources)
  ds.assign(toAssign = paste0(fit$coefficients[,1][1]), newobj = "intercept", datasources = datasources)

  if(length(vars) == 0){
    # Get RSS if the model was empty (pred ~ 1)
    ds.assign(toAssign = paste0("(intercept - ", data, "$", pred, ")^2"), newobj = "res2", datasources = datasources)
  } else if (length(vars) == 1){
    # Get RSS if the model had one variable (pred ~ var1)
    ds.assign(toAssign = paste0("c(0,", paste0(fit$coefficients[,1][-1], collapse = ","), ")"), newobj = "coeffs", datasources = datasources)
    ds.assign(toAssign = paste0("coeffs[as.numeric(", data, "$", vars[1], ")]"), newobj = "terms", datasources = datasources)
    ds.assign(toAssign = paste0("((intercept + terms) - ", data, "$", pred, ")^2"), newobj = "res2", datasources = datasources)
  } else if (length(vars) == 2){
    # Get RSS if the model had two variables (pred ~ var1 + var2)
    l1 <- length(ds.levels(x = paste0(data, "$", vars[1]), datasources = datasources)[[1]]$Levels) - 2
    ds.assign(toAssign = paste0("c(0,", paste0(as.numeric(fit$coefficients[,1][2:(2 + l1)]), collapse = ","), ")"), newobj = "coeffs1", datasources = datasources)
    ds.assign(toAssign = paste0("coeffs1[as.numeric(", data, "$", vars[1], ")]"), newobj = "terms1", datasources = datasources)

    if(!grepl(":", vars[2], fixed = TRUE)){
      l2 <- length(ds.levels(x = paste0(data, "$", vars[2]), datasources = datasources)[[1]]$Levels) - 2
      ds.assign(toAssign = paste0("c(0,", paste0(as.numeric(fit$coefficients[,1][(2 + l1 + 1):(2 + l1 + l2 + 1)]), collapse= ","), ")"), newobj = "coeffs2", datasources = datasources)
      ds.assign(toAssign = paste0("coeffs2[as.numeric(", data, "$", vars[2], ")]"), newobj = "terms2", datasources = datasources)
    } else {
      # This is the case where the second variable is the interaction term
      # ds.glm cannot be called with formulas pred ~ var1 + var1:var2 or pred ~ var1*var2 - var2 however
      # with this case I could make type III.
    }

    ds.assign(toAssign = paste0("((intercept + terms1 + terms2) - ", data, "$", pred, ")^2"), newobj = "res2", datasources = datasources)
  } else if (length(vars) == 3){
    # Get RSS if the model had two variables (pred ~ var1 * var2)
    l1 <- length(ds.levels(x = paste0(data, "$", vars[1]), datasources = datasources)[[1]]$Levels) - 2
    ds.assign(toAssign = paste0("c(0,", paste0(as.numeric(fit$coefficients[,1][2:(2 + l1)]), collapse = ","), ")"), newobj = "coeffs1", datasources = datasources)
    ds.assign(toAssign = paste0("coeffs1[as.numeric(", data, "$", vars[1], ")]"), newobj = "terms1", datasources = datasources)

    l2 <- length(ds.levels(x = paste0(data, "$", vars[2]), datasources = datasources)[[1]]$Levels) - 2
    ds.assign(toAssign = paste0("c(0,", paste0(as.numeric(fit$coefficients[,1][(2 + l1 + 1):(2 + l1 + l2 + 1)]), collapse= ","), ")"), newobj = "coeffs2", datasources = datasources)
    ds.assign(toAssign = paste0("coeffs2[as.numeric(", data, "$", vars[2], ")]"), newobj = "terms2", datasources = datasources)

    ds.assign(toAssign = paste0("c(0,", paste0(as.numeric(fit$coefficients[,1][-c(1:(2 + l1 + l2 + 1))]), collapse= ","), ")"), newobj = "coeffs3", datasources = datasources)
    ds.assign(toAssign = paste0("(as.numeric(", data, "$", vars[1], ") - 1) * (as.numeric(", data, "$", vars[2], ") - 1) + 1" ), newobj = "indices", datasources = datasources)
    ds.assign(toAssign = paste0("coeffs3[indices]"), newobj = "terms3", datasources = datasources)

    ds.assign(toAssign = paste0("((intercept + terms1 + terms2 + terms3) - ", data, "$", pred, ")^2"), newobj = "res2", datasources = datasources)
  }

  # This is to get the sum, there probably is a better way of doing it.
  temp <- ds.mean(x = "res2", datasources = datasources)$Mean.by.Study
  if (nrow(temp) == 1) {return (temp[1] * temp[3])}
  else {return ( sum(temp[,1] * temp[,3]))}
}
