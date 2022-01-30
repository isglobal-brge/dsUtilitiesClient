#' Title
#'
#' @param formula
#' @param data
#' @param type
#' @param datasources
#'
#' @return
#' @export
#'
#' @examples
ds.aov <- function(formula=NULL, data=NULL, type = "split", datasources = NULL){
  # Check input
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }

  # ensure datasources is a list of DSConnection-class
  if(!(is.list(datasources) && all(unlist(lapply(datasources, function(d) {methods::is(d,"DSConnection")}))))){
    stop("The 'datasources' were expected to be a list of DSConnection-class objects", call.=FALSE)
  }

  # verify that 'formula' was sets
  if(is.null(formula)){
    stop(" Please provide a valid formula!", call.=FALSE)
  }
  formula <- stats::as.formula(formula)  # Should you check if  the formula is valid for the data?

  # if the argument 'data' is set, check that the data frame is defined (i.e. exists) on the server site
  if(is.null(data)){
    stop("data=NULL. Please provide the name of a matrix or dataframe!", call.=FALSE)  # What is  the  .call for?
  }else{
    defined <- dsBaseClient:::isDefined(datasources, data)  # What  do you do with the result of this?
  }

  var <- formula[[3]]  # Need to change these for two-way
  pred <- formula[[2]]

  # Get means by subsets
  group.means <- list()
  overall.means <- list()
  sum.sq.errors <- list()
  for (i in 1:length(datasources)){
    levels <- ds.levels(x = paste0(data, "$", var), datasources = datasources[i])
    group.means.study <- list()
    sum.sq.errors.study <- 0
    for (lvl in levels[[1]]$Levels) {
      ds.make(paste0(data, "$", var, "==", "'", lvl, "'"), "bools")
      ds.asNumeric("bools", "bools_to_ones")

      ds.dataFrameSubset(
        df.name = data,
        newobj = "subtable",
        V1.name = "bools_to_ones",
        V2.name = "1",
        Boolean.operator = "==",
        keep.NAs = FALSE,
        datasources = datasources[i]
      )
      group.means.study[[lvl]] <- ds.mean(x = paste0("subtable$", pred), datasources = datasources[i])$Mean.by.Study
      overall.means[[i]] <- ds.mean(x = paste0(data, "$", pred), datasources = datasources[i])$Mean.by.Study

      ds.make(paste0("(subtable$", pred, "-", group.means.study[[lvl]][1],")^2"), "squared.diffs", datasources = datasources[i])
      temp <- ds.mean(x = "squared.diffs", datasources = datasources[i])$Mean.by.Study
      sum.sq.errors.study <- sum.sq.errors.study + temp[1] * temp[3]
    }
    group.means[[i]] <- group.means.study
    sum.sq.errors[[i]] <- sum.sq.errors.study
  }


  if (type == "split"){
    metrics.study <- list()
    for (i in 1:length(datasources)){
      group.means.study <- group.means[[i]]
      SSE <- sum.sq.errors[[i]]

      # SSE H0
      ds.make(paste0("(", data, "$", pred, "-", overall.means[[i]][1],")^2"), "squared.diffs", datasources = datasources[i])
      temp <- ds.mean(x = "squared.diffs", datasources = datasources[i])$Mean.by.Study
      SSE.H0 <- temp[1] * temp[3]

      df <- length(group.means.study) - 1
      error.df <- temp[3] - length(group.means.study)
      explained.SSE <- SSE.H0 - sum.sq.errors[[i]]

      explained.MSE <- explained.SSE / df
      MSE <- SSE / error.df
      F.val <- explained.MSE / MSE

      metrics <- data.frame(Df = c(df, error.df))
      metrics$'Sum Sq' <- c(explained.SSE, SSE)
      metrics$'Mean Sq' <- c(explained.MSE, MSE)
      metrics$'F value' <- c(F.val, NA)
      metrics$'Pr(>F)' <- c(pf(F.val, df, error.df, lower.tail = FALSE), NA)
      rownames(metrics) <- c(var, "Residuals")

      metrics.study[[i]] <- metrics
    }
  }


}
