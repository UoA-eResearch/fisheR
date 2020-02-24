#' Conversion function
#'
#' Conversion function attempts to force input data structure to correct format (i.e., numeric data). *NOTE, does not attempt to rearrange columns or format dates*. Outputs must still be manually checked.
#'
#' @param df A dataframe containing time information and variables of interest
#' @returns A dataframe with each column as numeric data


conversion_func <-  function(df) {
  cat("Current structure: \n") # Should be able to do this in a single cat call right?!
  cat(str(df), "\n")
  cat("Current class: \n")
  cat(class(df), "\n")

  if (class(df)[1] != "data.frame") {
    warning(
      "fisher() requires a dataframe (tibble format not suitable), attempting to converting to data.frame."
    )
    df <- as.data.frame(df)
    cat("Attempted to convert to data.frame. Check class: \n")
    cat(class(df), "\n")
    df
  }


  if (!all(lapply(df, class) %in% c("numeric", "integer", "factor", "character"))) {
    stop(
      "Unable to convert classes that are not factor or character. All columns must be numeric, check structure for date type objects. First column should be time and all following columns the variables of interest"
    )
  } else {
    if (any(lapply(df, class) %in% c("factor", "character"))) {
      warning("Class factor or character found, attempting to convert to numeric. \n")
    }
  }

  if (any(unlist(lapply(df, is.numeric)) == FALSE)) {
    df[] <- lapply(df, function(x) {
      if (is.factor(x)) {
        as.numeric(as.character(x))
      } else {
        if (is.numeric(x) == FALSE) {
          as.numeric(x)
        } else {
          x
        }
      }
    })
    df
  }
  cat("Succesfully attempted to convert data structure, be careful of result! \n")
  cat("new structure: \n")
  cat(str(df), "\n")
  return(df)
}

# ### examples:
# df <- data.frame(date = c("1990", "1980", "1970"),
#                  fac = as.factor(c(4,2,10)),
#                  var = rnorm(3), stringsAsFactors = F)
#
#
# df1 <- conversion_func(df) # Converts columns to numeric
# str(df1)
#
# date_df <- df
# date_df$date = as.Date(date_df$date, "%Y")
# conversion_func(date_df) # Will error as we have not attempted to convert dates
#
#
# df2 <- conversion_func(as_tibble(df)) # Requires dplyr. Will convert tibble to data.frame
