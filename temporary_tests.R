df <- data.frame(date = c("1990", "1980", "1970"),
                 fac = as.factor(c(4,2,10)),
                 var = rnorm(3), stringsAsFactors = F)
df
str(df)

as_tibble(df)


df <- apply(df, MARGIN = 2, as.numeric)
str(df)
df


lapply(df, as.numeric)



is_factor <- vapply(df, is.factor, logical(1))
df[is_factor] <- as.numeric(as.character(df[is_factor]))
df
str(df)
as.character(df$fac)





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
str(df)




if (any(unlist(lapply(df, is.numeric))) == FALSE) {
  warning("All columns must be numeric, will attempt to convert to numeric")
  # df[] <- lapply(df, function(x) {
  #   if (is.factor(x)) {
  #     as.numeric(as.character(x))
  #   } else {
  #     if (is.numeric(x) == FALSE) {
  #       as.numeric(x)
  #     } else {
  #       x
  #     }
  #   }
  # })
  # df
}
df
str(df)









if (is_tibble(df)) {
  warning("Tibble format not suited to fishR, converting to data.frame")
  df <- as.data.frame(df)
}
