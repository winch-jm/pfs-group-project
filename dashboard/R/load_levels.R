load_levels <- function(path = "data/leiden_levels.json") {
  jsonlite::fromJSON(path, simplifyDataFrame = TRUE)
}
