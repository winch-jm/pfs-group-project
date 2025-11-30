compute_layout <- function(g) {
  ggraph::create_layout(g, layout = "fr")
}