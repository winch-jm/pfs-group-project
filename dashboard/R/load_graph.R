load_graph <- function(path = "data/graph_edges.csv") {
  edges <- readr::read_csv(path)
  graph_from_data_frame(edges, directed = FALSE)
}