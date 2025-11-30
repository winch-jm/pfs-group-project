# app.R
library(shiny)
library(readr)
library(dplyr)
library(igraph)
library(ggraph)
library(plotly)

#----------------- load data & build graph (once) -----------------#

edges <- read_csv("../data/graph_edges_mcf7.csv",
                  col_types = cols(
                    src    = col_character(),
                    dist    = col_character(),
                    weight = col_double()
                  ))

# igraph uses 'name' as vertex ID by default; we'll treat src/dst as labels
g <- graph_from_data_frame(edges, directed = FALSE)
E(g)$weight <- edges$weight

metadata <- read.table("../data/siginfo_beta.txt",sep="\t",header=TRUE)
cp_info <- read.table("../data/compoundinfo_beta.txt",sep="\t",header=TRUE,fill = TRUE,
                      quote = "",
                      comment.char = "")
atc_info <- read_csv("../data/drug2atc.csv")

umap_coords <- read_csv("../data/MCF7_umap.csv")
# print(names(atc_info))

#print(names(cp_info))
metadata <- metadata %>%
  inner_join(cp_info, by = c("cmap_name","pert_id"))

metadata <- metadata %>%
  inner_join(umap_coords, by="sig_id")

metadata <- metadata %>%
  inner_join(atc_info, by="inchi_key")

partition <- readr::read_csv("../data/partition_mcf7.csv",
                             col_types = cols(
                               node      = col_character(),
                               node_id   = col_character(),
                               community = col_integer()
                             ))


# make community a factor for nicer colors
partition <- partition %>%
  dplyr::mutate(
    node      = as.character(node),
    node_id   = as.character(node_id),
    community = factor(community)
  )

# Merge partition + metadata
node_info <- partition %>%
  inner_join(metadata, by = c("node_id" = "sig_id"))


# match by vertex name (we used src/dst as character, so V(g)$name is "0","1",...)
V(g)$community  <- node_info$community[match(V(g)$name, node_info$node)]
V(g)$cell_iname <- node_info$cell_iname[match(V(g)$name, node_info$node)]
V(g)$pert_id    <- node_info$pert_id[match(V(g)$name, node_info$node)]
V(g)$tas        <- node_info$tas[match(V(g)$name, node_info$node)]
V(g)$atc_class  <- node_info$level1[match(V(g)$name, node_info$node)]
V(g)$atc_desc   <- node_info$level1_description[match(V(g)$name, node_info$node)]
V(g)$moa        <- node_info$moa[match(V(g)$name, node_info$node)]
V(g)$drug       <- node_info$cmap_name[match(V(g)$name, node_info$node)]

print(table(V(g)$atc_class))   # how many nodes per community?
print(vcount(g))            # total node count

# compute a layout once so it doesn't change between redraws
# (you can swap "fr" for "kk", "lgl", etc.)
graph_layout <- create_layout(g, layout = "fr")

# ---- graph layout as plain df for plotly ----
layout_df <- graph_layout %>%
  as.data.frame() %>%
  dplyr::rename(graph_x = x, graph_y = y) %>%
  dplyr::select(name, graph_x, graph_y)

node_info <- node_info %>%
  left_join(layout_df, by = c("node" = "name"))

print(nrow(graph_layout))
print(sum(is.na(graph_layout$x) | is.na(graph_layout$y))  )


# -------------------------------------------------------------------------


#----------------- Shiny app -----------------#

# ui <- fluidPage(
#   titlePanel("Graph from CSR (minimal view)"),
#   
#   sidebarLayout(
#     sidebarPanel(
#       selectInput(
#         "color_by",
#         "Color nodes by:",
#         choices = c(
#           "community",
#           "cell_iname",
#           "tas",
#           "atc_class",
#           "atc_desc",
#           "moa"
#         ),
#         selected = "community"
#       ),
#       checkboxInput("drop_no_atc", "Drop nodes without ATC label", value = FALSE),
#       sliderInput("edge_alpha", "Edge transparency:",
#                   min = 0.01, max = 0.8, value = 0.02, step = 0.01),
#       sliderInput("node_size", "Node size:",
#                   min = 0.1, max = 8, value = 1,step=0.025),
#       checkboxInput("show_labels", "Show node labels", value = FALSE),
#       checkboxInput("show_edges", "Show edges", value = FALSE),
#       selectInput(
#         "atc_level",
#         "ATC level:",
#         choices = c("Level 1" = "level1",
#                     "Level 2" = "level2",
#                     "Level 3" = "level3",
#                     "Level 4" = "level4"),
#         selected = "level1"
#       )
#     ),
#     mainPanel(
#       plotlyOutput("graph_plot", height = "800px"),   # <- plotly instead of visNetwork
#       br(),
#       h3("Community sizes"),
#       plotOutput("comm_bar", height = "300px"),
#       br(),
#       h3("TAS score distribution"),
#       plotOutput("tas_hist", height = "300px"),
#       br(),
#       h3("ATC Enrichment"),
#       plotOutput("atc_enrichment", height = "600px"),
#       br(),
#       h3("MoA Enrichment"),
#       plotOutput("moa_enrichment", height = "600px")
#     )
#   )
# )
ui <- fluidPage(
  titlePanel("Graph from CSR (minimal view)"),
  
  sidebarLayout(
    sidebarPanel(
      width=1,
      selectInput(
        "color_by",
        "Color nodes by:",
        choices = c(
          "community",
          "cell_iname",
          "tas",
          "atc_class",
          "atc_desc",
          "moa"
        ),
        selected = "community"
      ),
      checkboxInput("drop_no_atc", "Drop nodes without ATC label", value = FALSE),
      sliderInput("edge_alpha", "Edge transparency:",
                  min = 0.01, max = 0.8, value = 0.02, step = 0.01),
      sliderInput("node_size", "Node size:",
                  min = 0.1, max = 8, value = 1, step = 0.025),
      checkboxInput("show_labels", "Show node labels", value = FALSE),
      checkboxInput("show_edges", "Show edges", value = FALSE),
      selectInput(
        "atc_level",
        "ATC level:",
        choices = c("Level 1" = "level1",
                    "Level 2" = "level2",
                    "Level 3" = "level3",
                    "Level 4" = "level4"),
        selected = "level1"
      )
    ),
    
    mainPanel(
      # Row 1: big graph tile
      fluidRow(
        column(
          width = 10,
          h4("Graph layout"),
          plotlyOutput("graph_plot", height = "500px")
        ),
        column(
          width = 10,
          h4("UMAP"),
          plotlyOutput("umap_plot", height = "500px")
        )
      ),
      br(),
      
      # Row 2: community sizes + TAS histogram
      fluidRow(
        column(
          width = 6,
          h4("Community sizes"),
          plotOutput("comm_bar", height = "250px")
        ),
        column(
          width = 6,
          h4("TAS score distribution"),
          plotOutput("tas_distribution", height = "250px")
        )
      ),
      br(),
      
      # Row 3: ATC + MoA enrichment
      fluidRow(
        column(
          width = 12,
          h4("ATC Enrichment"),
          plotOutput("atc_enrichment", height = "800px")
        )),
      fluidRow(
        column(
          width = 12,
          h4("MoA Enrichment"),
          plotOutput("moa_enrichment", height = "800px")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  filtered_nodes <- reactive({
    if (!input$drop_no_atc) {
      return(node_info)
    }
    node_info %>%
      filter(
        !is.na(level1) | 
          !is.na(level2) | 
          !is.na(level3) | 
          !is.na(level4)
      )
  })
  
  
  # subset precomputed layout based on filtered nodes
  filtered_layout <- reactive({
    keep <- filtered_nodes()$node
    graph_layout[graph_layout$name %in% keep, ]
  })
  
  output$selected_table <- renderTable({
  ids <- selected_ids()
  if (!length(ids)) return(NULL)

  node_info %>%
    dplyr::filter(node %in% ids) %>%
    dplyr::select(
      node_id, cmap_name, cell_iname, tas,
      community, level1, moa
    ) %>%
    head(30)  # cap to avoid huge tables
})
  
  output$umap_plot <- renderPlotly({
    df <- node_info
    
    plot_ly(
      data   = df,
      x      = ~UMAP_1,
      y      = ~UMAP_2,
      type   = "scattergl",
      mode   = "markers",
      color  = ~.data[[input$color_by]],   # reuse your color_by control
      key    = ~node,                      # important: key used for linked brushing
      source = "umap",
      marker = list(
        size    = input$node_size,
        opacity = 0.8
      )
    ) %>%
      layout(
        dragmode = "lasso",
        xaxis = list(title = "UMAP 1"),
        yaxis = list(title = "UMAP 2")
      )
  })
  
  selected_ids <- reactive({
    # use lasso/box selection on the UMAP plot
    d <- plotly::event_data("plotly_selected", source = "umap")
    if (is.null(d) || nrow(d) == 0) {
      return(character(0))
    }
    unique(d$key)  # 'key' will be node (see plots below)
  })
  
  atc_data <- reactive({
    lvl <- input$atc_level   # "level1"..."level4"
    
    df <- filtered_nodes() %>%      # <- use filtered nodes if you want enrichment to respect drop_no_atc
      filter(!is.na(.data[[lvl]])) %>%
      select(community, atc = all_of(lvl))
    
    N <- nrow(df)
    
    enr <- df %>%
      count(community, atc, name = "k") %>%
      left_join(df %>% count(atc,       name = "K"), by = "atc") %>%
      left_join(df %>% count(community, name = "n"), by = "community") %>%
      mutate(
        p         = phyper(k - 1, K, N - K, n, lower.tail = FALSE),
        score     = -log10(p + 1e-12),
        score_sig = ifelse(p < 0.05, score, NA_real_)
      )
    
    order <- enr %>%
      group_by(community) %>%
      summarize(max_score = max(score, na.rm = TRUE)) %>%
      arrange(desc(max_score)) %>%
      pull(community)
    
    enr$community <- factor(enr$community, levels = order)
    enr$atc       <- factor(enr$atc)
    
    enr
  })

  # ---- zoomable graph with plotly ----
  output$graph_plot <- renderPlotly({
    df  <- node_info
    sel <- selected_ids()
    
    # base layer: all nodes
    p <- plot_ly(
      data  = df,
      x     = ~graph_x,
      y     = ~graph_y,
      type  = "scattergl",
      mode  = "markers",
      color = ~.data[[input$color_by]],
      key   = ~node,
      marker = list(
        size    = input$node_size,
        opacity = if (length(sel) == 0) 0.8 else 0.15
      ),
      showlegend = TRUE
    )
    
    # overlay highlight for selected nodes
    if (length(sel) > 0) {
      sel_df <- df[df$node %in% sel, ]
      
      p <- p %>%
        add_markers(
          data = sel_df,
          x    = ~graph_x,
          y    = ~graph_y,
          type = "scattergl",
          mode = "markers",
          marker = list(
            size    = input$node_size + 2,
            color   = "black",
            opacity = 1,
            line    = list(color = "white", width = 1)
          ),
          inherit = FALSE,
          showlegend = FALSE
        )
    }
    
    p %>%
      layout(
        xaxis = list(title = "", showticklabels = FALSE, zeroline = FALSE),
        yaxis = list(title = "", showticklabels = FALSE, zeroline = FALSE)
      )
  })
  # ------------------------------------
  
  
  output$comm_bar <- renderPlot({
    node_info %>%
      count(community) %>%
      ggplot(aes(x = community, y = n)) +
      geom_col() +
      labs(title = "Community Sizes",
           x = "Community", y = "# Nodes") +
      theme_minimal()
  })
  
  output$tas_distribution <- renderPlot({ggplot(node_info, aes(x = community, y = tas, fill = community)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
    theme_minimal() +
    labs(
      title = "TAS Distribution per Community",
      x = "Community",
      y = "TAS"
    ) +
    theme(legend.position = "none")
  })
  
  output$atc_enrichment <- renderPlot({
    
    enr <- atc_data()
    
    level_label <- switch(
      input$atc_level,
      "level1" = "ATC Level 1",
      "level2" = "ATC Level 2",
      "level3" = "ATC Level 3",
      "level4" = "ATC Level 4"
    )
    
    ggplot(enr, aes(x = community, y = atc, fill = score_sig)) +
      geom_tile() +
      scale_fill_viridis_c(
        option   = "C",
        na.value = "grey95",
        name     = "-log10 p"
      ) +
      labs(
        title = paste("ATC Enrichment by Community:", level_label),
        x = "Community",
        y = level_label
      ) +
      theme_minimal()
  })
  
  # ---- MoA enrichment (per community) ----
  moa_data <- reactive({
    # start from annotated nodes
    df <- node_info %>%
      filter(!is.na(moa), moa != "") %>%
      select(community, moa)
    
    N <- nrow(df)
    
    enr <- df %>%
      count(community, moa, name = "k") %>%
      left_join(df %>% count(moa,       name = "K"), by = "moa") %>%
      left_join(df %>% count(community, name = "n"), by = "community") %>%
      mutate(
        p     = phyper(k - 1, K, N - K, n, lower.tail = FALSE),
        score = -log10(p + 1e-12),
        fold  = (k / n) / (K / N)      # enrichment over background
      ) %>%
      filter(p < 0.05)                 # keep only significant MoAs
    
    # keep top 3 MoAs per community
    enr <- enr %>%
      group_by(community) %>%
      slice_max(order_by = score, n = 3, with_ties = FALSE) %>%
      ungroup()
    
    # order communities by max enrichment
    comm_order <- enr %>%
      group_by(community) %>%
      summarize(max_score = max(score, na.rm = TRUE)) %>%
      arrange(desc(max_score)) %>%
      pull(community)
    
    # order MoAs by mean enrichment across communities
    moa_order <- enr %>%
      group_by(moa) %>%
      summarize(mean_score = mean(score, na.rm = TRUE)) %>%
      arrange(desc(mean_score)) %>%
      pull(moa)
    
    enr$community <- factor(enr$community, levels = comm_order)
    enr$moa       <- factor(enr$moa,       levels = moa_order)
    
    enr
  })
  
  
  output$moa_enrichment <- renderPlot({
    enr <- moa_data()
    
    ggplot(enr, aes(x = community, y = moa)) +
      geom_point(
        aes(size = fold, color = score),
        alpha = 0.9
      ) +
      scale_color_viridis_c(
        option = "C",
        name   = "-log10 p"
      ) +
      scale_size(
        range = c(2, 8),
        name  = "Fold enrichment"
      ) +
      labs(
        title = "MoA Enrichment by Community (dotplot, p < 0.05)",
        x = "Community",
        y = "Mechanism of Action"
      ) +
      theme_minimal() +
      theme(
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
  })
  
}

shinyApp(ui, server)
