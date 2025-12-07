library(shiny)
# -------------------------------------------------------------------------
# UI
# -------------------------------------------------------------------------

ui <- fluidPage(
  titlePanel("Drug Association Network Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      width = 2,
      selectInput(
        "color_by",
        "Color nodes by:",
        choices = c(
          "community",
          "cell_iname",
          "tas",
          "atc_level1" = "level1",
          "atc_level2" = "level2",
          "atc_level3" = "level3",
          "atc_level4" = "level4",
          "moa",
          "pert_time",
          "pert_dose"
        ),
        selected = "community"
      ),
      # --- NEW FILTERS ---
      selectizeInput(
        "cell_filter",
        "Cell line:",
        choices  = sort(unique(drug_nodes$cell_iname)),
        selected = sort(unique(drug_nodes$cell_iname)),
        multiple = TRUE
      ),
      selectizeInput(
        "time_filter",
        "Perturbation time (h):",
        choices  = sort(unique(drug_nodes$pert_time)),
        selected = sort(unique(drug_nodes$pert_time)),
        multiple = TRUE
      ),
      selectizeInput(
        "dose_filter",
        "Perturbation dose:",
        choices  = sort(unique(drug_nodes$pert_dose)),
        selected = sort(unique(drug_nodes$pert_dose)),
        multiple = TRUE
      ),
      checkboxInput("drop_no_atc", "Drop nodes without ATC label", value = FALSE),
      sliderInput("edge_alpha", "Edge transparency:",
                  min = 0.0, max = 1.0, value = 0.7, step = 0.1),
      sliderInput("node_size", "Base node size:",
                  min = 1, max = 10, value = 5, step = 1),
      checkboxInput("show_edges", "Show drug-level edges", value = TRUE),
      selectInput(
        "atc_level",
        "ATC level:",
        choices = c("Level 1" = "level1",
                    "Level 2" = "level2",
                    "Level 3" = "level3",
                    "Level 4" = "level4"),
        selected = "level1"
      ),
      selectInput(
        "selected_comm",
        "Community for summary:",
        choices  = sort(unique(as.character(drug_nodes$community))),
        selected = sort(unique(as.character(drug_nodes$community)))[1]
      )
    ),
    
    mainPanel(
      fluidRow(
        column(
          width = 12,
          h4("Graph layout (drug nodes)"),
          plotlyOutput("graph_plot", height = "500px")
        )
      ),
      br(),
      fluidRow(
        column(
          width = 12,
          h4("UMAP (drug nodes)"),
          plotlyOutput("umap_plot", height = "500px")
        )
      ),
      br(),
      fluidRow(
        column(
          width = 6,
          h4("Community sizes (drug-level)"),
          plotOutput("comm_bar", height = "250px")
        ),
        column(
          width = 6,
          h4("TAS distribution per community (drug-level)"),
          plotOutput("tas_distribution", height = "250px")
        )
      ),
      # br(),
      # fluidRow(
      #   column(
      #     width = 12,
      #     h4("Community robustness: within vs between edge weights"),
      #     plotOutput("robustness_scatter", height = "400px")
      #   )
      # ),
      br(),
      fluidRow(
        column(
          width=12,
          h4("Community robustness dumbbell: within vs between edge weights"),
          plotOutput("robustness_dumbbell", height = "600px")
        )
      ),
      br(),
      # fluidRow(
      #   column(
      #     width = 12,
      #     h4("Within-community signature-level edge weights (all communities)"),
      #     plotOutput("edge_ridge_all", height = "600px")
      #   )
      # ),
      # br(),
      # fluidRow(
      #   column(
      #     width = 12,
      #     h4("Between-community signature-level edge weights (all communities)"),
      #     plotOutput("edge_ridge_between", height = "600px")
      #   )
      # ),
      
      # Row 2b: MoA vs TAS
      fluidRow(
        column(
          width = 12,
          h4("TAS by MoA (drug-level)"),
          plotOutput("moa_tas_plot", height = "500px")
        )
      ),
      br(),
      fluidRow(
        column(
          width = 12,
          h4("ATC Enrichment"),
          plotOutput("atc_enrichment", height = "600px")
        )
      ),
      fluidRow(
        column(
          width = 12,
          h4("MoA Enrichment"),
          plotOutput("moa_enrichment", height = "600px")
        )
      ),
      br(),
      fluidRow(
        column(
          width = 12,
          h4("MoA \u2194 ATC cross-ontology mapping"),
          plotOutput("moa_atc_alluvial", height = "500px")
        )
      ),
      fluidRow(
        column(
          width = 12,
          h4(textOutput("summary_title")),
          tableOutput("community_summary"),
          
          # overall community stats
          tableOutput("community_summary"),
          
          br(),
          h5("Top MoAs in Community"),
          tableOutput("community_top_moas"),
          
          br(),
          h5("Top ATC Classes in Community"),
          tableOutput("community_top_atc"),
          br(),
          h5("Drugs in this Community"),
          tableOutput("community_drugs")
        )
      ),
      br(),
      fluidRow(
        column(
          width = 12,
          h4("Community MoA entropy"),
          plotOutput("moa_entropy", height = "350px")
        )
      ),
      br(),
      fluidRow(
        column(
          width = 12,
          h4("Mean TAS vs Community Size (drug-level)"),
          plotOutput("tas_vs_size", height = "500px")
        )
      ),
      br()
    )
  )
)
