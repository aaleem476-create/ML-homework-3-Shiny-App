library(shiny)
library(ggplot2)
library(NMF)
library(uwot)
library(mclust)

# Read and prepare the TB NanoString data.
tb_data <- readRDS("TBnanostring.rds")
status <- factor(tb_data$TB_Status)
x <- tb_data[, setdiff(names(tb_data), "TB_Status"), drop = FALSE]
x[] <- lapply(x, as.numeric)
x <- as.matrix(x)

# Precompute the three two-dimensional representations.
pca_fit <- prcomp(x, center = TRUE, scale. = FALSE)
pca_xy <- as.data.frame(pca_fit$x[, 1:2])
names(pca_xy) <- c("Dim1", "Dim2")

set.seed(123)
umap_xy <- as.data.frame(uwot::umap(
  x,
  n_neighbors = 15,
  n_components = 2,
  metric = "euclidean"
))
names(umap_xy) <- c("Dim1", "Dim2")

# NMF requires non-negative values. Shift each gene by its minimum.
x_nmf <- sweep(x, 2, apply(x, 2, min, na.rm = TRUE), FUN = "-")
set.seed(123)
nmf_fit <- NMF::nmf(x_nmf, rank = 2, nrun = 10)
nmf_xy <- as.data.frame(NMF::basis(nmf_fit))
names(nmf_xy) <- c("Dim1", "Dim2")

# Store the precomputed coordinates in a named list.
embedding_list <- list(
  PCA = pca_xy,
  UMAP = umap_xy,
  NMF = nmf_xy
)

ui <- fluidPage(
  titlePanel("TB NanoString: Dimension Reduction and K-means"),
  sidebarLayout(
    sidebarPanel(
      selectInput(
        "method",
        "Dimension-reduction method",
        choices = c("PCA", "UMAP", "NMF"),
        selected = "PCA"
      ),
      numericInput(
        "clusters",
        "Number of clusters (k)",
        value = 2,
        min = 2,
        max = 10,
        step = 1
      )
    ),
    mainPanel(
      plotOutput("cluster_plot"),
      h4("Cluster by TB status"),
      tableOutput("cluster_table"),
      h4("Adjusted Rand index"),
      textOutput("ari_text")
    )
  )
)

server <- function(input, output, session) {
  
  # Select the already-computed two-dimensional coordinates.
  selectedData <- reactive({
    embedding_list[[input$method]]
  })
  
  # Recalculate k-means when the selected method or k changes.
  clusters <- reactive({
    set.seed(123)
    kmeans(selectedData(), centers = input$clusters, nstart = 25)
  })
  
  output$cluster_plot <- renderPlot({
    plot_data <- selectedData()
    plot_data$Cluster <- factor(clusters()$cluster)
    plot_data$TB_Status <- status
    
    center_data <- as.data.frame(clusters()$centers)
    names(center_data) <- c("Dim1", "Dim2")
    
    ggplot(
      plot_data,
      aes(
        x = Dim1,
        y = Dim2,
        color = Cluster,
        shape = TB_Status
      )
    ) +
      geom_point(size = 2.5, alpha = 0.8) +
      geom_point(
        data = center_data,
        aes(x = Dim1, y = Dim2),
        inherit.aes = FALSE,
        shape = 4,
        size = 6,
        stroke = 2,
        color = "black"
      ) +
      labs(
        title = paste(input$method, "with k-means clustering"),
        x = paste(input$method, "dimension 1"),
        y = paste(input$method, "dimension 2"),
        color = "K-means cluster",
        shape = "TB status"
      ) +
      theme_minimal()
  })
  
  output$cluster_table <- renderTable({
    table(
      Cluster = clusters()$cluster,
      TB_Status = status
    )
  }, rownames = TRUE)
  
  output$ari_text <- renderText({
    ari <- mclust::adjustedRandIndex(clusters()$cluster, status)
    paste("Adjusted Rand index:", round(ari, 3))
  })
}

shinyApp(ui = ui, server = server)
