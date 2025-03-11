# When subsetting a Xenium dataset in Seurat, the default subset function reduces the dataset by removing specific cells and their associated cell-by-transcript expression matrix. However, this operation does not modify or subset the spatial coordinates of individual transcripts stored in the Xenium assay. Thus, when plotting transcript-level data (e.g., using ImageFeaturePlot), the visualization still includes transcripts from the entire, original dataset, even though the cell-level data have been subsetted.

# To accurately visualize subsets at the transcript-coordinate level, additional spatial filtering steps (e.g., custom spatial subsetting functions like the image.subset function) must be performed, explicitly removing transcript coordinates outside of the defined cell subsets.

image.subset <- function(xenium.obj.subset, cell_boundaries, genes = NULL) {
  library(sf)
  library(dplyr)
  library(purrr)
  ### 1. Read and Process Cell Boundaries ###
  cell_boundaries <- read.csv(gzfile(cell_boundaries), stringsAsFactors = FALSE)
  # For each cell_id, construct a polygon ensuring it's closed
  polys <- cell_boundaries %>%
    group_by(cell_id) %>%
    group_map(~{
      # Extract coordinates as a matrix for the current cell group
      coords <- as.matrix(.x[, c("vertex_x", "vertex_y")])
      # If the polygon isn't closed, append the first coordinate to the end
      if (!all(coords[1,] == coords[nrow(coords),])) {
        coords <- rbind(coords, coords[1,])
      }
      st_polygon(list(coords))
    })
  # Extract unique cell_ids using distinct() instead of summarise()
  cell_ids <- cell_boundaries %>%
    distinct(cell_id) %>%
    pull(cell_id)
  # Build an sf object from the polygons and cell_ids
  cell_polygons_sf <- st_sf(cell_id = cell_ids, geometry = st_sfc(polys, crs = NA_character_))
  ### 2. Identify Which Genes to Filter ###
  available_genes <- names(xenium.obj.subset@images[["fov"]]@molecules$molecules)
  if (is.null(genes)) {
    genes <- available_genes
  } else {
    genes <- intersect(genes, available_genes)
  }
  ### 3. Process Genes Sequentially ###
  for (gene in genes) {
    gene_coords <- xenium.obj.subset@images[["fov"]]@molecules$molecules[[gene]]@coords
    if (nrow(gene_coords) == 0) next  # Skip if no transcripts detected
    # Convert gene coordinates to an sf points object
    points_sf <- st_as_sf(data.frame(gene_coords), coords = c("x", "y"), crs = NA)
    joined_data <- st_join(points_sf, cell_polygons_sf, left = FALSE) %>%
      mutate(x = st_coordinates(geometry)[,1],
             y = st_coordinates(geometry)[,2])
    # Now you can select x, y, and cell_id:
    points_subset <- joined_data %>%
      select(x, y, cell_id) %>%
      filter(cell_id %in% colnames(xenium.obj.subset))
    # Update the Xenium object with the filtered data
    coords_matrix <- st_coordinates(points_subset)
    colnames(coords_matrix) <- c("x", "y")
    xenium.obj.subset@images[["fov"]]@molecules$molecules[[gene]]@coords <- coords_matrix
  }
  return(xenium.obj.subset)
} # applicable for both all genes and apart of genes
