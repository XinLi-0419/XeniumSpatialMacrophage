library("sf")
library("purrr")
library("furrr")  # Parallel processing
library("Seurat")
library("ggplot2")
library("tidyr")
library("dplyr")
library("plyr")
library("grid")
library("gridExtra")
library("stringr")
library("ggrepel")

# home.dir <- "/Volumes/Jakubzick/20240823__201225__Xin_082324/"
# sample <- "output-XETG00256__0033690__Region_2__20240823__201234" #LPS # gene_expression_graphclust 10 and 24 for IMs # 26 not overlap with Pf4, C1qb, C1qc (checked with Explorer) # also not express C5ar1, which is different from the LPS ones
# resolutions <- "gene_expression_graphclust"

image.subset_by.gene <- function(xenium.obj.subset, cell_boundaries, genes = NULL) {
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




# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# image.subset <- function(xenium.obj.subset, cell_boundaries) {
#   library(sf)
#   library(dplyr)
#   library(purrr)
#   ### 1. Read and Process Cell Boundaries ###
#   cell_boundaries <- read.csv(gzfile(cell_boundaries), stringsAsFactors = FALSE)
#   # Convert boundary coordinates to polygons for each cell
#   cell_polygons_list <- cell_boundaries %>%
#     group_by(cell_id) %>%
#     summarise(geometry = list(st_polygon(list(
#       rbind(cbind(vertex_x, vertex_y), c(vertex_x[1], vertex_y[1])) # Ensure closure
#     )))) %>%
#     ungroup()
#   # cell_polygons_list <- cell_boundaries %>% # buggy
#   #   group_by(cell_id) %>%
#   #   summarise(geometry = list(st_polygon(list(cbind(vertex_x, vertex_y))))) %>%
#   #   ungroup()
#   # Convert to an sf object
#   cell_polygons_sf <- st_as_sf(cell_polygons_list)
#   ### 2. Identify Which Genes to Filter ###
#   gene_list <- names(xenium.obj.subset@images[["fov"]]@molecules$molecules)
#   ### 3. Loop Over Each Gene and Filter Spatial Data ###
#   for (gene in gene_list) {
#     # Extract gene transcript coordinates
#     gene_coords <- xenium.obj.subset@images[["fov"]]@molecules$molecules[[gene]]@coords
#     if (nrow(gene_coords) == 0) next  # Skip if no transcripts detected for this gene
#     # Convert gene coordinates to an sf points object
#     points_sf <- st_as_sf(data.frame(gene_coords), coords = c("x", "y"))
#     # Determine which points fall inside which cell polygon
#     point_in_poly <- st_within(points_sf, cell_polygons_sf)
#     # Attach a cell_id to each transcript
#     points_df <- data.frame(gene_coords) %>%
#       mutate(
#         point_id = row_number(),
#         cell_id = map_chr(point_in_poly, ~ if(length(.x) > 0) cell_polygons_sf$cell_id[.x[1]] else NA_character_)
#       )
#     # 4. Filter transcripts to only include those in the selected cluster
#     points_subset <- points_df %>% filter(cell_id %in% colnames(xenium.obj.subset))
#     # 5. Update the Xenium object with the filtered data
#     xenium.obj.subset@images[["fov"]]@molecules$molecules[[gene]]@coords <- as.matrix(points_subset[, c("x", "y")])
#   }
#   # Return the updated Xenium object
#   return(xenium.obj.subset)
# } # take a while # no parallel
# 
# image.subset <- function(xenium.obj.subset, cell_boundaries, parallel = TRUE) {
#   library(sf)
#   library(dplyr)
#   library(purrr)
#   library(furrr)  # Parallel processing
#   ### 1. Read and Process Cell Boundaries ###
#   cell_boundaries <- read.csv(gzfile(cell_boundaries), stringsAsFactors = FALSE)
#   # Convert boundary coordinates to polygons for each cell
#   cell_polygons_list <- cell_boundaries %>%
#     group_by(cell_id) %>%
#     summarise(geometry = list(st_polygon(list(cbind(vertex_x, vertex_y))))) %>%
#     ungroup()
#   # Convert to an sf object
#   cell_polygons_sf <- st_as_sf(cell_polygons_list)
#   ### 2. Identify Which Genes to Process ###
#   gene_list <- names(xenium.obj.subset@images[["fov"]]@molecules$molecules)
#   ### 3. Function to Process a Single Gene ###
#   process_gene <- function(gene) {
#     # Extract gene transcript coordinates
#     gene_coords <- xenium.obj.subset@images[["fov"]]@molecules$molecules[[gene]]@coords
#     if (nrow(gene_coords) == 0) return(NULL)  # Skip if no transcripts detected for this gene
#     # Convert gene coordinates to an sf points object
#     points_sf <- st_as_sf(data.frame(gene_coords), coords = c("x", "y"))
#     # Determine which points fall inside which cell polygon
#     point_in_poly <- st_within(points_sf, cell_polygons_sf)
#     # Attach a cell_id to each transcript
#     points_df <- data.frame(gene_coords) %>%
#       mutate(
#         point_id = row_number(),
#         cell_id = map_chr(point_in_poly, ~ if(length(.x) > 0) cell_polygons_sf$cell_id[.x[1]] else NA_character_)
#       )
#     # Filter transcripts to only include those in the selected cells
#     points_subset <- points_df %>% filter(cell_id %in% colnames(xenium.obj.subset))
#     # Return the filtered coordinates
#     return(as.matrix(points_subset[, c("x", "y")]))
#   }
#   ### 4. Run Processing in Parallel or Sequentially ###
#   if (parallel) {
#     plan(multisession)  # Enable parallel processing
#     gene_results <- future_map(gene_list, process_gene, .progress = TRUE)
#   } else {
#     gene_results <- map(gene_list, process_gene)
#   }
#   ### 5. Update the Xenium Object ###
#   names(gene_results) <- gene_list  # Assign names for lookup
#   for (gene in gene_list) {
#     if (!is.null(gene_results[[gene]])) {
#       xenium.obj.subset@images[["fov"]]@molecules$molecules[[gene]]@coords <- gene_results[[gene]]
#     }
#   }
#   # Return the updated Xenium object
#   return(xenium.obj.subset)
# } # parallel
# 
# image.subset_by.gene <- function(xenium.obj.subset, cell_boundaries, genes = NULL, parallel = TRUE) {
#   library(sf)
#   library(dplyr)
#   library(purrr)
#   library(furrr)  # Parallel processing
#   ### 1. Read and Process Cell Boundaries ###
#   cell_boundaries <- read.csv(gzfile(cell_boundaries), stringsAsFactors = FALSE)
#   # Convert boundary coordinates to polygons for each cell
#   cell_polygons_list <- cell_boundaries %>%
#     group_by(cell_id) %>%
#     summarise(geometry = list(st_polygon(list(cbind(vertex_x, vertex_y))))) %>%
#     ungroup()
#   # Convert to an sf object
#   cell_polygons_sf <- st_as_sf(cell_polygons_list)
#   ### 2. Identify Which Genes to Filter ###
#   available_genes <- names(xenium.obj.subset@images[["fov"]]@molecules$molecules)
#   # Use all genes if not specified
#   if (is.null(genes)) {
#     genes <- available_genes
#   } else {
#     # Ensure requested genes exist in the dataset
#     genes <- intersect(genes, available_genes)
#   }
#   ### 3. Function to Process a Single Gene ###
#   process_gene <- function(gene) {
#     # Extract gene transcript coordinates
#     gene_coords <- xenium.obj.subset@images[["fov"]]@molecules$molecules[[gene]]@coords
#     if (nrow(gene_coords) == 0) return(NULL)  # Skip if no transcripts detected for this gene
#     # Convert gene coordinates to an sf points object
#     points_sf <- st_as_sf(data.frame(gene_coords), coords = c("x", "y"))
#     # Determine which points fall inside which cell polygon
#     point_in_poly <- st_within(points_sf, cell_polygons_sf)
#     # Attach a cell_id to each transcript
#     points_df <- data.frame(gene_coords) %>%
#       mutate(
#         point_id = row_number(),
#         cell_id = map_chr(point_in_poly, ~ if(length(.x) > 0) cell_polygons_sf$cell_id[.x[1]] else NA_character_)
#       )
#     # Filter transcripts to only include those in the selected cluster
#     points_subset <- points_df %>% filter(cell_id %in% colnames(xenium.obj.subset))
#     # Return the filtered coordinates
#     return(as.matrix(points_subset[, c("x", "y")]))
#   }
#   ### 4. Run Processing in Parallel or Sequentially ###
#   if (parallel) {
#     plan(multisession)  # Enable parallel processing
#     gene_results <- future_map(genes, process_gene, .progress = TRUE)
#   } else {
#     gene_results <- map(genes, process_gene)
#   }
#   ### 5. Update the Xenium Object ###
#   names(gene_results) <- genes  # Assign names for lookup
#   for (gene in genes) {
#     xenium.obj.subset@images[["fov"]]@molecules$molecules[[gene]]@coords <- gene_results[[gene]]
#   }
#   # Return the updated Xenium object
#   return(xenium.obj.subset)
# } # Parallel processing takes too much memory and crushs the computer
# 
# image.subset_by.gene <- function(xenium.obj.subset, cell_boundaries, genes = NULL) { # Parallel processing takes too much memory and crushs the computer
#   library(sf)
#   library(dplyr)
#   library(purrr)
#   ### 1. Read and Process Cell Boundaries ###
#   cell_boundaries <- read.csv(gzfile(cell_boundaries), stringsAsFactors = FALSE)
#   # Convert boundary coordinates to polygons for each cell
#   cell_polygons_list <- cell_boundaries %>%
#     group_by(cell_id) %>%
#     summarise(geometry = list(st_polygon(list(
#       rbind(cbind(vertex_x, vertex_y), c(vertex_x[1], vertex_y[1])) # Ensure closure
#     )))) %>%
#     ungroup()
#   # Convert to an sf object
#   cell_polygons_sf <- st_as_sf(cell_polygons_list)
#   ### 2. Identify Which Genes to Filter ###
#   available_genes <- names(xenium.obj.subset@images[["fov"]]@molecules$molecules)
#   # Use all genes if not specified
#   if (is.null(genes)) {
#     genes <- available_genes
#   } else {
#     # Ensure requested genes exist in the dataset
#     genes <- intersect(genes, available_genes)
#   }
#   ### 3. Process Each Gene Sequentially ###
#   for (gene in genes) {
#     # Extract gene transcript coordinates
#     gene_coords <- xenium.obj.subset@images[["fov"]]@molecules$molecules[[gene]]@coords
#     if (nrow(gene_coords) == 0) next  # Skip if no transcripts detected for this gene
#     # Convert gene coordinates to an sf points object
#     points_sf <- st_as_sf(data.frame(gene_coords), coords = c("x", "y"))
#     # Determine which points fall inside which cell polygon
#     point_in_poly <- st_within(points_sf, cell_polygons_sf)
#     # Attach a cell_id to each transcript
#     points_df <- data.frame(gene_coords) %>%
#       mutate(
#         point_id = row_number(),
#         cell_id = map_chr(point_in_poly, ~ if(length(.x) > 0) cell_polygons_sf$cell_id[.x[1]] else NA_character_)
#       )
#     # Filter transcripts to only include those in the selected cluster
#     points_subset <- points_df %>% filter(cell_id %in% colnames(xenium.obj.subset))
#     # Update the Xenium object with the filtered data
#     xenium.obj.subset@images[["fov"]]@molecules$molecules[[gene]]@coords <- as.matrix(points_subset[, c("x", "y")])
#   }
#   # Return the updated Xenium object
#   return(xenium.obj.subset)
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# library(sf)
# library(dplyr)
# library(purrr)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# home.dir <- "/Volumes/Jakubzick/20240823__201225__Xin_082324/"
# sample <- "output-XETG00256__0033690__Region_2__20240823__201234" #LPS # gene_expression_graphclust 10 and 24 for IMs # 26 not overlap with Pf4, C1qb, C1qc (checked with Explorer) # also not express C5ar1, which is different from the LPS ones
# resolutions <- "gene_expression_graphclust"
# 
# mLu.Xenium.LPS.mLu.IMs <- image.subset_by.gene(xenium.obj.subset = mLu.Xenium.LPS.mLu.IMs, cell_boundaries = paste0(home.dir, sample, "/cell_boundaries.csv.gz"), genes = c("Ccl7", "Ccl2", "Ccl12", "Cxcl14", "Ccl5", "Ccl3", "Ccl4", "Cxcl1", "Cxcl2", "Cxcl3", "Ccl8", "Ccl6", "Ccl9", "Cxcl10", "Cxcl9", "Cxcl13", "Ccl24", "Cd3e", "Cd19")) # did not end up doing, always freeze
# 
# xenium.obj.subset = mLu.Xenium.LPS.mLu.IMs.Original
# cell_boundaries = paste0(home.dir, sample, "/cell_boundaries.csv.gz")
# genes = c("Ccl24")
# 
# 
# library(sf)
# library(dplyr)
# library(purrr)
# 
# ### 1. Read and Process Cell Boundaries ###
# cell_boundaries <- read.csv(gzfile(cell_boundaries), stringsAsFactors = FALSE)
# 
# # Convert boundary coordinates to polygons
# cell_polygons_sf <- cell_boundaries %>%
#   group_by(cell_id) %>%
#   summarise(geometry = st_sfc(st_polygon(list(
#     rbind(cbind(vertex_x, vertex_y), c(vertex_x[1], vertex_y[1])) # Ensure closure
#   ))), .groups = "drop") %>%
#   st_as_sf()  # Ensure sf object
# 
# 
# 
# 
# cell_polygons_sf <- cell_boundaries %>%
#   group_by(cell_id) %>%
#   summarise(geometry = list(st_polygon(list(
#     rbind(
#       cbind(vertex_x, vertex_y),   # Original boundary coordinates
#       c(vertex_x[1], vertex_y[1])  # Ensure the polygon closes
#     )
#   ))), .groups = "drop") %>%
#   st_as_sf()
# 
# 
# 
# cell_boundaries <- cell_boundaries[cell_boundaries$label_id %in% 1:10, ]
# 
# 
# cell_polygons_list <- cell_boundaries %>% # buggy
#   group_by(cell_id) %>%
#   summarise(geometry = list(st_polygon(list(cbind(vertex_x, vertex_y))))) %>%
#   ungroup()
# cell_polygons_sf <- st_as_sf(cell_polygons_list)
# 
# 
# 
# 
# library(sf)
# library(dplyr)
# library(dplyr)
# library(sf)
# 
# # For each cell_id, construct a polygon ensuring it's closed
# polys <- cell_boundaries %>%
#   group_by(cell_id) %>%
#   group_map(~{
#     # Extract coordinates as a matrix for the current cell group
#     coords <- as.matrix(.x[, c("vertex_x", "vertex_y")])
#     # If the polygon isn't closed, append the first coordinate to the end
#     if (!all(coords[1,] == coords[nrow(coords),])) {
#       coords <- rbind(coords, coords[1,])
#     }
#     st_polygon(list(coords))
#   })
# # Extract unique cell_ids using distinct() instead of summarise()
# cell_ids <- cell_boundaries %>%
#   distinct(cell_id) %>%
#   pull(cell_id)
# # Build an sf object from the polygons and cell_ids
# cell_polygons_sf <- st_sf(cell_id = cell_ids, geometry = st_sfc(polys, crs = NA_character_))
# 
# ### 2. Identify Which Genes to Filter ###
# available_genes <- names(xenium.obj.subset@images[["fov"]]@molecules$molecules)
# 
# if (is.null(genes)) {
#   genes <- available_genes
# } else {
#   genes <- intersect(genes, available_genes)
# }
# 
# ### 3. Process Genes Sequentially and Remove Empty Ones ###
# for (gene in genes) {
#   gene_coords <- xenium.obj.subset@images[["fov"]]@molecules$molecules[[gene]]@coords
#   if (nrow(gene_coords) == 0) {
#     xenium.obj.subset@images[["fov"]]@molecules$molecules[[gene]] <- NULL  # Remove empty genes
#     next
#   }
#   
#   # Convert gene coordinates to an sf points object
#   points_sf <- st_as_sf(data.frame(gene_coords), coords = c("x", "y"), crs = NA)
#   
#   # Use st_join() instead of st_within() to speed up processing
#   # joined_data <- st_join(points_sf, cell_polygons_sf, left = FALSE)
#   # 
#   # # Ensure cell_id column exists before filtering
#   # if (!"cell_id" %in% colnames(joined_data)) {
#   #   warning(paste("No cell_id found for gene:", gene))
#   #   xenium.obj.subset@images[["fov"]]@molecules$molecules[[gene]] <- NULL  # Remove genes with no valid cells
#   #   next
#   # }
#   # 
#   # # Keep only relevant columns and filter
#   # points_subset <- joined_data %>%
#   #   select(x, y, cell_id) %>%
#   #   filter(cell_id %in% colnames(xenium.obj.subset))
#   
#   
#   
#   
#   
#   
#   
#   # After performing st_join(), extract x and y coordinates:
#   joined_data <- st_join(points_sf, cell_polygons_sf, left = FALSE) %>%
#     mutate(x = st_coordinates(geometry)[,1],
#            y = st_coordinates(geometry)[,2])
#   
#   # Now you can select x, y, and cell_id:
#   points_subset <- joined_data %>%
#     select(x, y, cell_id) %>%
#     filter(cell_id %in% colnames(xenium.obj.subset))
#   
#   
#   
#   
#   
#   
#   
#   # If no remaining coordinates, remove the gene
#   if (nrow(points_subset) == 0) {
#     xenium.obj.subset@images[["fov"]]@molecules$molecules[[gene]] <- NULL  # Remove empty genes
#   } else {
#     # Update the Xenium object with the filtered data
#     xenium.obj.subset@images[["fov"]]@molecules$molecules[[gene]]@coords <- as.matrix(points_subset[, c("x", "y")])
#   }
# }
# 
# return(xenium.obj.subset)

