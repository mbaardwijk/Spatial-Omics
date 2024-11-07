library(SpatialExperiment)
library(Seurat)
library(SeuratData)
library(dplyr)
library(ggplot2)
library(cowplot)
library(grid)

rotate_image <- function(p, rot_angle) {
  gt <- ggplot_gtable(ggplot_build(p))
  panel_idx <- which(gt$layout$name == "panel")
  rot_vp <- viewport(angle = rot_angle)
  gt[["grobs"]][[panel_idx]] <- editGrob(gt[["grobs"]][[panel_idx]], vp = rot_vp)
  p_rot <- ggdraw() + draw_grob(gt)
  
  return(p_rot)
}

## Function (adapted from https://github.com/drighelli/SpatialExperiment/issues/115)
seurat_to_spe <- function(seu, sample_id, img_id, coords) {
  ## Convert to SCE
  sce <- Seurat::as.SingleCellExperiment(seu)
  
  ## Extract spatial coordinates
  spatialCoords <- as.matrix(
    coords[, c("imagecol", "imagerow")])
  
  ## Extract and process image data
  img <- SpatialExperiment::SpatialImage(
    x = as.raster(seu@images[[img_id]]@image))
  
  imgData <- DataFrame(
    sample_id = sample_id,
    image_id = img_id,
    data = I(list(img)),
    scaleFactor = seu@images[[img_id]]@scale.factors$lowres)
  
  # Convert to SpatialExperiment
  spe <- SpatialExperiment(
    assays = assays(sce),
    rowData = rowData(sce),
    colData = colData(sce),
    metadata = metadata(sce),
    reducedDims = reducedDims(sce),
    altExps = altExps(sce),
    sample_id = sample_id,
    spatialCoords = spatialCoords,
    imgData = imgData
  )
  # indicate all spots are on the tissue
  spe$in_tissue <- 1
  spe$sample_id <- sample_id
  # Return Spatial Experiment object
  spe
}

# Function to remove mitochondrial genes from a Seurat object
remove_mt_genes <- function(seu){
  counts <- GetAssayData(seu, assay = "Spatial")
  counts <- counts[-(grep(pattern = "^MT-", x = rownames(counts))),]
  seu <- subset(seu, features = rownames(counts))
  return(seu)
}
