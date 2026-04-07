load_xenium <- function(
    sample_id = "output-220149",
    data.folder = "data/Raw_data_files",
    mtx.dir = "cell_feature_matrix/matrix.mtx.gz",
    barcode.dir = "cell_feature_matrix/barcodes.tsv.gz",
    features.dir = "cell_feature_matrix/features.tsv.gz",
    cell.data = "cells.csv.gz",
    cell.data.barcodes = "cell_id",
    additional_metadata = c("Tissue"="Vagina", "Condition" = "HSV-"),
    coordinates = c("x_centroid", "y_centroid")
) {
  data.dir = file.path(data.folder, sample_id)
  library(Matrix)
  library(Seurat)
  library(tidyverse)

  mtx=readMM(file.path(data.dir, mtx.dir))
  barcodes=read.delim(file.path(data.dir, barcode.dir), header=F)
  features=read.delim(file.path(data.dir, features.dir), header=F)
  celldat=read.csv(file.path(data.dir, cell.data)) %>% column_to_rownames(cell.data.barcodes)

  colnames(mtx) <- barcodes$V1
  rownames(mtx) <- features$V2

  mtx <- mtx[features$V3 == 'Gene Expression',]
  features <- features[features$V3 == 'Gene Expression',]

  suppressWarnings(dat <- CreateSeuratObject(mtx, assay = "Xenium", meta.data = celldat, project = sample_id))

  if(!is_empty(additional_metadata)){
    for(i in 1:length(additional_metadata)){
      dat[[names(additional_metadata)[i]]] <- additional_metadata[i] %>% as.character()
    }
  }

  dat@assays$Xenium@meta.data <- features

  if(!is_empty(coordinates)){
    dat@reductions[["coordinates"]] <- CreateDimReducObject(as.matrix(data.frame(coord_1 = dat[[coordinates[1]]] %>% unlist(), coord_2=dat[[coordinates[2]]] %>% unlist(), row.names = colnames(dat))),
                                assay="Xenium",
                                key="coord_"
    )
  }
  return(dat)
}
dat <- load_xenium(
  sample_id = "output-220149",
  data.folder = "data/Raw_data_files",
  mtx.dir = "cell_feature_matrix/matrix.mtx.gz",
  barcode.dir = "cell_feature_matrix/barcodes.tsv.gz",
  features.dir = "cell_feature_matrix/features.tsv.gz",
  cell.data = "cells.csv.gz",
  cell.data.barcodes = "cell_id",
  additional_metadata = c("Tissue"="Vagina", "Condition" = "HSV-"),
  coordinates = c("x_centroid", "y_centroid"))
