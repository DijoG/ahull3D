![C++](https://img.shields.io/badge/C++-RcppCGAL-blue?style=for-the-badge&logo=cplusplus)
![R](https://img.shields.io/badge/-%2764?style=for-the-badge&logo=r&logoColor=grey)
![Topology](https://img.shields.io/badge/Topology-Alpha%20Hull%20-blue?style=for-the-badge)

# ahull3D: Fast 3D Alpha Hull with Label Propagation

Ultra-fast 3D alpha shape computation with label propagation from input points to hull vertices and faces. Optimized for tree segmentation workflows.

## Installation

```r
remotes::install_github("DijoG/ahull3D")
```

## Usage example

```r
library(lidR)
library(dplyr)
library(ahull3D)

# Read tree point cloud
trees <- readLAS("forest.las")
trees_filtered <- filter_poi(trees, treeid %in% c(1,2,3,4,5))

# Prepare data
lasdff <- trees_filtered@data[,c("X", "Y", "Z", "treeid")] %>%
  as.data.frame() %>%
  distinct(across(1:3), .keep_all = TRUE) %>%
  as.matrix()

# Compute alpha hull with tree labels ("treeid" propagation)
tictoc::tic()
a <- 
  ahull3D::ahull3d(
  points = lasdff[,1:3],
  input_truth = lasdff[,4], 
  alpha = 0.1
  )
tictoc::toc()  # ~25 seconds for hulling/meshing 1214385 points

# Access and check results
input_truth <- attr(a, "input_truth")    # a vector of n vertices
face_truth <- attr(a, "face_truth")      # a matrix of 3*n faces 


length(input_truth) == ncol(a$vb)        # should be TRUE
dim(face_truth) == dim(a$it)             # should be TRUE TRUE
```

## Feautures

  - 🚀 Fast: 25 seconds for 130k points
  - 🏷️ Label propagation: Input labels → hull vertices → faces
  - 📏 Correct dimensions: input_truth matches vertices, face_truth is 3×faces
  
  

  