#' @useDynLib ahull3D, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom rgl tmesh3d
#' @importFrom Rvcg vcgClean
NULL

#' Fast 3D Alpha Hull with Label Propagation
#' @param points Nx3 matrix of points
#' @param alpha Alpha value
#' @param input_truth Optional labels (point id)
#' @param volume Compute volume (default: FALSE)
#' @export
ahull3D <- function(points, alpha, input_truth = NULL, volume = FALSE) {
  # Input validation
  stopifnot(is.matrix(points), ncol(points) == 3L)
  stopifnot(alpha >= 0)
  
  if(is.null(input_truth)) {
    input_truth = rep(1, nrow(points))
    return_labels = FALSE
  } else {
    stopifnot(length(input_truth) == nrow(points))
    return_labels = TRUE
  }
  
  # Call C++ function
  result = .Call("_ahull3D_FAS_cpp_with_labels", 
                 t(points), alpha, input_truth, volume,
                 PACKAGE = "ahull3D")
  
  vertices = result$vertices
  vertex_labels = result$vertex_labels
  
  # Create mesh (original makeMesh logic)
  nvertices = ncol(vertices)
  if(nvertices == 0L) {
    message("The alpha shape is empty.")
    return(NULL)
  }
  
  mesh0 = tmesh3d(
    vertices    = vertices,
    indices     = matrix(1L:nvertices, nrow = 3L),
    homogeneous = FALSE
  )
  
  # Clean mesh
  mesh = vcgClean(mesh0, sel = c(0L, 7L), silent = TRUE)
  mesh[["remvert"]] = NULL
  mesh[["remface"]] = NULL
  
  # Add labels after deduplication
  if(return_labels) {
    # Map labels to deduplicated vertices
    if(requireNamespace("RANN", quietly = TRUE)) {
      cleaned_vertices = t(mesh$vb[1:3, ])
      orig_vertices = t(vertices[1:3, ])
      nn_result = RANN::nn2(orig_vertices, cleaned_vertices, k = 1)
      final_labels = vertex_labels[nn_result$nn.idx[, 1]]
    } else {
      # Simple fallback
      final_labels = vertex_labels[1:ncol(mesh$vb)]
    }
    
    attr(mesh, "input_truth") = final_labels
    attr(mesh, "face_truth") = matrix(final_labels[mesh$it], nrow = 3)
  }
  
  # Add volume if requested
  if(volume) {
    attr(mesh, "volume") = attr(vertices, "volume")
  }
  
  attr(mesh, "alpha") = alpha
  
  return(mesh)
}