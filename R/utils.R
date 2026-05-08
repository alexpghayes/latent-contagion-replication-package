left_padded_sequence <- function(x) {
  original <- withr::with_options(
    c(scipen = 999),
    as.character(x)
  )

  max_digits <- max(vapply(original, nchar, integer(1)))
  formatC(x, width = max_digits, format = "d", flag = "0")
}

ASE <- function(ufm) {
  s <- fastRG::svds(ufm)
  k <- ufm$k

  # if (k > 1) {
  #   S <- diag(sqrt(s$d * 2))
  # } else {
  #   S <- matrix(sqrt(s$d * 2))
  # }

  if (k > 1) {
    S <- diag(sqrt(s$d))
  } else {
    S <- matrix(sqrt(s$d))
  }

  US <- s$u %*% S
  colnames(US) <- paste0("US", 1:ncol(US))
  US
}

US <- function(A, rank, ...) {
  s <- RSpectra::svds(A, k = rank, nu = rank, nv = rank, ...)
  if (rank > 1) {
    us <- s$u %*% diag(sqrt(s$d))
  } else {
    us <- s$u %*% matrix(sqrt(s$d))
  }

  colnames(us) <- as.character(1:rank)
  us
}

# rotate Y to align it to X
align <- function(X, Y) {
  XY <- crossprod(X, Y)
  s <- svd(XY)
  rotation <- s$v %*% t(s$u)
  Y %*% rotation
}


#' Take network averages as if na.rm = TRUE, kinda
#'
#' NAs propagate partially -- if x is NA, Ax will be NA for that same entry. however, that NA value of x will
#' just be ignored by X's peers, such that x's NA entry does not propage to all it's neighbors
#'
#' @param A sparse adjacency matrix (symmetric)
#' @param x vector (with NA values that should be ignored, as in na.rm = TRUE, except network average edition)
#'
#' A methodological note about this choice is important, since it makes a substantial difference
#'
safe_network_average <- function(A, x) {
  na_indices <- which(is.na(x))

  A_safe <- as.matrix(A)
  A_safe[, na_indices] <- 0
  A_safe[na_indices, ] <- 0
  A_safe <- drop0(A_safe)

  G_safe <- rowScale(A_safe, 1 / rowSums(A_safe))

  out <- drop(G_safe %*% x)
  out[na_indices] <- NA
  out
}

safe_network_average_matrix <- function(A, X) {
  apply(X, 2, safe_network_average, A = as.matrix(A))
}
