add_gaussian_noise <- function(A, sigma = 0.1) {
  n <- nrow(A)
  noise <- matrix(rnorm(n * n, mean = 0, sd = sigma), nrow = n, ncol = n)
  noise <- (noise + t(noise)) / 2
  diag(noise) <- 0
  as.matrix(A) + noise
}

cap_degree <- function(A, max_degree = 5) {
  triplets <- summary(triu(A))

  edges <- data.frame(
    i = triplets$i,
    j = triplets$j,
    x = triplets$x
  )

  calc_degrees <- function(edge_df, num_nodes) {
    d_i <- tapply(edge_df$x, edge_df$i, sum)
    d_j <- tapply(edge_df$x, edge_df$j, sum)
    degs <- numeric(num_nodes)

    degs[as.integer(names(d_i))] <- degs[as.integer(names(d_i))] + d_i
    degs[as.integer(names(d_j))] <- degs[as.integer(names(d_j))] + d_j

    degs
  }

  num_nodes <- nrow(A)
  current_degrees <- calc_degrees(edges, num_nodes)

  while (any(current_degrees > max_degree)) {
    bad_nodes <- which(current_degrees > max_degree)

    target_node <- 1
    if (length(bad_nodes) > 1) {
      target_node <- sample(bad_nodes, 1)
    } else {
      target_node <- bad_nodes
    }

    relevant_indices <- which(edges$i == target_node | edges$j == target_node)

    if (length(relevant_indices) == 0) {
      next
    }

    edge_weights <- edges$x[relevant_indices]

    if (length(relevant_indices) == 1) {
      idx_to_decrement <- relevant_indices
    } else {
      idx_to_decrement <- sample(relevant_indices, 1, prob = edge_weights)
    }

    edges$x[idx_to_decrement] <- edges$x[idx_to_decrement] - 1

    node_a <- edges$i[idx_to_decrement]
    node_b <- edges$j[idx_to_decrement]

    current_degrees[node_a] <- current_degrees[node_a] - 1
    current_degrees[node_b] <- current_degrees[node_b] - 1

    if (edges$x[idx_to_decrement] <= 0) {
      edges <- edges[-idx_to_decrement, ]
    }
  }

  sparseMatrix(
    i = edges$i,
    j = edges$j,
    x = edges$x,
    dims = dim(A),
    symmetric = TRUE
  )
}

flip_edges <- function(A, flip_rate = 0.1) {
  triplets <- summary(A)
  edges <- triplets[triplets$i < triplets$j, ]

  n_edges <- nrow(edges)
  n_swaps_target <- floor(n_edges * flip_rate)

  swaps_done <- 0

  while (swaps_done < n_swaps_target) {
    idx <- sample(n_edges, 2)
    idx1 <- idx[1]
    idx2 <- idx[2]

    if (edges$x[idx1] == edges$x[idx2]) {
      next
    }

    u <- edges$i[idx1]
    v <- edges$j[idx1]
    x <- edges$i[idx2]
    y <- edges$j[idx2]

    if (runif(1) > 0.5) {
      new_v <- y
      new_y <- v
    } else {
      new_v <- x
      new_y <- u # This would imply swapping u and x effectively
      # To keep it simple in code, let's just swap j coordinates (target nodes)
      # But we must handle the case where u--y or x--v creates a self-loop
    }

    if (u == y || x == v) {
      next
    }

    edges$j[idx1] <- y
    edges$j[idx2] <- v

    if (edges$i[idx1] > edges$j[idx1]) {
      tmp <- edges$i[idx1]
      edges$i[idx1] <- edges$j[idx1]
      edges$j[idx1] <- tmp
    }

    if (edges$i[idx2] > edges$j[idx2]) {
      tmp <- edges$i[idx2]
      edges$i[idx2] <- edges$j[idx2]
      edges$j[idx2] <- tmp
    }

    swaps_done <- swaps_done + 1
  }

  sparseMatrix(
    i = edges$i,
    j = edges$j,
    x = edges$x,
    dims = dim(A),
    symmetric = TRUE
  )
}

# given a sparse adjacency matrix A, we want to set entries of A to NA
# via a missing at random assumption. this requires that we track zeroes
# and missingness separately, which is an enormous PITA
remove_edges <- function(A, missing_rate = 0.1) {
  A_upper <- triu(A, k = 1)

  mask <- rsparsematrix(
    nrow = nrow(A),
    ncol = ncol(A),
    symmetric = TRUE,
    density = missing_rate,
    rand.x = \(x) 1
  )

  mask_tbl <- mask |>
    triu() |>
    as("TsparseMatrix") |>
    summary() |>
    as_tibble()

  # explicit zeroes for all the non-missing zeroes
  explicitT <- A_upper |>
    as("TsparseMatrix") |>
    summary() |>
    as_tibble() |>
    complete(i, j) |>
    anti_join(mask_tbl, by = join_by(i, j)) |>
    mutate(
      x = replace_na(x, 0)
    )

  upper <- sparseMatrix(
    i = explicitT$i,
    j = explicitT$j,
    x = explicitT$x,
    dims = dim(A)
  )

  upper + t(upper)
}


#' Generate data with specified correlation to X
#' @param y Vector to correlate with
#' @param corr Target correlation value
#' @return Vector with specified correlation to y
complement <- function(y, corr = 0.5) {
  x <- rnorm(length(y))
  y.perp <- residuals(lm(x ~ y))
  corr * sd(y.perp) * y + y.perp * sd(y) * sqrt(1 - corr^2)
}

generate_traits <- function(U, corr, num_traits) {
  k <- ncol(U)

  if (length(corr) == 1) {
    corr <- rep(corr, num_traits)
  } else {
    # must specify corr either once as a scalar, or as
    # a vector for all traits all at once
    stopifnot(length(corr) == num_traits)
  }

  traits <- matrix(0, nrow = nrow(U), ncol = num_traits)

  for (trait in 1:num_traits) {
    # repeat traits if num_traits > k
    trait_index <- trait %% k
    if (trait_index == 0) {
      trait_index <- k
    }
    u <- U[, trait_index, drop = TRUE]
    traits[, trait] <- complement(u, corr = corr[trait])
  }

  traits
}

#' Ego-centric network sampling
#'
#' Randomly sample a fraction of nodes as "egos" and retain only edges incident
#' to at least one ego node. All other edges are set to zero. The returned
#' matrix has the same dimensions as the input.
#'
#' @param A Symmetric adjacency matrix (sparse or dense).
#' @param fraction Proportion of nodes to sample as egos (default 0.5).
#' @return A list with components:
#'   - `A_ego`: A sparse adjacency matrix where only edges incident to at least
#'     one sampled ego are retained.
#'   - `egos`: Integer vector of sampled ego node indices.
ego_sample <- function(A, fraction = 0.5) {
  n <- nrow(A)
  n_egos <- floor(n * fraction)
  egos <- sort(sample.int(n, size = n_egos))

  # Zero out rows AND columns for non-ego nodes
  non_egos <- setdiff(seq_len(n), egos)

  A_ego <- as(A, "dgCMatrix")
  A_ego[non_egos, non_egos] <- 0

  list(
    A_ego = drop0(A_ego),
    egos = egos
  )
}

aggregated_relational_data <- function(
  A,
  population,
  corr,
  num_traits = 5
) {
  s_pop <- svds(population$A_model)
  X <- s_pop$u %*% diag(sqrt(s_pop$d))

  traits <- generate_traits(X, corr = corr, num_traits = num_traits)
  Y <- A %*% traits

  list(
    Y = Y,
    W = traits
  )
}
