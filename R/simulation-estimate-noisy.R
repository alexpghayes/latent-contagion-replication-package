estimate_P_gaussian <- function(A_noisy, rank) {
  s <- irlba(A_noisy, nv = rank, nu = rank)
  s$u %*% diag(s$d, nrow = rank) %*% t(s$v)
}

estimate_P_capped <- function(A_capped, rank) {
  s <- irlba(A_capped, nv = rank, nu = rank)
  s$u %*% diag(s$d, nrow = rank) %*% t(s$v)
}

estimate_P_flipped <- function(A_flipped, rank) {
  s <- irlba(A_flipped, nv = rank, nu = rank)
  s$u %*% diag(s$d, nrow = rank) %*% t(s$v)
}

estimate_P_missing <- function(A_missing, rank) {
  # assumes A is a sparseMatrix where implicit zeroes
  # are missing, and all explicit zeroes are recorded as
  # zero-valued edge entries. this is computationally
  # inefficient, can return if it's too slow for our purposes
  mf <- adaptive_impute(
    A_missing,
    rank = rank,
    max_iter = 75,
    check_interval = 10
  )

  AT <- as(A_missing, "TsparseMatrix")

  # Create a mask of explicitly stored positions
  obs <- sparseMatrix(
    i = AT@i + 1,
    j = AT@j + 1,
    x = TRUE,
    dims = dim(AT)
  )

  missing <- as(1 - obs, "sparseMatrix")
  imputed <- fastadi:::masked_approximation(mf, missing)

  # note quite estimate P in this case, rather just imputing the missing
  # entries of A based on low-rank structure
  drop0(A_missing + imputed)
}

#' Estimate P from an ego-sampled partial network using the LE algorithm
#'
#' Implements Algorithm 1 + full matrix recovery from Chan & Li (2023),
#' "Fitting Low-rank Models on Egocentrically Sampled Partial Networks".
#'
#' @param A_ego Sparse adjacency matrix with non-ego × non-ego block zeroed out.
#' @param egos Integer vector of ego node indices.
#' @param rank Target rank K for the low-rank approximation.
#' @return A dense N × N estimated probability matrix P_hat with rows/columns
#'   in the original node ordering.
estimate_P_ego <- function(A_ego, egos, rank) {
  N <- nrow(A_ego)
  non_egos <- setdiff(seq_len(N), egos)
  n <- length(egos)

  # Extract blocks in ego/non-ego ordering
  A11 <- A_ego[egos, egos, drop = FALSE]
  A12 <- A_ego[egos, non_egos, drop = FALSE]

  # Step 1-2: Rank-K truncated SVD of A11 -> P_tilde_11
  s11 <- irlba(A11, nv = rank, nu = rank)
  # P_tilde_11 = U_K D_K V_K^T
  # P_tilde_11^+ = V_K D_K^{-1} U_K^T
  d_inv <- 1 / s11$d
  P_tilde_11_pinv <- s11$v %*% diag(d_inv, nrow = rank) %*% t(s11$u)

  # Step 3: P_hat_22 = A12^T P_tilde_11^+ A12
  P_hat_22 <- t(A12) %*% P_tilde_11_pinv %*% A12

  # Full matrix recovery (Section 2.4):
  # Rank-K SVD of A_obs = [A11, A12] for P_hat_11 and P_hat_12
  A_obs <- cbind(A11, A12)
  s_obs <- irlba(A_obs, nv = rank, nu = rank)
  P_tilde_obs <- s_obs$u %*% diag(s_obs$d, nrow = rank) %*% t(s_obs$v)

  P_hat_11 <- P_tilde_obs[, seq_len(n), drop = FALSE]
  # Symmetrize P_hat_11
  P_hat_11 <- (P_hat_11 + t(P_hat_11)) / 2
  P_hat_12 <- P_tilde_obs[, (n + 1):ncol(P_tilde_obs), drop = FALSE]
  P_hat_21 <- t(P_hat_12)

  # Assemble in ego/non-ego block order
  P_hat_blocked <- rbind(
    cbind(P_hat_11, P_hat_12),
    cbind(P_hat_21, P_hat_22)
  )

  # Reorder back to original node indices
  reorder <- order(c(egos, non_egos))
  P_hat <- P_hat_blocked[reorder, reorder]

  as.matrix(P_hat)
}

estimate_P_aggregated <- function(Y, W, rank) {
  s <- svd(Y)
  U <- s$u

  UtY <- crossprod(U, Y)
  UtW <- crossprod(U, W)

  Ntilde <- tcrossprod(UtY, UtW)
  Dtilde <- tcrossprod(UtW)

  Sigmatilde <- Ntilde %*% solve(Dtilde)

  Ptilde <- U %*% tcrossprod(Sigmatilde, U)
  Ptilde
}

estimate_P_none <- function(A, rank) {
  s <- irlba(A, nv = rank, nu = rank)
  s$u %*% diag(s$d, nrow = rank) %*% t(s$v)
}

corrupt_and_estimate <- function(A, population, noise_type, rank, ...) {
  noise_type <- match.arg(
    noise_type,
    c("none", "gaussian", "capped", "flipped", "missing", "ego", "aggregated")
  )

  result <- switch(
    noise_type,
    none = {
      P_hat <- estimate_P_none(A, rank)
      list(noisy_A = A, P_hat = P_hat, noise_type = "none")
    },
    gaussian = {
      A_noisy <- add_gaussian_noise(A, ...)
      P_hat <- estimate_P_gaussian(A_noisy, rank)
      list(noisy_A = A_noisy, P_hat = P_hat, noise_type = "gaussian")
    },
    capped = {
      A_noisy <- cap_degree(A, ...)
      P_hat <- estimate_P_capped(A_noisy, rank)
      list(noisy_A = A_noisy, P_hat = P_hat, noise_type = "capped")
    },
    flipped = {
      A_noisy <- flip_edges(A, ...)
      P_hat <- estimate_P_flipped(A_noisy, rank)
      list(noisy_A = A_noisy, P_hat = P_hat, noise_type = "flipped")
    },
    missing = {
      A_noisy <- remove_edges(A, ...)
      P_hat <- estimate_P_missing(A_noisy, rank)
      list(
        noisy_A = A_noisy,
        P_hat = P_hat,
        noise_type = "missing"
      )
    },
    ego = {
      ego_result <- ego_sample(A, ...)
      P_hat <- estimate_P_ego(ego_result$A_ego, ego_result$egos, rank)
      list(noisy_A = ego_result$A_ego, P_hat = P_hat, noise_type = "ego")
    },
    aggregated = {
      ard_data <- aggregated_relational_data(A, population, ...)
      P_hat <- estimate_P_aggregated(ard_data$Y, ard_data$W, rank)
      list(noisy_A = ard_data$Y, P_hat = P_hat, noise_type = "aggregated")
    }
  )

  result$rank <- rank
  result
}

Xhat_from_P <- function(P_hat, rank) {
  s <- irlba(P_hat, nv = rank, nu = rank)

  if (rank > 1) {
    Xhat <- s$u %*% diag(sqrt(s$d))
  } else {
    Xhat <- s$u * sqrt(s$d)
  }

  colnames(Xhat) <- paste0("Xhat", left_padded_sequence(1:rank))
  Xhat
}

compute_auxiliary_data_noisy <- function(tbl_graph, population, P_hat) {
  data <- tidygraph::as_tibble(tbl_graph)
  A <- as_adjacency_matrix(tbl_graph)
  deg <- degree(tbl_graph)
  DinvA <- Matrix::rowScale(A, 1 / deg)

  y_peer_no_x <- data$y_peer_no_x
  y_peer_x <- data$y_peer_x
  y_latent_no_x <- data$y_latent_no_x
  y_latent_x <- data$y_latent_x

  X <- ASE(population$A_model)
  colnames(X) <- paste0("Xhat", left_padded_sequence(seq_len(ncol(X))))
  X_reg <- X

  P <- tcrossprod(X)
  ED <- rowSums(P)
  EDinvP <- Matrix::rowScale(P, 1 / ED)
  diag(EDinvP) <- 0

  Xhat <- Xhat_from_P(P_hat, population$k)
  colnames(Xhat) <- paste0("Xhat", left_padded_sequence(seq_len(ncol(Xhat))))

  aligned <- align(X, Xhat)
  Xhat_reg <- aligned
  colnames(Xhat_reg) <- paste0(
    "Xhat",
    left_padded_sequence(seq_len(ncol(Xhat)))
  )

  Phat <- tcrossprod(Xhat)
  EDhat <- rowSums(Phat)
  EDinvPhat <- Matrix::rowScale(Phat, 1 / EDhat)
  diag(EDinvPhat) <- 0

  g_y_peer_no_x <- DinvA %*% y_peer_no_x
  gtilde_y_peer_no_x <- EDinvP %*% y_peer_no_x
  ghat_y_peer_no_x <- EDinvPhat %*% y_peer_no_x

  g_y_peer_x <- DinvA %*% y_peer_x
  gtilde_y_peer_x <- EDinvP %*% y_peer_x
  ghat_y_peer_x <- EDinvPhat %*% y_peer_x

  g_y_latent_no_x <- DinvA %*% y_latent_no_x
  gtilde_y_latent_no_x <- EDinvP %*% y_latent_no_x
  ghat_y_latent_no_x <- EDinvPhat %*% y_latent_no_x

  g_y_latent_x <- DinvA %*% y_latent_x
  gtilde_y_latent_x <- EDinvP %*% y_latent_x
  ghat_y_latent_x <- EDinvPhat %*% y_latent_x

  W <- as.matrix(population$W)

  Z_peer_no_x <- as.matrix(cbind(W, DinvA %*% W, DinvA %*% DinvA %*% W))
  Ztilde_peer_no_x <- as.matrix(cbind(W, EDinvP %*% W, EDinvP %*% EDinvP %*% W))
  Zhat_peer_no_x <- as.matrix(cbind(
    W,
    EDinvPhat %*% W,
    EDinvPhat %*% EDinvPhat %*% W
  ))

  colnames(Z_peer_no_x) <- paste0(
    "Z_peer_no_x",
    left_padded_sequence(seq_len(ncol(Z_peer_no_x)))
  )
  colnames(Ztilde_peer_no_x) <- paste0(
    "Ztilde_peer_no_x",
    left_padded_sequence(seq_len(ncol(Ztilde_peer_no_x)))
  )
  colnames(Zhat_peer_no_x) <- paste0(
    "Zhat_peer_no_x",
    left_padded_sequence(seq_len(ncol(Zhat_peer_no_x)))
  )

  Z_peer_x <- as.matrix(cbind(
    W,
    X_reg,
    DinvA %*% W,
    DinvA %*% X_reg,
    DinvA %*% DinvA %*% W,
    DinvA %*% DinvA %*% X_reg
  ))
  Ztilde_peer_x <- as.matrix(cbind(
    W,
    X_reg,
    EDinvP %*% W,
    EDinvP %*% X_reg,
    EDinvP %*% EDinvP %*% W,
    EDinvP %*% EDinvP %*% X_reg
  ))
  Zhat_peer_x <- as.matrix(cbind(
    W,
    X_reg,
    EDinvPhat %*% W,
    EDinvPhat %*% X_reg,
    EDinvPhat %*% EDinvPhat %*% W,
    EDinvPhat %*% EDinvPhat %*% X_reg
  ))

  colnames(Z_peer_x) <- paste0(
    "Z_peer_x",
    left_padded_sequence(seq_len(ncol(Z_peer_x)))
  )
  colnames(Ztilde_peer_x) <- paste0(
    "Ztilde_peer_x",
    left_padded_sequence(seq_len(ncol(Ztilde_peer_x)))
  )
  colnames(Zhat_peer_x) <- paste0(
    "Zhat_peer_x",
    left_padded_sequence(seq_len(ncol(Zhat_peer_x)))
  )

  Z_latent_no_x <- as.matrix(cbind(W, DinvA %*% W, DinvA %*% DinvA %*% W))
  Ztilde_latent_no_x <- as.matrix(cbind(
    W,
    EDinvP %*% W,
    EDinvP %*% EDinvP %*% W
  ))
  Zhat_latent_no_x <- as.matrix(cbind(
    W,
    EDinvPhat %*% W,
    EDinvPhat %*% EDinvPhat %*% W
  ))

  colnames(Z_latent_no_x) <- paste0(
    "Z_latent_no_x",
    left_padded_sequence(seq_len(ncol(Z_latent_no_x)))
  )
  colnames(Ztilde_latent_no_x) <- paste0(
    "Ztilde_latent_no_x",
    left_padded_sequence(seq_len(ncol(Ztilde_latent_no_x)))
  )
  colnames(Zhat_latent_no_x) <- paste0(
    "Zhat_latent_no_x",
    left_padded_sequence(seq_len(ncol(Zhat_latent_no_x)))
  )

  Z_latent_x <- as.matrix(cbind(
    W,
    X_reg,
    DinvA %*% W,
    DinvA %*% X_reg,
    DinvA %*% DinvA %*% W,
    DinvA %*% DinvA %*% X_reg
  ))
  Ztilde_latent_x <- as.matrix(cbind(
    W,
    X_reg,
    EDinvP %*% W,
    EDinvP %*% X_reg,
    EDinvP %*% EDinvP %*% W,
    EDinvP %*% EDinvP %*% X_reg
  ))
  Zhat_latent_x <- as.matrix(cbind(
    W,
    X_reg,
    EDinvPhat %*% W,
    EDinvPhat %*% X_reg,
    EDinvPhat %*% EDinvPhat %*% W,
    EDinvPhat %*% EDinvPhat %*% X_reg
  ))

  colnames(Z_latent_x) <- paste0(
    "Z_latent_x",
    left_padded_sequence(seq_len(ncol(Z_latent_x)))
  )
  colnames(Ztilde_latent_x) <- paste0(
    "Ztilde_latent_x",
    left_padded_sequence(seq_len(ncol(Ztilde_latent_x)))
  )
  colnames(Zhat_latent_x) <- paste0(
    "Zhat_latent_x",
    left_padded_sequence(seq_len(ncol(Zhat_latent_x)))
  )

  auxiliary_df <- tibble(
    g_y_peer_no_x = as.numeric(g_y_peer_no_x),
    gtilde_y_peer_no_x = as.numeric(gtilde_y_peer_no_x),
    ghat_y_peer_no_x = as.numeric(ghat_y_peer_no_x),

    g_y_peer_x = as.numeric(g_y_peer_x),
    gtilde_y_peer_x = as.numeric(gtilde_y_peer_x),
    ghat_y_peer_x = as.numeric(ghat_y_peer_x),

    g_y_latent_no_x = as.numeric(g_y_latent_no_x),
    gtilde_y_latent_no_x = as.numeric(gtilde_y_latent_no_x),
    ghat_y_latent_no_x = as.numeric(ghat_y_latent_no_x),

    g_y_latent_x = as.numeric(g_y_latent_x),
    gtilde_y_latent_x = as.numeric(gtilde_y_latent_x),
    ghat_y_latent_x = as.numeric(ghat_y_latent_x)
  ) |>
    bind_cols(
      as_tibble(Z_peer_no_x),
      as_tibble(Ztilde_peer_no_x),
      as_tibble(Zhat_peer_no_x),
      as_tibble(Z_peer_x),
      as_tibble(Ztilde_peer_x),
      as_tibble(Zhat_peer_x),
      as_tibble(Z_latent_no_x),
      as_tibble(Ztilde_latent_no_x),
      as_tibble(Zhat_latent_no_x),
      as_tibble(Z_latent_x),
      as_tibble(Ztilde_latent_x),
      as_tibble(Zhat_latent_x),
      as_tibble(Xhat_reg)
    )

  # sometimes instruments have columns of all NA and this breaks things downstream
  # i should probably investigate why this happens but first we're trying a hacky solution to
  # see if it fixes things
  auxiliary_df |>
    select(where(\(x) !all(is.na(x))))
}


get_noisy_estimates <- function(
  tbl_graph,
  population,
  noise_type,
  noise_param,
  noise_param_name,
  expected_degree
) {
  A <- as_adjacency_matrix(tbl_graph)

  result <- do.call(
    corrupt_and_estimate,
    c(
      list(
        A = A,
        population = population,
        noise_type = noise_type,
        rank = population$k
      ),
      noise_param
    )
  )
  P_hat <- result$P_hat

  data <- tidygraph::as_tibble(tbl_graph)
  auxiliary_data <- compute_auxiliary_data_noisy(tbl_graph, population, P_hat)
  data <- bind_cols(data, auxiliary_data)

  truth_no_x <- c(population$beta_w, population$rho)
  truth_x <- c(population$beta_w, population$beta_x, population$rho)

  estimates <- bind_rows(
    estimate_helper(data, "peer", "x", "ghat", truth_x),
    estimate_helper(data, "latent", "x", "ghat", truth_x)
  )

  estimates$n <- population$n
  estimates$rank <- population$k
  estimates$noise_type <- noise_type
  estimates$noise_param <- noise_param
  estimates$noise_param_name <- noise_param_name
  estimates$expected_degree <- expected_degree
  estimates
}
