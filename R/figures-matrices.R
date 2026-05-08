dense_to_long_tidy <- function(dense) {
  dense |>
    as.matrix() |>
    as_tibble(rownames = "row") |>
    gather(col, value, -row) |>
    mutate_all(as.numeric)
}

make_matrix_figures <- function() {
  set.seed(26)

  k <- 5

  B <- matrix(0.1, k, k)
  diag(B) <- 0.8

  theta <- runif(100, min = 1, max = 10)

  model <- dcsbm(
    theta = theta,
    B = B,
    expected_density = 0.2,
    sort_nodes = TRUE
  )
  model

  P <- expectation(model)

  # note that color scales aren't going to match across these guys

  P_plot <- P |>
    dense_to_long_tidy() |>
    ggplot(aes(x = col, y = row, fill = value)) +
    geom_raster() +
    scale_y_reverse() +
    scale_fill_gradient2(
      high = "black",
      limits = c(0, 5),
      oob = oob_squish
    ) +
    theme_void() +
    theme(
      legend.position = "none"
    )

  P_path <- here("figures", "matrices", "P.png")

  ggsave(
    P_path,
    height = 5,
    width = 5,
    dpi = 500,
    create.dir = TRUE,
    bg = "white"
  )

  A <- sample_sparse(model)

  rownames(A) <- rownames(P)
  colnames(A) <- colnames(P)

  A_plot <- A |>
    as.matrix() |>
    dense_to_long_tidy() |>
    ggplot(aes(x = col, y = row, fill = value)) +
    geom_raster() +
    scale_y_reverse() +
    scale_fill_gradient2(
      high = "black",
      limits = c(0, 5),
      oob = oob_squish
    ) +
    theme_void() +
    theme(
      legend.position = "none"
    )

  A_path <- here("figures", "matrices", "A.png")

  ggsave(
    A_path,
    height = 5,
    width = 5,
    dpi = 500,
    create.dir = TRUE,
    bg = "white"
  )

  s <- irlba(A, k = k)
  Phat <- tcrossprod(s$u, rowScale(s$v, s$d))
  rownames(Phat) <- rownames(P)
  colnames(Phat) <- colnames(P)

  Phat_plot <- Phat |>
    dense_to_long_tidy() |>
    ggplot(aes(x = col, y = row, fill = value)) +
    geom_raster() +
    scale_y_reverse() +
    scale_fill_gradient2(high = "black") +
    theme_void() +
    theme(
      legend.position = "none"
    )

  Phat_path <- here("figures", "matrices", "Phat.png")

  ggsave(
    Phat_path,
    height = 5,
    width = 5,
    dpi = 500,
    create.dir = TRUE,
    bg = "white"
  )

  noise <- matrix(
    rnorm(prod(dim(A)), sd = 0.2),
    nrow = nrow(A),
    ncol = ncol(A)
  )

  A_dense <- as.matrix(A)
  A_noise <- A_dense + noise
  rownames(A_noise) <- rownames(P)
  colnames(A_noise) <- colnames(P)

  A_noise_plot <- A_noise |>
    dense_to_long_tidy() |>
    ggplot(aes(x = col, y = row, fill = value)) +
    geom_raster() +
    scale_y_reverse() +
    scale_fill_gradient2(
      high = "black",
      limits = c(0, 5),
      oob = oob_squish
    ) +
    theme_void() +
    theme(
      legend.position = "none"
    )

  A_noise_plot

  A_noise_path <- here("figures", "matrices", "A_noise.png")

  ggsave(
    A_noise_path,
    height = 5,
    width = 5,
    dpi = 500,
    create.dir = TRUE,
    bg = "white"
  )

  missing_indicator <- rsparsematrix(
    nrow = nrow(A),
    ncol = ncol(A),
    density = 0.1,
    rand.x = \(x) NA
  )
  A_missing <- A_dense + as.matrix(missing_indicator)

  A_missing_plot <- A_missing |>
    dense_to_long_tidy() |>
    ggplot(aes(x = col, y = row, fill = value)) +
    geom_raster() +
    scale_y_reverse() +
    scale_fill_gradient2(
      high = "black",
      limits = c(0, 5),
      oob = oob_squish,
      na.value = "firebrick"
    ) +
    theme_void() +
    theme(
      legend.position = "none"
    )

  A_missing_path <- here("figures", "matrices", "A_missing.png")

  ggsave(
    A_missing_path,
    height = 5,
    width = 5,
    dpi = 500,
    create.dir = TRUE,
    bg = "white"
  )

  s_noise <- irlba(A_noise, k = k)
  Phat_noise <- tcrossprod(s_noise$u, rowScale(s_noise$v, s_noise$d))
  rownames(Phat_noise) <- rownames(P)
  colnames(Phat_noise) <- colnames(P)

  Phat_noise_plot <- Phat_noise |>
    dense_to_long_tidy() |>
    ggplot(aes(x = col, y = row, fill = value)) +
    geom_raster() +
    scale_y_reverse() +
    scale_fill_gradient2(high = "black") +
    theme_void() +
    theme(
      legend.position = "none"
    )

  Phat_noise_path <- here("figures", "matrices", "Phat_noise.png")

  ggsave(
    Phat_noise_path,
    height = 5,
    width = 5,
    dpi = 500,
    create.dir = TRUE,
    bg = "white"
  )

  idx <- which(!is.na(A_missing), arr.ind = TRUE)
  vals <- A_missing[!is.na(A_missing)]

  # explicit zeroes retained, implicit zeroes denote missing values
  # inefficient but that's fine
  A_missing_sparse <- sparseMatrix(
    i = idx[, 1],
    j = idx[, 2],
    x = vals,
    dims = dim(A_missing),
    giveCsparse = TRUE
  )

  s_missing <- adaptive_impute(A_missing_sparse, rank = k)

  Phat_missing <- tcrossprod(s_missing$u, rowScale(s_missing$v, s_missing$d))
  rownames(Phat_missing) <- rownames(P)
  colnames(Phat_missing) <- colnames(P)

  Phat_missing_plot <- Phat_missing |>
    dense_to_long_tidy() |>
    ggplot(aes(x = col, y = row, fill = value)) +
    geom_raster() +
    scale_y_reverse() +
    scale_fill_gradient2(high = "black") +
    theme_void() +
    theme(
      legend.position = "none"
    )

  Phat_missing_path <- here("figures", "matrices", "Phat_missing.png")

  ggsave(
    Phat_missing_path,
    height = 5,
    width = 5,
    dpi = 500,
    create.dir = TRUE,
    bg = "white"
  )

  facet_data_helper <- function(
    X,
    dgp = c("true", "noisy", "missing"),
    matrix = c("observed", "population", "estimate")
  ) {
    dgp <- arg_match(dgp)
    matrix <- arg_match(matrix)

    X <- as.matrix(X)
    tidied <- dense_to_long_tidy(X)
    tidied$dgp <- dgp
    tidied$matrix <- matrix
    tidied
  }

  all_matrices_long <- bind_rows(
    facet_data_helper(A, "true", "observed"),
    facet_data_helper(A_noise, "noisy", "observed"),
    facet_data_helper(A_missing, "missing", "observed"),
    facet_data_helper(P, "true", "population"),
    facet_data_helper(P, "noisy", "population"),
    facet_data_helper(P, "missing", "population"),
    facet_data_helper(Phat, "true", "estimate"),
    facet_data_helper(Phat_noise, "noisy", "estimate"),
    facet_data_helper(Phat_missing, "missing", "estimate"),
  ) |>
    mutate(
      dgp = fct_relevel(as.factor(dgp), "true", "noisy", "missing"),
      matrix = fct_relevel(
        as.factor(matrix),
        "observed",
        "population",
        "estimate"
      )
    )

  all_matrices_plot <- all_matrices_long |>
    ggplot(aes(x = col, y = row, fill = value)) +
    geom_raster() +
    scale_y_reverse() +
    scale_fill_gradient2(
      high = "black",
      # limits = c(0, 5),
      oob = oob_squish,
      na.value = "steelblue"
    ) +
    theme_void(14) +
    theme(
      legend.position = "none",
      strip.text.y = element_text(angle = -90),
      text = element_text(family = "TeX Gyre Pagella")
    ) +
    facet_grid(
      rows = vars(matrix),
      cols = vars(dgp),
      labeller = labeller(
        dgp = c(
          true = "Exact network",
          noisy = "Noisy network",
          missing = "Partial network"
        ),
        matrix = c(
          observed = "Observed data",
          population = "Population structure",
          estimate = "Estimated structure"
        )
      )
    )

  all_matrices_path <- here("figures", "matrices", "all_matrices.png")

  ggsave(
    all_matrices_path,
    height = 6.5,
    width = 6.5,
    dpi = 500,
    create.dir = TRUE,
    bg = "white"
  )

  c(
    A_path,
    A_noise_path,
    A_missing_path,
    P_path,
    Phat_path,
    Phat_noise_path,
    Phat_missing_path,
    all_matrices_path
  )
}
