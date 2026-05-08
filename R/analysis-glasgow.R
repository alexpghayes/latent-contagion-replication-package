clean_glasgow_network <- function(graph) {
  graph |>
    activate(edges) |>
    filter(friendship != "Structurally missing") |>
    activate(nodes) |>
    mutate(
      in_degree = centrality_degree(mode = "in"),
      out_degree = centrality_degree(mode = "out")
    ) |>
    filter(in_degree > 0 | out_degree > 0) |>
    mutate(
      tobacco_dimaria = as.numeric(tobacco_int > 1),
      alcohol_dimaria = as.numeric(alcohol_int > 2),
      cannabis_dimaria = as.numeric(cannabis_int > 2)
    )
}

uncentered_vif <- function(object, ...) {
  V <- summary(object)$cov.unscaled
  Vi <- crossprod(model.matrix(object))
  nam <- names(coef(object))
  v1 <- diag(V)
  v2 <- diag(Vi)
  structure(v1 * v2, names = nam)
}

tidy_plus <- function(fit, homophily) {
  tidied <- try(
    # covariance calculation can fail at low ranks for numerical reasons,
    # in which case we'll just punt on anything but point estimation
    tidy(fit, conf.int = TRUE, conf.level = 0.95),
    silent = TRUE
  )

  if (inherits(tidied, "try-error")) {
    tidied <- enframe(coef(fit), name = "term", value = "estimate")
  }

  vifs <- try(
    uncentered_vif(fit),
    silent = TRUE
  )

  if (!inherits(vifs, "try-error")) {
    vif_df <- enframe(vifs, name = "term", value = "vif")
    tidied <- tidied |>
      left_join(vif_df, by = join_by(term))
  }

  coef_names <- names(coef(fit))
  Xhat_names <- coef_names[str_detect(coef_names, "Xhat")]

  if (homophily) {
    xhat_anova_pvalue <- linearHypothesis(
      fit,
      Xhat_names,
      rep(0, length(Xhat_names))
    ) |>
      clean_names() |>
      filter(!is.na(pr_f)) |>
      pull(pr_f)
  } else {
    xhat_anova_pvalue <- NA
  }

  tidied |>
    mutate(
      xhat_anova_pvalue = xhat_anova_pvalue
    )
}

get_glasgow_estimates <- function(graph) {
  covariates <- c(
    "sex_fct",
    "age",
    "leisure_church",
    "money"
  )

  outcome_chr <- "tobacco_dimaria"
  outcome_sym <- rlang::sym(outcome_chr)

  graph <- graph |>
    activate(nodes) |>
    mutate(outcome := !!outcome_sym) |>
    select(outcome, one_of(covariates), -!!outcome_sym)

  A <- sign(as_adjacency_matrix(graph))

  nodes <- graph |>
    activate(nodes) |>
    as_tibble() |>
    mutate(
      peer_avg_outcome = safe_network_average(A, outcome)
    )

  ranks <- 2:25

  mf_svd <- irlba(A, max(ranks))
  Xhat_max <- mf_svd$u %*% diag(sqrt(mf_svd$d))
  colnames(Xhat_max) <- paste0("Xhat", 1:max(ranks))

  tsls_peer_at_rank <- function(rank = rank, homophily = TRUE) {
    Xhat <- Xhat_max[, 1:rank, drop = FALSE]
    Xhat_tbl <- as.data.frame(Xhat) |>
      as_tibble()

    features <- covariates

    if (homophily) {
      features <- c(features, colnames(Xhat_tbl))
    }

    # t4 has lots of missing data, beware
    data <- nodes |>
      bind_cols(Xhat_tbl) |>
      select(one_of(features), outcome, peer_avg_outcome)

    W4 <- model.matrix.lm(
      outcome ~ . - peer_avg_outcome,
      data = data,
      na.action = "na.pass" # keep, let downstream modeling functions drop NAs
    )

    # don't use intercept as an instrument (sanity check this choice)
    W4 <- W4[, -1]

    G4W4 <- safe_network_average_matrix(A, W4)
    G4G4W4 <- safe_network_average_matrix(A, G4W4)

    Z4 <- cbind(W4, G4W4, G4G4W4)
    Z4 <- as.matrix(Z4)

    tsls <- ivreg(
      outcome ~ . | Z4,
      data = data
    )

    tsls |>
      tidy_plus(homophily) |>
      mutate(
        rank = rank,
        estimator = "tsls",
        network = "peer",
        homophily = homophily
      )
  }

  # has numerical issues just like the OLS latent function, need to
  # investigate what is happening there
  tsls_latent_at_rank <- function(rank = rank, homophily = TRUE) {
    Xhat <- Xhat_max[, 1:rank, drop = FALSE]
    Xhat_tbl <- as.data.frame(Xhat) |>
      as_tibble()

    D <- diag(mf_svd$d[1:rank])

    if (rank == 1) {
      # diag() behaves differently for length 1 and length 1+ inputs, smh
      D <- mf_svd$d[1]
    }

    Phat <- mf_svd$u[, 1:rank, drop = FALSE] %*%
      D %*%
      t(mf_svd$v[, 1:rank, drop = FALSE])

    features <- covariates

    if (homophily) {
      features <- c(features, colnames(Xhat_tbl))
    }

    # t4 has lots of missing data, beware
    data <- nodes |>
      bind_cols(Xhat_tbl) |>
      select(one_of(features), outcome) |>
      # construct the latent contagion effect, but abuse notation call it the same
      # thing we called the peer contagion effect so we don't have to re-write code
      # that already works for the peer contagione effect
      mutate(
        peer_avg_outcome = safe_network_average(Phat, outcome)
      )

    W4 <- model.matrix.lm(
      outcome ~ . - peer_avg_outcome,
      data = data,
      na.action = "na.pass" # keep, let downstream modeling functions drop NAs
    )

    # don't use intercept as an instrument (sanity check this choice)
    W4 <- W4[, -1]

    G4W4 <- safe_network_average_matrix(Phat, W4)
    G4G4W4 <- safe_network_average_matrix(Phat, G4W4)

    Z4 <- cbind(W4, G4W4, G4G4W4)
    Z4 <- as.matrix(Z4)

    colnames(Z4) <- paste0("Z", 1:ncol(Z4))
    covariates <- setdiff(colnames(data), "outcome")
    instruments <- colnames(Z4)

    merged <- bind_cols(
      data,
      as.data.frame(Z4)
    )

    formula <- glue(
      "outcome ~ {paste0(covariates, collapse = ' + ')} | {paste0(instruments, collapse = ' + ')}"
    )

    tsls <- ivreg(
      as.formula(formula),
      data = merged
    )

    # can fail due to numerical issues at low ranks -- Ghat and Xhat are colinear at rank 1 definitionally,
    # and can have issues at other low ranks depending on the precise structure of Xhat
    tidied <- tidy_plus(tsls, homophily)

    mutate(
      tidied,
      rank = rank,
      estimator = "tsls",
      network = "latent",
      homophily = homophily
    )
  }

  estimates_homophily <- bind_rows(
    map_dfr(ranks, tsls_peer_at_rank),
    map_dfr(ranks, tsls_latent_at_rank)
  )

  estimates_no_homophily <- bind_rows(
    map_dfr(ranks, tsls_peer_at_rank, homophily = FALSE),
    map_dfr(ranks, tsls_latent_at_rank, homophily = FALSE)
  )

  bind_rows(
    estimates_homophily,
    estimates_no_homophily
  )
}

plot_glasgow_estimates <- function(estimates) {
  estimates |>
    filter(
      term == "peer_avg_outcome",
      rank > 2
    ) |>
    ggplot() +
    aes(
      x = rank,
      ymin = conf.low,
      y = estimate,
      ymax = conf.high,
      fill = homophily,
      color = homophily
    ) +
    geom_line() +
    geom_hline(yintercept = -1, color = "darkgray", linetype = "dashed") +
    geom_hline(yintercept = 0, color = "darkgray", linetype = "dashed") +
    geom_hline(yintercept = 1, color = "darkgray", linetype = "dashed") +
    facet_grid(
      cols = vars(network),
      labeller = labeller(
        network = c(
          peer = "Peer contagion",
          latent = "Latent contagion"
        )
      )
    ) +
    geom_ribbon(alpha = 0.4) +
    scale_color_viridis_d(begin = 0.15, end = 0.85) +
    scale_fill_viridis_d(begin = 0.15, end = 0.85) +
    theme_minimal(11) +
    theme(text = element_text(family = "TeX Gyre Pagella")) +
    labs(
      y = "Estimated contagion coefficient",
      x = "Embedding dimension",
      color = "X included",
      fill = "X included",
    )

  anxiety_estimates_path <- here("figures", "glasgow", "estimates.png")

  ggsave(
    anxiety_estimates_path,
    height = 6.5 * 9 / 16,
    width = 6.5,
    dpi = 300,
    # device = cairo_pdf,
    create.dir = TRUE
  )

  estimates |>
    filter(
      term == "peer_avg_outcome",
      !is.na(xhat_anova_pvalue)
    ) |>
    mutate(across(estimator, toupper)) |>
    ggplot() +
    aes(
      x = rank,
      y = xhat_anova_pvalue
    ) +
    geom_line() +
    facet_grid(
      # rows = vars(estimator),
      cols = vars(network),
      labeller = labeller(
        network = c(
          peer = "Peer contagion",
          latent = "Latent contagion"
        ),
        homophily = c(
          "TRUE" = "Adjusting for X",
          "FALSE" = "Not adjusting for X"
        ),
        estimator = c(
          ols = "OLS",
          tsls = "TSLS"
        )
      )
    ) +
    theme_minimal(15) +
    theme(text = element_text(family = "TeX Gyre Pagella")) +
    labs(
      title = expression(paste(
        "Omnibus testing ",
        beta[x] == 0,
        " and ",
        theta[x] == 0
      )),
      y = "P-value",
      x = "Embedding dimension"
    )

  omnibus_pvalue_path <- here("figures", "glasgow", "omnibus_pvalue.png")

  ggsave(
    omnibus_pvalue_path,
    height = 6.5 * 9 / 16,
    width = 6.5,
    dpi = 300,
    # device = cairo_pdf,
    create.dir = TRUE
  )

  c(anxiety_estimates_path, omnibus_pvalue_path)
}
