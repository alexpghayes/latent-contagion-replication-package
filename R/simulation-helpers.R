estimate_helper <- function(
  data,
  contagion = c("peer", "latent"),
  homophily = c("x", "no_x"),
  estimator = c("g", "gtilde", "ghat"),
  truth
) {
  contagion <- match.arg(contagion)
  homophily <- match.arg(homophily)
  estimator <- match.arg(estimator)

  outcome <- glue("y_{contagion}_{homophily}")

  contagion_term <- glue("{estimator}_{outcome}")

  if (homophily == "no_x") {
    covariates_vec <- data |>
      select(
        matches("C[0-9]*")
      ) |>
      colnames()
  } else {
    covariates_vec <- data |>
      select(
        matches("C[0-9]*"),
        matches("Xhat[0-9]*")
      ) |>
      colnames()
  }

  # covariate names
  covariates <- paste0(covariates_vec, collapse = " + ")

  # instrument names
  estimator_short <- str_remove(estimator, "g")
  instruments_base <- glue("Z{estimator_short}_{contagion}_{homophily}")

  instruments_vec <- data |>
    select(matches(instruments_base)) |>
    colnames()

  instruments <- paste0(instruments_vec, collapse = " + ")

  formula <- glue(
    "{outcome} ~ {covariates} + {contagion_term} + 0 | {instruments} + 0"
  )
  tsls_fit <- ivreg(formula, data = data)

  tidied <- tsls_fit |>
    tidy(conf.int = TRUE, conf.level = 0.8) |>
    mutate(
      estimator_type = "tsls",
      term_clean = if_else(str_detect(term, "_y_"), "rho", term)
    )

  tidied$contagion <- contagion
  tidied$homophily <- homophily
  tidied$estimator <- estimator
  tidied$model <- outcome

  truth_df <- enframe(truth, "term_clean", "truth")

  left_join(tidied, truth_df, by = join_by(term_clean))
}

compute_auxiliary_data <- function(tbl_graph, population) {
  # Extract outcomes and graph structure
  data <- tidygraph::as_tibble(tbl_graph)
  A <- as_adjacency_matrix(tbl_graph)
  deg <- degree(tbl_graph)
  DinvA <- rowScale(A, 1 / deg)

  # Extract outcomes
  y_peer_no_x <- data$y_peer_no_x
  y_peer_x <- data$y_peer_x
  y_latent_no_x <- data$y_latent_no_x
  y_latent_x <- data$y_latent_x

  # Get population ASE for oracle estimator
  X <- ASE(population$A_model)
  colnames(X) <- paste0("Xhat", left_padded_sequence(seq_len(ncol(X))))
  X_reg <- X

  # Compute expected adjacency matrix for oracle estimator
  P <- tcrossprod(X)
  ED <- rowSums(P)
  EDinvP <- rowScale(P, 1 / ED)
  diag(EDinvP) <- 0

  # Compute estimated ASE
  Xhat <- US(A, rank = population$k)
  colnames(Xhat) <- paste0("Xhat", left_padded_sequence(seq_len(ncol(Xhat))))

  # Align Xhat to X
  aligned <- align(X, Xhat)
  Xhat_reg <- aligned
  colnames(Xhat_reg) <- paste0(
    "Xhat",
    left_padded_sequence(seq_len(ncol(Xhat)))
  )

  # Compute Ghat = D^{-1} Xhat Xhat^T
  Phat <- tcrossprod(Xhat)
  EDhat <- rowSums(Phat)
  EDinvPhat <- rowScale(Phat, 1 / EDhat)
  diag(EDinvPhat) <- 0

  # Compute contagion terms: g_y, gtilde_y, ghat_y for each outcome
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

  # Compute instruments for TSLS
  # For a given G, we want instruments: X, W, GX, GW, G^2W, G^2X
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

  # Combine all auxiliary data into a data frame
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

get_estimates <- function(tbl_graph, population, expected_degree) {
  data <- tidygraph::as_tibble(tbl_graph)

  auxiliary_data <- compute_auxiliary_data(tbl_graph, population)
  data <- bind_cols(data, auxiliary_data)

  truth_no_x <- c(population$beta_w, population$rho)
  truth_x <- c(population$beta_w, population$beta_x, population$rho)

  estimates <- bind_rows(
    estimate_helper(data, "peer", "no_x", "g", truth_no_x),
    estimate_helper(data, "peer", "no_x", "gtilde", truth_no_x),
    estimate_helper(data, "peer", "no_x", "ghat", truth_no_x),

    estimate_helper(data, "latent", "no_x", "g", truth_no_x),
    estimate_helper(data, "latent", "no_x", "gtilde", truth_no_x),
    estimate_helper(data, "latent", "no_x", "ghat", truth_no_x),

    estimate_helper(data, "peer", "x", "g", truth_x),
    estimate_helper(data, "peer", "x", "gtilde", truth_x),
    estimate_helper(data, "peer", "x", "ghat", truth_x),

    estimate_helper(data, "latent", "x", "g", truth_x),
    estimate_helper(data, "latent", "x", "gtilde", truth_x),
    estimate_helper(data, "latent", "x", "ghat", truth_x)
  )

  estimates$n <- population$n
  estimates$rank <- population$k
  estimates$expected_degree <- expected_degree
  estimates
}

plot_estimates <- function(estimates) {
  if (!dir.exists(here("figures/simulations"))) {
    dir.create(here("figures/simulations"))
  }

  summarized <- estimates |>
    mutate(
      bias = estimate - truth,
      covered = conf.low <= truth & truth <= conf.high
    ) |>
    summarize(
      mean_squared_error = mean(bias^2),
      mean_vif = mean(vif),
      coverage = mean(covered),
      .by = c(
        term_clean,
        n,
        rank,
        model,
        estimator,
        estimator_type,
        expected_degree
      )
    )

  path1 <- character(0)

  for (edeg in c("n^1/4", "n^1/2", "n^3/4")) {
    plot1 <- summarized |>
      filter(
        expected_degree == edeg,
        !str_detect(model, "no_x"),
        estimator != "gtilde"
      ) |>
      mutate(
        full_estimator = paste(estimator, estimator_type, sep = "_"),
        term_clean = case_match(
          term_clean,
          "C1" ~ "Z[1]",
          "C2" ~ "Z[2]",
          "C3" ~ "Z[3]",
          "rho" ~ "GY",
          "Xhat1" ~ "X[1]",
          "Xhat2" ~ "X[2]",
          "Xhat3" ~ "X[3]",
          "Xhat4" ~ "X[4]",
          "Xhat5" ~ "X[5]",
        )
      ) |>
      ggplot(aes(
        n,
        mean_squared_error,
        color = term_clean
      )) +
      geom_point() +
      geom_line() +
      scale_x_log10(labels = label_log(digits = 2)) +
      scale_y_log10(labels = label_log(digits = 2)) +
      scale_color_viridis_d(labels = label_parse()) +
      facet_grid(
        cols = vars(full_estimator),
        rows = vars(model),
        labeller = labeller(
          full_estimator = c(
            g_ols = "Estimator: Peer OLS",
            g_tsls = "Estimator: Peer TSLS",
            ghat_ols = "Estimator: Latent OLS",
            ghat_tsls = "Estimator: Latent TSLS"
          ),
          model = c(
            y_latent_x = "Model: Latent contagion",
            y_peer_x = "Model: Peer contagion"
          )
        )
      ) +
      labs(
        y = "Mean squared error (log scale)",
        x = "Number of nodes (log scale)",
        color = "Coef"
      ) +
      theme_minimal(12) +
      theme(text = element_text(family = "TeX Gyre Pagella"))

    if (edeg == "n^1/2") {
      edeg_nice <- "n12"
    } else if (edeg == "n^1/4") {
      edeg_nice <- "n14"
    } else if (edeg == "n^3/4") {
      edeg_nice <- "n34"
    }

    path <- here(
      "figures",
      "simulations",
      glue("mean_squared_error-{edeg_nice}.png")
    )

    ggsave(
      path,
      plot = plot1,
      width = 6.7,
      height = 6.7 * 9 / 16,
      dpi = 500,
      bg = "white"
    )

    path1 <- c(path1, path)
  }

  plot2 <- summarized |>
    filter(
      expected_degree == "n^3/4",
      !str_detect(model, "no_x"),
      estimator != "gtilde"
    ) |>
    mutate(
      term_clean = case_match(
        term_clean,
        "C1" ~ "Z[1]",
        "C2" ~ "Z[2]",
        "C3" ~ "Z[3]",
        "rho" ~ "GY",
        "Xhat1" ~ "X[1]",
        "Xhat2" ~ "X[2]",
        "Xhat3" ~ "X[3]",
        "Xhat4" ~ "X[4]",
        "Xhat5" ~ "X[5]",
      )
    ) |>
    ggplot(aes(n, coverage, color = term_clean)) +
    geom_point() +
    geom_line() +
    scale_x_log10(labels = label_log(digits = 2)) +
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8)) +
    scale_color_viridis_d(labels = label_parse()) +
    facet_grid(
      cols = vars(estimator),
      rows = vars(model),
      labeller = labeller(
        estimator = c(
          g = "Peer Estimator",
          ghat = "Latent Estimator"
        ),
        model = c(
          y_latent_x = "Latent contagion",
          y_peer_x = "Peer contagion"
        )
      ),
    ) +
    labs(
      y = "Coverage rate",
      x = "Number of nodes (log scale)"
    ) +
    theme_minimal(16)

  path2 <- here("figures", "simulations", "coverage-n34.png")

  ggsave(
    path2,
    plot = plot2,
    width = 12,
    height = 8,
    dpi = 300,
    bg = "white"
  )

  plot_estimates <- function(term_clean = "rho") {
    plot4 <- estimates |>
      filter(
        term_clean == {{ term_clean }},
        expected_degree == "n^3/4",
        !str_detect(model, "no_x")
      ) |>
      ggplot(aes(n, estimate, color = estimator)) +
      geom_jitter(position = position_dodge2(width = 1.5)) +
      geom_hline(yintercept = 0.2, linetype = "dashed", color = "darkgray") +
      scale_x_log10(labels = label_log(digits = 2)) +
      scale_color_viridis_d() +
      facet_grid(
        cols = vars(estimator, estimator_type),
        rows = vars(model),
        labeller = labeller(.rows = label_value, .cols = label_both),
        scales = "free"
      ) +
      labs(
        y = "Estimate",
        x = "Number of nodes (log scale)",
        title = "Rho point estimates",
        subtitle = "dashed darkgray line at true value of rho (not intercept!)"
      ) +
      theme_bw(16)

    path4 <- here(
      "figures",
      "simulations",
      glue("estimates-{term_clean}-n34.png")
    )

    ggsave(
      path4,
      plot = plot4,
      width = 8,
      height = 8,
      dpi = 300,
      bg = "white"
    )

    path4
  }

  terms <- c("rho", "Xhat1", "Xhat2", "Xhat3", "Xhat4", "Xhat5")
  term_paths <- map_chr(terms, plot_estimates)

  c(path1, path2, term_paths)
}
