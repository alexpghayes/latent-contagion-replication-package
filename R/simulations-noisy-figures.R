plot_noisy_estimates <- function(estimates) {
  paths <- c()

  has_noise_levels <- "noise_param_value" %in% names(estimates)

  group_vars <- c(
    "term_clean",
    "n",
    "rank",
    "model",
    "estimator",
    "estimator_type",
    "noise_type",
    "homophily",
    "expected_degree"
  )
  if (has_noise_levels) {
    group_vars <- c(group_vars, "noise_param_value")
  }

  summarized <- estimates |>
    mutate(
      bias = estimate - truth,
      covered = conf.low <= truth & truth <= conf.high
    ) |>
    summarize(
      mean_squared_error = mean(bias^2),
      mean_bias = mean(bias),
      sd_estimate = sd(estimate),
      mean_vif = mean(vif),
      coverage = mean(covered),
      .by = all_of(group_vars)
    )

  label_coef <- function(x) {
    case_match(
      x,
      "C1" ~ "W[1]",
      "C2" ~ "W[2]",
      "C3" ~ "W[3]",
      "rho" ~ "rho",
      "Xhat1" ~ "X[1]",
      "Xhat2" ~ "X[2]",
      "Xhat3" ~ "X[3]",
      "Xhat4" ~ "X[4]",
      "Xhat5" ~ "X[5]",
      .default = x
    )
  }

  noise_labels <- c(
    none = "True A (baseline)",
    gaussian = "Gaussian noise",
    capped = "Degree capped",
    flipped = "Edges flipped",
    missing = "Edges missing",
    ego = "Ego-centric",
    aggregated = "Aggregated (Y=AW)"
  )

  if (has_noise_levels) {
    middle_params <- estimates |>
      filter(noise_type != "none") |>
      group_by(noise_type) |>
      summarize(
        middle_param = sort(unique(noise_param_value))[ceiling(
          length(unique(noise_param_value)) / 2
        )],
        .groups = "drop"
      )

    estimates_main <- estimates |>
      left_join(middle_params, by = "noise_type") |>
      filter(noise_type == "none" | noise_param_value == middle_param) |>
      select(-middle_param)

    summarized_main <- summarized |>
      left_join(middle_params, by = "noise_type") |>
      filter(noise_type == "none" | noise_param_value == middle_param) |>
      select(-middle_param)
  } else {
    estimates_main <- estimates
    summarized_main <- summarized
  }

  for (est_type in c("ols", "tsls")) {
    for (model_type in c("peer", "latent")) {
      model_name <- paste0("y_", model_type, "_x")

      plot_data <- summarized_main |>
        filter(
          estimator_type == est_type,
          model == model_name
        ) |>
        mutate(term_label = label_coef(term_clean))

      if (nrow(plot_data) == 0) {
        next
      }

      p <- plot_data |>
        ggplot(aes(
          x = n,
          y = mean_squared_error,
          color = term_label
        )) +
        geom_point(size = 1.5, alpha = 0.8) +
        geom_line(linewidth = 0.6) +
        scale_x_log10(labels = label_log(digits = 2)) +
        scale_y_log10(labels = label_log(digits = 2)) +
        scale_color_viridis_d(labels = label_parse()) +
        facet_wrap(
          vars(noise_type),
          ncol = 3,
          labeller = as_labeller(noise_labels)
        ) +
        labs(
          x = "Number of nodes (log scale)",
          y = "Mean squared error (log scale)",
          color = "Coefficient",
          title = paste0(
            "MSE by coefficient: ",
            toupper(est_type),
            " estimator"
          ),
          subtitle = paste0(
            str_to_title(model_type),
            " contagion model (with X)"
          )
        ) +
        theme_minimal(11) +
        theme(
          legend.position = "right",
          plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold")
        )

      filename <- paste0("mse_all_coef_", est_type, "_", model_type, "_x.png")
      path <- here(
        "figures",
        "simulations",
        "noisy",
        "mse_coef",
        filename
      )
      ggsave(
        path,
        plot = p,
        width = 10,
        height = 7,
        dpi = 300,
        create.dir = TRUE
      )
      paths <- c(paths, path)
    }
  }

  for (est_type in c("ols", "tsls")) {
    for (model_type in c("peer", "latent")) {
      model_name <- paste0("y_", model_type, "_no_x")

      plot_data <- summarized_main |>
        filter(
          estimator_type == est_type,
          model == model_name
        ) |>
        mutate(term_label = label_coef(term_clean))

      if (nrow(plot_data) == 0) {
        next
      }

      p <- plot_data |>
        ggplot(aes(
          x = n,
          y = mean_squared_error,
          color = term_label
        )) +
        geom_point(size = 1.5, alpha = 0.8) +
        geom_line(linewidth = 0.6) +
        scale_x_log10(labels = label_log(digits = 2)) +
        scale_y_log10(labels = label_log(digits = 2)) +
        scale_color_viridis_d(labels = label_parse()) +
        facet_wrap(
          vars(noise_type),
          ncol = 3,
          labeller = as_labeller(noise_labels)
        ) +
        labs(
          x = "Number of nodes (log scale)",
          y = "Mean squared error (log scale)",
          color = "Coefficient",
          title = paste0(
            "MSE by coefficient: ",
            toupper(est_type),
            " estimator"
          ),
          subtitle = paste0(
            str_to_title(model_type),
            " contagion model (no X)"
          )
        ) +
        theme_minimal(11) +
        theme(
          legend.position = "right",
          plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold")
        )

      filename <- paste0(
        "mse_all_coef_",
        est_type,
        "_",
        model_type,
        "_no_x.png"
      )
      path <- here(
        "figures",
        "simulations",
        "noisy",
        "mse_coef",
        filename
      )
      ggsave(
        path,
        plot = p,
        width = 10,
        height = 7,
        dpi = 300,
        create.dir = TRUE
      )
      paths <- c(paths, path)
    }
  }

  coverage_data <- summarized_main |>
    filter(term_clean == "rho") |>
    mutate(
      model_label = case_match(
        model,
        "y_peer_x" ~ "Peer (X)",
        "y_peer_no_x" ~ "Peer (no X)",
        "y_latent_x" ~ "Latent (X)",
        "y_latent_no_x" ~ "Latent (no X)"
      )
    )

  for (est_type in c("ols", "tsls")) {
    plot_data <- coverage_data |>
      filter(estimator_type == est_type, !str_detect(model, "no_x"))

    if (nrow(plot_data) == 0) {
      next
    }

    p <- plot_data |>
      ggplot(aes(
        x = factor(n),
        y = noise_type,
        fill = coverage
      )) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(
        aes(label = sprintf("%.0f%%", coverage * 100)),
        size = 2.5
      ) +
      scale_fill_gradient2(
        low = "#d73027",
        mid = "#ffffbf",
        high = "#1a9850",
        midpoint = 0.8,
        limits = c(0, 1),
        labels = percent
      ) +
      scale_y_discrete(labels = noise_labels) +
      facet_wrap(vars(model_label), ncol = 2) +
      labs(
        x = "Number of nodes",
        y = NULL,
        fill = "Coverage",
        title = paste0("80% CI coverage for rho: ", toupper(est_type)),
        subtitle = "Green = good coverage (~80%), Red = poor coverage"
      ) +
      theme_minimal(11) +
      theme(
        legend.position = "right",
        plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        panel.grid = element_blank()
      )

    filename <- paste0("coverage_heatmap_rho_", est_type, ".png")
    path <- here("figures", "simulations", "noisy", "coverage", filename)
    ggsave(
      path,
      plot = p,
      width = 10,
      height = 6,
      dpi = 300,
      create.dir = TRUE
    )
    paths <- c(paths, path)
  }

  for (edeg in c("n^3/4")) {
    ribbon_data <- estimates_main |>
      filter(
        term_clean == "rho",
        expected_degree == edeg
      ) |>
      mutate(sq_error = (estimate - truth)^2) |>
      summarize(
        mse = mean(sq_error),
        .by = c(n, noise_type, estimator_type, model)
      ) |>
      mutate(
        model_label = case_match(
          model,
          "y_peer_x" ~ "Peer contagion",
          "y_peer_no_x" ~ "Peer (no X)",
          "y_latent_x" ~ "Latent contagion",
          "y_latent_no_x" ~ "Latent (no X)"
        ),
        estimator_label = paste("Latent", toupper(estimator_type))
      )

    p <- ribbon_data |>
      ggplot() +
      aes(
        x = n,
        y = mse,
        color = noise_type
      ) +
      geom_line() +
      geom_point() +
      scale_x_log10(labels = label_log(digits = 2)) +
      scale_y_log10(labels = label_log(digits = 2)) +
      scale_color_viridis_d(labels = noise_labels) +
      facet_grid(
        rows = vars(estimator_label),
        cols = vars(model_label)
      ) +
      labs(
        x = "Number of nodes (log scale)",
        y = "Mean squared error (log scale)",
        color = "Network measurement"
      ) +
      theme_minimal(12) +
      theme(
        text = element_text(family = "TeX Gyre Pagella"),
        legend.position = "bottom"
      )

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
      "noisy",
      glue("mse_rho_all-{edeg_nice}.png")
    )

    ggsave(
      path,
      plot = p,
      width = 6.7,
      height = 6.7 * 9 / 16,
      dpi = 500,
      create.dir = TRUE
    )
    paths <- c(paths, path)
  }
  paths
}
