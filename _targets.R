library(targets)
library(crew)
library(dplyr)
library(purrr)

data(glasgow, package = "latentnetmediate")

tar_option_set(
  packages = c(
    "broom",
    "car",
    "dplyr",
    "fastadi",
    "fastRG",
    "forcats",
    "glue",
    "ggplot2",
    "here",
    "ids",
    "igraph",
    "irlba",
    "ivreg",
    "janitor",
    "Matrix",
    "methods",
    "purrr",
    "readr",
    "rlang",
    "scales",
    "stringr",
    "tibble",
    "tidyr",
    "tidygraph",
    "withr"
  ),
  controller = crew_controller_local(workers = 12),
  format = "qs",
  memory = "transient",
  storage = "worker",
  retrieval = "worker",
  garbage_collection = TRUE,

  trust_timestamps = TRUE
)

tar_source()

list(
  tar_target(chunk_size, 10),
  tar_target(num_chunks, 10),

  tar_target(n, c(100, 163, 264, 430, 698, 1135, 1845)), #, 3000)),
  tar_target(rank, c(5)),
  tar_target(expected_degree, c("n^1/2", "n^1/4", "n^3/4")),

  tar_target(chunk_indices, 1:num_chunks),

  tar_target(
    params,
    tibble::tibble(
      n = n,
      rank = rank,
      expected_degree = expected_degree
    ),
    pattern = cross(n, rank, expected_degree)
  ),

  tar_target(
    population,
    purrr::map(
      1:chunk_size,
      ~ model_contagion_basic(
        n = params$n,
        k = params$rank,
        expected_degree = params$expected_degree
      )
    ),
    iteration = "list",
    pattern = cross(chunk_indices, map(params))
  ),

  tar_target(
    tbl_graph,
    purrr::map(population, sample_tidygraph_contagion),
    iteration = "list",
    pattern = map(population)
  ),

  tar_target(
    estimate,
    purrr::map2(tbl_graph, population, get_estimates, params$expected_degree),
    iteration = "list",
    pattern = map(tbl_graph, population, cross(chunk_indices, map(params)))
  ),

  tar_target(
    combined_estimates,
    dplyr::bind_rows(estimate, .id = "rep_within_chunk"),
    pattern = map(estimate)
  ),

  tar_target(
    basic_estimate_plots,
    plot_estimates(combined_estimates),
    format = "file"
  ),

  tar_target(
    noise_config,
    tibble::tribble(
      ~noise_type  , ~noise_param , ~noise_param_name ,
      # "none"       ,  0           , "no noise"        ,
      "gaussian"   ,  0.1         , "sigma"           ,
      "capped"     , 20           , "max_degree"      ,
      "flipped"    ,  0.15        , "flip_rate"       ,
      "missing"    ,  0.3         , "missing_rate"    ,
      "ego"        ,  0.5         , "fraction"        ,
      "aggregated" ,  0.8         , "corr"            ,
    )
  ),

  tar_target(
    noisy_estimate,
    map2(
      tbl_graph,
      population,
      get_noisy_estimates,
      noise_type = noise_config$noise_type,
      noise_param = noise_config$noise_param,
      noise_param_name = noise_config$noise_param_name,
      expected_degree = params$expected_degree
    ),
    pattern = cross(
      slice(
        map(tbl_graph, population, cross(chunk_indices, map(params))),
        c(3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54)
      ),
      map(noise_config)
    ),
    iteration = "list"
  ),

  tar_target(
    combined_noisy_estimates,
    bind_rows(noisy_estimate, .id = "rep_within_chunk"),
    pattern = map(noisy_estimate)
  ),

  tar_target(
    noisy_estimate_plots,
    plot_noisy_estimates(combined_noisy_estimates),
    format = "file"
  ),

  tar_target(
    concentration_figures,
    make_matrix_figures(),
    format = "file"
  ),

  tar_target(
    glasgow_clean,
    map(glasgow, clean_glasgow_network)
  ),

  tar_target(
    glasgow_estimates,
    get_glasgow_estimates(glasgow_clean[[1]])
  ),

  tar_target(
    glasgow_figures,
    plot_glasgow_estimates(glasgow_estimates),
    format = "file"
  )
)
