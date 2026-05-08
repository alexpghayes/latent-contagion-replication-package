#' @examples
#'
#' set.seed(26)
#'
#' b <- model_contagion_basic(n = 200, k = 5)
#'
#' graph <- sample_tidygraph_contagion(b)
#' graph
#'
#' A <- as_adjacency_matrix(graph)
#'
#' data <- as_tibble(graph)
#' truth <- c(b$beta_w, b$beta_x, b$rho)
#'
#' estimate_helper(data, "peer", "x", "g", truth)
#'
#' graph |>
#'   as_tibble() |>
#'   select(contains("y"), -contains("latent")) |>
#'   View()
#'
#' plot_tidygraph_matrix(graph)
#' plot_tidygraph_matrix_rownormalized(graph)
#'
#' hist(igraph::degree(graph))
#'
#'
model_contagion_basic <- function(
  n,
  k = 5,
  dim_c = 3,
  rho = 0.2,
  expected_degree = c("n^1/2", "n^1/4", "n^3/4")
) {
  expected_degree <- match.arg(expected_degree)

  if (expected_degree == "n^1/2") {
    expected_degree <- n^(1 / 2)
  } else if (expected_degree == "n^1/4") {
    expected_degree <- n^(1 / 4)
  } else if (expected_degree == "n^3/4") {
    expected_degree <- n^(3 / 4)
  }

  B <- matrix(0.05, nrow = k, ncol = k)
  diag(B) <- 0.5

  theta <- rexp(n, 1 / 4) + 1

  A_model <- dcsbm(
    theta = theta, # heterogeneous degrees
    B = B,
    allow_self_loops = FALSE,
    poisson_edges = TRUE,
    expected_degree = expected_degree
  )

  # A_model <- dcsbm(
  #   theta = rexp(n) + 1, # rep(1, n), # rexp(n),  # heterogeneous degrees
  #   B = B,
  #   pi = pi,
  #   allow_self_loops = FALSE,
  #   poisson_edges = FALSE
  # )
  #
  # desired_min_degree <- 4 * n^0.6 # max(4 * sqrt(n), 0.3 * n)
  # min_expected_degree <- min(expected_degrees(A_model))
  # rescale_by <- desired_min_degree / min_expected_degree
  # A_model$S <- A_model$S * rescale_by

  # trt column expected later on even though this is a bit silly
  # trt <- stats::rbinom(n, size = 1, prob = 0.5)

  C <- matrix(
    stats::rnorm(n * dim_c, mean = 1),
    nrow = n,
    ncol = dim_c
  )

  W <- C
  W[, 1] <- 1
  colnames(W) <- paste0("C", 1:dim_c)

  mod <- model_contagion(
    A_model = A_model,
    W = W,
    rho = rho,
    sigma = 1,
    subclass = "basic"
  )

  mod$id <- ids::uuid()
  mod
}
