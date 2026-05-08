new_contagion <- function(
  n,
  k,
  W,
  beta_w,
  beta_x,
  A_model,
  rho,
  sigma,
  ...,
  subclass = character()
) {
  rlang::check_dots_unnamed()

  model <- list(
    n = n,
    k = k,
    W = W,
    beta_w = beta_w,
    beta_x = beta_x,
    A_model = A_model,
    rho = rho,
    sigma = sigma,
    model_name = subclass,
    ...
  )

  class(model) <- c(subclass, "contagion")
  model
}

model_contagion <- function(
  A_model,
  W,
  ...,
  beta_w = NULL,
  beta_x = NULL,
  rho = 1,
  sigma = 1
) {
  # this class is a bit of hack. the idea is that, we when sample from it,
  # we actually sample from four different, but closely related models, all
  # at once. the modeling options come from a 2x2 grid:
  #
  #  (peer contagion vs latent contagion) x (X in regression, X not in regression)
  #
  # notation: coefficients in peer contagion models called beta, coefficients
  #   in latent contagion models called gamma. also, let W = [1 T C] be
  #   observed nodal covariates, A be the adjacency matrix, d_i be ith node
  #   degree, P be the expected value of A, Y be nodal outcomes, and
  #   varepsilon be i.i.d. errors. let D be a diagonal matrix of node degrees
  #
  # peer contagion models:
  #
  #   Y = W betaw + D^{-1} A Y betay + varepsilon
  #   Y = W betaw + X betax + D^{-1} A Y betay + varepsilon
  #
  # latent contagion models:
  #
  #   Y = W gammaw + E[D|X]^{-1} E[A|X] Y gammay + varepsilon
  #   Y = W gammaw + X gammax + E[D|X]^{-1} E[A|X] Y gammay + varepsilon
  #
  # all of these outcome models can be rearranged into "reduced form"
  #
  # peer contagion models:
  #
  #   Y = (I - betay D^{-1} A)^{-1} (W betaw + varepsilon)
  #   Y = (I - betay D^{-1} A)^{-1} (W betaw + X betax + varepsilon)
  #
  # latent contagion models:
  #
  #   Y = (I - gammay E[D|X]^{-1} E[A|X])^{-1} (W gammaw + varepsilon)
  #   Y = (I - gammay E[D|X]^{-1} E[A|X])^{-1} (W gammaw + X gammax + varepsilon)
  #
  # we will set betay = gammay = rho, and betax = gammax, and betaw = gammaw
  # and sigma will correspond to the variance of N(0, sigma^2) varepsilon

  # NOTE: this is *not* A_model$X. there is a notational clash here. we use X
  #   to mean the population ASE. fastRG uses X to mean the varimax X matrix.
  #   not the same!

  dim_w <- ncol(W)
  dim_x <- A_model$k

  if (is.null(beta_w)) {
    beta_w <- rep(5, dim_w) # stats::rnorm(dim_w, mean = 2, sd = 0.1)
  } else {
    stopifnot(length(beta_w) == dim_w)
  }

  if (is.null(beta_x)) {
    # going to drop one column of X for identifiability in the regression
    # only
    beta_x <- rep(2, dim_x) # stats::rnorm(dim_x , mean = 2, sd = 0.1)
    names_before_dropping_one <- paste0("Xhat", left_padded_sequence(1:dim_x))
    names(beta_x) <- names_before_dropping_one #e[-1]
  } else {
    stopifnot(length(beta_x) == dim_x)
  }

  names(beta_w) <- colnames(W)
  names(rho) <- "rho"

  new_contagion(
    n = nrow(W),
    k = A_model$k,
    W = W,
    beta_w = beta_w,
    beta_x = beta_x,
    A_model = A_model,
    rho = rho,
    sigma = sigma,
    ...
  )
}

sample_tidygraph_contagion <- function(model, ...) {
  # model$W should already have appropriate column names
  W_df <- tibble::as_tibble(as.matrix(model$W))

  graph <- sample_tidygraph(
    model$A_model
  ) |>
    tidygraph::arrange(as.numeric(name)) |>
    tidygraph::mutate(!!!W_df)

  varepsilon <- stats::rnorm(model$n, sd = model$sigma)
  X <- ASE(model$A_model)
  colnames(X) <- paste0("Xhat", left_padded_sequence(1:ncol(X)))

  # identifiability hack
  X_reg <- X #[, -1, drop = FALSE]

  # a little weird to do this here rather than in model_contagion(), but whatever
  # names(model$beta_x) <- colnames(X_reg)

  I <- Matrix::Diagonal(model$n, 1)

  #### peer models

  A <- igraph::as_adjacency_matrix(graph)
  deg <- igraph::degree(graph)
  DinvA <- Matrix::rowScale(A, 1 / deg)

  # a little weird to do this here rather than in model_contagion(), but whatever
  names(model$beta_x) <- colnames(X)

  # fingers crossed this doesn't start exploding. this will break if there
  # are identification issues. possible wrap in a tryCatch to alert if this
  # happens deep in a simulation loop
  y_peer_no_x <- solve(
    I - model$rho * DinvA,
    model$W %*% model$beta_w + varepsilon
  )
  y_peer_x <- solve(
    I - model$rho * DinvA,
    model$W %*% model$beta_w + X_reg %*% model$beta_x + varepsilon
  )

  ##### latent models

  # this will be memory intensive because these matrices are going to be dense
  # if this turns out to be prohibitive, we'll want to compute the inverse
  # of the pre-multiplication matrix more efficiently, probably by leverage
  # the svds(A_model) idea. that is, we can compute the pseudoinverse
  # from the svd of I - rho E[D]^{-1} P, and we can make a truncated svd
  # of that very fast

  P <- tcrossprod(X)
  ED <- rowSums(P)
  EDinvP <- Matrix::rowScale(P, 1 / ED)
  diag(EDinvP) <- 0

  # fingers crossed this doesn't start exploding. this will break if there
  # are identification issues. possible wrap in a tryCatch to alert if this
  # happens deep in a simulation loop
  y_latent_no_x <- solve(
    I - model$rho * EDinvP,
    model$W %*% model$beta_w + varepsilon
  )
  y_latent_x <- solve(
    I - model$rho * EDinvP,
    model$W %*% model$beta_w + X_reg %*% model$beta_x + varepsilon
  )

  # Add only the outcome variables to the graph
  # Auxiliary data for estimators (g_y_*, gtilde_y_*, ghat_y_*, Z_*, etc.)
  # will be computed in the estimation logic
  graph |>
    tidygraph::activate(nodes) |>
    tidygraph::mutate(
      y_peer_no_x = drop(y_peer_no_x),
      y_peer_x = drop(y_peer_x),
      y_latent_no_x = drop(y_latent_no_x),
      y_latent_x = drop(y_latent_x)
    )
}
