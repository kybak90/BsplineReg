b_spline = function(x, knots, degree, i) {
  # Base case: 0th degree (piecewise constant)
  if (degree == 0) {
    return(ifelse(knots[i] <= x & x < knots[i + 1], 1, 0))
  }

  # Recursive case: degree > 0
  B_i_d1 = b_spline(x, knots, degree - 1, i)
  B_i1_d1 = b_spline(x, knots, degree - 1, i + 1)

  denom1 = knots[i + degree] - knots[i]
  denom2 = knots[i + degree + 1] - knots[i + 1]

  term1 = if (denom1 == 0) 0 else
    ((x - knots[i]) / denom1) * B_i_d1
  term2 = if (denom2 == 0) 0 else
    ((knots[i + degree + 1] - x) / denom2) * B_i1_d1

  return(term1 + term2)
}

create_design_matrix = function(x_values, knots, degree)
{
  n = length(x_values)  # Number of data points
  num_basis = length(knots) - degree - 1  # Number of basis functions
  design_matrix = matrix(0, nrow = n, ncol = num_basis)  # Initialize matrix

  for (j in 1:num_basis)
    for (i in 1:n)
      design_matrix[i, j] = b_spline(x_values[i], knots, degree, j)

  return(design_matrix)
}

add_boundary_knots = function(x, interior_knots, degree = 3, tiny = 1e-5)
{
  knots = c(rep(min(x) - tiny, degree + 1), interior_knots, rep(max(x) + tiny, degree + 1))
  return(knots)
}

#' Generating interior knots sequence for given data points
#'
#' This function computes a knot sequence based on the quantiles of the data points
#'
#' @param x Numeric vector representing data points
#' @param dimension Numeric value indicating the number of basis functions
#' @param degree Numeric value representing the spline degree
#'
#' @return Numeric vector with interior knots sequence
#'
#' @export
#'
#' @examples
#' x = runif(50, 0, 1)
#' knots = knots_quantile(x, 5, 3)
#' print(knots)
knots_quantile = function(x, dimension, degree = 3)
{
  dimension = max(dimension, degree + 1)
  number_interior_knots = dimension - degree - 1
  if (number_interior_knots > 0)
    probs = (1 : number_interior_knots) / (number_interior_knots + 1)
  else
    probs = NULL

  interior_knots = quantile(x, probs, type = 1)
  return(interior_knots)
}

#' Fitting regression spline estimator
#'
#' This function computes the coefficients spline regression estimator.
#'
#' @param x_values,y_values Numeric vectors
#' @param interior_knots Numeric vector representing interior knots sequence. It can be obtained from \code{\link{knots_quantile}} function based on the given x points.
#' @param degree Numeric value indicating the spline degree
#'
#' @return List containing the spline coefficients estimates, knots sequence, spline degree
#'
#' @examples
#' set.seed(923)
#' n = 30
#' x_values = sort(runif(n, 0, 1))
#' y_values = sin(2 * pi * x_values) + cos(4 * pi * x_values) + rnorm(n, sd = 0.2)
#' #spline degree specification
#' degree = 3
#' # knot generation
#' num_interior_knots = 5
#' interior_knots = knots_quantile(x_values, num_interior_knots)
#' # model fitting
#' model = fit_spline(x_values, y_values, interior_knots, degree)
#' print(model)
#'
#' @export
fit_spline = function(x_values, y_values, interior_knots, degree)
{
  knots = add_boundary_knots(x_values, interior_knots, degree)
  G = create_design_matrix(x_values, knots, degree)
  beta = solve(t(G) %*% G) %*% t(G) %*% y_values
  return(list(beta = beta, knots = knots, degree = degree))
}

#' Predicting the values of y for given x based on the spline estimator
#'
#' This function computes the predicted values of y for given x based on the spline estimator.
#'
#' @param model List object obtained from the \code{\link{fit_spline}} function
#'
#' @param new_x Numeric vector representing a grid of evaluation points
#'
#' @return Numeric vector with predicted values at new_x
#'
#' @export
predict_spline = function(model, new_x)
{
  G_new = create_design_matrix(new_x, model$knots, model$degree)
  y_pred = G_new %*% model$beta
  return(y_pred)
}

#' Plotting spline estimator with scatter plots
#'
#' This function provides a plot of the spline estimator with data points
#'
#' @param x_values,y_values Numeric vector
#' @param model List object obtained from the \code{\link{fit_spline}} function
#' @param grid_x Numeric vector with a grid of evaluation points
#'
#' @export
#' @import ggplot2
plot_spline = function(x_values, y_values, model, grid_x)
{
  y_pred = predict_spline(model, grid_x)
  data_plot = data.frame(x = x_values, y = y_values)
  spline_plot = data.frame(x = grid_x, y = y_pred)

  ggplot() +
    geom_point(data = data_plot, aes(x, y), color = "black") +
    geom_line(data = spline_plot, aes(x, y), color = "blue", linewidth = 1.2) +
    labs(title = "Fitted B-spline Regression", x = "x", y = "y") +
    xlim(c(min(x_values), max(x_values))) +
    theme_minimal()
}
