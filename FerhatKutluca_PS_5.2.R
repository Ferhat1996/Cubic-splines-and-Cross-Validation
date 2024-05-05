#*******************************Exercise_2*************************************#

#******************************Exercise_2.i************************************#

# # Define a B-spline basis function based on Lecture 8's formulation
b_spline_base <- function(x, knots, r) {
  n <- length(knots)
  a <- knots[n] - knots[1]                              # Assuming equal spaced knots                         
  bs <- vector("numeric", length(x))
  
  for (j in 0:(r-1)) {
    bs <- bs + (-1)^j * choose(r-1, j) * (pmax(0, x - knots[j+1]))^(r-1)
  }
  
  bs / (a^(r-1) * factorial(r-1))
}

# Generate design matrix for B-splines
create_design_matrix <- function(X, knots, r) {
  sapply(knots, function(knot) b_spline_base(X, c(rep(min(X), r-1), knot, rep(max(X), r-1)), r))
}

# Estimate conditional expectation using custom B-splines and calculate leverage scores
estimate_conditional_expectation_with_leverage <- function(X, Y, num_knots, r) {
  knots <- quantile(X, probs = seq(0, 1, length.out = num_knots))
  B <- create_design_matrix(X, knots, r)
  fit <- lm(Y ~ B - 1)
  leverage_scores <- hatvalues(fit)
  return(list(fit = fit, leverage = leverage_scores))
}

#******************************Exercise_2.ii***********************************#

# Generate sample data
set.seed(123)
n <- 1000                                               # n can be changed to observe in different sample sizes                                         
X <- runif(n, 0, 1)                                     # X ∼ Unif[0, 1]
epsilon <- runif(n, -1, 1)                              # ε ∼ Unif[−1, 1]

# Define the true Conditional Expectation Function (CEF)
true_CEF <- function(X) { sin(2 * pi * X) }             # freely chosen m(·)
Y <- true_CEF(X) + epsilon                              # Y = m(X)+ε


# Summarize response and true CEF
summary(Y)
summary(true_CEF(X))

# Determine optimal number of subintervals for B-splines via cross-validation
cv_errors <- sapply(1:30, function(k) {
  result <- estimate_conditional_expectation_with_leverage(X, Y, k, 4) # r = 4 for cubic splines
  leverage_scores <- result$leverage
  residuals <- Y - predict(result$fit)
  sum((residuals / (1 - leverage_scores))^2)
})

optimal_k <- which.min(cv_errors)
print(paste("Optimal knot number: ", optimal_k))

# Fit model using optimal number of knots
fit_bspline <- estimate_conditional_expectation_with_leverage(X, Y, optimal_k, 4)$fit

# Predict using the fitted model
predicted_values <- predict(fit_bspline, newdata = data.frame(B = create_design_matrix(X, quantile(X, probs = seq(0, 1, length.out = optimal_k)), 4)))
summary(predicted_values)


# Fit a smoothing spline for comparison
fit_smoothing_spline <- smooth.spline(X, Y, df = 16)

# Set up the plot area
par(mfrow = c(1, 2))

# Plot 1: Display Fitted Values with Points
plot(X, Y, col = 'blue', pch = 16, xlab = 'X', ylab = 'Y', main = "Fitted Values with Points")
points(X, predicted_values, col = 'red', pch = 16)
legend("topright", legend = c("Real Data Points", "Cubic Spline Points"), col = c("blue", "red"), pch = 16)
mtext(side = 4, line = 1, paste("Knot number =", optimal_k, "and", "n =", n))

# Plot 2: Display lines of Smoothing Cubic Spline Fit vs True CEF
plot(NA, NA, xlim=range(X), ylim=range(c(Y, predict(fit_smoothing_spline, X))), 
     xlab = 'X', ylab = 'Y', main="Smoothing Cubic Spline Fİt vs True CEF")
lines(fit_smoothing_spline, col="red", lwd=2)
curve(true_CEF, add = TRUE, col = "blue", lwd = 2)
legend("topright", legend = c("Smoothing Cubic Spline Fit", "True CEF"), col = c("red", "blue"), lwd = 2)
mtext(side = 4, line = 1, paste("Knot number =", optimal_k, "and", "n =", n))


# Plot 3: Combine both plots into one
par(mfrow = c(1, 1))
plot(X, Y, col = 'blue', pch = 16, xlab = 'X', ylab = 'Y', main = "Cubic Spline Points vs Smoothing Cubic Spline Fit vs True CEF")
points(X, predicted_values, col = 'red', pch = 16)
lines(fit_smoothing_spline, col = "green", lwd = 2)
curve(true_CEF, add = TRUE, col = "black", lwd = 2)  # Add true CEF
legend("topright", legend = c("Cubic Spline Points", "Smoothing Cubic Spline Fit", "True CEF"), col = c("red", "green", "black"), pch = c(16, NA, NA), lwd = 2)
mtext(side = 4, line = 1, paste("Knot number =", optimal_k, "and", "n =", n))

