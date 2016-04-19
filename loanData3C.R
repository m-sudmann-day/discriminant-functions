
# genData.R
#
# Matthew Sudmann-Day
# Barcelona GSE Data Science
#
# Generates a dataset containing three classifications of loan applications:
# Approved, Denied, and Undecided.  It separates these classifications with
# two separate discriminant functions and generates a graph showing the data
# points and boundary lines.  The first boundary line splits denied loans
# from all others because, presumably, it is better to deny a loan that
# should be approved than it is to approve a loan that should be denied.
# Second, the remaining points are split with another boundary that separates
# approved from undecided loans.  Note that this boundary line does not
# extend beyond the first boundary line so as to avoid creating ambiguous
# regions.
#
# Uses R packages: mvtnorm, ggplot2
#
# Public function: run.analysis()
#

library(mvtnorm)
library(ggplot2)

# .gen.data() - a function for creating synthetic data for a bivariate
#               normal distribution with given label, count, mean, standard
#               deviation, and correlation
#
# label - the label or classification to associate with the data to be created
# count - the number of data points to create
# mu - the mean of the data as a vector with two elements
# sd - the standard deviation of the data as a vector with two elements
# rho - the correlation between the variables
.gen.data <- function(label, count, mu, sd, rho)
{
  cov <- rho * sd[1] * sd[2]
  sigma <- matrix(c(sd[1]^2, cov, cov, sd[2]^2), 2, 2)
  data <- as.data.frame(rmvnorm(n=count, mean=mu, sigma=sigma))
  data$group = label

  colnames(data) <- c("x", "y", "group")
  return(data)
}

# .calculate.boundary() - a function to calculate the end coordinates of a
#                         boundary line
#
# data - the data frame containing the binary classified observations
# x1 - the x-coordinate of one end of the boundary line
# x2 - the x-coordinate of the other end of the boundary line
# label - the label to associate with the boundary
#
# The x-coordinates must be passed in because the data parameter may
# not contain all points that the user sees, thereby creating too short
# of a boundary line for cosmetic purposes.
#
# returns a list containing:
#   intercept - the intercept of the boundary line on the y-axis
#   slope - the slope of the boundary line
#   df - a two-row data.frame prepared in a friendly manner for use
#        in ggplot's geom_line function
#   fit - a vector of fitted values (0 or 1) based on the linear model used
.calculate.boundary <- function(training.data, data.to.fit, group)
{
  # Perform a linear regression, explicitly defining an intercept column.
  model <- lm(target ~ x + y + 1, training.data)
  
  # Extract weight  s and bias (intercept) from the results of the regression.
  weight.x <- coefficients(model)['x']
  weight.y <- coefficients(model)['y']
  bias <- coefficients(model)[1]
  
  # Determine an intercept and slope that would be appropriate for display
  # on a plot chart that shows the two underlying variables on the two axes.
  intercept <- (0.5 - bias) / weight.y
  slope <- -1 * weight.x / weight.y
  
  # Generate two end-points for a boundary line.
  x1 <- min(data.to.fit$x)
  x2 <- max(data.to.fit$x) 
  y1 = intercept + x1 * slope
  y2 = intercept + x2 * slope
  
  # Return the endpoints and the label of the line as a data frame.
  df <- data.frame(group = rep(group, 2))
  df$x <- c(x1, x2)
  df$y <- c(y1, y2)

  # Use our model to predict the classifications of points within the data
  # we want to fit.  
  fit <- predict(model, newdata=data.to.fit)

  return(list(intercept=intercept, slope=slope, df=df, fit=fit))
}

# .calculate.intersection.of.two.boundary.lines() - a function to
# calculate the intersection of two boundary lines as produced by
# the function .calculate.boundary().
#
# line1 - the first boundary line
# line2 - the second boundary line

# returns a list representing the intersection's coordinates (x, y)
.calculate.intersection.of.two.boundary.lines <- function(line1, line2)
{
  m1 <- line1$slope
  m2 <- line2$slope
  b1 <- line1$intercept
  b2 <- line2$intercept
  x <- (b2-b1)/(m1-m2)
  y <- b1 + m1*x
  
  return(list(x=x, y=y))
}

# run.analysis() - a function to analyze classified loan data, produce
#                  a linear model from, and then apply decisions based on
#                  that linear model.
#
# Create classification data for three categories of loan data:
#   Denied, Approved, Undecided.
#
# Generate a boundary line between Denied loans and all others.
# Generate a boundary line between Approved and Undecided loans, but not
# extended into the Denied region to avoid producing an ambiguous region.
run.analysis <- function()
{
  #Create the data points.
  data <- rbind(
    .gen.data("Denied", 50, c(6,150), c(3,20), 0.5),
    .gen.data("Approved", 50, c(10,100), c(2,30), 0.7),
    .gen.data("Undecided", 50, c(5,80), c(2,25), 0.8))
    
  # Set a target of 1 for all denied loans, then calculate the boundary line, and
  # store the fitted results back into the data set.
  data$target <- ifelse(data$group == "Denied", 1, 0)
  denial.boundary <- .calculate.boundary(data, data, "1 Denied vs. All Others")
  data$deny.fit <- denial.boundary$fit
  data$not.deny.fit <- 1 - denial.boundary$fit

  # Make a copy of the data, including only the loans that were not denied.  
  not.denied <- data[data$group != "Denied",]
  
  # Set a target of 1 for all approved loans, then calculate the boundary line, and
  # store the fitted results back into the data set.
  not.denied$target <- ifelse(not.denied$group == "Approved", 1, 0)
  approval.boundary <- .calculate.boundary(not.denied, data, "2 Approved vs. Undecided")
  data$approve.fit <- approval.boundary$fit
  data$undecided.fit <- 1 - approval.boundary$fit
  
  # Set the overall decision based on the two binary fit results.
  data$decision <- ifelse(data$deny.fit > 0.5, "Denied", ifelse(data$approve.fit > 0.5, "Approved", "Undecided"))
  data$correctness <- ifelse(data$group == data$decision, "Correct", "Incorrect")

  # Calculate the intersection of the two boundary lines.
  intersection.of.boundaries <- .calculate.intersection.of.two.boundary.lines(denial.boundary, approval.boundary)
  
  # We do not want the approval boundary line to cross the denial boundary line, so determine
  # whether we will replace the min or the max from the approval boundary line with the
  # intersection point between the two lines.
  if (approval.boundary$slope > 0) {
    replace.row <- which.max(approval.boundary$df$x)
  } else {
    replace.row <- which.min(approval.boundary$df$x)
  }
  
  # Set the intersection point as the end of the approval boundary line.
  approval.boundary$df$x[replace.row] <- intersection.of.boundaries$x
  approval.boundary$df$y[replace.row] <- intersection.of.boundaries$y
  
  # Avoid a ggplot error by giving a fit category to the lines even though it is unused.
  approval.boundary$df$decision <- "Undecided"
  denial.boundary$df$decision <- "Undecided"
  approval.boundary$df$Correct <- "Correct"
  denial.boundary$df$Correct <- "Correct"
  
  # Eliminate a temporary column that we do not want to return.
  data$target <- NULL

  return(list(data=data, denial.boundary=denial.boundary$df, approval.boundary=approval.boundary$df))
}

# Set a random seed.
set.seed(1221)

# Run the analysis.
analysis <- run.analysis()

# Extract the final result including original status, deny fit, approve fit, and decision.
result.df <- analysis$data

# Write the data frame to a file.
write.csv(result.df, "predictions.csv", row.names=FALSE)

# Render and save the data as a scatter plot with boundary lines.
plot <- ggplot(data=analysis$data, aes(x=x, y=y, colour=group, shape=decision))
plot <- plot + geom_point(aes(alpha=correctness))
plot <- plot + scale_alpha_manual(name="Correctness (Opacity)", values=c("Correct"=0.5, "Incorrect"=1))
plot <- plot + scale_shape_discrete(name="Decision (Shape)")
plot <- plot + scale_colour_discrete(name="Training (Color)")
plot <- plot + geom_line(data=analysis$denial.boundary)
plot <- plot + geom_line(data=analysis$approval.boundary)
plot <- plot + xlab("Solvency") + ylab("PI Ratio")
plot <- plot + ggtitle("Loan Applications\n(Training Classification vs. Model Decision)")
plot <- plot + coord_cartesian(xlim=c(min(analysis$data$x), max(analysis$data$x)), ylim=c(min(analysis$data$y), max(analysis$data$y)))
ggsave(filename="discFunction3C.pdf", width=8, height=6, units="in", dpi=100, plot=plot)
plot
