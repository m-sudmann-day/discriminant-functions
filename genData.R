
# genData.R
#
# Matthew Sudmann-Day
# Barcelona GSE Data Science
#
# Generates a dataset that could be troublesome for discriminant functions
# and/or linear probability models to correctly identify binary classified
# data due to the presence of high-leverage outliers.  The dataset is
# created in such a way as to be easy for humans to visually identify a
# 100% correct boundary despite the model's difficulty in doing so.
#
# Uses R packages: mvtnorm, ggplot2
#
# Public functions: generate.dataset(), generate.troublesome.dataset()
#

library(mvtnorm)
library(ggplot2)

# .generatesubset() - a function to generate synthetic normally-distributed
#                     multivariate data with two variables (x, y)
#
# count - the number of observations to generate
# mean.x - the mean of the first variable
# mean.y - the mean of the second variable
# sd.x - the standard deviation of the first variable
# sd.y - the standard deviation of the second variable
# rho - the correlation between the two variables
# value - the binary classification (0/1) to which all rows should be associated
# label - the string classification to which all rows should be associated
#
# returns a data frame containing the generated data
.generate.subset <- function(count, mean.x, mean.y, sd.x, sd.y, rho, value, label)
{
  # Calculate the covariance and the variance/covariance matrix for the
  # two variables.
  cov <- rho * sd.x * sd.y
  vc <- matrix(c(sd.x ^ 2, cov, cov, sd.y ^2), 2, 2)
  
  # Use rmvnorm to generate multivariate random numbers.
  data <- rmvnorm(count, c(mean.x, mean.y), vc)
  
  # Associate the binary and string classifications with the generated numeric
  # values and return in a data frame.
  result <- data.frame(x=data[,1], y=data[,2], value=value, group=label)
  return(result)
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
# returns a two-row data.frame, prepared in a friendly manner for use
# in ggplot's geom_line function
.calculate.boundary <- function(data, x1, x2, label)
{
  # Perform a linear regression, explicitly defining an intercept column.
  fit <- lm(value ~ x + y + 1, data)
  
  # Extract weight  s and bias (intercept) from the results of the regression.
  weight.x <- coefficients(fit)['x']
  weight.y <- coefficients(fit)['y']
  bias <- coefficients(fit)[1]
  
  # Determine an intercept and slope that would be appropriate for display
  # on a plot chart that shows the two underlying variables on the two axes.
  intercept <- (0.5 - bias) / weight.y
  slope <- -1 * weight.x / weight.y
  
  # Generate two end-points for a boundary line.
  y1 = intercept + x1 * slope
  y2 = intercept + x2 * slope
  
  # Return the endpoints and the label of the line as a data frame.
  df <- data.frame(group=rep(label, 2))
  df$x <- c(x1, x2)
  df$y <- c(y1, y2)
  df$label <- label
  return(df)
}

# generate.dataset() - a function to generate a dataset containing binary
#                      classified data, the two groups being generated from
#                      two different distributions
#
# normal.count - the number of normal points
# extreme.count - the number of extreme points (if normal points far outnumber
#                                               extreme points, then the extreme
#                                               points are probably outliers)
# label1 - the label for the first group of points
# label2 - the label for the second group of points
# save.data - indicates that the function should write the generated data out
#             to a file (default:TRUE)
# save.plot - indicates that the function should write the generated data out
#             to a plot (default:TRUE)
#
# This function expects the caller to set a working directory.
#
# returns nothing
generate.dataset <- function(normal.count, extreme.count, label1, label2,
                             save.data=TRUE, save.plot=TRUE)
{
  # Determine the number of normal and extreme values to generate from each group
  # of the two classifications.
  normal.count.group1 <- round(normal.count / 2)
  normal.count.group2 <- normal.count - normal.count.group1
  extreme.count.group1 <- round(extreme.count / 2)
  extreme.count.group2 <- extreme.count - extreme.count.group1
  
  # Generate each subset of data separately.
  normals1 <- .generate.subset (normal.count.group1, 1, 12, 1.5, 1, 0.8, 0, label1)
  normals2 <- .generate.subset (normal.count.group2, 2, 7, 1.5, 1, 0.9, 1, label2)
  extremes1 <- .generate.subset (extreme.count.group1, 3, 35, 1, 1, -0.3, 0, label1)
  extremes2 <- .generate.subset (extreme.count.group2, -1, -15, 1, 1, -0.3, 1, label2)
  
  # Create data sets with and without the extreme points.
  data.without.extremes <- rbind(normals1, normals2)
  data.with.extremes <- rbind(data.without.extremes, extremes1, extremes2)

  # Determine the begin and end points on the x-dimension for calculating boundary lines.
  x.min <- min(data.with.extremes$x)
  x.max <- max(data.with.extremes$x)
  
  # Calculate boundary lines that the linear probability model puts between the two
  # classifications, first without extremes then with extremes.
  boundary1 <- .calculate.boundary(data.without.extremes, x.min, x.max, "Boundary Without Extremes")
  boundary2 <- .calculate.boundary(data.with.extremes, x.min, x.max, "Boundary With Extremes")
  
  # If the caller wants the data saved to file, do so.
  if (save.data)
  {
    write.csv(data.with.extremes, "dataset.csv", row.names=FALSE)
  }
  
  # If the caller wants a graphical plot saved to file, do so.
  if (save.plot)
  {
    plot <- ggplot(data.with.extremes, aes(x=x, y=y, color=group, group=group))
    plot <- plot + ggtitle("Classification of Binary Data by LPM")
    plot <- plot + labs(x="Var. X", y="Var. Y")
    plot <- plot + geom_point()
    plot <- plot + geom_line(data=boundary1)
    plot <- plot + geom_line(data=boundary2)

    ggsave(filename="dataPlot.pdf", width=8, height=6, units="in", dpi=100, plot=plot)
  }
}

# generate.troublesome.dataset() - a function that calls generate.dataset() in a fashion
#                                  that ensures that a linear probability model will do a
#                                  poor job of classifying the data
#
# save.data - indicates that the function should write the generated data out
#             to a file (default:TRUE)
# save.plot - indicates that the function should write the generated data out
#             to a plot (default:TRUE)
#
# This function expects the caller to set a working directory.
#
# returns nothing
generate.troublesome.dataset <- function(save.data=TRUE, save.plot=TRUE)
{
  set.seed(11111)

  generate.dataset(60, 2, "Group 1", "Group 2", save.data, save.plot)
}

# run generate.troublesome.dataset() to create plot and csv file.
generate.troublesome.dataset()
