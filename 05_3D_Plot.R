## Load the plot3D library
library("plot3D")

## Import the matrix
A <- read.csv("MCP.csv")

## Isolate the factors from the matrix
x <- A$B.cells
y <- A$Endothelial.cells
z <- A$Neutrophils

# Fit a linear model to the data
fit <- lm(z ~ x + y)

# Predict values on a regular xy grid
grid.lines <- 26
x.pred <- seq(min(x), max(x), length.out = grid.lines)
y.pred <- seq(min(y), max(y), length.out = grid.lines)
xy <- expand.grid(x = x.pred, y = y.pred)
z.pred <- matrix(predict(fit, newdata = xy), 
                 nrow = grid.lines, ncol = grid.lines)

# Predicted z-values, fitted points for drop lines to the surface
fitpoints <- predict(fit) 

# Plot the result
scatter3D(x = x, y = y, z = z, 
          pch = 16, cex = 0.5, theta = 45, phi = -10, ticktype = "detailed",
          xlab = "Endothelial Cells", ylab = "B Cells", zlab = "Neutrophils", 
          clab = "Neutrophils", 
          surf = list(x = x.pred, y = y.pred, z = z.pred, facets = NA, fit = fitpoints),
          colkey = list(length = 0.8, width = 0.4),
          bty = "g",
          main = "MCP-Counter")