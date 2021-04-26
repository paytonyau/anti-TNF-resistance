require(plot3D)
A <- read.csv("MCP.csv")
x <- A$B.cells
y <- A$Endothelial.cells
z <- A$Neutrophils

# linear fit
fit <- lm(z ~ x + y)

# predict values on regular xy grid
grid.lines = 26
x.pred <- seq(min(x), max(x), length.out = grid.lines)
y.pred <- seq(min(y), max(y), length.out = grid.lines)
xy <- expand.grid( x = x.pred, y = y.pred)
z.pred <- matrix(predict(fit, newdata = xy), 
                 nrow = grid.lines, ncol = grid.lines)

# predicted z-values, fitted points for droplines to surface
fitpoints <- predict(fit) 

scatter3D(z = z, x = x , y = y, 
          pch = 16, cex = 0.5, 
          theta = 45, phi = -10, ticktype = "detailed",
          xlab = "Endothelial Cells", ylab = "B Cells", zlab = "Neutrophils", 
          clab = "Neutrophils", 
          surf = list(x = x.pred, y = y.pred, z = z.pred, facets = NA, fit = fitpoints),
          colkey = list(length = 0.8, width = 0.4),
          bty = "g",
          main = "xCELL")