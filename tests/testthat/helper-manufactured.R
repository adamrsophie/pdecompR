u_sin <- function(x, y) sin(pi * x) * sin(pi * y)
f_sin <- function(x, y) 2 * pi^2 * sin(pi * x) * sin(pi * y)
g_sin <- function(x, y) cbind(pi * cos(pi * x) * sin(pi * y), pi * sin(pi * x) * cos(pi * y))
