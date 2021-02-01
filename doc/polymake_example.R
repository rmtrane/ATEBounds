## -----------------------------------------------------------------------------
library(ACEBounds, quietly = TRUE)

vertices <- create_vertices(
  n_z_levels = 2, 
  data_format = "bivariate",
  x_monotone = FALSE,
  y_monotone = FALSE,
  include_alpha = TRUE
)

vertices

## ----include = FALSE----------------------------------------------------------
file.remove("example")

## -----------------------------------------------------------------------------
polymake_script <- create_script(vertices = vertices, 
                                 output_file_name = "example",
                                 polymake_script_name = "polymake_script",
                                 overwrite = TRUE,
                                 output_folder = ".")

polymake_script

## -----------------------------------------------------------------------------
system(paste("/Applications/polymake.app/Contents/MacOS/polymake.run --script", polymake_script$polymake_script))

## -----------------------------------------------------------------------------
from_polymake <- read_polymake_results(
  output_file = polymake_script$output_file, 
  data_format = "bivariate"
)

print(from_polymake, n = Inf)

## ----eval = FALSE-------------------------------------------------------------
#  ## Data type
#  for (data_format in c("bivariate", "trivariate")){
#    ## Number of levels for IV
#    for (n_z_levels in 2:4){
#      ## Monotonicity of X?
#      for (x_mono in c(FALSE, TRUE)){
#        ## Monotonicity of Y?
#        for (y_mono in c(FALSE, TRUE)){
#          ## Create vertices
#          vertices <- create_vertices(n_z_levels, x_monotone = x_mono, y_monotone = y_mono,
#                                      data_format = data_format, include_alpha = TRUE)
#  
#          ## Create script
#          script <- create_script(vertices,
#                                  output_file_name = paste("n_z_levels", n_z_levels,
#                                                           "x_monotone", x_mono,
#                                                           "y_monotone", y_mono,
#                                                           sep = "-"),
#                                  overwrite = FALSE,
#                                  output_folder = here::here(paste(data_format, "bound_matrices", sep = "_")))
#  
#          ## Run polymake using script created if script is new (i.e. error = 0)
#          if(script$error < 1){
#            system(paste("/Applications/polymake.app/Contents/MacOS/polymake.run --script", script$polymake_script))
#            system(paste("rm", script$polymake_script))
#          }
#        }
#      }
#    }
#  }

