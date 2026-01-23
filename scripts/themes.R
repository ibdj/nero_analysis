# install the pseudo-package from this Github repository

library(devtools)
devtools::install_github("max-alletsee/rstudio-themes")

library(rstudiothemes) # ... then load the library

# example 1: bulk-install all light themes
install_rstudio_themes(theme = "all_dark")

# example 2: install two specific light themes
install_rstudio_themes(theme = c("Ayu Light", "Github {rsthemes}"))
