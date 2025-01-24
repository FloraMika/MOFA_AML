#### #### Requierements
#
# S. Worthington, ipak.R, (2012), GitHub repository, https://gist.github.com/stevenworthington/3178163

# ipak function: install and load multiple R packages.
# check to see if packages are installed. Install them if they are not, then load them into the R session.

#!/usr/bin/env Rscript

ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg))
        install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
}

# usage
packages <- c("sva","dplyr","pca3d", "affy", "randomForest", "RColorBrewer", "Rtsne")
ipak(packages)

sessionInfo()
