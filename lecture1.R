if(!"tinytex" %in% rownames(installed.packages())) install.packages("tinytex")

tinytex::install_tinytex()

install.packages("cli")

#Install devtools package if necessary
if(!"devtools" %in% rownames(installed.packages())) install.packages("devtools")

library("devtools")


# Install the stable development verions from GitHub
devtools::install_github("crsh/papaja")
library('papaja')
