## try http:// if https:// URLs are not supported
install.packages(c("caTools", "snow"),repos="http://cran.us.r-project.org")
#source("https://bioconductor.org/biocLite.R")
source("http://bioconductor.org/biocLite.R")
biocLite(ask=FALSE)
biocLite("Rsamtools",ask=FALSE)
