# ChIP-seq processing pipeline
An R package for anlaysis of ChIP-seq and other functional sequencing data.
Please see [package homepage](http://compbio.med.harvard.edu/Supplements/ChIP-seq/) for details.

## Requirements
A unix-flavored OS with R and [Boost C++](http://www.boost.org/) installed.

## Installation
You can use modtools to install SPP:
```
require(devtools)
devtools::install_github('hms-dbmi/spp', build_vignettes = FALSE)
``` 
Alternatively, download a .tar.gz containing the [latest release](https://github.com/hms-dbmi/spp/releases) and use the standard R installation command, e.g.:
```
R CMD INSTALL spp_1.13.tar.gz
```

Note: if Boost libraries are installed in non-standard location, please specify the location in a BOOST_ROOT environment variable prior to running the installation.
```
export BOOST_ROOT=/path/boost_1_58_0/
```
