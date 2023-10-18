---
title: "IMCthreshViewer"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{IMCthreshViewer}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

### Installation

To install the latest version from the github repository, use:


```r
if(!require(devtools)){
  install.packages("devtools") # If not already installed
}
devtools::install_github("LAChpx/IMCthreshViewer")
```


### Usage

After installing the package, use the following code to run the app:


```r
library(IMCthreshViewer)
shiny_IMCthreshViewer()
```
