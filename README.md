
<img src="vignettes/LogoCIM.png" align="right" alt="" width="150" />

IMCthreshViewer
===============

### IMCthreshViewer

Package to visualize the tiffs obtained in IMC but also to establish marker positivity thresholds from fcs obtained in post-segmentation.


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


### Help

<img src="vignettes/ImageAcc1.png" align="center" alt="" />
<img src="vignettes/ImageAcc2.png" align="center" alt="" />
