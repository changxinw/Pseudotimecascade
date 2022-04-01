# STIP
STIP: State Transition Inference Prediction

### Depends
- R (>= 3.5.0)
- ComplexHeatmap 
- VGAM
- circlize
- dplyr
- ggplot2
- clusterProfiler

### Installation

```{r }
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", "clusterProfiler")

BiocManager::install(c("ComplexHeatmap", ""))

install.packages(c("VGAM", "circlize", "dplyr", "ggplot2"))

if (!require("devtools"))
  install.packages("devtools")
devtools::install_github("changxinw/STIP")
```

### Example
Please refer to example folder for a Seurat example.