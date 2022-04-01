# STIP
STIP: State Transition Inference Prediction

### Depends
- R (>= 3.5.0)
- ComplexHeatmap 
- VGAM
- circlize
- dplyr
- grDevices

### Installation
```{r }
if (!require("devtools"))
  install.packages("devtools")
devtools::install_github("changxinw/STIP")
```

### Example
``` {r}
library(STIP)
data(tcell)
fit_data <- fitData(tcell)
gl = c("Rrp15", "Celf2", "Pcmtd1", "Plxdc2")
p <- PreprocessSTIP(fit_data, gl)
pdf("STIP_example_package.pdf", width = 7, height=10)
ComplexHeatmap::draw(p)
dev.off()
```