## Reproducibility Files for Papers Using `duqling`

This repository contains scripts and supporting materials to reproduce results from papers that use the `duqling` framework.

---

### All Emulators are Wrong, Many are Useful, and Some are More Useful Than Others

**Paper:**  
[All Emulators are Wrong, Many are Useful, and Some are More Useful Than Others](https://arxiv.org/pdf/2512.09060)

**Citation:**  
Rumsey, K.N., Gibson, G.C., Francom, D. and Morris, R., 2025. *All Emulators are Wrong, Many are Useful, and Some are More Useful Than Others: A Reproducible Comparison of Computer Model Surrogates.* arXiv preprint arXiv:2512.09060.

**Main script:**  
`R/khaos/khaos_paper_main.R`

**Required packages:**
```r
- library(duqling) # devtools::install_github("knrumsey/duqling")
- library(tidyverse)
- library(dbscan)
- library(ggrepel)
```

---

### Bayesian Adaptive Polynomial Chaos

**Paper:**  
[Bayesian Adaptive Polynomial Chaos Expansions](https://onlinelibrary.wiley.com/doi/pdf/10.1002/sta4.70151)

**Citation:**  
Rumsey, K.N., Francom, D., Gibson, G., Tucker, J.D., and Huerta, G. (2026).  
*Bayesian Adaptive Polynomial Chaos Expansions*. Stat, 15(1), e70151.

**Main script:**  
`R/khaos/khaos_paper_main.R`

**Required packages:**
```r
library(duqling)   # devtools::install_github("knrumsey/duqling")
library(khaos)     # devtools::install_github("knrumsey/khaos")
library(tidyverse)
library(BART)
library(laGP)
library(ggplot2)
library(tidyr)
library(dplyr)
library(stargazer)
```
