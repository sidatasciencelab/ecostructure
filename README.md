# ecostructure
ecostructure: grade-of-membership clustering and visualization for ecology
=============================================================
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip) 
[![CircleCI](https://circleci.com/gh/sidatasciencelab/ecostructure.svg?style=svg&circle-token=c025eb22e0d0a3d4efae2b44ee70bcd114f7fc59)](https://app.circleci.com/pipelines/github/sidatasciencelab/ecostructure)
[![DOI](https://zenodo.org/badge/45345.svg)](https://zenodo.org/badge/latestdoi/####)
[![Follow](https://img.shields.io/twitter/follow/SIDataScience.svg?style=social&style=flat&logo=twitter&label=Follow&color=blue)](https://twitter.com/SIDataScience)

---
### :warning: **Version 0.99.1 of ecostructure (github.com/kkdey/ecostructure) is deprecated** :warning:  
Please make sure you are using the code and docs maintained by the Smithsonian Data Science Lab (github.com/sidataciencelab/ecostructure). 

This package contains functions for fitting [STRUCTURE](http://www.genetics.org/content/155/2/945) type models to ecological data, both at local and global geographic scales, together with GIS based visualizations of the fitted models. These *grade-of-membership* models can be used to assess the local representation of large regional biotas, their degree of intermixing in local assemblages, and their rate of turnover across geographic space owing to environmental or climatic turnover. ecostructure makes use of advances in clustering algorithms, first from the package CountClust in 0.99.1 and now leveraging its successor fastTopics in 2.0. 

---

### Installation

Install **ecostructure** following the instructions below.

```R
remotes::install_github("stephenslab/fastTopics")
remotes::install_github("linxihui/NNLM")
remotes::install_github("sidatasciencelab/ecostructure")
```
**ecostructure** requires access to the "gfortran" library. Mac OS X
users may encounter the error "library not found for -lgfortran" when
installing. To fix this error, please follow the instructions at this
[link](https://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks--lgfortran-and--lquadmath-error/). Apple Silicon users should instead follow the directions found [here](https://github.com/bioinfoDZ/RISC/issues/10). 

Then load **ecostructure**

```R
library(ecostructure)
```

---

### Visualizations

Some examples of visualizations produced using our **ecostructure** package

<img src="bin/ecostructure.2.001.jpeg" alt="misc" align = "middle">

If you want to try **ecostructure** and replicate figures like this, please check our tutorial [here](https://sidatasciencelab.github.io/ecostructure/).

### Citation

If you are using **ecostructure** or our code, please cite our papers:

White, Alexander E. and Dey, Kushal K. and Mohan, Dhananjai and Stephens, Matthew and Price, Trevor D. *Regional influences on community structure across the tropical-temperate divide*. Nature Communications. 2019. 10 (1). 2646. [10.1038/s41467-019-10253-6](https://doi.org/10.1038/s41467-019-10253-6)

White, Alexander E. and Dey, Kushal K. and and Stephens, Matthew and Price, Trevor D. *Dispersal syndromes drive the formation of biogeographical regions, illustrated by the case of Wallace’s Line*. Global Ecology and Biogeography. 2021.  [10.1111/geb.13250](https://doi.org/10.1111/geb.13250)
### Authors 
 
- [Alex White](https://sidatasciencelab.github.io/people/leads/alex.html)
- [Kushal K Dey](https://www.deylab.net/)
- [Matthew Stephens](http://stephenslab.uchicago.edu/)
- [Trevor Price](https://pondside.uchicago.edu/ecol-evol/people/price.html)

---

### Questions?

For any queries or concerns related to the software, you can open an issue [here](https://github.com/sidatasciencelab/ecostructure/issues). 

The methods used to develop this framework have advanced since our initial release and may yet evolve in future iterations of this package. The Smithsonian Data Science Lab continues to collaborate with the Dey Lab (Sloan Kettering), the Price Lab (U Chicago) and the Stephens Lab (U Chicago) to develop this work.

You are welcome to contribute to **ecostructure** by forking this repo.









