# ecostructure
### An R package for clustering and visualization of structure in ecological data

This package contains functions for fitting [STRUCTURE](http://www.genetics.org/content/155/2/945) type models to ecological data, both at local and global geographic scales, together with GIS based visualizations of the fitted models. These *grade-of-membership* models can be used to assess the local representation of large regional biotas, their degree of intermixing in local assemblages, and their rate of turnover across geographic space owing to environmental or climatic turnover. ecostructure makes use of advances in clustering alrogrithms, first from the package CountClust in 1.0 and now leveraging its successor fastTopics in 2.0. 

### Authors 
 
- [Alex White](https://sidatasciencelab.github.io/)
- [Kushal K Dey](https://www.deylab.net/)

This package grew out of a collaboration between Kushal Dey and Alex White while they were graduate students at the University of Chicago. Their advisors played a large role in developing these methods: 

- [Matthew Stephens](http://stephenslab.uchicago.edu/)
- [Trevor Price](https://pondside.uchicago.edu/ecol-evol/people/price.html)


### Citation

If you are using **ecostructure** or our code, please cite our papers:

White, Alexander E. and Dey, Kushal K. and Mohan, Dhananjai and Stephens, Matthew and Price, Trevor D. *Regional influences on community structure across the tropical-temperate divide*. Nature Communications. 2019. 10 (1). 2646. 10.1038/s41467-019-10253-6.

White, Alexander E. and Dey, Kushal K. and and Stephens, Matthew and Price, Trevor D. *Dispersal syndromes drive the formation of biogeographical regions, illustrated by the case of Wallace’s Line*. Global Ecology and Biogeography. 2021.  https://doi.org/10.1111/geb.13250


### Installation

Install **ecostructure** following the instructions below.

```R
install.packages(devtools)
install.packages("sf")
devtools::install_github("kkdey/methClust")
devtools::install_github("kkdey/maptpx"). # this is an updated version of CRAN package maptpx
devtools::install_github("kkdey/CountClust")
devtools::install_github("sidatsciencelab/ecostructure")
```
Then load **ecostructure**

```R
library(ecostructure)
```

### Demo

Some examples of produced using our **ecostructure** package

<img src="bin/ecostructure.2.001.jpeg" alt="misc" align = "middle">

If you want to try **ecostructure** and replicate figures like this, please check our tutorial [here](https://kkdey.github.io/ecostructure/).


### Questions?

For any queries or concerns related to the software, you can open an issue [here](https://github.com/kkdey/ecostructure/issues). Or you can contact 
us. Kushal Dey - *kshldey@gmail.com* or Alex White -
*aewhite100@gmail.com*.

We thank Prof. Trevor D. Price and Dhananjai Mohan for data collection for our
accompanying paper. We thank our mentors Prof. Matthew Stephens and 
Prof. Trevor Price for helpful suggestions and discussions. 

You are also welcome to contribute to **ecostructure** by forking this repo.









