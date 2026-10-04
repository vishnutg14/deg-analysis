# Analysis of GSE116583 obtained from NCBI

This dataset was obtained from [this paper](https://doi.org/10.1165/rcmb.2017-0430TR) titled *A Beginner’s Guide to Analysis of RNA Sequencing Data*. More details of it can be found in the [references](#references) below.

This time I wanted to use just R and not touch Python. In case you want to follow, you could just copy the `renv.lock` file to your system and run these commands so that it installs required libraries (Keep your global environment clean)
```
install.packages("renv") #if you don't have renv installed
renv::restore()
```

# References
1. Koch CM, Chiu SF, Akbarpour M, Bharat A, Ridge KM, Bartom ET, Winter DR. A Beginner's Guide to Analysis of RNA Sequencing Data. Am J Respir Cell Mol Biol. 2018 Aug;59(2):145-157. doi: 10.1165/rcmb.2017-0430TR. PMID: 29624415; PMCID: PMC6096346. 
