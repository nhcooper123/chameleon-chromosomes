# chameleon-chromosomes
# Microchromosome fusions underpin convergent evolution of chameleon karyotypes

Author(s): [Natalie Cooper](mailto:natalie.cooper.@nhm.ac.uk)

This repository contains all the code and some data used in the [paper](TO ADD link). 

To cite the paper: 
> Marcello Mezzasalma, Jeffrey W. Streicher, Fabio M. Guarino, Marc E.H. Jones, Simon P. Loader, Gaetano Odierna, and Natalie Cooper 2023. Microchromosome fusions underpin convergent evolution of chameleon karyotypes. Evolution. In press.

To cite this repo: 
> Natalie Cooper. 2023. GitHub: nhcooper123/chameleon-chromosomes: code for the paper. Zenodo. DOI: TO ADD.

![alt text](https://github.com/nhcooper123/chameleon-chromosomes/raw/master/manuscript/figures/ChromoSSE_plot.png)

## Data
All raw and cleaned data are available from the [NHM Data Portal](https://doi.org/10.5519/rmaopg7d).

For reproducibility purposes, if you want to rerun the data wrangling steps in script `analyses\01-data-wrangling`, download `CHROMO_STATES_2020-10-16.csv` and place it into a `raw-data/` folder along with the dataset of [Meiri 2018](https://onlinelibrary.wiley.com/doi/full/10.1111/geb.12773).

The cleaned datasets, along rest of the ancillary data required to run all analyses and produce all figures is within the `data` folder.

* *Bayesian-tree-species.tre* - tree needed for the main analyses in the paper
* *chromosome-data-species.csv* - data needed for the main analyses in the paper
* *Bayesian-tree.tre* - tree for the all OTUs analyses in the supplementary materials
* *chromosome-data.csv* - data needed for the all OTUs analyses in the supplementary materials
* *chromevol-simulations-outputs-best.csv* - outputs from the "best" model of chromosome evolution used in the main text. The equivalent for the all OTUs analyses is in the folder `data/OTUs`.
* *null-simulations-outputs.csv* - outputs from the null simulations using the best model of chromosome evolution used in the main text. The equivalent for the all OTUs analyses is in the folder `data/OTUs`.

If you use the cleaned data please cite as follows: 
> Marcello Mezzasalma, Jeffrey W. Streicher, Fabio M. Guarino, Marc E.H. Jones, Simon P. Loader, Gaetano Odierna, and Natalie Cooper 2022. CHROMREP [Data set]. Natural History Museum. https://doi.org/10.5519/rmaopg7d.

The life-history data all comes from Meiri 2018 so please cite this paper if you use it.

> Meiri, S. Traits of lizards of the world: Variation around a successful evolutionary design. Global Ecol Biogeogr. 2018; 27: 1168– 1172

## ChromEvol analyses (`chromEvol/` folder)

The majority of analyses run in R, but we also used ChromEvol, a stand alone package available from http://chromevol.tau.ac.il/. The associated paper for the package is:

> Lior Glick, Itay Mayrose, ChromEvol: Assessing the Pattern of Chromosome Number Evolution and the Inference of Polyploidy along a Phylogeny, Molecular Biology and Evolution, Volume 31, Issue 7, July 2014, Pages 1914–1922, https://doi.org/10.1093/molbev/msu122   

To run these analyses you will need to download ChromEvol and save it in the `chromEvol/` folder. This unzips into a folder called `chromEvol_source/` that contains the package as an executable called `chromEvol`. I have not uploaded this to the repo because it's not my package. 

ChromEvol analyses require a dataset, a tree and a control file. These are within the `chromEvol/` folder preceded by `data_for_chromevol` or `tree_for_chroevol` or with the word `control` somewhere in their name respectively. Several versions exist for the different sensitivity analyses. These analyses save their outputs into the `chromEvol/results` folder, in folders named based on the names of the control files. Full details of how to build control files is in the ChromEvol manual at http://chromevol.tau.ac.il/. We provide code for creating appropriate datasets and trees in `analyses/02-`

Running chromEvol analyses is simple. Just open a Terminal in the `chromEvol` folder and run:

`chromEvol_source/chromEvol name-of-control-file`

You can run all analyses from our paper using the shell script: `run_all_models.sh`. 

The `chromEvol/` folder also includes a folder called `root_freq`. This has the three different root frequency definitions used in the different ChromEvol analyses (see text).

Finally the `OTUs` folder within `chromEvol/` contain the same things but for analyses using all OTUs (found in the supplementary materials of paper).

## ChromoSSE analyses (`ChromoSSE/` folder)

The majority of analyses run in R, but we also used [ChromoSSE](https://revbayes.github.io/tutorials/chromo/), which runs within [RevBayes](https://revbayes.github.io/). The associated paper for ChromoSSE is:

> Freyman W.A., Höhna S. 2018. Cladogenetic and anagenetic models of chromosome number evolution: a Bayesian model averaging approach. Systematic Biology. 67:1995–215.

And for RevBayes:

> Höhna, Landis, Heath, Boussau, Lartillot, Moore, Huelsenbeck, Ronquist. 2016. RevBayes: Bayesian phylogenetic inference using graphical models and an interactive model-specification language. Systematic Biology, 65:726-736.

To run these analyses you will need to download RevBayes. This unzips into a folder called `rb` that contains the package as an executable called `chromEvol`. I have not uploaded this to the repo because it's not my package. 

ChromoSSE analyses require a dataset, a tree and a control file. These are within the `chromEvol/` folder preceded by `data_for_chromevol` or `tree_for_chroevol` or with the word `control` somewhere in their name respectively. Several versions exist for the different sensitivity analyses. These analyses save their outputs into the `chromEvol/results` folder, in folders named based on the names of the control files. Full details of how to build control files is in the ChromEvol manual at http://chromevol.tau.ac.il/. We provide code for creating appropriate datasets and trees in `analyses/02-`

To run the ChromoSSE analyses, open `rb` and a Terminal window should open. Next set the working directory so RevBayes knows where to look for the data i.e. `setwd("YOUR_FOLDER/chameleon-chromosomes/ChromoSSE")`. Decide on the number of iterations and then run the code to set this: `iterations = 10000`. Then run your RevBayes control file (files ending in `.Rev`) using source, e.g.

`source("ChromoSSE_exclude_n18.Rev")`

The `ChromoSSE/` folder includes RevBayes control files to replicate all models in the chromEvol model set.

**These analyses take a long time to run! One model is approximately two weeks of run time.**

-------
## R analyses
The analysis code is divided into `.Rmd` files that run the analyses and plot the figures for each section of the paper/supplementary materials. 

Note that several of these scripts use an R package written by N.Cusimano. You can download this here https://www.en.sysbot.bio.lmu.de/people/employees/cusimano/use_r/. The easiest way to use this is to download the package and save the `ChromEvol.R` file. You can then install these functions in R using `source`. In my repo I kept this within the `chromEvol` folder so used the code `source("chromEvol/ChromEvol.R")` to install the functions.

1. **00-fix-tip-labels-bayesian-tree.R** changes the tip labels of the Bayesian tree to match the dataset.
1. **01-data-wrangling.Rmd** wrangles the data into a clean format for analyses.
2. **02A-data-for-chromosome-evolution-models_main.Rmd** creates data and tree files for the ChromEvol analyses.
3. **02B-data-for-chromosome-evolution-models_OTUs.Rmd** as above but for all OTUs analyses.
4. **03A-visualising-data-on-phylogenies_main.Rmd** visualises chromosome properties onto a phylogeny.
5. **03B-visualising-data-on-phylogenies_OTUs.Rmd** as above but for all OTUs analyses.
6. **04A-chromosome-evolution-models.Rmd** extract the outputs from the ChromEvol models, but only for the two main models of interest (constant rates and linear rates).
7. **04A-supp-chromosome-evolution-models-all_main.Rmd** extract the outputs from the ChromEvol models for all models.
8. **04B-chromosome-evolution-models_OTUs.Rmd** as 4A above but for all OTUs analyses.
9. **04B-supp-chromosome-evolution-models-all_OTUs.Rmd** as 4A-supp above but for all OTUs analyses.
10. **05A-micro-macro-correlates-analyses_main.Rmd** GLMs of haploid numbers of chromosomes against number of micro- and macro-chromosomes, numbers of micro- vs numbers of micro-chromosomes, ITS versus haploid numbers of chromosomes, and number of micro- and macro-chromosomes. Plus figures.
11. **05B-micro-macro-correlates-analyses_OTUs.Rmd** as above but for all OTUs analyses.
12. **06A-run-or-collate-simulations_main.Rmd** collates simulations outputs from ChromEvol and runs null simulations of haploid chromosome numbers.
13. **06B-run-or-collate-simulations_OTUs.Rmd** as above but for all OTUs analyses.
14. **06C-figures-simulation-models_main.Rmd** creates figures using simulations.
15. **06D-figures-simulation-models_OTUs.Rmd** as above but for all OTUs analyses.
16. **07A-taxon-species-pairs_main.Rmd** creates plots of species chromosome number distance versus phylogenetic distance.
17. **07B-taxon-species-pairs_OTUs.Rmd** as above but for all OTUs analyses.
18. **08A-ecological-correlates-figures.Rmd** creates figures of haploid chromosome number and a range of ecological, life-history and geographical variable for each species. FYI there are no all OTUs equivalents of these analyses as data are only available at the species level.
19. **08B-ecological-correlates-analyses.Rmd** MCMCglmm models of haploid chromosome number and a range of ecological, life-history and geographical variable for each species. FYI there are no all OTUs equivalents of these analyses as data are only available at the species level.
20. **09A-morphology-chromosome-no_main.Rmd** correlations of different aspects of chromosome morphology with haploid number of chromosomes.
21. **09B-morphology-chromosome-no_OTUs.Rmd** as above but for all OTUs analyses.
22. **10A-data-for-chromoSSE-models_main.Rmd** Extracts data for the ChromoSSE models in RevBayes.
23. **11A-ChromoSSE-assess-convergence-get-parameters_main.Rmd**. Checks the outputs of the ChromoSSE models for convergence (using traceplots and ESS).
24. **12A-ChromoSEE-plotting.R**. Plots results from ChromoSSE models to create Figure 7 in the main text.
25. **13A-ChromoSEE-comparisons.R**. Plots correlated evolution plots to create supplemental figure.

-------
## Other folders

* `/outputs` contains the figures and tables
* `/manuscript` contains the manuscript materials in LaTeX format

-------
## Session Info
For reproducibility purposes, here is the output of `devtools::session_info()` used to perform the analyses in the publication.

    ─ Session info ─────────────────────────────────────────────────────────────────────────────
    setting  value
    version  R version 4.2.0 (2022-04-22)
    os       macOS Big Sur 11.4
    system   x86_64, darwin17.0
    ui       RStudio
    language (EN)
    collate  en_US.UTF-8
    ctype    en_US.UTF-8
    tz       Europe/Dublin
    date     2023-04-19
    rstudio  2023.03.0+386 Cherry Blossom (desktop)
    pandoc   2.19.2 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/ (via rmarkdown)

   ─ Packages ─────────────────────────────────────────────────────────────────────────────────
   package           * version    date (UTC) lib source
   ape               * 5.7-1      2023-03-13 [1] CRAN (R 4.2.0)
   aplot               0.1.9      2022-11-24 [1] CRAN (R 4.2.0)
   backports           1.4.1      2021-12-13 [1] CRAN (R 4.2.0)
   bit                 4.0.5      2022-11-15 [1] CRAN (R 4.2.0)
   bit64               4.0.5      2020-08-30 [1] CRAN (R 4.2.0)
   broom             * 1.0.3      2023-01-25 [1] CRAN (R 4.2.0)
   cachem              1.0.6      2021-08-19 [1] CRAN (R 4.2.0)
   callr               3.7.3      2022-11-02 [1] CRAN (R 4.2.0)
   cli                 3.6.0      2023-01-09 [1] CRAN (R 4.2.0)
   clusterGeneration   1.3.7      2020-12-15 [1] CRAN (R 4.2.0)
   coda              * 0.19-4     2020-09-30 [1] CRAN (R 4.2.0)
   codetools           0.2-18     2020-11-04 [1] CRAN (R 4.2.0)
   colorspace          2.1-0      2023-01-23 [1] CRAN (R 4.2.0)
   combinat            0.0-8      2012-10-29 [1] CRAN (R 4.2.0)
   corpcor             1.6.10     2021-09-16 [1] CRAN (R 4.2.0)
   crayon              1.5.2      2022-09-29 [1] CRAN (R 4.2.0)
   cubature            2.0.4.4    2022-03-22 [1] CRAN (R 4.2.0)
   DBI                 1.1.3      2022-06-18 [1] CRAN (R 4.2.0)
   devtools            2.4.5      2022-10-11 [1] CRAN (R 4.2.0)
   digest              0.6.31     2022-12-11 [1] CRAN (R 4.2.0)
   doParallel          1.0.17     2022-02-07 [1] CRAN (R 4.2.0)
   dplyr             * 1.1.0      2023-01-29 [1] CRAN (R 4.2.0)
   ellipsis            0.3.2      2021-04-29 [1] CRAN (R 4.2.0)
   evaluate            0.20       2023-01-17 [1] CRAN (R 4.2.0)
   expm                0.999-7    2023-01-09 [1] CRAN (R 4.2.0)
   fansi               1.0.4      2023-01-22 [1] CRAN (R 4.2.0)
   fastmap             1.1.0      2021-01-25 [1] CRAN (R 4.2.0)
   fastmatch           1.1-3      2021-07-23 [1] CRAN (R 4.2.0)
   forcats           * 1.0.0      2023-01-29 [1] CRAN (R 4.2.0)
   foreach             1.5.2      2022-02-02 [1] CRAN (R 4.2.0)
   fs                  1.6.1      2023-02-06 [1] CRAN (R 4.2.0)
   generics            0.1.3      2022-07-05 [1] CRAN (R 4.2.0)
   ggfortify         * 0.4.14     2022-01-03 [1] CRAN (R 4.2.0)
   ggfun               0.0.9      2022-11-21 [1] CRAN (R 4.2.0)
   ggnewscale        * 0.4.8      2022-10-06 [1] CRAN (R 4.2.0)
   ggplot2           * 3.4.1      2023-02-10 [1] CRAN (R 4.2.0)
   ggplotify           0.1.0      2021-09-02 [1] CRAN (R 4.2.0)
   ggtree            * 3.6.2      2022-11-10 [1] Bioconductor
   glue                1.6.2      2022-02-24 [1] CRAN (R 4.2.0)
   gridExtra           2.3        2017-09-09 [1] CRAN (R 4.2.0)
   gridGraphics        0.5-1      2020-12-13 [1] CRAN (R 4.2.0)
   gtable              0.3.1      2022-09-01 [1] CRAN (R 4.2.0)
   here              * 1.0.1      2020-12-13 [1] CRAN (R 4.2.0)
   hms                 1.1.2      2022-08-19 [1] CRAN (R 4.2.0)
   htmltools           0.5.4      2022-12-07 [1] CRAN (R 4.2.0)
   htmlwidgets         1.5.4      2021-09-08 [1] CRAN (R 4.2.0)
   httpuv              1.6.6      2022-09-08 [1] CRAN (R 4.2.0)
   igraph              1.3.5      2022-09-22 [1] CRAN (R 4.2.0)
   iterators           1.0.14     2022-02-05 [1] CRAN (R 4.2.0)
   jsonlite            1.8.4      2022-12-06 [1] CRAN (R 4.2.0)
   knitr               1.42       2023-01-25 [1] CRAN (R 4.2.0)
   later               1.3.0      2021-08-18 [1] CRAN (R 4.2.0)
   lattice             0.20-45    2021-09-22 [1] CRAN (R 4.2.0)
   lazyeval            0.2.2      2019-03-15 [1] CRAN (R 4.2.0)
   lifecycle           1.0.3      2022-10-07 [1] CRAN (R 4.2.0)
   lubridate         * 1.9.2      2023-02-10 [1] CRAN (R 4.2.0)
   magrittr            2.0.3      2022-03-30 [1] CRAN (R 4.2.0)
   maps              * 3.4.1      2022-10-30 [1] CRAN (R 4.2.0)
   MASS                7.3-56     2022-03-23 [1] CRAN (R 4.2.0)
   Matrix            * 1.4-1      2022-03-23 [1] CRAN (R 4.2.0)
   MCMCglmm          * 2.34       2022-06-21 [1] CRAN (R 4.2.0)
   memoise             2.0.1      2021-11-26 [1] CRAN (R 4.2.0)
   mime                0.12       2021-09-28 [1] CRAN (R 4.2.0)
   miniUI              0.1.1.1    2018-05-18 [1] CRAN (R 4.2.0)
   mnormt              2.1.1      2022-09-26 [1] CRAN (R 4.2.0)
   munsell             0.5.0      2018-06-12 [1] CRAN (R 4.2.0)
   nlme                3.1-157    2022-03-25 [1] CRAN (R 4.2.0)
   numDeriv            2016.8-1.1 2019-06-06 [1] CRAN (R 4.2.0)
   optimParallel       1.0-2      2021-02-11 [1] CRAN (R 4.2.0)
   patchwork         * 1.1.2      2022-08-19 [1] CRAN (R 4.2.0)
   phangorn          * 2.11.1     2023-01-23 [1] CRAN (R 4.2.0)
   phytools          * 1.5-1      2023-02-19 [1] CRAN (R 4.2.0)
   pillar              1.8.1      2022-08-19 [1] CRAN (R 4.2.0)
   pkgbuild            1.4.0      2022-11-27 [1] CRAN (R 4.2.0)
   pkgconfig           2.0.3      2019-09-22 [1] CRAN (R 4.2.0)
   pkgload             1.3.2      2022-11-16 [1] CRAN (R 4.2.0)
   plotrix             3.8-2      2021-09-08 [1] CRAN (R 4.2.0)
   prettyunits         1.1.1      2020-01-24 [1] CRAN (R 4.2.0)
   processx            3.8.0      2022-10-26 [1] CRAN (R 4.2.0)
   profvis             0.3.7      2020-11-02 [1] CRAN (R 4.2.0)
   promises            1.2.0.1    2021-02-11 [1] CRAN (R 4.2.0)
   ps                  1.7.2      2022-10-26 [1] CRAN (R 4.2.0)
   purrr             * 1.0.1      2023-01-10 [1] CRAN (R 4.2.0)
   quadprog            1.5-8      2019-11-20 [1] CRAN (R 4.2.0)
   R6                  2.5.1      2021-08-19 [1] CRAN (R 4.2.0)
   Rcpp                1.0.10     2023-01-22 [1] CRAN (R 4.2.0)
   readr             * 2.1.4      2023-02-10 [1] CRAN (R 4.2.0)
   remotes             2.4.2      2021-11-30 [1] CRAN (R 4.2.0)
   RevGadgets        * 1.1.0      2023-02-03 [1] Github (revbayes/RevGadgets@6f8491b)
   rlang               1.0.6      2022-09-24 [1] CRAN (R 4.2.0)
   rmarkdown           2.20       2023-01-19 [1] CRAN (R 4.2.0)
   rprojroot           2.0.3      2022-04-02 [1] CRAN (R 4.2.0)
   rstudioapi          0.14       2022-08-22 [1] CRAN (R 4.2.0)
   scales              1.2.1      2022-08-20 [1] CRAN (R 4.2.0)
   scatterplot3d       0.3-42     2022-09-08 [1] CRAN (R 4.2.0)
   sessioninfo         1.2.2      2021-12-06 [1] CRAN (R 4.2.0)
   shiny               1.7.3      2022-10-25 [1] CRAN (R 4.2.0)
   stringi             1.7.12     2023-01-11 [1] CRAN (R 4.2.0)
   stringr           * 1.5.0      2022-12-02 [1] CRAN (R 4.2.0)
   tensorA             0.36.2     2020-11-19 [1] CRAN (R 4.2.0)
   tibble            * 3.1.8      2022-07-22 [1] CRAN (R 4.2.0)
   tidyr             * 1.3.0      2023-01-24 [1] CRAN (R 4.2.0)
   tidyselect          1.2.0      2022-10-10 [1] CRAN (R 4.2.0)
   tidytree            0.4.2      2022-12-18 [1] CRAN (R 4.2.0)
   tidyverse         * 2.0.0      2023-02-22 [1] CRAN (R 4.2.0)
   timechange          0.2.0      2023-01-11 [1] CRAN (R 4.2.0)
   treeio              1.22.0     2022-11-01 [1] Bioconductor
   tzdb                0.3.0      2022-03-28 [1] CRAN (R 4.2.0)
   urlchecker          1.0.1      2021-11-30 [1] CRAN (R 4.2.0)
   usethis             2.1.6      2022-05-25 [1] CRAN (R 4.2.0)
   utf8                1.2.3      2023-01-31 [1] CRAN (R 4.2.0)
   vctrs               0.5.2      2023-01-23 [1] CRAN (R 4.2.0)
   vroom               1.6.1      2023-01-22 [1] CRAN (R 4.2.0)
   withr               2.5.0      2022-03-03 [1] CRAN (R 4.2.0)
   xfun                0.37       2023-01-31 [1] CRAN (R 4.2.0)
   xtable              1.8-4      2019-04-21 [1] CRAN (R 4.2.0)
   yaml                2.3.7      2023-01-23 [1] CRAN (R 4.2.0)
   yulab.utils         0.0.6      2022-12-20 [1] CRAN (R 4.2.0)

   [1] /Library/Frameworks/R.framework/Versions/4.2/Resources/library

## Checkpoint for reproducibility
To rerun all the code with packages as they existed on CRAN at time of our analyses we recommend using the `checkpoint` package, and running this code prior to the analysis:

```{r}
checkpoint("2023-04-05")
```

