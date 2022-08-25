# chameleon-chromosomes
# Microchromosome fusions underpin convergent evolution of chameleon karyotypes

Author(s): [Natalie Cooper](mailto:natalie.cooper.@nhm.ac.uk)

THIS IS A WORK IN PROGRESS!

This repository contains all the code and some data used in the [paper](TO ADD link). 

To cite the paper: 
> Marcello Mezzasalma, Jeffrey W. Streicher, Fabio M. Guarino, Marc E.H. Jones, Simon P. Loader, Gaetano Odierna, and Natalie Cooper 2022. Microchromosome fusions underpin convergent evolution of chameleon karyotypes. TBC.

To cite this repo: 
> Natalie Cooper. 2022. GitHub: nhcooper123/chameleon-chromosomes: code for the paper. Zenodo. DOI: TO ADD.

![alt text](https://github.com/nhcooper123/chameleon-chromosomes/raw/master/manuscript/figures/new-tree-plot4.png)

## Data
All raw and cleaned data are available from the [NHM Data Portal](TO ADD).

For reproducibility purposes, if you want to rerun the data wrangling steps in script `analyses\01-data-wrangling`, download `CHROMO_STATES_2020-10-16.csv` and place it into a `raw-data/` folder along with the dataset of Meiri 2018: https://onlinelibrary.wiley.com/doi/full/10.1111/geb.12773.

The cleaned datasets, along rest of the ancillary data required to run all analyses and produce all figures is within the `data` folder.

* *ML-tree-species.tre* - tree needed for the main analyses in the paper
* *chromosome-data-species.csv* - data needed for the main analyses in the paper
* *ML-tree.tre* - tree for the all OTUs analyses in the supplementary materials
* *chromosome-data.csv* - data needed for the main analyses in the paper
* *ML-tree.tre* - tree for the all OTUs analyses in the supplementary materials
* *chromevol-simulations-outputs-best.csv* - outputs from the "best" model of chromosome evolution used in the main text. The equivalent for the all OTUs analyses is in the folder `data/OTUs`.
* *null-simulations-outputs.csv* - outputs from the null simulations using the best model of chromosome evolution used in the main text. The equivalent for the all OTUs analyses is in the folder `data/OTUs`.

If you use the cleaned data please cite as follows: 
> Marcello Mezzasalma, Jeffrey W. Streicher, Fabio M. Guarino, Marc E.H. Jones, Simon P. Loader, Gaetano Odierna, and Natalie Cooper 2022. CHROMREP [Data set]. Natural History Museum. https://doi.org/10.5519/rmaopg7d.

The life-history data all comes from Meiri 2018 so please cite this paper if you use it.

> Meiri, S. Traits of lizards of the world: Variation around a successful evolutionary design. Global Ecol Biogeogr. 2018; 27: 1168– 1172

## ChromEvol analyses (`ChromEvol/` folder)

The majority of analyses run in R, but we also used ChromEvol, a stand alone package available from http://chromevol.tau.ac.il/. The associated paper for the package is:

> Lior Glick, Itay Mayrose, ChromEvol: Assessing the Pattern of Chromosome Number Evolution and the Inference of Polyploidy along a Phylogeny, Molecular Biology and Evolution, Volume 31, Issue 7, July 2014, Pages 1914–1922, https://doi.org/10.1093/molbev/msu122   

To run these analyses you will need to download ChromEvol and save it in the `chromEvol/` folder. This unzips into a folder called `chromEvol_source-current/` that contains the package as an executable called `chromEvol`. I have not uploaded this to the repo because it's not my package. 

ChromEvol analyses require a dataset, a tree and a control file. These are within the `chromEvol/` folder preceded by `data_for_chromevol` or `tree_for_chroevol` or with the word `control` somewhere in their name respectively. Several versions exist for the different sensitivity analyses. These analyses save their outputs into the `chromEvol/results` folder, in folders named based on the names of the control files. Full details of how to build control files is in the ChromEvol manual at http://chromevol.tau.ac.il/. We provide code for creating appropriate datasets and trees in `analyses/02-`

Running chromEvol analyses is simple. Just open a Terminal in the `chromEvol` folder and run:

`chromEvol_source-current/chromEvol name-of-control-file`

You can run all analyses from our paper using the shell script: `run_all_models.sh`. 

The `ChromEvol/` folder also includes a folder called `root_freq`. This has the three different root frequency definitions used in the different ChromEvol analyses (see text).

Finally the `OTUs` folder within `ChromEvol/` contain the same things but for analyses using all OTUs (found in the supplementary materials of paper).

-------
## R analyses
The analysis code is divided into `.Rmd` files that run the analyses and plot the figures for each section of the paper/supplementary materials. 

Note that several of these scripts use an R package written by N.Cusimano. You can download this here https://www.en.sysbot.bio.lmu.de/people/employees/cusimano/use_r/. The easiest way to use this is to download the package and save the `ChromEvol.R` file. You can then install these functions in R using `source`. In my repo I kept this within the `chromEvol` folder so used the code `source("chromEvol/ChromEvol.R")` to install the functions.

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
22. **F01-figure-tree-chromosomes.R** - creates figure 1 from ChromEvol outputs.
23. **figure-tree-chromosomes_main.R** - creates several versions of figure 1 from ChromEvol outputs.

-------
## Other folders

* `/outputs` contains the figures and tables
* `/manuscript` contains the manuscript materials in LaTeX format

-------
## Session Info
For reproducibility purposes, here is the output of `devtools::session_info()` used to perform the analyses in the publication.

## Checkpoint for reproducibility
To rerun all the code with packages as they existed on CRAN at time of our analyses we recommend using the `checkpoint` package, and running this code prior to the analysis:

```{r}
checkpoint("2022-09-10")
```

