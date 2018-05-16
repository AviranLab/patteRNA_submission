# patteRNA Analysis Scripts

Analysis tools for patteRNA data.


## What Is It?

This repository contains all code and score files necessary to produce figures included in our submission to Genes.

## Getting Started

The following files must be unzipped (`gunzip`) before starting:
* `PARS/GM12878.pars.gz`
* `PARS/4SL_d75/scores.txt.gz`
* `PARS/5SL_d75/scores.txt.gz`

These tools use Python v3. You must also have the following dependencies installed in order to run these tools:

* NumPy
* matplotlib
* SciPy

For more advanced users, the use of a virtual python environment is recommended (see [venv](https://docs.python.org/3/library/venv.html#module-venv) and [pyenv](https://github.com/pyenv/pyenv)).

## Data

Output score files from patteRNA are held in two folders, `genome_scores` and `scores`. Scores from the whole-genome searches of HIV-1 are kept in `genome_scores`, and the names of score files follow a convention to indicate the type of search. Files containing "full-length" are results when scoring the full 232-nt RRE motif in the data. Files containing "short" are results when searching for the 59-nt SL III/SL IV region. The tag `sc` indicates that sequence constraints were enforced. The rest of the score files are in the `scores` folder and their names indicate their contents. Finally, raw HIV RRE data are included as necessary to produce *in silico* mixtures.

## Scripts

Within the `analysis` folder are a collection of scripts to produce each relevant figure from the main text and supplementary figures. Most scripts are self-contained; just run the script to compile the data and produce figures. The exception is an additional script used to generate *in silico* mixtures of RRE, which generates SHAPE files instead of figures.

## Auxillary Files

Other folders in this repository contain files helpful during analysis. These folders are `motif_files`, which contains lists of motifs used in certain searches, and `sample_nulls`, which contains a temporary null distribution file used to demonstrate the normalization process for Figure 2. Two local libraries of useful methods for our analysis is also included (`bundlelib` and `misclib`).
