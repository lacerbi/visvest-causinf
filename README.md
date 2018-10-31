# Bayesian Comparison of Explicit and Implicit Causal Inference Strategies in Multisensory Heading Perception

This repository accompanies the manuscript by Acerbi et al. (2018), published in *PLoS Computational Biology* [[1](https://github.com/lacerbi/visvest-causinf/blob/master/README.md#reference)].
It includes human subjects' data and the code used for fitting and comparing the models reported in the paper.

## The repository

The original working name of the repository was `VestBMS` (standing for *Bayesian Model Selection of Vestibular experiment*), which appears as prefix of many functions and folders. 

The structure of the repository is as follows: 

- *Analytics*: Functions to analyze data and results of model fits.
- *ModelWork*: Core model fitting functions (to be used in conjunction with the [ModelWork](https://github.com/lacerbi/ModelWork) toolbox).
- *PlotFunctions*: Functions for plotting results.
- *data*: Raw data for each subject, in `.mat` form.
- *scripts*: Bash scripts used to submit jobs on the computer cluster (using [SLURM](https://slurm.schedmd.com/)).
- *supplement*: Files used to analyze data in [de Winkel et al. (2017)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0169676), for a control analysis reported in the Supporting Information.

Please contact me at <luigi.acerbi@unige.ch> if you have any questions.

## Reference

1. Acerbi\*, L., Dokka\*, K., Angelaki, D. E. & Ma, W. J. (2018). Bayesian Comparison of Explicit and Implicit Causal Inference Strategies in Multisensory Heading Perception, *PLoS Computational Biology* 14(7): e1006110. (\*equal contribution; [link](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006110))

### License

This code is released under the terms of the [MIT License](https://github.com/lacerbi/visvest-causinf/blob/master/LICENSE.txt).

