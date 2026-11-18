# seafar: Sparse Exploratory Approximate Factor Analysis for High-dimensional Data in R <img src="man/figures/logo.png" align="right" height="139" alt="" />

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.20326842-blue.svg)](https://doi.org/10.5281/zenodo.20326842)

The package contains functions to perform a Sparse Exploratory Approximate Factor Analysis (SEAFA). These methods are implementations and extensions of the Regularized ESEM method ([Le et al. (2026)](https://doi.org/10.3758/s13428-026-02960-y)) that uses approximate factor model with regularization (cardinality constraint and the LASSO penalty) to obtain simple structure for the loading matrix. A vignette was created as a tutorial paper for interested users.  

Please cite the package as shown in `CITATION.cff`.

## Installation instructions

This package requires `PCAtools` which was not available on [CRAN](https://cran.r-project.org/) at the time of writing. This can lead to problems when installing `seafar` from GitHub. We therefore recommend installing `PCAtools` from [Bioconductor](https://cran.r-project.org/web/packages/BiocManager/vignettes/BiocManager.html) before installing `seafar`.

```r
install.packages("BiocManager")
BiocManager::install("PCAtools")
```

Currently the only way to install the package with the vignette is to build it from source. This might take a few minutes depending on your system.

```r
install.packages("remotes")
remotes::install_github(
  "trale97/seafar@*release",
  build_vignettes = TRUE,
  dependencies = TRUE
)
```

You can also skip the vignette for a quicker installation:

```r
install.packages("remotes")
remotes::install_github("trale97/seafar@*release")
```

## Getting started

The best way to get started with `seafar` is to follow the tutorial paper we published for the package. You can access it in R like this:

```r
library(seafar)
browseVignettes("seafar")
```

## Reproducibility statement

The tutorial paper can be recreated from source by following the steps below:

1. Make sure you have `seafar` installed.

2. Install [Quarto](https://quarto.org/docs/get-started/).

3. If you don't have a LaTeX distribution yet, we recommend TinyTeX:

```r
install.packages("tinytex")
tinytex::install_tinytex()
```

4. Navigate into `vignettes/`

```bash
cd vignettes
```

5. Render the document by running the following command:

```bash
quarto render SEAFAR_tutorial.qmd --to apaquarto-pdf
```

## Development guidelines

Contributions to the development of `seafar` are welcome. The author of the package will review all pull request. In this section, you can find the preferred development workflows for this package.

### Workflow for contribution

1. Create an issue on the `seafar` [GitHub issues page](https://github.com/trale97/seafar/issues) to describe you intended contribution.
2. Create a [branch for the issue](https://docs.github.com/en/issues/tracking-your-work-with-issues/using-issues/creating-a-branch-for-an-issue).
3. Create a [fork](https://docs.github.com/en/pull-requests/how-tos/work-with-forks/fork-a-repo#forking-a-repository) of `seafar` or sync your existing fork with the upstream repository
4. Make changes to the fork
5. Run a check with `devtools::check()`
7. Review any potential errors or warnings that are highlighted by the `devtools` report
8. Commit changes to the fork
9. Make sure the fork is up-to-date with the upstream branch (resolve any conflicts yourself)
10. Create a merge request to the branch of the issue on the upstream repository.

### Recommended workflow for submission of package to CRAN

We follow the main workflow highlighted in the [R Packages](https://r-pkgs.org/release.html#double-r-cmd-checking) book.

#### Run checks

1. [ ] `devtools::check(remote = TRUE, manual = TRUE)` on the computer you are using to develop the package, with the current version of R. It includes some extra checks usually turned off to make day-to-day development faster.
2. [ ] `devtools::check_win_devel(email = your@email.com)` on the computer you are using to develop the package. This will send the package to 

#### Update release documentation

1. [ ] Make sure you have update the version number in both [DESCRIPTION](DESCRIPTION) and [CITATION](CITATION.cff).
2. [ ] Update the [NEWS.md](NEWS.md) or create one via `usethis::use_news_md()`.
3. [ ] `usethis::use_cran_comments()` to add a [cran-comments.md](cran-comments.md) to communicate to CRAN the output of the checks we have run.

#### Compile

In a terminal, compile the package running the following code in the parent directory of the cloned repository on your computer:

```bash
R CMD build seafar
```

This will create a seafar_x.y.z.tar.gz in the parent directory

#### Check again

Run a final check on seafar_x.y.z.tar.gz

```bash
R CMD check --as-cran seafar_x.y.z.tar.gz
```

#### Upload

Upload the seafar_x.y.z.tar.gz without errors and warnings to the [CRAN Submission Portal](https://cran.r-project.org/submit.html)

## Terms of use

This package is licensed under the MIT License.

The autism dataset is licensed under the terms specified in `data/README.md`.

The [APA Quarto template](https://github.com/wjschne/apaquarto) by W. J. Schneider is licensed under [CC0-1.0](https://creativecommons.org/publicdomain/zero/1.0/).
