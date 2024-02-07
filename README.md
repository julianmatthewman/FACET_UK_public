# Analysis code and codelists for the UK part of the study "Association of different prescribing patterns for oral corticosteroids with fracture preventive care among older adults: Population-based cohort studies in the UK and Canada"

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7950694.svg)](https://doi.org/10.5281/zenodo.7950694)

Link to published article: [https://doi.org/10.1001%2Fjamadermatol.2023.2495](https://doi.org/10.1001%2Fjamadermatol.2023.2495)

## How to run

1.  Open the R console and call `renv::restore()` to install the required R packages.
2.  Provide a `paths.R` file in the `paths/` folder (see below for using dummy data)
3.  Call `tar_make()` to run the pipeline.

## How to inspect targets

-   Call `tar_read(target)` to retrieve a specified target.
-   Call `tar_visnetwork(targets_only = TRUE)` to visualise the pipeline.

## File structure

| File                      | Purpose                                                                                                        |
|-----------------------|-------------------------------------------------|
| [\_targets.R](_targets.R) | Declares the [`targets`](https://docs.ropensci.org/targets) pipeline. See `tar_script()` for details.          |
| [R/](R/)                  | Contains R scripts with functions to be used in the pipeline.                                                  |
| [codelists/](codelists/)  | Contains all codelists.                                                                                        |
| [renv.lock](renv.lock)    | The [`renv`](https://rstudio.github.io/renv/articles/renv.html) lock file that specifies all package versions. |

## Run with dummy data

To run analyses with dummy data, provide a `paths.R` file in the `paths/` folder, with the following contents:

```         
path_in <- "dummy_data/"
path_j_drive <- "dummy_data/"
path_z_drive <- "dummy_data/"
path_z_drive_linked <- "dummy_data/"
path_common_dosages <- "dummy_data/Lookups_2022_01/common_dosages.txt"
```
