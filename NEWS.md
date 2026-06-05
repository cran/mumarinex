# mumarinex 2.0

## > New function

* Added the `mum_GUI` function

## > Minor changes

* Added a `log` argument (default `TRUE`)to the `mumarinex`, `decomplot` and `diagnostic_tool` functions, allowing users to disable log-transformation prior to CBCS sub-index computation for communities without overabundant taxa.

* Implemented calculation of relative abundance decrease and increase percentages for CBCS within the `diagnostic_tool` function (columns `Relative_D` and `Relative_I`, respectively).

* Updated the `Simulated_data` dataset to ensure consistency with example data used in the reference publication.

* Updated function examples to match the new structure of the `Simulated_data` format.

* Revised sub-index nomenclature to ensure consistency with terminology updated following reviewer feedback.

## > Bug correction

* Fixed a graphical rendering issue in the `decomplot` function related to incorrect plot window scaling.

# mumarinex 1.0

* Initial CRAN submission.
