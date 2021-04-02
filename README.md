## Welcome

Files in this GitHub repository can be used for analyses of
longitudinal biomarkers of the Maastr**ICC**ht cohort.

## Description

* The file _COVID19\_DATABEW\_BIOM.R_ is used to read in data from
  several files to be combined into a single data file, which is saved
  as _d\_long_

* Next, _COVID19\_ANALYSES\_BIOM.R_ is used to perform linear
  mixed-effects regression after loading the previously made dataset
  _d\_long_. The crude regression models are saved for making figures.

* The file _COVID19\_FIGURES\_BIOM.R_ plots scatter plots with
  regression lines based on the crude models saved in the previous
  steps.
