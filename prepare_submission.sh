#! /bin/bash

Rscript -e "library(knitr); knit('report.Rnw')"
jupyter nbconvert --to latex models.ipynb 
pdftk report.pdf models.pdf cat output EDA_modelling.pdf
zip -r predictions.zip predictions
