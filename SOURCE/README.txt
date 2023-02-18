1) build the cuda code using the following code
nvcc -O2 InterOptCuda.cu -o InterOptCuda -std=c++11

2) run RunAll_Experiments.R
3) run BuildResults.R

The imputation is done using the Imputation.R script
The imputed data are already provided in the DATA directory.

Required R libraries:

library('matrixStats')

library('edgeR')
library('miRBaseConverter')
library('GEOquery')

library('MASS')
library('parallel')
library('CovTools')

library('pheatmap')
library('ggplot2')
library('reshape')
library('scales')
library('cowplot')
library('ggsignif')
library('xlsx')
library('ggsci')