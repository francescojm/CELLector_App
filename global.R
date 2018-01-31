library(shiny)
library(shinythemes)
library(stringr)
library(plotrix) #polar.plot(), pie3D()
library(devtools)
library(RColorBrewer)
library(CELLector)

#install_github("sjp/grImport2")
library(grImport2)

## Loading Primary Tumour Binary Event Matrices
data(CELLector.PrimTum.BEMs)
data(CELLector.CFEs)

## Deriving available TCGA labels
TCGALabels<-names(CELLector.PrimTum.BEMs)

tumours<-CELLector.PrimTum.BEMs$COREAD
features<-CELLector.CFEs

