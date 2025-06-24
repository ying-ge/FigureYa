## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load-package, message=FALSE-----------------------------------------
library(DealGPL570)

## ---- message=FALSE, eval=FALSE------------------------------------------
#  GEOquery::getGEOSuppFiles("GSE104683", makeDirectory = FALSE, baseDir = tempdir())
#  file <- list.files(path = tempdir(), pattern = "GSE104683_RAW.tar", full.names = TRUE)
#  file
#  result <- DealGPL570(file = file)
#  result[1:5, 1:3]

## ---- eval=FALSE---------------------------------------------------------
#  result <- DealGPL570(file = file, type = "probeID")
#  result[1:5, 1:2]

## ---- eval=FALSE---------------------------------------------------------
#  result <- DealGPL570(file = file, type = "geneSymbol")
#  result[1:5, 1:2]

