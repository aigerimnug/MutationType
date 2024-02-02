## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----style, echo = FALSE, results = 'asis'------------------------------------
library(BiocStyle)

## ----echo = FALSE-------------------------------------------------------------
library(knitr)

## ----setup, results="hide", include = FALSE-----------------------------------
library(MutationType)

## -----------------------------------------------------------------------------
library("VariantAnnotation")
library("BSgenome")
library("BSgenome.Hsapiens.UCSC.hg19")

## -----------------------------------------------------------------------------
vcffile <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
vcf <- readVcf(vcffile, "hg19")

## -----------------------------------------------------------------------------
determine_mutation_types(vcf, Hsapiens, 3)

## -----------------------------------------------------------------------------
count_SNVmutations(vcf, Hsapiens)

