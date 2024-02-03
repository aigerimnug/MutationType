library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)

test_that("values for context length are correct", {
  vcffile <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
  vcf <- readVcf(vcffile, "hg19")
  
  expect_error(MutationsType(vcf, Hsapiens, '3'))
  expect_equal(MutationsType(vcf, Hsapiens, 1.5), character(0))
})

