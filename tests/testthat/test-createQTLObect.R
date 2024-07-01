library(methods)

# Load test data and create a mock eQTLObject
data(testGene)
data(testSNP)

test_that("createQTLObject handles matrix input correctly", {
  eqtl <- createQTLObject(snpMatrix = testSNP, genedata = testGene)
  expect_true(inherits(eqtl, "eQTLObject"))
  expect_equal(nrow(eqtl@rawData$rawExpMat), 100)
  expect_equal(ncol(eqtl@rawData$rawExpMat), 10)
  expect_equal(nrow(eqtl@rawData$snpMat), 100)
  expect_equal(ncol(eqtl@rawData$snpMat), 10)
})

test_that("createQTLObject handles Seurat object input correctly", {
  # Assuming you have a Seurat object or a way to mock one
  seuratMock <- list(
    counts = testGene,
    meta.data = data.frame(celltype = rep("TypeA", 10))
  )
  eqtl <- createQTLObject(snpMatrix = testSNP, genedata = seuratMock)
  expect_true(inherits(eqtl, "eQTLObject"))
  expect_equal(nrow(eqtl@rawData$rawExpMat), 100)
  expect_equal(ncol(eqtl@rawData$rawExpMat), 10)
  expect_equal(nrow(eqtl@rawData$snpMat), 100)
  expect_equal(ncol(eqtl@rawData$snpMat), 10)
  expect_equal(nrow(eqtl@groupBy), 10)
})

test_that("createQTLObject handles incorrect input gracefully", {
  expect_error(createQTLObject(snpMatrix = testSNP, genedata = "invalid_input"))
  # Add more expect_error tests for other cases as needed
})
