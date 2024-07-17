data(testGene)
data(testSNP)
eqtl <- createQTLObject(snpMatrix = testSNP,
                        genedata = testGene,
                        biClassify = FALSE,
                        species = 'human',
                        group = NULL)

eqtl <- normalizeGene(eqtl, method = "logNormalize")

eqtl <- filterGeneSNP(eQTLObject = eqtl,
                      snp.number.of.cells.percent = 2,
                      expression.min = 0,
                      expression.number.of.cells.percent = 2)

eqtl <- callQTL(eQTLObject = eqtl,
                gene_ids = NULL,
                downstream = NULL,
                upstream = NULL,
                p.adjust.method = "bonferroni",
                useModel = "zinb",
                p.adjust.Threshold = 0.05,
                logfc.threshold = 0.1)


test_that("visualizeQTL function behaves as expected", {

  # test QTLplot result
  plot1 <- visualizeQTL(eqtl, SNPid = "1:632647", Geneid = "RPS27", plottype = "QTLplot")
  expect_true(is.list(plot1))

  # test violin result
  plot2 <- visualizeQTL(eqtl, SNPid = "1:632647", Geneid = "RPS27", plottype = "violin")
  expect_true(is.list(plot2))

  # test boxplot result
  plot3 <- visualizeQTL(eqtl, SNPid = "1:632647", Geneid = "RPS27", plottype = "boxplot")
  expect_true(is.list(plot3))

  # test histplot result
  plot4 <- visualizeQTL(eqtl, SNPid = "1:632647", Geneid = "RPS27", plottype = "histplot")
  expect_true(is.list(plot4))

  # when specific group
  data(testSeurat)
  eqtl <- createQTLObject(snpMatrix = testSNP,
                          genedata = testSeurat,
                          biClassify = FALSE,
                          species = 'human',
                          group = "celltype")
  eqtl <- normalizeGene(eqtl, method = "logNormalize")
  eqtl <- filterGeneSNP(eQTLObject = eqtl,
                        snp.number.of.cells.percent = 2,
                        expression.min = 0,
                        expression.number.of.cells.percent = 2)
  eqtl <- callQTL(eQTLObject = eqtl,
                  useModel = "poisson")
  plot5 <- visualizeQTL(eqtl, SNPid = "1:632647", Geneid = "RPS27", groupName = "GMP", plottype = "QTLplot")
  expect_true(is.list(plot5))

  # test invalid plot types
  expect_error(visualizeQTL(eqtl, SNPid = "1:632647", Geneid = "RPS27", plottype = "invalid_type"),
               "Invalid plottype,
         Please choose from 'QTLplot', 'violin' , 'boxplot' or 'histplot'.")
})
