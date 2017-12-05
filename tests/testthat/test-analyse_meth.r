context("Analysing methylome")

test_that("Sun exposure data could be analysed", {
  cols = as.numeric(sunexp_design$sex) * 2
  names(cols) = rownames(sunexp_design)

  gene = genes[1,]
  res1 = analyse_meth(gene, sunexp_data, sunexp_platform, cols, PLOT=TRUE)
  layout(matrix(1:6,2))
  res1 = analyse_meth(genes[1,], sunexp_data, sunexp_platform)  
  res1 = analyse_meth(genes[2,], sunexp_data, sunexp_platform)  
  res1 = analyse_meth(genes[3,], sunexp_data, sunexp_platform)  
  res1 = analyse_meth(genes[4,], sunexp_data, sunexp_platform)  
  res1 = analyse_meth(genes[5,], sunexp_data, sunexp_platform)  
  res1 = analyse_meth(genes[6,], sunexp_data, sunexp_platform)  
})
