context("miniplot")

test_that("miniplot works as expected", {
  expect_error(miniplot())

  tf <- tempfile(fileext = ".png")
  png(tf)
  expect_null(miniplot(1, 2))
  dev.off()

  unlink(tf)
})
