#Function file: montecarlo.R

# A hand-built cn_result so the expected p-values are known exactly.
# Simulated maxima: 2 x 10, 3 x 5, 4 x 4, 5 x 1  (nsim = 20)
.fake_cn <- function(alpha = 0.05) {
  sims <- c(rep(2, 10), rep(3, 5), rep(4, 4), 5)
  list(critical_n = stats::quantile(sims, 1 - alpha),
       ci_lower = NA_real_, ci_upper = NA_real_,
       simulated_values = sims, alpha = alpha)
}

# Observed counts: A = 6, B = 5, C = 4, D = 1
.sigs <- list(
  s1 = c("A", "B", "C", "D"), s2 = c("A", "B", "C"), s3 = c("A", "B", "C"),
  s4 = c("A", "B", "C"), s5 = c("A", "B"), s6 = "A"
)

test_that("p_empirical is (1 + hits) / (nsim + 1) and never zero", {
  res <- significantTaxa(.sigs, .fake_cn(), plot = FALSE)
  p <- setNames(res$p_empirical, as.character(res$taxon))
  expect_equal(p[["A"]], 1 / 21)          # no simulation reaches 6
  expect_equal(p[["B"]], (1 + 1) / 21)    # one simulation reaches 5
  expect_equal(p[["C"]], (1 + 5) / 21)
  expect_equal(p[["D"]], (1 + 20) / 21)
  expect_true(all(res$p_empirical > 0))
})

test_that("significant agrees with p_empirical <= alpha", {
  res <- significantTaxa(.sigs, .fake_cn(alpha = 0.05), plot = FALSE)
  expect_identical(res$significant, res$p_empirical <= 0.05)
  expect_identical(as.character(res$taxon[res$significant]), "A")

  # a looser alpha picks up B as well
  res10 <- significantTaxa(.sigs, .fake_cn(alpha = 0.10), plot = FALSE)
  expect_setequal(as.character(res10$taxon[res10$significant]), c("A", "B"))
})

test_that("cn_result without alpha falls back to 0.05", {
  cn <- .fake_cn()
  cn$alpha <- NULL
  res <- significantTaxa(.sigs, cn, plot = FALSE)
  expect_identical(res$significant, res$p_empirical <= 0.05)
})
