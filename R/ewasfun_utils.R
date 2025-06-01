#' @title EWAS Model Computation Utilities
#' @description Model fitting functions for EWAS analysis using linear, mixed, or Cox models.
#' These are designed to be used inside parallel loops with minimal memory footprint.
#'
#' @param cg A row vector representing one CpG site's beta values across samples.
#' @param ff A model formula object (e.g., cpg ~ var1 + var2).
#' @param cov A data.frame containing the covariates for the model.
#' @param facnum (Only for lm/lmer) Number of factor levels for exposure variable.
#'
#' @return A numeric vector with model coefficients, standard errors, and p-values.
#'
#' @export
ewasfun_lm <- function(cg, ff, cov, facnum) {
  cov$cpg <- as.vector(t(cg))
  res <- tryCatch({
    out <- summary(stats::lm(ff, data = cov))
    unlist(lapply(2:facnum, function(i) out$coefficients[i, c(1, 2, 4)]))
  }, error = function(e) {
    rep(NA_real_, 3 * (facnum - 1))
  })
  return(res)
}

#' @rdname ewasfun_lm
#' @export
ewasfun_lmer <- function(cg, ff, cov, facnum) {
  cov$cpg <- as.vector(t(cg))
  res <- tryCatch({
    out <- summary(lmerTest::lmer(ff, data = cov))
    unlist(lapply(2:facnum, function(i) out$coefficients[i, c(1, 2, 5)]))
  }, error = function(e) {
    rep(NA_real_, 3 * (facnum - 1))
  })
  return(res)
}

#' @rdname ewasfun_lm
#' @export
ewasfun_cox <- function(cg, ff, cov) {
  cov$cpg <- as.vector(t(cg))
  res <- tryCatch({
    out <- summary(survival::coxph(ff, data = cov))
    c(as.vector(out$conf.int[1, c(1, 3, 4)]), out$coefficients[1, 5])
  }, error = function(e) {
    rep(NA_real_, 4)
  })
  return(res)
}
