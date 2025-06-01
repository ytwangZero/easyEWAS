# File: R/ewasfun_utils.R

#' @title EWAS Model Computation Utilities
#' @description Functions to perform EWAS model fitting for lm, lmer, and cox

# Linear model version of ewasfun
#' @export
ewasfun_lm <- function(cg, ff, cov, facnum) {
  cov$cpg <- as.vector(t(cg))
  out <- summary(lm(ff, data = cov))
  unlist(lapply(2:facnum, function(i) out$coefficients[i, c(1, 2, 4)]))
}

# Linear mixed-effect model version
#' @export
ewasfun_lmer <- function(cg, ff, cov, facnum) {
  cov$cpg <- as.vector(t(cg))
  out <- summary(lmer(ff, data = cov))
  unlist(lapply(2:facnum, function(i) out$coefficients[i, c(1, 2, 5)]))
}

# Cox proportional hazard model version
#' @export
ewasfun_cox <- function(cg, ff, cov) {
  cov$cpg <- as.vector(t(cg))
  out <- summary(coxph(ff, data = cov))
  c(as.vector(out$conf.int[1, c(1, 3, 4)]), out$coefficients[1, 5])
}
