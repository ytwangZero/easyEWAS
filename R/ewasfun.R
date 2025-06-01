# file: R/ewasfun.R

# Separate, standalone model functions (no closure pollution)
ewasfun_lm <- function(cg, formula, covdata, facnum) {
  covdata$cpg <- as.vector(t(cg))
  out <- summary(lm(formula, data = covdata))
  unlist(lapply(2:facnum, function(i) out$coefficients[i, c(1, 2, 4)]))
}

ewasfun_lmer <- function(cg, formula, covdata, facnum) {
  covdata$cpg <- as.vector(t(cg))
  out <- summary(lmer(formula, data = covdata))
  unlist(lapply(2:facnum, function(i) out$coefficients[i, c(1, 2, 5)]))
}

ewasfun_cox <- function(cg, formula, covdata) {
  covdata$cpg <- as.vector(t(cg))
  out <- summary(coxph(formula, data = covdata))
  c(as.vector(out$conf.int[1, c(1, 3, 4)]), out$coefficients[1, 5])
}
