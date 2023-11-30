kappa.fleiss <- function(ratings, exact = FALSE, detail = FALSE) {
  ratings <- as.matrix(na.omit(ratings))
  ns <- nrow(ratings)
  nr <- ncol(ratings)
  lev <- sort(unique(as.vector(ratings)))
  
  # Build frequency table for each subject
  ttab <- t(apply(ratings, 1, function(row) table(factor(row, levels = lev))))
  
  # Agreement proportion
  agreeP <- sum(rowSums(ttab^2) - nr) / (nr * (nr - 1) * ns)
  
  # Chance agreement proportion
  if (!exact) {
    method <- "Fleiss' Kappa for m Raters"
    chanceP <- sum(colSums(ttab)^2) / (ns * nr)^2
  } else {
    method <- "Fleiss' Kappa for m Raters (exact value)"
    rtab <- t(sapply(1:nr, function(i) table(factor(ratings[, i], levels = lev))))
    rtab <- rtab / ns
    
    prop_sum_squared <- sum(colSums(ttab)^2) / (ns * nr)^2
    var_across_raters <- sum(apply(rtab, 2, var) * (nr - 1) / nr) / (nr - 1)
    
    chanceP <- prop_sum_squared - var_across_raters
  }
  
  # Fleiss' Kappa for m raters
  value <- (agreeP - chanceP) / (1 - chanceP)
  
  # Standard error & test statistic
  if (!exact) {
    pj <- colSums(ttab) / (ns * nr)
    qj <- 1 - pj
    varkappa <- (2 / (sum(pj * qj)^2 * (ns * nr * (nr - 1)))) * (sum(pj * qj)^2 - sum(pj * qj * (qj - pj)))
    SEkappa <- sqrt(varkappa)
    u <- value / SEkappa
    p.value <- 2 * (1 - pnorm(abs(u)))
    
    if (detail) {
      pjk <- (colSums(ttab^2) - ns * nr * pj) / (ns * nr * (nr - 1) * pj)
      kappaK <- (pjk - pj) / (1 - pj)
      varkappaK <- 2 / (ns * nr * (nr - 1))
      SEkappaK <- sqrt(varkappaK)
      uK <- kappaK / SEkappaK
      p.valueK <- 2 * (1 - pnorm(abs(uK)))
      tableK <- round(cbind(kappaK, uK, p.valueK), 3)
      rownames(tableK) <- lev
      colnames(tableK) <- c("Kappa", "z", "p.value")
    }
  }
  
  # Returned information
  rval <- list(method = method, subjects = ns, raters = nr, value = value)
  
  if (!exact) {
    if (detail) {
      rval$detail <- tableK
    }
    rval <- c(rval, stat.name = "z", statistic = u, p.value = p.value)
  }
  
  return(rval)
}

