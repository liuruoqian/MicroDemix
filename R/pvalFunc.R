#' MD_pvalue
#'
#' This function generates a simulation based p-value of testing differential abundance
#' under the null hypothesis that p.r = p.o.Variability in the results is anticipated due to the
#' randomness embedded in this function.
#' @param p.r A vector for the relative abundance of corresponding
#' taxa in group 1 (e.g. the microbial relative abundance in rectum).
#' @param p.o A vector for the relative abundance of corresponding
#' taxa in group 2 (e.g. the microbial relative abundance in other GI locations).
#' @param N.r total counts across all samples in group 1 (i.e. the sum of sequence
#' depths of all samples in group 1).
#' @param N total counts across all samples in group 2 (i.e. the sum of sequence
#' depths of all samples in group 2).
#' @return A scalar in [0, 1].
#' @examples
#' p1 <- rep(0.1, 10)
#' p2 <- c(rep(0.09, 5), rep(0.11, 5))
#' MD_pvalue(p1, p2, 5000, 5000);
#' @export


MD_pvalue <- function(p.r, p.o, N.r, N){
  G <- length(p.o)
  total_counts.n <- N.r
  total_counts.c <- N

  p1.hat <- p.r

  Sig <- matrix(NA, G, G)
  for (i in 1:G){
    for (j in 1:G){
      Sig[i,j] <- -p1.hat[i]*p1.hat[j]
    }
  }

  Sig <- Sig + diag(p1.hat,G,G)

  samp1 <- mvrnorm(10000, mu = rep(0, G-1), Sigma = Sig[1:(G-1),1:(G-1)])
  samp2 <- mvrnorm(10000, mu = rep(0, G-1), Sigma = Sig[1:(G-1),1:(G-1)])


  samp1 <- samp1/sqrt((total_counts.n))
  samp2 <- samp2/sqrt((total_counts.c))


  samp_1 <- cbind(samp1, -rowSums(samp1))
  samp_2 <- cbind(samp2, -rowSums(samp2))

  # samp11 <- samp_1[,4:8]
  # samp21 <- samp_2[,4:8]


  SS_func <- function(v1,v2){
    return (sum(v1^2)+sum(v2^2)-2*sum(v1*v2))
  }

  null_samp <- mapply(SS_func, v1 = as.list(data.frame(t(samp_1))), v2=as.list(data.frame(t(samp_2))))

  t <- sum((p.r-p.o)^2)
  return (mean(null_samp >= t))
}
