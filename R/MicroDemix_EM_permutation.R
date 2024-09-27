#' MicroDemix_EM_permutation
#'
#' This function performs a permutation based test of differential abundance
#' under the null hypothesis that p.r = p.o. Variability in the results is anticipated due to the
#' randomness embedded in this function.
#' @import optimx
#' @import BB
#' @import MCMCpack
#' @param data1 stool (mixture) microbiome data set with taxa in
#'            columns and samples in rows.
#' @param data2 rectum (reference) microbiome data set with taxa in
#'            columns and samples in rows. The taxa in data1 and data2 should
#'            be matching.
#' @param cova1 covariates related to samples in data1 that have a matching order
#' with data1. 'cova1' must be a matrix with covariates in columns and observations
#' in rows.
#' @param cova2 covariates related to samples in data2 that have a matching order
#' with data2. 'cova2' must be a matrix with covariates in columns and observations
#' in rows.
#' @param maxiter The maximum number of iterations taken by the EM algorithm. Defaults to 10.
#' @param B The number of permutations.
#' @return The p-value. A scalar in [0, 1].
#' @examples
#' MicroDemix_EM_permutation(yig, yig.n, x, x.n, maxiter=5, B=10);
#' @export


MicroDemix_EM_permutation <- function(data1, data2, cova1, cova2, maxiter=10, B){

  loglik <- function(pi.s, y.vec, x.vec, para, pr){

    G <- length(y.vec)
    p <- length(x.vec)

    ### para[1:G]: u_g: po_g = exp(u_g)/(\sum exp(u_g)) ----
    po <- para[1:G]
    beta <- para[(G+1):(G+p+1)]

    gamma <- para[(G+p+2):(G+2*(p+1))]

    ps <- pr*pi.s + (1 - pi.s)*po # stool sample compositions ----
    loglik1 <- dmultinom(x = y.vec, prob=ps, log=T)
    #print(loglik1)

    mu_i <- exp( beta[1] + sum(beta[-1]*x.vec) )/( 1+ exp( beta[1] + sum(beta[-1]*x.vec) )) ## glm: transformation ----
    phi_i <- exp( gamma[1] + sum(gamma[-1]*x.vec) )


    a_1i <- phi_i*mu_i
    a_2i <- phi_i*(1 - mu_i)



    loglik2 <- dbeta( x=pi.s, shape1=a_1i, shape2 = a_2i, log=T )
    #print(loglik2)

    loglik_sum <- loglik1+loglik2
    return(loglik_sum)

  }

  loglik_all <- function( pi, y,  x, para, pr){

    n <- dim(y)[1]
    G <- dim(y)[2]
    p <- dim(x)[2]

    po <- para[1:G]
    beta <- para[(G+1):(G+p+1)]
    #print(beta)
    gamma <- para[(G+p+2):(G+2*(p+1))]
    #print(gamma)

    loglik_sum <- 0

    for(i in 1:n){

      ps <- pr*pi[i] + (1 - pi[i])*po # stool sample compositions ----
      loglik1 <- dmultinom(x = y[i,], prob=ps, log=T)
      #print(loglik1)

      mu_i <- exp( beta[1] + sum(beta[-1]*x[i,]) )/( 1+ exp( beta[1] + sum(beta[-1]*x[i,]) )) ## glm: transformation ----
      phi_i <- exp( gamma[1] + sum(gamma[-1]*x[i,]) )
      a_1i <- phi_i*mu_i
      a_2i <- phi_i*(1 - mu_i)

      loglik2 <- dbeta( x=pi[i], shape1=a_1i, shape2 = a_2i, log=T )
      #print(loglik2)

      loglik_sum <- loglik_sum+loglik1+loglik2

    }

    return(loglik_sum)
  }

  E_loglik <- function(MCMC_pi, y, x, para, pr){

    n <- dim(y)[1]

    loglik_temp <- 0
    for(i in 1:n){

      #MCMC_pi.i <- as.matrix(MCMC_pi[i,])
      MCMC_pi.i <- MCMC_pi[i,]
      loglik_temp <- loglik_temp + mean(sapply(MCMC_pi.i, loglik, y.vec = y[i,], x.vec=x[i,], para=para, pr=pr))

    }

    return(-loglik_temp)

  }

  obj_EM <- function(MCMC_pi, y, x, para, pr, Const){

    temp1 <- E_loglik(MCMC_pi, y, x, para, pr) ## negative loglik
    temp2 <- Const*(sum(para[1:G]) - 1)^2
    return(temp1 + temp2)
  }

  EM_MicroDeMix <- function(y, x, para, pr, eps=1e-4, maxit=10, Const=10000){
    options(warn=-1)

    n <- dim(y)[1]
    G <- dim(y)[2]
    p <- dim(x)[2]

    para_old <- para
    diff<-1 ## difference between para_old & para_new
    iter<-1
    while(diff > eps & iter<=maxit){
      print(paste("EM iteration:",iter))

      #### Step1: sampling ----
      MCMCsample <- NULL

      for(i in 1:n){

        #if(pi.ind[i] == F){

        mcmc_loglik_try <- try(MCMCpack::MCMCmetrop1R(
          loglik,
          theta.init=0.5,
          burnin = 500,
          mcmc = 2000, ### sample size of the MCMC sample ----
          thin = 1,
          tune = 1,
          verbose = 0,
          seed = NA,
          logfun = TRUE,
          force.samp = FALSE,
          V = NULL,
          #optim.method = "L-BFGS-B",
          optim.method = "Nealder-Mead",
          optim.lower = 0.05,
          optim.upper = 0.95,
          optim.control = list(fnscale = -1, trace = 0, REPORT = 10, maxit = 500),
          y.vec = y[i,], ## stool data
          x.vec = x[i,], ## covariate data
          para = para_old,
          pr = pr ## rectum proportion estimated from the reference data
        ),silent = TRUE)

        #MCMCsample <- rbind(MCMCsample, t( mcmc_loglik_try) )

        if(class( mcmc_loglik_try) != "try-error"){
          MCMCsample <- rbind(MCMCsample, t( mcmc_loglik_try) )
        }
        else{
          beta <- para_old[(G+1):(G+p)]
          gamma <- para_old[(G+p+1):(G+2*p)]

          mu_i <- exp( beta[1] + sum(beta[-1]*x[i,]) )/( 1+ exp( beta[1] + sum(beta[-1]*x[i,]) )) ## glm: transformation ----
          phi_i <- exp( gamma[1] + sum(gamma[-1]*x[i,]) )
          a_1i <- phi_i*mu_i
          a_2i <- phi_i*(1 - mu_i)
          sample_beta <- rbeta(2000, a_1i, a_2i)
          sample_beta[sample_beta > 0.95] <- 0.95
          sample_beta[sample_beta < 0.05] <- 0.05
          MCMCsample <- rbind(MCMCsample, sample_beta )
        }

        #}


      }



      #### Step2: E-step & M-step -----
      ### use optimx: do this separately for po and other parameters ----

      lower.vec <- c(rep(0, G), rep(-Inf, 2*(p+1)))
      upper.vec <- c(rep(1,G), rep(Inf, 2*(p+1)))


      ## try together: Nelder-Mead ----
      #t2 <- optim(par=para_old, fn=E_loglik, gr=E_loglik_dev, method=c("BFGS"), control=list(maxit=10), MCMC_pi=MCMCsample, y=y, x=x, pr)
      ## 2-step update is ok ---
      # t2.bfgs <- optimx(par= para_old, fn= obj_EM, gr= NULL,
      #                   lower=lower.vec, upper=upper.vec,
      #                   method=c("L-BFGS-B"), control=list(maxit=2), MCMC_pi=MCMCsample, y=y, x=x, pr=pr, Const=10000)
      #

      t2.spg <- BB::spg(par=para_old, fn=obj_EM, gr=NULL, method=3, lower=lower.vec, upper=upper.vec,
                    project=NULL, projectArgs=NULL,
                    control=list(maxit=2), quiet=FALSE, alertConvergence=FALSE, MCMC_pi=MCMCsample, y=y, x=x, pr=pr, Const=10000)


      ## if not work, will try to BFGS with derivatives ----

      #para_new <- unlist(t2.bfgs[1:length(para)])
      para_new <- t2.spg$par[1:length(para)]
      diff <- sum(abs(para_new - para_old))
      print(diff)

      para_old <- para_new

      iter <- iter+1
    }

    return(para_old)

  }

  ###################################################################


  yig_n.type1 <- data1
  yig.n <- data2
  x <- cova1
  x.n <- cova2

  s <- nrow(yig_n.type1)
  s.r <- nrow(yig.n)

  G <- dim(yig_n.type1)[2]
  d <- dim(x)[2]

  p1.hat <- colSums(yig.n)/sum(yig.n)

  p.mle <- colSums(yig_n.type1)/sum(yig_n.type1)

  para_init <- c(p.mle, rep(0,2*(d+1))) ## reasonable starting value

  fit.EM <- EM_MicroDeMix(y=yig_n.type1, x=x, para=para_init, pr=p1.hat, maxit=maxiter)

  pc.est <- fit.EM[1:G]
  pc.est <- pc.est/sum(pc.est)

  t <- sum((p1.hat - pc.est)^2)

  data <- rbind(yig.n, yig_n.type1)
  x.all <- rbind(x.n, x)
  s.all <- nrow(data)

  t.v <- c()
  for (i in 1:B){
    print(paste("permutation iter:",i))
    set.seed(100+i)

    ind <- sample(1:s.all)

    data.new <- data[ind, ]
    x.new <- as.matrix(x.all[ind, ])

    yig.n_perm <- data.new[1:s.r, ]
    yig_perm <- data.new[(s.r+1):s.all, ]

    x_perm <- as.matrix(x.new[(s.r+1):s.all, ])

    p1.perm <- colSums(yig.n_perm)/sum(yig.n_perm)
    pmle.perm <- colSums(yig_perm)/sum(yig_perm)

    para_init <- c(pmle.perm, rep(0,2*(d+1))) ## reasonable starting value

    fit.perm <- EM_MicroDeMix(y=yig_perm, x=x_perm, para=para_init, pr=p1.perm, maxit=maxiter)

    pc.perm <- fit.perm[1:G]
    pc.perm <- pc.perm/sum(pc.perm)

    t.perm <- sum((pc.perm-p1.perm)^2)
    t.v <- append(t.v, t.perm)

  }

  pval <- mean(t.v >= t)
  return (pval)
}


