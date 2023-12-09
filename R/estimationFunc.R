
#' MD_estimate
#'
#' This function estimate the microbial relative abundance in stool (mixture),
#' rectum (reference) and other GI locations.
#' @param data1 stool (mixture) microbiome data set with taxa in
#'            columns and samples in rows.
#' @param data2 rectum (reference) microbiome data set with taxa in
#'            columns and samples in rows. The taxa in data1 and data2 should
#'            be matching.
#' @param cova covariates related to samples in data1 that have a matching order
#' with data1. 'cova' must be a matrix with covariates in columns and observations
#' in rows.
#' @return A list containing the following components:
#' \item{pr.est}{the estimated relative abundance for the rectum (reference) microbiome}
#' \item{ps.est}{the estimated relative abundance for the stool (mixture) microbiome}
#' \item{po.est}{the estimated relative abundance for the microbiome in other GI locations}
#' @examples
#' MD_estimate(yig, yig.n, x);
#' @export


MD_estimate <- function(data1, data2, cova){
  # source("R/baseFunc.R")

  func_mu <- function(xi, b0_mu, b_mu){
    return (exp(b0_mu + t(xi)%*%b_mu)/(1+exp(b0_mu + t(xi)%*%b_mu)))
  }

  func_phi <- function(xi, b0_phi, b_phi){
    return (exp(b0_phi + t(xi)%*%b_phi))
  }

  samp_multi <- function(p){
    return (rmultinom(1, counts, prob = p))
  }

  func_d <- function(yi,pi, p1,u_v,a1,a2){

    # u_G <- -sum(u)
    # u_v <- c(u, u_G) # u vector with G elements

    pc <- exp(u_v)/sum(exp(u_v))
    #d <- prod((pi*p1 + (1-pi)*pc)^yi) * (pi^(a1-1)*(1-pi)^(a2-1))/beta(a1,a2)
    prob2 <- pi*p1 + (1-pi)*pc
    d <- dmultinom(yi, prob=prob2) * (pi^(a1-1)*(1-pi)^(a2-1))/beta(a1,a2)
    return (d)
  }

  # integrate density over pi for each i
  # change p1=p1 to p1=p1.hat on 6/12/2022
  func_inte <- function(u_inte,yig,a1_v,a2_v){
    inte <- cubintegrate(f=func_d, lower=0, upper=1,
                         yi=yig, p1=p1.hat, u=u_inte, a1=a1_v, a2=a2_v)$integral
    return(inte)
  }


  u_to_pc <- function(u_vec){
    return (exp(u_vec)/sum(exp(u_vec)))
  }


  optFunc_beta <- function(theta, u_em, yig, x1){
    u_para <- u_em
    b0_mu1 <- theta[1]
    b_mu1 <- theta[(2):(1+d)]
    b0_phi1 <- theta[2+d]
    b_phi1 <- theta[(3+d):length(theta)]

    mu_v <- apply(x1,1,func_mu, b0_mu=b0_mu1, b_mu=b_mu1)
    phi_v <- apply(x1,1,func_phi, b0_phi=b0_phi1, b_phi=b_phi1)

    a1_v <- mu_v*phi_v
    a2_v <- (1-mu_v)*phi_v

    u_mtx <- matrix(rep(u_para,s.c),G,s.c) # change G-1 to G on 1/7/2023

    # return negative loglikelihood
    return (-sum(log(mapply(func_inte, as.list(data.frame(u_mtx)),
                            as.list(data.frame(t(yig))), as.list(a1_v),as.list(a2_v)))))
  }

  optFunc_pc <- function(u, theta, yig, x1){
    u_para <- u
    b0_mu1 <- theta[1]
    b_mu1 <- theta[(2):(1+d)]
    b0_phi1 <- theta[2+d]
    b_phi1 <- theta[(3+d):length(theta)]

    mu_v <- apply(x1,1,func_mu, b0_mu=b0_mu1, b_mu=b_mu1)
    phi_v <- apply(x1,1,func_phi, b0_phi=b0_phi1, b_phi=b_phi1)

    a1_v <- mu_v*phi_v
    a2_v <- (1-mu_v)*phi_v

    s <- nrow(yig)  # revised on 8/25
    u_mtx <- matrix(rep(u_para,s),G, s)

    # return negative loglikelihood
    return (-sum(log(mapply(func_inte, as.list(data.frame(u_mtx)),
                            as.list(data.frame(t(yig))), as.list(a1_v),as.list(a2_v)))))
  }

  # -----------------------------------------------------------
  yig <- data1
  yig.n <- data2
  x<- cova

  p1.hat <- colSums(yig.n)/sum(yig.n)
  p.mle <- colSums(yig)/sum(yig)

  G <- ncol(yig)

  model <- function(x) {
    F1 <- numeric(G)
    for(i in 1:(G-1)){
      F1[i] <- exp(x[i])/sum(exp(x)) - p.mle[i]
    }
    F1[G] <- sum(x)
    c(F1 = F1)
  }

  u.hat <- multiroot(f = model, start = rep(0,G))$root


  s.n <- nrow(yig.n)
  s.c <- nrow(yig)

  d <- dim(x)[2]

  # u_para <- u.hat
  #
  # mu_v <- apply(x,1,func_mu, b0_mu=0, b_mu=c(0,0))
  #
  # phi_v <- apply(x,1,func_phi, b0_phi=0, b_phi=c(0,0))
  #
  # a1_v <- mu_v*phi_v
  # a2_v <- (1-mu_v)*phi_v
  #
  # u_mtx <- matrix(rep(u_para,s.c),G,s.c)
  # p1_mtx <- matrix(rep(p1.hat,s.c),G,s.c)
  #
  #
  # which(log(mapply(func_inte, as.list(data.frame(u_mtx)),
  #                  as.list(data.frame(t(yig))), as.list(a1_v),as.list(a2_v))) == -Inf)

  b <- suppressWarnings(optimx(par= u.hat,fn=optFunc_pc, theta = rep(0, (2+2*d)), yig=yig, x1=x,
               method = "Nelder-Mead",
              control = list(maxit=100)))

  b1 <- b
  LL1 <- b1$value
  LL_v <- rep(NA,31)
  LL_v[1] <- LL1

  b_init <- rep(0, (2+2*d))
  u_init <- u.hat
  for (iter in 1:30){


    b2 <- suppressWarnings(optimx(par= b_init,fn=optFunc_beta, u_em= unlist(b1[1:G]), yig=yig,
                 x1=x, method = "Nelder-Mead", control = list(maxit=100)))
    b_init <- unlist(b2[1:(2+2*d)])


    b1 <- suppressWarnings(optimx(par= u_init, fn=optFunc_pc, theta = unlist(b2[1:(2+2*d)]), yig=yig, x1=x,
                 method = "Nelder-Mead",
                 control = list(maxit=100)))
    u_init <- unlist(b1[1:G])

    if (abs(LL1 - b1$value) <= 0.001) break

    LL1 <- b1$value
    LL_v[iter+1] <- LL1
  }

  pc.hat <- u_to_pc(unlist(b1[1:G]))
  results <- list()
  results$pr.est <- p1.hat
  results$ps.est <- p.mle
  results$po.est <- pc.hat

  return(results)


}

#' RA_plot
#'
#' visualization of microbial relative abundance (proportions) at three different locations.
#' @param p.r relative abundance at location 1 (e.g. rectum)
#' @param p.s relative abundance at location 2 (e.g. stool)
#' @param p.o relative abundance at location 3 (e.g. other GI locations).
#' The taxa in p.r, p.s, and p.o should be matching.
#' @param G the number of taxa in p.r
#' @param taxon the names of taxa in p.r
#' @return A stack plot
#' @examples
#' p1 <- rep(0.2, 5)
#' p2 <- c(0.1, 0.2, 0.25, 0.25, 0.2)
#' p3 <- c(rep(0.1, 3),rep(0.35, 2))
#' RA_plot(p1,p2,p3,5,c("A","B","C","D","E"));
#' @export

RA_plot <- function(p.r, p.s, p.o, G, taxon){

  loc_phy <- rep(c(1,2,3),each=G)  # x Axis
  phy <- c(p.r, p.s, p.o)
  Taxon <- rep(taxon,3)
  data_phy <- data.frame(loc_phy, phy, Taxon)

  stack_phy <- ggplot(data_phy, aes(x=loc_phy, y=phy, fill=Taxon)) +
    geom_area()

  stack_phy <- stack_phy + scale_x_continuous(breaks = c(1,2,3),
                                              labels=c('proportion1',  'proportion2','proportion3'))+
    # ggtitle("Phylum-level Analysis") +
    theme(plot.title = element_text( size=12),
          axis.title.x=element_blank(),
          axis.title.y = element_text(size=10),
          legend.text = element_text(size=10))+
    labs(y = "Relative Abundance")


  return (stack_phy)

}
