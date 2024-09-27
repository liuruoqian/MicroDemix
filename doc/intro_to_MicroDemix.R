## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  warning = FALSE,
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(MicroDemix)

## -----------------------------------------------------------------------------
library(optimx)
library(cubature)
library(rootSolve)
library(MASS)
library(ggplot2)
library(BB)
library(MCMCpack)

## -----------------------------------------------------------------------------
data(yig)
head(yig)

## -----------------------------------------------------------------------------
data(yig.n)
head(yig.n)

## -----------------------------------------------------------------------------
data(x)
dim(x)

## -----------------------------------------------------------------------------
data(x.n)
dim(x.n)

## -----------------------------------------------------------------------------
est1 <- MD_estimate(yig, yig.n, x)
est1

## ----fig.width=7.5, fig.height=4.5--------------------------------------------
p.r <- c(0, 0, 0.0195415822, 0.52043881, 0.0002028398, 0, 0.3991845842, 0, 0.0606321839)
p.s <- c(0.02713150778, 0.00075219743, 0.00073529412, 0.43284820825, 0.00010311021, 0.00008958756, 0.53113421231, 0.00650608519, 0.00069979716)
p.o <- c(0.0571234686474, 0.0014037570334, 0.0000033854869, 0.122693737055, 0.0000030544905, 0.0000021475416, 0.8050093482991, 0.0137601784493, 0.0000009229968)
families <- c("Clostridiaceae", "Enterococcaceae",    "Erysipelotrichaceae", "Lachnospiraceae", "Peptococcaceae",  "Planococcaceae",  "Ruminococcaceae",    "Turicibacteraceae", "Veillonellaceae")
RA_plot(p.r, p.s, p.o, 9, families)

## -----------------------------------------------------------------------------
MD_pvalue(est1$pr.est, est1$po.est, sum(yig.n), sum(yig))

## -----------------------------------------------------------------------------
p1 <- rep(0.1, 10)
p2 <- c(rep(0.09, 5), rep(0.11, 5))
MD_pvalue(p1, p2, 5000, 5000)

## ----eval=FALSE---------------------------------------------------------------
#  system.time(est2 <- MicroDemix_EM(yig, yig.n, x, maxiter=5))
#  #> [1] "EM iteration: 1"
#  #> iter:  0  f-value:  4704.532  pgrad:  281.5245
#  #> [1] 0.002658093
#  #> [1] "EM iteration: 2"
#  #> iter:  0  f-value:  4701.743  pgrad:  277.8827
#  #> [1] 0.0008061163
#  #> [1] "EM iteration: 3"
#  #> iter:  0  f-value:  4701.414  pgrad:  276.3458
#  #> [1] 0.0002402273
#  #> [1] "EM iteration: 4"
#  #> iter:  0  f-value:  4701.446  pgrad:  276.3473
#  #> [1] 0.0004213928
#  #> [1] "EM iteration: 5"
#  #> iter:  0  f-value:  4701.244  pgrad:  274.4289
#  #> [1] 0.0001306647
#  #>    user  system elapsed
#  #>  531.33    0.98  940.91
#  
#  est2
#  #> $pr.est
#  #>     Taxon1     Taxon2     Taxon3     Taxon4     Taxon5     Taxon6     Taxon7
#  #> 0.03272567 0.04125933 0.05265133 0.06706233 0.08539033 0.10519300 0.13534633
#  #>     Taxon8     Taxon9    Taxon10
#  #> 0.17478833 0.22439333 0.08119000
#  #>
#  #> $ps.est
#  #>     Taxon1     Taxon2     Taxon3     Taxon4     Taxon5     Taxon6     Taxon7
#  #> 0.03213167 0.04096100 0.05220233 0.06671833 0.08465967 0.10574300 0.13572900
#  #>     Taxon8     Taxon9    Taxon10
#  #> 0.17467567 0.22533033 0.08184900
#  #>
#  #> $po.est
#  #>     Taxon1     Taxon2     Taxon3     Taxon4     Taxon5     Taxon6     Taxon7
#  #> 0.03191630 0.04084895 0.05203789 0.06658997 0.08438937 0.10594348 0.13587368
#  #>     Taxon8     Taxon9    Taxon10
#  #> 0.17463619 0.22567804 0.08208611

## ----eval=FALSE---------------------------------------------------------------
#  MicroDemix_EM_permutation(yig, yig.n, x, x.n, maxiter=2, B=10)
#  
#  # [1] "EM iteration: 1"
#  # iter:  0  f-value:  4704.429  pgrad:  281.1108
#  # [1] 0.002657345
#  # [1] "EM iteration: 2"
#  # iter:  0  f-value:  4701.638  pgrad:  277.4302
#  # [1] 0.0008550884
#  # [1] "permutation iter: 1"
#  # [1] "EM iteration: 1"
#  # iter:  0  f-value:  4688.334  pgrad:  279.2411
#  # [1] 0.0008837178
#  # [1] "EM iteration: 2"
#  # iter:  0  f-value:  4688.151  pgrad:  279.6902
#  # [1] 0.0003699496
#  # [1] "permutation iter: 2"
#  # [1] "EM iteration: 1"
#  # iter:  0  f-value:  4680.651  pgrad:  280.0751
#  # [1] 0.0006625398
#  # [1] "EM iteration: 2"
#  # iter:  0  f-value:  4680.5  pgrad:  279.7522
#  # [1] 0.0002500549
#  # [1] "permutation iter: 3"
#  # [1] "EM iteration: 1"
#  # iter:  0  f-value:  4704.84  pgrad:  280.6383
#  # [1] 0.0008771262
#  # [1] "EM iteration: 2"
#  # iter:  0  f-value:  4704.655  pgrad:  280.7079
#  # [1] 0.0003814814
#  # [1] "permutation iter: 4"
#  # [1] "EM iteration: 1"
#  # iter:  0  f-value:  4695.358  pgrad:  280.5077
#  # [1] 0.001261839
#  # [1] "EM iteration: 2"
#  # iter:  0  f-value:  4694.589  pgrad:  280.0775
#  # [1] 0.0004829247
#  # [1] "permutation iter: 5"
#  # [1] "EM iteration: 1"
#  # iter:  0  f-value:  4714.963  pgrad:  280.0942
#  # [1] 0.0006955507
#  # [1] "EM iteration: 2"
#  # iter:  0  f-value:  4714.625  pgrad:  278.8059
#  # [1] 0.0002040188
#  # [1] "permutation iter: 6"
#  # [1] "EM iteration: 1"
#  # iter:  0  f-value:  4692.693  pgrad:  279.8199
#  # [1] 0.0007145038
#  # [1] "EM iteration: 2"
#  # iter:  0  f-value:  4692.438  pgrad:  279.3542
#  # [1] 0.0004163259
#  # [1] "permutation iter: 7"
#  # [1] "EM iteration: 1"
#  # iter:  0  f-value:  4701.11  pgrad:  281.4275
#  # [1] 0.001500689
#  # [1] "EM iteration: 2"
#  # iter:  0  f-value:  4699.981  pgrad:  278.8331
#  # [1] 0.0002896197
#  # [1] "permutation iter: 8"
#  # [1] "EM iteration: 1"
#  # iter:  0  f-value:  4688.652  pgrad:  278.8074
#  # [1] 0.001133446
#  # [1] "EM iteration: 2"
#  # iter:  0  f-value:  4688.138  pgrad:  280.7433
#  # [1] 0.0005703671
#  # [1] "permutation iter: 9"
#  # [1] "EM iteration: 1"
#  # iter:  0  f-value:  4708.529  pgrad:  280.8886
#  # [1] 0.001269645
#  # [1] "EM iteration: 2"
#  # iter:  0  f-value:  4707.73  pgrad:  278.8487
#  # [1] 0.0002213754
#  # [1] "permutation iter: 10"
#  # [1] "EM iteration: 1"
#  # iter:  0  f-value:  4712.357  pgrad:  280.8193
#  # [1] 0.0005794702
#  # [1] "EM iteration: 2"
#  # iter:  0  f-value:  4712.134  pgrad:  279.7022
#  # [1] 0.0001730983
#  # [1] 0  #pvalue

