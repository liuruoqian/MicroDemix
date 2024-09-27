# MicroDemix
R package for Micro-Demix, a method that decomposes the stool microbiome at compositional level to understand the
heterogeneity of the gut microbiome across various GI locations. In this section, we introduce the installation process of 'MicroDemix' and demonstrate the features of this R package using example datasets. 
## Installation
### Install package "MicroDemix"
```
library(devtools)
install_github("liuruoqian/MicroDemix", build_vignettes = TRUE)
library(MicroDemix)
```
### Load the following R packages
```
library(optimx)
library(cubature)
library(rootSolve)
library(MASS)
library(ggplot2)
library(BB)
library(MCMCpack)
```

## Example datasets in "MicroDemix"
The example datasets are from the simulation study with 10 taxa and 100 samples.
### 'yig': A stool (mixture) microbiome dataset
‘yig’ is an example dataset for stool (mixture) microbiome data with 100 samples in rows and 10 taxa in columns. We can access this dataset by running the following:
```
data(yig)
head(yig)
#>      Taxon1 Taxon2 Taxon3 Taxon4 Taxon5 Taxon6 Taxon7 Taxon8 Taxon9 Taxon10
#> [1,]    985   1245   1571   2030   2517   3192   4099   5201   6707    2453
#> [2,]    944   1234   1549   2083   2498   3096   4194   5201   6776    2425
#> [3,]    940   1218   1599   1953   2539   3210   4116   5211   6718    2496
#> [4,]    949   1277   1561   2029   2522   3161   4006   5316   6695    2484
#> [5,]    979   1235   1560   2039   2559   3202   4029   5284   6660    2453
#> [6,]    969   1218   1535   1994   2543   3221   3933   5277   6734    2576
```
### 'yig.n': A rectum (reference) microbiome dataset
‘yig.n’ is an example dataset for rectum (reference) microbiome data with 100 samples in rows and 10 taxa in columns. We can access this dataset by running the following:
```
data(yig.n)
head(yig.n)
#>      Taxon1 Taxon2 Taxon3 Taxon4 Taxon5 Taxon6 Taxon7 Taxon8 Taxon9 Taxon10
#> [1,]    999   1219   1678   2040   2583   3094   4031   5285   6649    2422
#> [2,]    980   1222   1587   1975   2488   3026   4022   5407   6805    2488
#> [3,]    964   1224   1555   2016   2623   3191   4066   5340   6566    2455
#> [4,]    978   1201   1626   2041   2670   2994   4066   5321   6614    2489
#> [5,]    992   1242   1558   2066   2538   3126   4137   5204   6697    2440
#> [6,]    928   1231   1581   2090   2568   3247   4058   5276   6625    2396
```
### 'x': A covariates dataset related to samples in data 'yig' 
‘x’ is a dataset contains covariates corresponding to samples in data ‘yig’. We can access this dataset by running the following:
```
data(x)
dim(x)
#> [1] 100   1
```
### ‘x.n’: A covariates dataset related to samples in data ‘yig.n’
‘x.n’ is a dataset contains covariates corresponding to samples in data ‘yig.n’. We can access this dataset by running the following:
```
data(x.n)
dim(x.n)
#> [1] 100   1
```
## Functions and data analysis with "MicroDemix"
### Function 1: MD_estimate(data1, data2, cova)
This function estimates the microbial relative abundance in stool (mixture), rectum (reference) and other GI locations. For example: 
```
est1 <- MD_estimate(yig, yig.n, x)
est1
#> $pr.est
#>     Taxon1     Taxon2     Taxon3     Taxon4     Taxon5     Taxon6     Taxon7 
#> 0.03272567 0.04125933 0.05265133 0.06706233 0.08539033 0.10519300 0.13534633 
#>     Taxon8     Taxon9    Taxon10 
#> 0.17478833 0.22439333 0.08119000 
#> 
#> $ps.est
#>     Taxon1     Taxon2     Taxon3     Taxon4     Taxon5     Taxon6     Taxon7 
#> 0.03213167 0.04096100 0.05220233 0.06671833 0.08465967 0.10574300 0.13572900 
#>     Taxon8     Taxon9    Taxon10 
#> 0.17467567 0.22533033 0.08184900 
#> 
#> $po.est
#>         p1         p2         p3         p4         p5         p6         p7 
#> 0.03158059 0.04068860 0.05170581 0.06653634 0.08392623 0.10626569 0.13600713 
#>         p8         p9        p10 
#> 0.17470683 0.22604926 0.08253354
```
### Function 2: RA_plot(p.r, p.s, p.o, G, taxon)
We can visualize microbial relative abundances (proportions) at three different locations using the ‘RA_plot()’ function. 
The following example replicates one of the plots from the paper, analyzing the phylum Firmicutes at the family level.
```
p.r <- c(0, 0, 0.0195415822, 0.52043881, 0.0002028398, 0, 0.3991845842, 0, 0.0606321839)
p.s <- c(0.02713150778, 0.00075219743, 0.00073529412, 0.43284820825, 0.00010311021, 0.00008958756, 0.53113421231, 0.00650608519, 0.00069979716)
p.o <- c(0.0571234686474, 0.0014037570334, 0.0000033854869, 0.122693737055, 0.0000030544905, 0.0000021475416, 0.8050093482991, 0.0137601784493, 0.0000009229968)
families <- c("Clostridiaceae", "Enterococcaceae",    "Erysipelotrichaceae", "Lachnospiraceae", "Peptococcaceae",  "Planococcaceae",  "Ruminococcaceae",    "Turicibacteraceae", "Veillonellaceae")
RA_plot(p.r, p.s, p.o, 9, families)
```
![image](https://github.com/liuruoqian/MicroDemix/assets/25522356/7e0c72a6-8cb1-4873-bbe1-e2f804d20282)
### Function 3: MD_pvalue(p.r, p.o, N.r, N)
We obtain a simulation-based p-value for testing differential abundance under the null hypothesis that p.r = p.o, where p.r and p.o are the microbial relative abundance 
in the rectum and other GI locations, respectively. 
```
MD_pvalue(est1$pr.est, est1$po.est, sum(yig.n), sum(yig))
[1] 0
```
We can obtain a p-value under other circumstances, for example:
```
p1 <- rep(0.1, 10)
p2 <- c(rep(0.09, 5), rep(0.11, 5))
MD_pvalue(p1, p2, 5000, 5000)
[1] 0.0031
```
### Function 4: MicroDemix_EM(data1, data2, cova, maxiter=10)
This function using the EM algorithm to estimate the microbial relative abundance in stool (mixture), rectum (reference) and other GI locations. 
The EM algorithm is more computationally intensive due to its relatively slow convergence rate and the repeated use of the Metropolis-Hastings (MH) algorithm in the E-step. In the following example, each EM iteration takes approximately 2 minutes on a local machine. However, the estimation accuracy improves significantly after just a few iterations.
```
system.time(est2 <- MicroDemix_EM(yig, yig.n, x, maxiter=5))
#> [1] "EM iteration: 1"
#> iter:  0  f-value:  4704.532  pgrad:  281.5245 
#> [1] 0.002658093
#> [1] "EM iteration: 2"
#> iter:  0  f-value:  4701.743  pgrad:  277.8827 
#> [1] 0.0008061163
#> [1] "EM iteration: 3"
#> iter:  0  f-value:  4701.414  pgrad:  276.3458 
#> [1] 0.0002402273
#> [1] "EM iteration: 4"
#> iter:  0  f-value:  4701.446  pgrad:  276.3473 
#> [1] 0.0004213928
#> [1] "EM iteration: 5"
#> iter:  0  f-value:  4701.244  pgrad:  274.4289 
#> [1] 0.0001306647
#>    user  system elapsed 
#>  531.33    0.98  940.91

est2
#> $pr.est
#>     Taxon1     Taxon2     Taxon3     Taxon4     Taxon5     Taxon6     Taxon7 
#> 0.03272567 0.04125933 0.05265133 0.06706233 0.08539033 0.10519300 0.13534633 
#>     Taxon8     Taxon9    Taxon10 
#> 0.17478833 0.22439333 0.08119000 
#> 
#> $ps.est
#>     Taxon1     Taxon2     Taxon3     Taxon4     Taxon5     Taxon6     Taxon7 
#> 0.03213167 0.04096100 0.05220233 0.06671833 0.08465967 0.10574300 0.13572900 
#>     Taxon8     Taxon9    Taxon10 
#> 0.17467567 0.22533033 0.08184900 
#> 
#> $po.est
#>     Taxon1     Taxon2     Taxon3     Taxon4     Taxon5     Taxon6     Taxon7 
#> 0.03191630 0.04084895 0.05203789 0.06658997 0.08438937 0.10594348 0.13587368 
#>     Taxon8     Taxon9    Taxon10 
#> 0.17463619 0.22567804 0.08208611
```
### Function 5: MicroDemix_EM_permutation(data1, data2, cova1, cova2, maxiter=10, B)
This function performs a permutation based test of differential abundance under the null hypothesis that p.r = p.o. 
The permutation test requires multiple runs of the EM algorithm. Therefore, we recommend performing this procedure using parallel computing or a high-performance machine.
```
MicroDemix_EM_permutation(yig, yig.n, x, x.n, maxiter=2, B=10)

# [1] "EM iteration: 1"
# iter:  0  f-value:  4704.429  pgrad:  281.1108 
# [1] 0.002657345
# [1] "EM iteration: 2"
# iter:  0  f-value:  4701.638  pgrad:  277.4302 
# [1] 0.0008550884
# [1] "permutation iter: 1"
# [1] "EM iteration: 1"
# iter:  0  f-value:  4688.334  pgrad:  279.2411 
# [1] 0.0008837178
# [1] "EM iteration: 2"
# iter:  0  f-value:  4688.151  pgrad:  279.6902 
# [1] 0.0003699496
# [1] "permutation iter: 2"
# [1] "EM iteration: 1"
# iter:  0  f-value:  4680.651  pgrad:  280.0751 
# [1] 0.0006625398
# [1] "EM iteration: 2"
# iter:  0  f-value:  4680.5  pgrad:  279.7522 
# [1] 0.0002500549
# [1] "permutation iter: 3"
# [1] "EM iteration: 1"
# iter:  0  f-value:  4704.84  pgrad:  280.6383 
# [1] 0.0008771262
# [1] "EM iteration: 2"
# iter:  0  f-value:  4704.655  pgrad:  280.7079 
# [1] 0.0003814814
# [1] "permutation iter: 4"
# [1] "EM iteration: 1"
# iter:  0  f-value:  4695.358  pgrad:  280.5077 
# [1] 0.001261839
# [1] "EM iteration: 2"
# iter:  0  f-value:  4694.589  pgrad:  280.0775 
# [1] 0.0004829247
# [1] "permutation iter: 5"
# [1] "EM iteration: 1"
# iter:  0  f-value:  4714.963  pgrad:  280.0942 
# [1] 0.0006955507
# [1] "EM iteration: 2"
# iter:  0  f-value:  4714.625  pgrad:  278.8059 
# [1] 0.0002040188
# [1] "permutation iter: 6"
# [1] "EM iteration: 1"
# iter:  0  f-value:  4692.693  pgrad:  279.8199 
# [1] 0.0007145038
# [1] "EM iteration: 2"
# iter:  0  f-value:  4692.438  pgrad:  279.3542 
# [1] 0.0004163259
# [1] "permutation iter: 7"
# [1] "EM iteration: 1"
# iter:  0  f-value:  4701.11  pgrad:  281.4275 
# [1] 0.001500689
# [1] "EM iteration: 2"
# iter:  0  f-value:  4699.981  pgrad:  278.8331 
# [1] 0.0002896197
# [1] "permutation iter: 8"
# [1] "EM iteration: 1"
# iter:  0  f-value:  4688.652  pgrad:  278.8074 
# [1] 0.001133446
# [1] "EM iteration: 2"
# iter:  0  f-value:  4688.138  pgrad:  280.7433 
# [1] 0.0005703671
# [1] "permutation iter: 9"
# [1] "EM iteration: 1"
# iter:  0  f-value:  4708.529  pgrad:  280.8886 
# [1] 0.001269645
# [1] "EM iteration: 2"
# iter:  0  f-value:  4707.73  pgrad:  278.8487 
# [1] 0.0002213754
# [1] "permutation iter: 10"
# [1] "EM iteration: 1"
# iter:  0  f-value:  4712.357  pgrad:  280.8193 
# [1] 0.0005794702
# [1] "EM iteration: 2"
# iter:  0  f-value:  4712.134  pgrad:  279.7022 
# [1] 0.0001730983
# [1] 0  #pvalue
```

