# MicroDemix
R package for Micro-Demix, a method that decomposes the stool microbiome at compositional level to understand the
heterogeneity of the gut microbiome across various GI locations. In this section, we introduce the installation process of 'MicroDemix' and demonstrate the features of this R package using example datasets. 
## Installation
### Install and load the following R packages
```
library(optimx)
library(cubature)
library(rootSolve)
library(MASS)
library(ggplot2)
```
### Install package "MicroDemix"
```
library(devtools)
install_github("liuruoqian/MicroDemix")
library(MicroDemix)
```
## Example datasets in "MicroDemix"

### 'yig': A stool (mixture) microbiome dataset
'yig' is an example dataset for stool (mixture) microbiome data with 12 samples in rows and 9 taxa (families) in columns. 
We can access this dataset by running the following:
```
data(yig)
# show the column (taxa) names 
colnames(yig)
[1] "Clostridiaceae"      "Enterococcaceae"     "Erysipelotrichaceae" "Lachnospiraceae"     "Peptococcaceae"     
[6] "Planococcaceae"      "Ruminococcaceae"     "Turicibacteraceae"   "Veillonellaceae"   
```
### 'yig.n': A rectum (reference) microbiome dataset
'yig.n' is an example dataset for rectum (reference) microbiome data with 30 samples in rows and 9 taxa (families) in columns. 
We can access this dataset by running the following:
```
data(yig.n)
# show the column (taxa) names 
colnames(yig.n)
[1] "Clostridiaceae"      "Enterococcaceae"     "Erysipelotrichaceae" "Lachnospiraceae"     "Peptococcaceae"     
[6] "Planococcaceae"      "Ruminococcaceae"     "Turicibacteraceae"   "Veillonellaceae"   
```
### 'x': A covariates dataset related to samples in data 'yig' 
'x' is a dataset contains information on consent_age and gender corresponding to samples in data 'yig'.
We can access this dataset by running the following:
```
data(x)
colnames(x)
[1] "consent_age" "gender" 
```
## Functions and data analysis with "MicroDemix"
### Function 1: MD_estimate(data1, data2, cova)
This function estimates the microbial relative abundance in stool (mixture), rectum (reference) and other GI locations. For example: 
```
est1 <- MD_estimate(yig, yig.n, x)
est1
$pr.est
     Clostridiaceae     Enterococcaceae Erysipelotrichaceae     Lachnospiraceae      Peptococcaceae 
       0.0000000000        0.0000000000        0.0195415822        0.5204388100        0.0002028398 
     Planococcaceae     Ruminococcaceae   Turicibacteraceae     Veillonellaceae 
       0.0000000000        0.3991845842        0.0000000000        0.0606321839 

$ps.est
     Clostridiaceae     Enterococcaceae Erysipelotrichaceae     Lachnospiraceae      Peptococcaceae 
       2.713151e-02        7.521974e-04        7.352941e-04        4.328482e-01        1.031102e-04 
     Planococcaceae     Ruminococcaceae   Turicibacteraceae     Veillonellaceae 
       8.958756e-05        5.311342e-01        6.506085e-03        6.997972e-04 

$po.est
          p1           p2           p3           p4           p5           p6           p7 
5.712347e-02 1.403757e-03 3.385487e-06 1.226937e-01 3.054491e-06 2.147542e-06 8.050093e-01 
          p8           p9 
1.376018e-02 9.229968e-07 
```
### Function 2: RA_plot(p.r, p.s, p.o, G, taxon)
We can visualize microbial relative abundance (proportions) at three different locations using RA_plot(). For example:
```
RA_plot(est1$pr.est, est1$ps.est, est1$po.est, 9, colnames(yig))
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
