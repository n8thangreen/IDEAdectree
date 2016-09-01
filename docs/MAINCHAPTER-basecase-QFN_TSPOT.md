# IDEA: Main Chapter of HTA Report Cost-Effectiveness Comparison
Nathan Green  
4 February 2016  

<!-- README.md is generated from README.Rmd. Please edit that file -->





```r
# source("../../../analysis scripts/IDEA/alt-YAML_Binomial_dectrees/indiv-dectree-sampling.R")
```



```r
library(IDEAdectree)
library(BCEA)
library(ggplot2)

# load("C:/Users/ngreen1/Dropbox/TB/IDEA/R/packages/IDEAdectree/data/TBdata_clinical_cleaned.RData")
load("../data/TBdata_clinical_cleaned.RData")
load("../data/COSTdistns_allerror.RData")
load("../data/senspec_env.RData")
load("../data/drug_dose-cost.RData")

## sensitivities and specificities from IDEA lab data
attach(senspec.env)

dat <- list()

yearindays <- 365
WTP <- c(20000, 30000)/yearindays
```


```r
##TODO##
## doesnt work at the moment...!!


dat1 <- IDEAdectree.simple(data=data, name.ruleout = "TSPOT",
                           SPEC = TSPOT.noIndet.spec.mean, SENS = TSPOT.noIndet.sens.mean, 
                           SPECvar = TSPOT.noIndet.spec.var, SENSvar = TSPOT.noIndet.sens.var, wholecohortstats = TRUE)
dat2 <- IDEAdectree.simple(data=data, name.ruleout = "TSPOT",
                           SPEC = TSPOT.HIV.noIndet.spec.mean, SENS = TSPOT.HIV.noIndet.sens.mean, 
                           SPECvar = TSPOT.HIV.noIndet.spec.var, SENSvar = TSPOT.HIV.noIndet.sens.var, wholecohortstats = TRUE)
dat3 <- IDEAdectree.simple(data=data, name.ruleout = "TSPOT",
                           SPEC = TSPOT.HIVneg.noIndet.spec.mean, SENS = TSPOT.HIVneg.noIndet.sens.mean, 
                           SPECvar = TSPOT.HIVneg.noIndet.spec.var, SENSvar = TSPOT.HIVneg.noIndet.sens.var, wholecohortstats = TRUE)

dat$e <- cbind(dat1$e, dat2$e[,2], dat3$e[,2])
dat$c <- cbind(dat1$c, dat2$c[,2], dat3$c[,2])

intlabels <- c("Current",
               "Enhanced TSPOT: All, no Indets", "Enhanced TSPOT: HIV positive, no Indets", "Enhanced TSPOT: HIV negative, no Indets")

my.plot.bcea(dat1, dat2, dat3, wtp=WTP, intlabels = intlabels, wtpNEG = "Y")
my.plot.bcea(dat1, dat2, dat3, wtp=WTP*yearindays, intlabels = c("All","HIV positive","HIV negative"), wtpNEG = "Y", SCALEcosts=FALSE, SCALEdays=FALSE, labelLong=FALSE, TITLE="(a)", N=769)

my.plot.bcea(dat1, dat2, dat3, wtp=WTP, intlabels = intlabels, contour=TRUE, wtpNEG = "Y")
my.plot.bcea(dat1, dat2, dat3, wtp=WTP*yearindays, intlabels = c("All","HIV positive","HIV negative"), contour=TRUE, wtpNEG = "Y", SCALEcosts=FALSE, SCALEdays=FALSE, labelLong=FALSE, TITLE="(a)", N=769)

my.plot.bcea(dat1, dat2, dat3, wtp=WTP, intlabels = intlabels, contour=TRUE, LEVELS=0.5, wtpNEG = "Y")
my.plot.bcea(dat1, dat2, dat3, wtp=WTP*yearindays, intlabels = c("All","HIV positive","HIV negative"), contour=TRUE, LEVELS=0.5, wtpNEG = "Y", SCALEcosts=FALSE, SCALEdays=FALSE, labelLong=FALSE, TITLE="(a)", N=769)

my.plot.ceac(dat1, dat2, dat3, intlabels = intlabels)
```



```r
dat1 <- IDEAdectree.simple(data=data, name.ruleout = "TSPOT",
                           SPEC = TSPOT.noIndet.spec.mean, SENS = TSPOT.noIndet.sens.mean, 
                           SPECvar = TSPOT.noIndet.spec.var, SENSvar = TSPOT.noIndet.sens.var)
```

```
## Loading required package: triangle
```

```
## Loading required package: assertive
```

```r
dat2 <- IDEAdectree.simple(data=data, name.ruleout = "TSPOT",
                           SPEC = TSPOT.HIV.noIndet.spec.mean, SENS = TSPOT.HIV.noIndet.sens.mean, 
                           SPECvar = TSPOT.HIV.noIndet.spec.var, SENSvar = TSPOT.HIV.noIndet.sens.var)
dat3 <- IDEAdectree.simple(data=data, name.ruleout = "TSPOT",
                           SPEC = TSPOT.HIVneg.noIndet.spec.mean, SENS = TSPOT.HIVneg.noIndet.sens.mean, 
                           SPECvar = TSPOT.HIVneg.noIndet.spec.var, SENSvar = TSPOT.HIVneg.noIndet.sens.var)

dat$e <- cbind(dat1$e, dat2$e[,2], dat3$e[,2])
dat$c <- cbind(dat1$c, dat2$c[,2], dat3$c[,2])

intlabels <- c("Current",
               "Enhanced TSPOT: All, no Indets", "Enhanced TSPOT: HIV positive, no Indets", "Enhanced TSPOT: HIV negative, no Indets")

# m <- bcea(e=dat$e, c=-dat$c, ref=1, interventions = intlabels)
# contour2(m, wtp=WTP, graph = "ggplot2", ICER.size=2, pos=c(0.1,0.9))+#, xlim=c(-5,5), ylim=c(-200,200)) +
#   ggtitle("") #+ geom_abline(intercept = 0, slope = WTP)

my.plot.bcea(dat1, dat2, dat3, wtp=WTP, intlabels = intlabels, wtpNEG = "Y")
```

```
## Loading required package: MASS
```

```
## Loading required package: car
```

![](MAINCHAPTER-QFN_TSPOT_files/unnamed-chunk-5-1.png)

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
```

```r
my.plot.bcea(dat1, dat2, dat3, wtp=WTP*yearindays, intlabels = c("All","HIV positive","HIV negative"), wtpNEG = "Y", SCALEcosts=FALSE, SCALEdays=FALSE, labelLong=FALSE, TITLE="(a)", N=769)
```

![](MAINCHAPTER-QFN_TSPOT_files/unnamed-chunk-5-2.png)

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
```

```r
my.plot.bcea(dat1, dat2, dat3, wtp=WTP, intlabels = intlabels, contour=TRUE, wtpNEG = "Y")
```

![](MAINCHAPTER-QFN_TSPOT_files/unnamed-chunk-5-3.png)

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
```

```r
my.plot.bcea(dat1, dat2, dat3, wtp=WTP*yearindays, intlabels = c("All","HIV positive","HIV negative"), contour=TRUE, wtpNEG = "Y", SCALEcosts=FALSE, SCALEdays=FALSE, labelLong=FALSE, TITLE="(a)", N=769)
```

![](MAINCHAPTER-QFN_TSPOT_files/unnamed-chunk-5-4.png)

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
```

```r
my.plot.bcea(dat1, dat2, dat3, wtp=WTP, intlabels = intlabels, contour=TRUE, LEVELS=0.5, wtpNEG = "Y")
```

![](MAINCHAPTER-QFN_TSPOT_files/unnamed-chunk-5-5.png)

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
```

```r
my.plot.bcea(dat1, dat2, dat3, wtp=WTP*yearindays, intlabels = c("All","HIV positive","HIV negative"), contour=TRUE, LEVELS=0.5, wtpNEG = "Y", SCALEcosts=FALSE, SCALEdays=FALSE, labelLong=FALSE, TITLE="(a)", N=769)
```

![](MAINCHAPTER-QFN_TSPOT_files/unnamed-chunk-5-6.png)

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
```

```r
my.plot.ceac(dat1, dat2, dat3, intlabels = c("All","HIV positive","HIV negative"), labelLong=FALSE, TITLE="(a)")
```

![](MAINCHAPTER-QFN_TSPOT_files/unnamed-chunk-5-7.png)

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
```

```r
my.plot.ceac(dat1, dat2, dat3, intlabels = c("All","HIV positive","HIV negative"), labelLong=FALSE, TITLE="(a)", SCALEcosts=FALSE, SCALEdays=FALSE, N=769)
```

![](MAINCHAPTER-QFN_TSPOT_files/unnamed-chunk-5-8.png)

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
```

```r
my.plot.ceac(dat1, dat2, dat3, intlabels = c("All","HIV positive","HIV negative"), labelLong=FALSE, TITLE="(a)", SCALEcosts=FALSE, SCALEdays=FALSE, N=769, CI=TRUE)
```

![](MAINCHAPTER-QFN_TSPOT_files/unnamed-chunk-5-9.png)

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
```

```r
## in years (not days)
# m <- bcea(e=dat$e/365, c=-dat$c, ref=1, interventions = intlabels)
# contour2(m, wtp=20000, graph = "ggplot2", ICER.size=2, pos=c(0.1,0.9))
# 
# ceac.plot(m)
```



```r
dat1 <- IDEAdectree.simple(data=data, name.ruleout = "TSPOT",
                           SPEC = TSPOT.Indet.spec.mean, SENS = TSPOT.Indet.sens.mean, 
                           SPECvar = TSPOT.Indet.spec.var, SENSvar = TSPOT.Indet.sens.var)
dat2 <- IDEAdectree.simple(data=data, name.ruleout = "TSPOT",
                           SPEC = TSPOT.HIV.Indet.spec.mean, SENS = TSPOT.HIV.Indet.sens.mean, 
                           SPECvar = TSPOT.HIV.Indet.spec.var, SENSvar = TSPOT.HIV.Indet.sens.var)
dat3 <- IDEAdectree.simple(data=data, name.ruleout = "TSPOT",
                           SPEC = TSPOT.HIVneg.Indet.spec.mean, SENS = TSPOT.HIVneg.Indet.sens.mean, 
                           SPECvar = TSPOT.HIVneg.Indet.spec.var, SENSvar = TSPOT.HIVneg.Indet.sens.var)

dat$e <- cbind(dat1$e, dat2$e[,2], dat3$e[,2])
dat$c <- cbind(dat1$c, dat2$c[,2], dat3$c[,2])

intlabels <- c("Current", "Enhanced TSPOT: All, with Indets", "Enhanced TSPOT: HIV positive, with Indets", "Enhanced TSPOT: HIV negative, with Indets")

# m <- bcea(e=dat$e, c=-dat$c, ref=1, interventions = intlabels)
# contour2(m, wtp=WTP, graph = "ggplot2", ICER.size=2, pos=c(0.1,0.9))+#, xlim=c(-5,5), ylim=c(-200,200)) +
#   ggtitle("") #+ geom_abline(intercept = 0, slope = WTP)

my.plot.bcea(dat1, dat2, dat3, wtp=WTP, intlabels = intlabels, wtpNEG = "Y")
```

![](MAINCHAPTER-QFN_TSPOT_files/unnamed-chunk-6-1.png)

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
```

```r
my.plot.bcea(dat1, dat2, dat3, wtp=WTP*yearindays, intlabels = c("All","HIV positive","HIV negative"), wtpNEG = "Y", SCALEcosts=FALSE, SCALEdays=FALSE, labelLong=FALSE, TITLE="(b)", N=769)
```

![](MAINCHAPTER-QFN_TSPOT_files/unnamed-chunk-6-2.png)

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
```

```r
my.plot.bcea(dat1, dat2, dat3, wtp=WTP, intlabels = intlabels, contour=TRUE, wtpNEG = "Y")
```

![](MAINCHAPTER-QFN_TSPOT_files/unnamed-chunk-6-3.png)

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
```

```r
my.plot.bcea(dat1, dat2, dat3, wtp=WTP*yearindays, intlabels = c("All","HIV positive","HIV negative"), contour=TRUE, wtpNEG = "Y", SCALEcosts=FALSE, SCALEdays=FALSE, labelLong=FALSE, TITLE="(b)", N=769)
```

![](MAINCHAPTER-QFN_TSPOT_files/unnamed-chunk-6-4.png)

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
```

```r
my.plot.bcea(dat1, dat2, dat3, wtp=WTP, intlabels = intlabels, contour=TRUE, LEVELS=0.5, wtpNEG = "Y")
```

![](MAINCHAPTER-QFN_TSPOT_files/unnamed-chunk-6-5.png)

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
```

```r
my.plot.bcea(dat1, dat2, dat3, wtp=WTP*yearindays, intlabels = c("All","HIV positive","HIV negative"), contour=TRUE, LEVELS=0.5, wtpNEG = "Y", SCALEcosts=FALSE, SCALEdays=FALSE, labelLong=FALSE, TITLE="(b)", N=769)
```

![](MAINCHAPTER-QFN_TSPOT_files/unnamed-chunk-6-6.png)

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
```

```r
my.plot.ceac(dat1, dat2, dat3, intlabels = c("All","HIV positive","HIV negative"), labelLong=FALSE, TITLE="(b)")
```

![](MAINCHAPTER-QFN_TSPOT_files/unnamed-chunk-6-7.png)

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
```

```r
my.plot.ceac(dat1, dat2, dat3, intlabels = c("All","HIV positive","HIV negative"), labelLong=FALSE, TITLE="(b)", SCALEcosts=FALSE, SCALEdays=FALSE, N=769)
```

![](MAINCHAPTER-QFN_TSPOT_files/unnamed-chunk-6-8.png)

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
```

```r
my.plot.ceac(dat1, dat2, dat3, intlabels = c("All","HIV positive","HIV negative"), labelLong=FALSE, TITLE="(b)", SCALEcosts=FALSE, SCALEdays=FALSE, N=769, CI=TRUE)
```

![](MAINCHAPTER-QFN_TSPOT_files/unnamed-chunk-6-9.png)

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
```

```r
## in years (not days)
# m <- bcea(e=dat$e/365, c=-dat$c, ref=1, interventions = intlabels)
# contour2(m, wtp=20000, graph = "ggplot2", ICER.size=2, pos=c(0.1,0.9))
# 
# ceac.plot(m)
```




```r
dat1 <- IDEAdectree.simple(data=data, name.ruleout = "QFN",
                           SPEC = QFN.noIndet.spec.mean, SENS = QFN.noIndet.sens.mean, 
                           SPECvar = QFN.noIndet.spec.var, SENSvar = QFN.noIndet.sens.var)
dat2 <- IDEAdectree.simple(data=data, name.ruleout = "QFN",
                           SPEC = QFN.HIV.noIndet.spec.mean, SENS = QFN.HIV.noIndet.sens.mean, 
                           SPECvar = QFN.HIV.noIndet.spec.var, SENSvar = QFN.HIV.noIndet.sens.var)
dat3 <- IDEAdectree.simple(data=data, name.ruleout = "QFN",
                           SPEC = QFN.HIVneg.noIndet.spec.mean, SENS = QFN.HIVneg.noIndet.sens.mean, 
                           SPECvar = QFN.HIVneg.noIndet.spec.var, SENSvar = QFN.HIVneg.noIndet.sens.var)

dat$e <- cbind(dat1$e, dat2$e[,2], dat3$e[,2])
dat$c <- cbind(dat1$c, dat2$c[,2], dat3$c[,2])

intlabels <- c("Current",
               "Enhanced QFN: All, no Indets", "Enhanced QFN: HIV positive, no Indets", "Enhanced QFN: HIV negative, no Indets")

# m <- bcea(e=dat$e, c=-dat$c, ref=1, interventions = intlabels)
# contour2(m, wtp=WTP, graph = "ggplot2", ICER.size=2, pos=c(0.1,0.9))+#, xlim=c(-5,5), ylim=c(-200,200)) +
#   ggtitle("") #+ geom_abline(intercept = 0, slope = WTP)

my.plot.bcea(dat1, dat2, dat3, wtp=WTP, intlabels = intlabels, wtpNEG = "Y")
```

![](MAINCHAPTER-QFN_TSPOT_files/unnamed-chunk-7-1.png)

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
```

```r
my.plot.bcea(dat1, dat2, dat3, wtp=WTP*yearindays, intlabels = c("All","HIV positive","HIV negative"), wtpNEG = "Y", SCALEcosts=FALSE, SCALEdays=FALSE, labelLong=FALSE, TITLE="(c)", N=769)
```

![](MAINCHAPTER-QFN_TSPOT_files/unnamed-chunk-7-2.png)

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
```

```r
my.plot.bcea(dat1, dat2, dat3, wtp=WTP, intlabels = intlabels, contour=TRUE, wtpNEG = "Y")
```

![](MAINCHAPTER-QFN_TSPOT_files/unnamed-chunk-7-3.png)

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
```

```r
my.plot.bcea(dat1, dat2, dat3, wtp=WTP*yearindays, intlabels = c("All","HIV positive","HIV negative"), contour=TRUE, wtpNEG = "Y", SCALEcosts=FALSE, SCALEdays=FALSE, labelLong=FALSE, TITLE="(c)", N=769)
```

![](MAINCHAPTER-QFN_TSPOT_files/unnamed-chunk-7-4.png)

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
```

```r
my.plot.bcea(dat1, dat2, dat3, wtp=WTP, intlabels = intlabels, contour=TRUE, LEVELS=0.5, wtpNEG = "Y")
```

![](MAINCHAPTER-QFN_TSPOT_files/unnamed-chunk-7-5.png)

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
```

```r
my.plot.bcea(dat1, dat2, dat3, wtp=WTP*yearindays, intlabels = c("All","HIV positive","HIV negative"), contour=TRUE, LEVELS=0.5, wtpNEG = "Y", SCALEcosts=FALSE, SCALEdays=FALSE, labelLong=FALSE, TITLE="(c)", N=769)
```

![](MAINCHAPTER-QFN_TSPOT_files/unnamed-chunk-7-6.png)

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
```

```r
my.plot.ceac(dat1, dat2, dat3, intlabels = c("All","HIV positive","HIV negative"), labelLong=FALSE, TITLE="(c)")
```

![](MAINCHAPTER-QFN_TSPOT_files/unnamed-chunk-7-7.png)

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
```

```r
my.plot.ceac(dat1, dat2, dat3, intlabels = c("All","HIV positive","HIV negative"), labelLong=FALSE, TITLE="(c)", SCALEcosts=FALSE, SCALEdays=FALSE, N=769)
```

![](MAINCHAPTER-QFN_TSPOT_files/unnamed-chunk-7-8.png)

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
```

```r
my.plot.ceac(dat1, dat2, dat3, intlabels = c("All","HIV positive","HIV negative"), labelLong=FALSE, TITLE="(c)", SCALEcosts=FALSE, SCALEdays=FALSE, N=769, CI=TRUE)
```

![](MAINCHAPTER-QFN_TSPOT_files/unnamed-chunk-7-9.png)

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
```

```r
## in years (not days)
# m <- bcea(e=dat$e/365, c=-dat$c, ref=1, interventions = intlabels)
# contour2(m, wtp=20000, graph = "ggplot2", ICER.size=2, pos=c(0.1,0.9))
# 
# ceac.plot(m)
```



```r
dat1 <- IDEAdectree.simple(data=data, name.ruleout = "QFN",
                           SPEC = QFN.Indet.spec.mean, SENS = QFN.Indet.sens.mean, 
                           SPECvar = QFN.Indet.spec.var, SENSvar = QFN.Indet.sens.var)
dat2 <- IDEAdectree.simple(data=data, name.ruleout = "QFN",
                           SPEC = QFN.HIV.Indet.spec.mean, SENS = QFN.HIV.Indet.sens.mean, 
                           SPECvar = QFN.HIV.Indet.spec.var, SENSvar = QFN.HIV.Indet.sens.var)
dat3 <- IDEAdectree.simple(data=data, name.ruleout = "QFN",
                           SPEC = QFN.HIVneg.Indet.spec.mean, SENS = QFN.HIVneg.Indet.sens.mean, 
                           SPECvar = QFN.HIVneg.Indet.spec.var, SENSvar = QFN.HIVneg.Indet.sens.var)

dat$e <- cbind(dat1$e, dat2$e[,2], dat3$e[,2])
dat$c <- cbind(dat1$c, dat2$c[,2], dat3$c[,2])

intlabels <- c("Current", "Enhanced QFN: All, with Indets", "Enhanced QFN: HIV positive, with Indets", "Enhanced QFN: HIV negative, with Indets")

# m <- bcea(e=dat$e, c=-dat$c, ref=1, interventions = intlabels)
# contour2(m, wtp=WTP, graph = "ggplot2", ICER.size=2, pos=c(0.1,0.9))+#, xlim=c(-5,5), ylim=c(-200,200)) +
#   ggtitle("") #+ geom_abline(intercept = 0, slope = WTP)

my.plot.bcea(dat1, dat2, dat3, wtp=WTP, intlabels = intlabels, wtpNEG = "Y")
```

![](MAINCHAPTER-QFN_TSPOT_files/unnamed-chunk-8-1.png)

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
```

```r
my.plot.bcea(dat1, dat2, dat3, wtp=WTP*yearindays, intlabels = c("All","HIV positive","HIV negative"), wtpNEG = "Y", SCALEcosts=FALSE, SCALEdays=FALSE, labelLong=FALSE, TITLE="(d)", N=769)
```

![](MAINCHAPTER-QFN_TSPOT_files/unnamed-chunk-8-2.png)

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
```

```r
my.plot.bcea(dat1, dat2, dat3, wtp=WTP, intlabels = intlabels, contour=TRUE, wtpNEG = "Y")
```

![](MAINCHAPTER-QFN_TSPOT_files/unnamed-chunk-8-3.png)

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
```

```r
my.plot.bcea(dat1, dat2, dat3, wtp=WTP*yearindays, intlabels = c("All","HIV positive","HIV negative"), contour=TRUE, wtpNEG = "Y", SCALEcosts=FALSE, SCALEdays=FALSE, labelLong=FALSE, TITLE="(d)", N=769)
```

![](MAINCHAPTER-QFN_TSPOT_files/unnamed-chunk-8-4.png)

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
```

```r
my.plot.bcea(dat1, dat2, dat3, wtp=WTP, intlabels = intlabels, contour=TRUE, LEVELS=0.5, wtpNEG = "Y")
```

![](MAINCHAPTER-QFN_TSPOT_files/unnamed-chunk-8-5.png)

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
```

```r
my.plot.bcea(dat1, dat2, dat3, wtp=WTP*yearindays, intlabels = c("All","HIV positive","HIV negative"), contour=TRUE, LEVELS=0.5, wtpNEG = "Y", SCALEcosts=FALSE, SCALEdays=FALSE, labelLong=FALSE, TITLE="(d)", N=769)
```

![](MAINCHAPTER-QFN_TSPOT_files/unnamed-chunk-8-6.png)

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
```

```r
my.plot.ceac(dat1, dat2, dat3, intlabels = c("All","HIV positive","HIV negative"), labelLong=FALSE, TITLE="(d)")
```

![](MAINCHAPTER-QFN_TSPOT_files/unnamed-chunk-8-7.png)

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
```

```r
my.plot.ceac(dat1, dat2, dat3, intlabels = c("All","HIV positive","HIV negative"), labelLong=FALSE, TITLE="(d)", SCALEcosts=FALSE, SCALEdays=FALSE, N=769)
```

![](MAINCHAPTER-QFN_TSPOT_files/unnamed-chunk-8-8.png)

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
```

```r
my.plot.ceac(dat1, dat2, dat3, intlabels = c("All","HIV positive","HIV negative"), labelLong=FALSE, TITLE="(d)", SCALEcosts=FALSE, SCALEdays=FALSE, N=769, CI=TRUE)
```

![](MAINCHAPTER-QFN_TSPOT_files/unnamed-chunk-8-9.png)

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
```

```r
## in years (not days)
# m <- bcea(e=dat$e/365, c=-dat$c, ref=1, interventions = intlabels)
# contour2(m, wtp=20000, graph = "ggplot2", ICER.size=2, pos=c(0.1,0.9))
# 
# ceac.plot(m)
```



```r
detach(senspec.env)
```


