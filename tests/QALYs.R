
dat1 <- IDEAdectree.simple(data=data, utility=0.01)
dat2 <- IDEAdectree.simple(data=data, utility=0.5)
dat3 <- IDEAdectree.simple(data=data, utility=0.99)

dat$e <- cbind(dat1$e, dat2$e[,2], dat3$e[,2])
dat$c <- cbind(dat1$c, dat2$c[,2], dat3$c[,2])

intlabels <- c("Current",
               "utility=0", "utility=0.5", "utility=1")

my.plot.bcea(dat1, dat2, dat3, wtp=WTP, intlabels = intlabels, wtpNEG = "Y")
my.plot.bcea(dat1, dat2, dat3, wtp=WTP, intlabels = intlabels, contour=TRUE, wtpNEG = "Y")
