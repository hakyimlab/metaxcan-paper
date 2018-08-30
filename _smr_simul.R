library(ggplot2)

chisqzy = rchisq(10e4,df=1,ncp=0)
chisqzx = rchisq(10e4,df=1,ncp=29)
instruments = which(chisqzx > qchisq(5e-8, df = 1, lower.tail=F))
chisqzy = chisqzy[instruments]
chisqzx = chisqzx[instruments]
smrstat = chisqzy * chisqzx / (chisqzy + chisqzx)

kk = rchisq(length(instruments),df=1,ncp=0)
qqplot(kk,smrstat);abline(0,1)

N <- 1000
simul <- c()

for (i in 1:N) {
  chisqzy = rchisq(10e4,df=1,ncp=0)
  chisqzx = rchisq(10e4,df=1,ncp=29)
  instruments = which(chisqzx > qchisq(5e-8, df = 1, lower.tail=F))
  chisqzy = chisqzy[instruments]
  chisqzx = chisqzx[instruments]
  smrstat = chisqzy * chisqzx / (chisqzy + chisqzx)
  simul <- c(simul, mean(smrstat))
}

data <- data.frame(mean=simul)
ggplot(data=data) + geom_density(aes(x=simul)) + scale_x_continuous(limits=c(0,1), breaks= c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) 
