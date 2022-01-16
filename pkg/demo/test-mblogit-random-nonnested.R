library(mclogit)
set.seed(534)

nwithin  <- 100
nbetween1 <- 20
nbetween2 <- 20

nbetween <- nbetween1*nbetween2

a1 <- -1
a2 <- 1

x <- seq(from=-2,to=2,length=nwithin)

x <- rep(x,nbetween)

u1_1 <- rnorm(nbetween1,sd=1)
u2_1 <- rnorm(nbetween1,sd=1)

u1_2 <- rnorm(nbetween2,sd=1)
u2_2 <- rnorm(nbetween2,sd=1)


g1 <- rep(1:nbetween1,each=nwithin*nbetween2)
g2 <- rep(1:nbetween2,each=nwithin,nbetween1)

eta1 <-  1*x + a1 + u1_1[g1] + u1_2[g2]
eta2 <- -1*x + a2 + u2_1[g1] + u2_2[g2]

exp.eta1 <- exp(eta1)
exp.eta2 <- exp(eta2)
sum.exp.eta <- 1 + exp.eta1 + exp.eta2

pi2 <- exp.eta1/sum.exp.eta
pi3 <- exp.eta2/sum.exp.eta
pi1 <- 1 - pi2 - pi3
pi <- cbind(pi1,pi2,pi3)

y <- sapply(1:length(x),
            function(i)sample.int(n=3,size=1,prob=pi[i,]))
y <- factor(y,labels=letters[1:3])

plot(y~x)


(fem <- mblogit(y~x))

(mxm_x <- mblogit(y~x,
                random=list(~1|g1,~1|g2),
                estimator="REML"
                ))
summary(mxm_x)

(mxm <- mblogit(y~x,
                random=~1|g1,
                estimator="REML"
                ))
summary(mxm)

pred_x <- predict(mxm_x,type="response")

pred_1 <- predict(mxm,type="response")

plot(pred_x,type="l")

plot(x,pred_x[,1],type="l")

plot(x,pred_x[,2],type="l")

plot(x,pred_x[,3],type="l")


plot(x,pi1,type="l")

plot(x,pi2,type="l")

plot(x,pi3,type="l")


plot(pi1,pred_x[,1],type="l")

plot(pi2,pred_x[,2],type="l")

plot(pi3,pred_x[,3],type="l")


epred_x <- predict(mxm_x)

plot(eta1,epred_x[,1],type="l")
abline(a=0,b=1,col="red")

plot(eta2,epred_x[,2],type="l")
abline(a=0,b=1,col="red")

Bmxm_x <- mclogit:::reff(mxm_x)

c(u1_1=mean(u1_1),
  u1_1_hat=mean(Bmxm_x[[1]][,1]))

plot(u1_1-mean(u1_1),Bmxm_x[[1]][,1])
abline(a=0,b=1)

plot(u2_1-mean(u2_1),Bmxm_x[[1]][,2])
abline(a=0,b=1)

plot(u1_2-mean(u1_2),Bmxm_x[[2]][,1])
abline(a=0,b=1)

plot(u2_2-mean(u2_2),Bmxm_x[[2]][,2])
abline(a=0,b=1)

(mxm_i <- mblogit(y~x,
                  random=~1+x|g1
                  ))

f <- sample(1:2,size=length(x),replace=TRUE)

(mxm_ii <- mblogit(y~x*f,
                  random=~1+x|g1
                  ))

summary(mxm_ii)
