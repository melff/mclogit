library(MASS) # For 'housing' data
library(nnet)
library(memisc)

(house.mult<- multinom(Sat ~ Infl + Type + Cont, weights = Freq,
                       data = housing))


(house.mblogit <- mblogit(Sat ~ Infl + Type + Cont, weights = Freq,
                         data = housing))

summary(house.mult)

summary(house.mblogit)

mtable(house.mblogit)
