Transport <- local({

  require(stats)
  set.seed(1)

  # Spatial positions of suburbs are
  # normally distributed

  N <- 10
  unit.len <- 10
  suburbs <- cbind(x=rnorm(n=N),y=rnorm(n=N))
  workpop.suburb <- rpois(n=N,lambda=200)


  # There are a couple of large firms in the area
  firms <- rbind(
          "Hollibarton"=c(-1,0),
          "Mikrobrainz"=c(1,1),
          "Whalemurd"=c(-1,1),
          "G. Gecko Mergers and Aquisitions"=c(1,-1),
          "Huddsucker Industries"=c(0,0)
          )

  colnames(firms) <- c("x","y")

  # There are also some bus stations scattered in
  # the area ...
  n.bus.stations <- as.integer(N/1.3)
  bus.stations <- cbind(  x=rnorm(n=n.bus.stations),
                          y=rnorm(n=n.bus.stations))

  # and some train stations

  n.train.stations <- as.integer(N/1.4)
  train.stations <- cbind(x=rnorm(n=n.train.stations),
                          y=rnorm(n=n.train.stations))

  # Oil is expensive these days (2020!).
  # Cars and car maintenance costs are high, too.
  car.cost <- 10
  bus.cost <- 4
  train.cost <- 3.7

  # Every firm has a bus station and a train station
  # and a parking place.
  # Therefore the distance of the nearest bus station
  # or train station, respectively, counts here.

  dist.car <- rep(0,N)
  dist.train <- sqrt(
          outer(suburbs[,1],train.stations[,1],"-")^2 +
          outer(suburbs[,2],train.stations[,2],"-")^2
          )
  dist.train <- apply(dist.train,1,min)

  dist.bus <- sqrt(
          outer(suburbs[,1],bus.stations[,1],"-")^2 +
          outer(suburbs[,2],bus.stations[,2],"-")^2
          )
  dist.bus <- apply(dist.bus,1,min)

  Transport <- expand.grid(
        transport=factor(1:3,labels=c("bus","car","train")),
        suburb=1:N)


  distance <- rep(NA,nrow(Transport))
  cost <- rep(NA,nrow(Transport))
  transport <- Transport$transport
  suburb <- Transport$suburb

  distance[transport=="bus"] <- dist.bus*unit.len
  distance[transport=="car"] <- dist.car*unit.len
  distance[transport=="train"] <- dist.train*unit.len

  cost[transport=="bus"] <- bus.cost
  cost[transport=="car"] <- car.cost
  cost[transport=="train"] <- train.cost

  working <- workpop.suburb[suburb]

  coef.true <- c(distance=-1.5,cost=-1)

  expeta.true <- exp(cbind(distance,cost) %*% coef.true)
  sumexpeta.true <- tapply(expeta.true,suburb,sum)
  prop.true <- expeta.true/c(sumexpeta.true[suburb])

  resp <- c(sapply(
                          split(data.frame(prop.true,working),
                                  suburb
                          ),
                  function(x) rmultinom(n=1,
                                  size=unique(x$working),
                                  prob=x$prop.true)
                  ))

  data.frame(transport,suburb,distance,cost,working,prop.true,resp)
})