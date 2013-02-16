electors <- local({

  require(stats)
  set.seed(1)

  unit.len=.5

  green <- c(
                  econ.left=.8,
                  welfare=0,
                  auth=-1.5
                  )

  labor <- c(
                  econ.left=2,
                  welfare=3,
                  auth=0
                  )

  communist <- c(
                  econ.left=5,
                  welfare=2,
                  auth=1
                  )

  conservative <- c(econ.left=-2,
                  welfare=-2,
                  auth=2
                  )

  liberal <- c(
                  econ.left=-2.5,
                  welfare=0,
                  auth=-1.5
                  )

  right.wing.populist <- c(
                  econ.left=0,
                  welfare=.2,
                  auth=4
                  )

  parties <- rbind(
          green,
          labor,
          communist,
          conservative,
          liberal,
          right.wing.populist
          )*unit.len


  working <- c(
                  econ.left=2,
                  welfare=5,
                  auth=1
                  )*unit.len

  new.middle <- c(
                  econ.left=-1,
                  welfare=1,
                  auth=-1
                  )

  old.middle <-c(
                  econ.left=-2.5,
                  welfare=0,
                  auth=2
                  )

  classes <- rbind(
          working,
          new.middle,
          old.middle
          )*unit.len

  T <- 25
  constant <- rep(1,nrow(classes)*T)

  util.econ <- -outer(
                  constant*classes[,"econ.left"],
                  parties[,"econ.left"],
                  "-"
                  )^2

  util.auth <- -outer(
                  constant*classes[,"auth"],
                  parties[,"auth"],
                  "-"
                  )^2

  util.welfare <- outer(
                  constant*classes[,"welfare"],
                  parties[,"welfare"]
                  )



  utils <- util.econ+util.auth+util.welfare

  rand <- matrix(rnorm(n=ncol(utils)*T,sd=1),nrow=T,ncol=ncol(utils))
  rand <- rand[rep(seq_len(nrow(rand)),each=nrow(classes)),]


  exp.eta <- exp(utils+rand)

  probs <- exp.eta/apply(exp.eta,1,sum)

  rownames(probs) <- rep(rownames(classes),T)

#   browser()

  the.sample<-c(
          working=500,
          new.middle=500,
          old.middle=500
          )

  the.tab <- structure(t(sapply(rownames(probs),
          function(nm) rmultinom(
                  n=1,
                  size=the.sample[nm],
                  prob=probs[nm,]
                  )
          )),
          dimnames=structure(dimnames(probs),names=c("class","party")),
          class="table"
          )

  electors <- within(as.data.frame(t(the.tab)),{
          time <- rep(1:T,each=nrow(classes)*nrow(parties))
          time <- time/T
          class <- factor(class,
                  levels=c(
                          "working",
                          "new.middle",
                          "old.middle"
                          )
                  )
          party <- factor(party,
                  levels=c(
                          "communist",
                          "labor",
                          "green",
                          "liberal",
                          "conservative",
                          "right.wing.populist"
                          )
                  )
  })

  parties <- as.data.frame(parties)
  parties <- parties[as.character(electors$party),]
  electors <- cbind(electors,parties)
  rownames(electors) <- NULL

  electors
})