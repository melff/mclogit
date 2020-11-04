\name{electors}
\alias{electors}
\title{Class, Party Position, and Electoral Choice}
\description{This is an artificial data set on
electoral choice as influenced by class and party positions.
}
\usage{data(electors)}
\format{A data frame containing the following variables:

\describe{
	\item{class}{class position of voters}
	\item{party}{party that runs for election}
	\item{Freq}{freqency by which each party list is chosen by members of each class}
        \item{time}{time variable, runs from zero to one}
	\item{econ.left}{economic-policy "leftness" of each party}
	\item{welfare}{emphasis of welfare expansion of each party}
	\item{auth}{position on authoritarian issues}
}
}
\examples{
data(electors)

summary(mclogit(
  cbind(Freq,interaction(time,class))~econ.left+welfare+auth,
  data=electors))

summary(mclogit(
  cbind(Freq,interaction(time,class))~econ.left/class+welfare/class+auth/class,
  data=electors))

\donttest{# This takes a bit longer.
summary(mclogit(
  cbind(Freq,interaction(time,class))~econ.left/class+welfare/class+auth/class,
  random=~1|party.time,
  data=within(electors,party.time<-interaction(party,time))))

summary(mclogit(
  cbind(Freq,interaction(time,class))~econ.left/(class*time)+welfare/class+auth/class,
  random=~1|party.time,
  data=within(electors,{
        party.time <-interaction(party,time)
        econ.left.sq <- (econ.left-mean(econ.left))^2
        })))
}
}

\keyword{datasets}
