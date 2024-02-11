\name{strategicplayers-package}
\alias{strategicplayers-package}
\alias{strategicplayers}
\docType{package}
\title{
\packageTitle{strategicplayers}
}
\description{
\packageDescription{strategicplayers}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{strategicplayers}
\packageIndices{strategicplayers}
use the sp function to get a list of strategic players indicies


}
\author{
\packageAuthor{strategicplayers}

Maintainer: \packageMaintainer{strategicplayers}
}
\references{
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6959850/
}

\keyword{ package }

\examples{
 
 
#I am commenting this all out so that the package won't require sna any more :)

#require(sna)

#generate a bernoulli random network on 20 nodes
#network<-rgraph(20, tprob=.2)          

#get the geodesic distances of the network
#geo<-geodist(network)[2]$gdist  

#defining the target group
#targets<-1:10                        

#defining the avoidance group
#avoids<-11:14                         

#defining the theta parameter
#theta<-.8                             

#find sp set of size 4
#spset<-sp(4, geo, targets, avoids, theta, n.loops=100)  
#spset

#calculates distance metric for spset
#distance(geo, targets, avoids, theta, spset)           

#plot the network with the strategic player set highlighted in yellow

#colors<-rep("white", 20)
#colors[targets]<-"green"
#colors[avoids]<-"red"
#colors[spset]<-"yellow"
#par(mar=c(1,1,1,1))
#gplot(network, vertex.col=colors, 
#usearrows=FALSE, edge.col="grey",
#vertex.border="grey", vertex.cex=1.7, pad=0, label=1:dim(network)[1])


}



