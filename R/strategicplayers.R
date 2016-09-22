

######Define distance function
#' distance
#' 
#' Takes in the geodesic distances, targets, avoiders, a parameter that prioritizes avoiding vs targetting, and the current players and returns the strategic players distance metric
#' 
#' @param gd a matrix of geodesic distances for the network of interest
#' @param targets a vector of indicies of the people you want to spread the intervention to
#' @param avoiders a vector of indicies of the people you don't want to spread the intervention to
#' @param theta a number between 0 and 1 which weights the distance metric, 1 only prioritizes closeness to targets, 0 only prioritizes maximizing distance from avoiders
#' @param players the indicies of people who you have chosen for the intervention (a subset of targets)
#' @return returns the distance metric for strategic players, which we want to maximize
#' @export
distance<-function(gd, targets, avoiders, theta, players){
  #gd= a matrix of geodesic distances for the network you are interested in
  #targets= the peole who you want to spread the intervention to
  #avoiders= the people you don't want to spread the intervention to
  #players=the people who you have chosen to get the intervention
  
  #returns the distance metric
  
  diag(gd)<-1
  gd.target<-as.matrix(gd[targets, players], ncol=length(players))
  if (dim(gd.target)[1]<2){d.vec.target<-min(gd.target)}
  if (dim(gd.target)[1]>1)d.vec.target<-apply(gd.target, 1, min)
  
  gd.avoiders<-as.matrix(gd[avoiders, players], ncol=length(players))
  if (dim(gd.avoiders)[1]<2){d.vec.avoiders<-min(gd.avoiders)}
  if (dim(gd.avoiders)[1]>1){d.vec.avoiders<-apply(gd.avoiders, 1, min)}
  
  return(theta*mean(1/d.vec.target)-(1-theta)*mean(1/d.vec.avoiders))  
}



######Define strategic players function
#' sp
#' 
#' Takes in the number of intervention subjects you wish to identify, geodesic distances, targets, avoiders, and a parameter that prioritizes avoiding vs targetting, and returns the indecies of the strategic players
#' 
#' @param n.players the number of intervention subjects you wish to identify
#' @param gd a matrix of geodesic distances for the network of interest
#' @param targets a vector of indicies of the people you want to spread the intervention to
#' @param avoiders a vector of indicies of the people you don't want to spread the intervention to
#' @param theta a number between 0 and 1 which weights the distance metric, 1 only prioritizes closeness to targets, 0 only prioritizes maximizing distance from avoiders.  Any number between 0 and 1 will be a compromise of these two goals.
#' @param n.loops the number of loops to run, the more loops you run the more likely you are to identify the optimal set of strategic players
#' @return returns the indicies for strategic players
#' @export
sp<-function(n.players, gd, targets, avoiders, theta=.5, n.loops=1000){
  #n.players= number of strategic players
  #gd= geodesic distance matrix of network you are looking at
  #target= indices of people who you want to spread to who you could pick as players
  #avoiders= people you don't want to spread the intervention to
  #theta= tuning parameter that lets you decide whether to prioritize diffusion or containment
  #n.loops tells how many loops we should do.  For the actual interveion, do many (at least 200?)
  
  #returns the best strategic player group found and the strategic player group
  
  if(length(setdiff(targets, avoiders)) < length(targets)) stop("targets cannot also be avoiders")
  if(n.players>=length(targets)) stop("choose n.players to be less than the number of targets")
  if(theta<0 |theta>1) stop("theta must be between 0 and 1 (inclusive)")
  if (n.players<1) stop("n.players must be at least one")
  if(length(dim(gd))!=2) stop("please provide a square matrix of geodesic differences for gd")
  if(dim(gd)[1]!=dim(gd)[2]) stop("please provide a square matrix of geodesic differences for gd")
  if(max(targets)>dim(gd)[1] | min(targets)<1) stop("Incorrect target identification specification")
  if(max(avoiders)>dim(gd)[1] | min(avoiders)<1) stop("Incorrect avoider identification specification")
  
  targets<-unique(targets)
  avoiders<-unique(avoiders)
  player.set<-matrix(rep(0, (n.players+1)*n.loops), nrow=n.loops)
  for (j in 1: n.loops){
    diff=.5
    players<-sample(targets, n.players)
    f=distance(gd, targets, avoiders, theta, players) 
    player.set[j,]<-c(f, players)
    loops=0
    stopit<-FALSE
    while(stopit==FALSE & loops<500){
      loops=loops+1
      potential.players<-rep(targets[!targets%in%players], length(players))
      potential.targets<-rep(players, each=length(targets[!targets%in%players]))
      f.prime.vec<-rep(-1, length(potential.players))
      for (i in 1: length(potential.players)){
        u<-potential.players[i]
        v<-potential.targets[i]
        new.players.now<-c(u, players[players!=v])
        f.prime.vec[i]<-distance(gd, targets, avoiders, theta, new.players.now)
      }
      max.f.prime<-max(f.prime.vec, na.rm=T)
      diff=max.f.prime-f 
      if (diff<=0 ){
        stopit<-TRUE
      }
      if (diff>0){
        swap<-(1:i)[f.prime.vec==max.f.prime ]
        swap<-swap[1]
        new.players<-c(potential.players[swap], players[players!=potential.targets[swap]])
        players<-new.players
        f=max.f.prime
        player.set[j,]<-c(f, players)
      }
    }
  }
  theset<-((1:n.loops)[player.set[,1]==max(player.set[,1])])[1]
  players.max<-player.set[theset, 2:(n.players+1)]
  return(players.max)
}



