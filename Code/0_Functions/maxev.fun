maxev<-function(sim, obs)
{
  k <- 0.05
  i <- 1
  eval <- data.frame(0)
  while(k < 1) 
  {
    a <- table(sim >= k, obs)[4]
    b <- table(sim >= k, obs)[2]
    c <- table(sim >= k, obs)[3]
    d <- table(sim >= k, obs)[1]
    N <- a + b + c + d
    eval[i,1] <- k
    eval[i,2] <- a
    eval[i,3] <- b
    eval[i,4] <- c
    eval[i,5] <- d
    eval[i,6] <- (a + c) / N
    eval[i,7] <- (b + d) / N
    eval[i,8] <- (a + d) / N
    eval[i,9] <- a / (a + c)
    eval[i,10] <- d / (b + d)
    eval[i,11] <- b / (b + d)
    eval[i,12] <- c / (a + c)
    eval[i,13] <- a / (a + b)
    eval[i,14] <- d / (c + d)
    eval[i,15] <- (b + c) / N
    eval[i,16] <- (a * d) / (c * b)
    eval[i,17] <- kappa2(table(sim >= k, obs))
    eval[i,18] <- (- (a*log(a)) - (b*log(b)) - (c*log(c)) - (d*log(d)) + ((a+b)*log(a+b)) + ((c+d)*log(c+d))) / ((N*log(N)) - ((a+c)*log(a+c)) - ((b+d)*log(b+d)))
    names(eval)<-c("cut","tp","fp","fn","tn","p","odp","ccr","se","sp","fpr","fnr","ppp","npp","mr","or","k","nmi")
    k <- k + 0.05
    i <- i + 1
  }
  return(eval)
}
