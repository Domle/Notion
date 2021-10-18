roc.plot<-function(z1,z2)
{
	eva <- data.frame(0)
	eva[1,1] <- 0
	eva[1,2] <- 1
	eva[1,3] <- 0
	names(eva)<-c("thrsh","se","sp")
	k <- 0.01
	i <- 2
	while(k < 0.99) {
		a <- table(z1 >= k, z2)[4]
		b <- table(z1 >= k, z2)[2]
		c <- table(z1 >= k, z2)[3]
		d <- table(z1 >= k, z2)[1]
		eva[i,"thrsh"] <- k
		eva[i,"se"] <- ifelse(is.na(a / (a + c)), 0, a / (a + c))
		eva[i,"sp"] <- d / (b + d)
		k <- k + 0.01
		i <- i + 1
	}
	eva[i,"thrsh"] <- 1
	eva[i,"se"] <- 0
	eva[i,"sp"] <- 1
	return(list(data = eva, auc = - roc.plot.auc(eva)))
}
