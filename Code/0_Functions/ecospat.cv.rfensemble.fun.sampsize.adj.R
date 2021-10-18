### ecospat.cv.rf function copied from : https://github.com/cran/ecospat/blob/master/R/ecospat.cv.R, accessed: 10.09.2016, since R 3.2.2 is not fully supporting the ecospat.package;

# ## Modified version (see comment)

# Help page
# rf.obj = any calibrated randomForest object with a binomial error distribution
# data.cv = A dataframe object containing the calibration data set with the same names for response and predictor variables 
# WHAT DOES THIS MEAN??
# -> It means some twitching is needed: data frame needs to have observations (the response) in first column, and the used predictors in the rest of the columns, and row ID needs to be according to length of data-frame

# For testing
# rf.obj <- mod_full
# data.cv <- data.frame("obs" = dat[ , c("obs")], dat[ , which(names(dat)%in%vrs)])
# K <- 8
# cv.lim <- 1
# jack.knife <- F
# ecospat.cv.rf (rf.obj, data.cv, K = 6, jack.knife = F)


# ECOSPAT CROSS-EVAL RANDOM FOREST FUNCTION
ecospat.cv.rf.adj <- function(rf.obj,data.cv, K=10, cv.lim=10, jack.knife = F){

	# CONTROL IF THE NUMBER OF OBSERVATIONS IS SUFFICIENT
	
	n <- nrow(data.cv)
	if ((K > n) || (K <= 1)) {stop("K outside allowable range")}
	n.0 <- length(rf.obj$y[rf.obj$y==0])
	n.1 <- length(rf.obj$y[rf.obj$y==1])

   # IF JACKNIVE IS ENABLED OR MIN PRESENCES VERY SMALL
	if (jack.knife == F && n.1 < cv.lim || jack.knife == F && n.1 < K || jack.knife == T){
		
	df.res <- data.frame(id=NA,predicted=NA)
	K <- nrow(data.cv)
	cat("K has been set to",K,"(leave-one-out CV is enabled!)","\n",append = F)
	
	for (i in 1:K){
	rf.tmp.cal <- randomForest(x=data.cv[-i,2:ncol(data.cv)], y=as.factor(data.cv[-i,1]), ntree=1000, importance=TRUE)
	rf.tmp.eval <- predict(rf.tmp.cal, data.cv[i,2:ncol(data.cv)], type="prob") 
	if (i == 1)
	{
	vect.id <- i
	vect.predicted <- as.double(rf.tmp.eval[,2])
	} else if (i > 1)
	{ 
	vect.id <- append(vect.id, i, after=length(vect.id))
	vect.predicted <- append(vect.predicted, as.double(rf.tmp.eval[,2]), after=length(vect.predicted))
	} 
	}
	
	df.tmp.res <- data.frame(cbind(id=round(as.numeric(vect.id),digits=0),predictions=as.numeric(vect.predicted))[order(vect.id),])
	df.res <- data.frame(id=df.tmp.res[,1],obs=as.vector(rf.obj$y),predictions=df.tmp.res[,2])

	# IF JACKNIFE IS DISABLED 
	} else {

	id.0 <- as.vector(row.names(data.cv)[data.cv[,1]==0], mode = "numeric") # 
	K.0 <- K
	K.lst.0 <- round(K)
	# some internal check for absences?
	kvals.0 <- unique(round(n.0/(1:floor(n.0/2))))
	temp.0 <- abs(kvals.0 - K.lst.0)
	if (!any(temp.0 == 0)) 
	{
	K.lst.0 <- kvals.0[temp.0 == min(temp.0)][1]
	}
	if (K.lst.0 != K.0) 
	{ 
	cat("K.0 has been set to", K.lst.0,"\n",append=F)
	}

	id.1 <- as.vector(row.names(data.cv)[data.cv[,1]==1], mode = "numeric")
	K.1 <- K
	K.lst.1 <- round(K)
	 # some internal check for presences?
	kvals.1 <- unique(round(n.1/(1:floor(n.1/2))))
	temp.1 <- abs(kvals.1 - K.lst.1)
	if (!any(temp.1 == 0)) 
	{ 
	K.lst.1 <- kvals.1[temp.1 == min(temp.1)][1]
	}
	if (K.lst.1 != K.1) 
	{ 
	cat("K.1 has been set to", K.lst.1,"\n",append=F)
	} 
	if (K.lst.0 != K.lst.1)
	{
	cat("P/A stratifications have not the same values","\n",append=F)
	min.K <- min(K.lst.0,K.lst.1)
	K.lst.0 <- min.K
	K.lst.1 <- min.K
	}
	cat("K has been finally set to",K.lst.0,"\n",append=F)
	K.lst <- K.lst.0

	# P/A STRATIFICATION 
	# ABSENCE STRATIFICATION 
	f.0 <- ceiling(n.0/K.lst.0)
	s.0 <- sample(rep(1:K.lst.0,f.0), n.0)
	n.s.0 <- table(s.0)
	ms.0 <- max(s.0)

	# PRESENCE STRATIFICATION 
	f.1 <- ceiling(n.1/K.lst.1)
	s.1 <- sample(rep(1:K.lst.1,f.1), n.1)
	n.s.1 <- table(s.1)
	ms.1 <- max(s.1)

	# RESPONSE PREDICTION: Commend DR. : EACH P/A POINT IS EXACTLY ONCE PREDICTED
	df.res <- data.frame(id=NA,predicted=NA)
	for (i in 1:K.lst){
		j.out <- id.0[(s.0 == i)]
		j.out <- append(j.out, id.1[(s.1 == i)], after=length(j.out)) 
		j.out <- sort(j.out)
		j.in <- id.0[(s.0 != i)]
		j.in <- append(j.in, id.1[(s.1 != i)], after=length(j.in)) 
		j.in <- sort(j.in)

		# Comment DR: Calibration data (all except the fraction in question, the "j.in" data)
		set.seed(71)
		rf.tmp.cal <- randomForest(x=data.cv[j.in,2:ncol(data.cv)], y=as.factor(data.cv[j.in,1]),  importance=T, nodesize=1, ntree=3000, sampsize = c(length(which (data.cv[j.in,1] ==1)), length(which (data.cv[j.in,1] ==1)))  )    # Compared to default code (downloaded from internet address) I add nodesized=T, ntree=3000 (instead of 1000), and sampsize = c(... ) arguments.
		rf.tmp.val <- predict(rf.tmp.cal, data.cv[j.out,2:ncol(data.cv)], type="prob")

		if (i == 1)
		{
		vect.id <- j.out
		vect.predicted <- as.vector(rf.tmp.val[,2],mode="numeric")
		} else if (i > 1)
		{ 
		vect.id <- append(vect.id, j.out, after=length(vect.id))
		vect.predicted <- append(vect.predicted, as.vector(rf.tmp.val[,2]), after=length(vect.predicted))
		}
	# Comment DR: Close loop across split-runs	 
	}

df.tmp.res <- data.frame(cbind(id=round(as.numeric(vect.id),digits=0),predictions=as.numeric(vect.predicted))[order(vect.id),])
df.res <- data.frame(id=df.tmp.res[,1],obs=as.vector(rf.obj$y),predictions=df.tmp.res[,2])

}

df.res[df.res[,3]>1,3] <- 1
df.res[df.res[,3]<0,3] <- 0
return(df.res)

}
