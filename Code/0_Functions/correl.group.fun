################################################################################################################################################
### Function written by Robin Engler.
### correl.group() will (hopefully) group the variables that are given as input according to their correlation value, and a given threshold
### 
### Xmat = an input matrix/data frame that contains all the variables to be tested in columns. There should be at least 3 variables in the matrix
###        and no value should be equal to NA (so make sure to clean your data before using the function, e.g. with "na.omit()")
### Threshold = The threshold above which two variables become grouped together. E.g. setting Threshold=0.7 will group together variables that
###             have a correlation >= 0.7
### cor.Method = The method to be used for computing correlation. See the help files of the cor() function (type "?cor") for help on this parameter.
### hclust.Method = the agglomeration method to be used. The default is "average" and this is the one I recomend to use.
###                 This should be (an unambiguous abbreviation of) one of "ward", "single", "complete", "average", "mcquitty", "median" or "centroid". See the help files of the hclust() function (type "?hclust") for help on this parameter.
### Absvalue = if set to TRUE, then the absolute values of correlation are considered instead of the correlation values.
### symetricToNA = if set to TRUE, the output correlation matrix will have only NA in its upper part (since the matrix is symetric no information is lost).
### PlotFileName = Enter the name to be given to the plot of the correlation-based tree graphic. If this parameter is left to "", then no plot is saved.
### GraphicFileType = "JPEG" will save the correlation-based tree plot as a JPEG file. "EMF" will save it as a windows enhanced metafile.
### OptionalTextForGraph = Any text that you would like to be added to the output graphs.
###



correl.group <- function(Xmat, Threshold=0.7, cor.Method="pearson", hclust.Method="average", Absvalue=T, symetricToNA=F, PlotFileName="", GraphicFileType="JPEG", OptionalTextForGraph=""){
	
	### Verify user input
	if(is.null(ncol(Xmat))){cat("INPUT ERROR: 'Xmat' must be a matrix with 3 or more columns \n"); return(NA)}
	if(ncol(Xmat)< 3){cat("INPUT ERROR: 'Xmat' must be a matrix with 3 or more columns \n"); return(NA)}
	if(all(cor.Method!=c("pearson","kendall","spearman"))){cat("INPUT ERROR: The 'cor.Method' argument must be either 'pearson','spearman' or 'kendall'. \n"); return(NA)}
	if(Absvalue==T) if(Threshold<=0 | Threshold>=1){cat("INPUT ERROR: The value of 'Threshold' must be in the range ]0:1[ \n"); return(NA)}
	if(Absvalue==F) if(Threshold<=-1 | Threshold>=1){cat("INPUT ERROR: The value of 'Threshold' must be in the range ]-1:1[ \n"); return(NA)}
	if(require(stats)==F){cat("INPUT ERROR: The library 'stats' is needed to run this function. Please instal 'stats' and try again. \n"); return(NA)}
	if(any(is.na(Xmat))){cat("INPUT ERROR: Your input matrix should not contain any NA values. You may e.g. remove them using 'na.omit()'. \n"); return(NA)}
	
	### Compute correlation between all variables (columns of Xmat). See how the variables cluster using hclust().
	if(Absvalue) Cor <- abs(cor(Xmat, method=cor.Method)) else Cor <- cor(Xmat, method=cor.Method)
	Tree <- hclust(as.dist(1-Cor), method=hclust.Method)
	### Note: hclust() objects contain a number of things, among others: -> Tree$height: The list of height at which the tree branches.
	###                                                                  -> Tree$merge: A matrix indicating how the variables group on the tree.
	
	### Group together the variables that have a correlation >= Threshold. This is easily done by "cutting" the tree
	Gr <- cutree(Tree,h=1-Threshold)
	### make sure the order of variables in "Gr" is the same as in the "Cor" matrix.
	if(any(names(Gr)!=colnames(Cor))){
		matchVal <- match(colnames(Cor),names(Gr))  # returns: which element in second argument matches the element in the first argument.
		Gr <- Gr[matchVal]
	}
	if(any(names(Gr)!=colnames(Cor))) cat("ERROR in data match. This should not happen \n")
	GrList <- list()
	GrNb <- max(Gr)
	for(J in 1:GrNb) GrList[[J]] <- names(Gr)[which(Gr==J)]
	names(GrList) <- paste("Group", 1:GrNb, sep="")
	
	### Compute average correlation between groups.
	BetweenGrCor <- matrix(NA, nrow=GrNb, ncol=GrNb, dimnames=list(names(GrList),names(GrList)))
	WithinGrCor <- rep(NA, GrNb)
	names(WithinGrCor) <- names(GrList)
	for(J in 1:GrNb){
		# 1. compute average correlation between groups.
		for(K in 1:GrNb) if(J>K) BetweenGrCor[J,K] <- round(mean(Cor[which(Gr==J),which(Gr==K)]), 3)
		# 2. compute average correlation within the group.
		Gr1 <- which(Gr==J)
		if(length(Gr1)>1) WithinGrCor[J] <- round((sum(Cor[Gr1,Gr1])-length(Gr1))/(length(Gr1)^2-length(Gr1)), 3)
		rm(Gr1)
	}
	
	### Add between group and within group correlation to function output
	GrList[[GrNb+1]] <- WithinGrCor
	GrList[[GrNb+2]] <- BetweenGrCor
	names(GrList)[GrNb+1] <- "MeanWithinGroupCorrelation"
	names(GrList)[GrNb+2] <- "MeanBetweenGroupCorrelation"
	### Add correlation values to ouput (if "symetricToNA=TRUE" we set all symetrical values to NA)
	if(symetricToNA) for(J in 1:nrow(Cor)) for(K in 1:nrow(Cor)) if(K>=J) Cor[J,K] <- NA
	GrList[[GrNb+3]] <- round(Cor,4)
	if(Absvalue) names(GrList)[GrNb+3] <- "AbsoluteValueCorrelations" else names(GrList)[GrNb+3] <- "Correlations"
	
	
	### If the user asked to save a plot of the Tree, then we do as she/he wishes...
	if(PlotFileName!=""){
		if(GraphicFileType=="JPEG") jpeg(filename=paste(PlotFileName, ".jpg", sep=""), h=800, w=800)
		if(GraphicFileType=="EMF") win.metafile(filename=paste(PlotFileName, ".emf", sep=""), width=8, height=8)
			Main <- "Correlation-based Tree"
			if(OptionalTextForGraph!="") Main <- paste(Main, "\n", OptionalTextForGraph, sep="")
			plot(Tree, ylab="1 - abs(Correlation)", main=Main, ylim=c(0,1))
			abline(h=1-Threshold, col="red", lty=2)
			if(Threshold<0.9) abline(h=0.9-Threshold, col="blue", lty=2)
		dev.off()
	}
	
	### Release memory, return function value...
	rm(Cor,Tree,WithinGrCor,BetweenGrCor,GrNb,Main)
	return(GrList)
}
################################################################################################################################################
################################################################################################################################################


