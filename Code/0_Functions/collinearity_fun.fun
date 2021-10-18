
# Function: 	collinearity()

# Description:  Tests a set of predictors for collinearities
#		and removes collinear predictors from a set of explanatory variables.
		
# Usage: 	collinearity(x)

# Arguments: 	x: a matrix containing the predictors, which will be checked for collinearity.

# Output:	A character vector containing the name of the selected predictors, 
#		among which none (or only a weak) collinearity exists 
#		(i.e. predictors with strong collinearity are excluded).

# Details: 	In the literature you can find different definitions for the condition index.
#		Here the condition index is defined as ka_j = d_1 / d_j,  for j = 1, 2, 3, ..., k. 
#		dj is the jth ordered singular value of x.
#		The largest condition index ka_k = d_1 / d_k is defined as the condition number.
#		The condition number is also defined as (+ sqrt(landa_max / landa_min)), 
#		where landa_max and landa_min are the maximal and minimal Eigenvalues of xTx respectively
#		(xT is the transpose of x).
#		A "large" condition number indicates a (near) linear dependency.
#		A condition number < 30 may indicate weak collinearity
#		A condition number > 30 indicates strong collinearity (> 100 extrene collinearity).




# References
# D. Belsley 1991: Conditioning Diagnostics: Collinearity and Weak Data in Regression, New York: John Wiley & Sons. 
# Ali S. Hadi 1996: Matrix Algebra as a tool, Wadsworth Publishing Company, Belmont.



 
# Date: September 2005
# Author: K. Steinmann






collinearity = function(predictor_matrix) {
 
alle.namen = dimnames(predictor_matrix)[[2]]
predictor_matrix = scale(predictor_matrix, center = FALSE, scale = TRUE)               

cn = numeric()
kicked.parameter = character()
dimension = dim(predictor_matrix)
d = dimension[[2]]


for (i in (1: (d-1)))  {									
		

	predictor_matrix.t = t(predictor_matrix)
	A = predictor_matrix.t %*% predictor_matrix
	eigen.object = eigen(A)
	dimension_a = dim(A)
	p = dimension_a[1]
	koeff.max = max(abs(eigen.object$vectors[,p]))					
	ind.koeff.max = which(abs(eigen.object$vectors[,p]) == koeff.max)		    
											
	namen.vector = dimnames(A)[[1]]							
	kicked.parameter[i] = namen.vector[ind.koeff.max]
	cn[i] = sqrt(abs(eigen.object$value[1])) / sqrt(abs(eigen.object$value[p]))			
											
	predictor_matrix = predictor_matrix[, - ind.koeff.max]
	
	
	}

common.names.ind = match(kicked.parameter, alle.namen, nomatch = NA, incomparables = FALSE)
additional.name = alle.namen[- common.names.ind]
kicked.parameter = c(kicked.parameter, additional.name)

condition.number = cn < 30											
selected.parameters = kicked.parameter[condition.number]
return(selected.parameters)

}




################################################################################################################################################

# Do not run

# Example

x1 = c(0.322, 0.012, 0.631, 0.316, 0, 0.631)
x2 = c(0.5, 0.5, 0.5, 0, 0, 0.5)
x3 = c(-0.5, 0.5, 0, 0.5, 0.5, 0)
x4 = c(0, -0.5, 0.5, 0.5, 0, 0.5)

x.df = data.frame(x1 = x1, x2=x2, x3=x3, x4=x4)
x.m = as.matrix(x.df)

variable.selection = collinearity(x.m)
variable.selection


# check with panel.cor function


## put (absolute) correlations on the upper panels,
## with size proportional to the correlations.
panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex * r)
}



pairs(x.df, lower.panel=panel.smooth, upper.panel=panel.cor)

