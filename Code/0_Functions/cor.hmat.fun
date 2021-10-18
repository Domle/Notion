# Written by : Niklaus E. Zimmermann, Ph.D.
#              Swiss Federal Research Institute WSL
#              Zuercherstrasse 111
#              CH-8903 Birmensdorf
#              nez@wsl.ch

# This function calculates a correlation semi-matrix 
# based on the Pearson correlation coefficient, The
# upper half above the diagonal is filled with NA vals.

cor.hmat <- function(objct)
{
   tmp<-dim(objct)
   k<-tmp[2]
   cor.matrix<-matrix(ncol=k,nrow=k)

   for (i in 1:k){
      for (j in 1:k){
         tmp<-cor.test(objct[,i], objct[,j], method="pearson")
         suff="   "
         if(tmp$p.value<0.100) suff=".  "
         if(tmp$p.value<0.050) suff="* "
         if(tmp$p.value<0.010) suff="** "
         if(tmp$p.value<0.001) suff="***"
         cor.matrix[i,j]<-paste(as.character(round(tmp$estimate,3)),suff,sep="")
      }
   }

   for (i in 2:k){
      for (j in 1:(i-1)){
         if (j>0){cor.matrix[j,i]<-""}
         
      }
   }

   tmp<-data.frame(cor.matrix)
   rownames(tmp)<-colnames(objct)
   colnames(tmp)<-colnames(objct)
   return(tmp)
}
