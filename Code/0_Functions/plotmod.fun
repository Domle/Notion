plotmod<-function(pres.obj,prob.obj,spp)
#Function PLOTMOD
#
#Author:  Niklaus E. Zimmermann
#Address: Swiss Federal Research Institute WSL, CH-8903 Birmensdorf, Switzerland
#Email:   niklaus.zimmermann@wsl.ch
#
#Date:    2/23/2005
#
#Purpose: This function plots the proportion of observed occurrences per 10%-class
#         of predicted occurrence, including 95% confidence intervals,
#
#NOTE:    Needs to be cleaned further. Is only reduced in the output, but still
#         calculates lots of unused v ariables.

plotmod<-function(pres.obj,prob.obj,spp)
{
   pltdata<-matrix(nrow=10,ncol=12)
   colnames(pltdata)<-c("Class","ProbAvg","ClsProp","CLlCls","CLuCls","Ncls","SumCls","CumProp","-Stdev","+Stdev","Ncum","SumCum")

   tmp0<-round((prob.obj$fitted+0.05),digits=1)
   tmp1<-split(prob.obj$fitted,tmp0)
   tmp2<-sapply(tmp1,mean)                        #mean probability per 10% class
   tmp3<-sapply(tmp1,length)                      #Df = class length-1; Here, it's length only, Df is adjusted below
   tmp4<-split(pres.obj,tmp0)                     #pres./abs. vector list per 10%-probability class
   tmp5<-sapply(tmp4,mean)                        #proportion of occurrence per probability class
   tmp6<-sapply(tmp4,sum)                         #sum of occurrence per probability class
   tmp7<-sqrt(sapply(tmp4,var))                   #st.dev
   tmp8<-qt(0.975,tmp3)*tmp7/sqrt(tmp3)           #Confidence Interval
   
   for (i in 1:10)
   {                                              #most measures below are for cumulative proportion plots
      tmpA<-(pres.obj[prob.obj$fitted<=(i/10)])   #vector of pres/abs up to probability class
      tmpB<-length(tmpA)                          #length of cumulative occ. vector
      tmpC<-sum(tmpA)                             #cumulative sum of occurrence
      tmpD<-mean(tmpA)                            #cumulative proportion of occurrence
      tmpE<-sqrt(var(tmpA))                       #st.dev
      tmpF<-qt(0.975,(tmpB-1))*tmpE/tmpB          #confidence interval for cumulative proportion plot
      pltdata[i,1]<-i/10                          # Class   
      pltdata[i,2]<-tmp2[i]                       # ProbAvg 
      pltdata[i,3]<-tmp5[i]                       # ClsProp 
      pltdata[i,4]<-tmp5[i]-tmp8[i]               # CLlCls  
      pltdata[i,5]<-tmp5[i]+tmp8[i]               # CLuCls  
      pltdata[i,6]<-tmp3[i]                       # Ncls    
      pltdata[i,7]<-tmp6[i]                       # SumCls  
      pltdata[i,8]<-tmpC/sum(pres.obj)            # CumProp 
      pltdata[i,9]<-pltdata[i,8]-tmpE             # -st.dev  
      pltdata[i,10]<-pltdata[i,8]+tmpE            # +st.dev  
      pltdata[i,11]<-tmpB                         # Ncum    
      pltdata[i,12]<-tmpC                         # SumCum  
   }
   
   plot(pltdata[,1],pltdata[,3],pch=16,xlim=c(0,1),ylim=c(0,1),xlab="Predicted Probability of Occurrence",
     ylab="Fraction of Occurrence per Class",main=spp,cex.lab=1.6,cex.axis=1.2,cex.main=2)
   segments((pltdata[,1]-.01),pltdata[,4],(pltdata[,1]+.01),pltdata[,4],lwd=1.5)
   segments((pltdata[,1]-.01),pltdata[,5],(pltdata[,1]+.01),pltdata[,5],lwd=1.5)
   segments(pltdata[,1],pltdata[,4],pltdata[,1],pltdata[,5],lwd=1.5)
   abline(0,1)
   text(pltdata[,1],-.017,as.character(pltdata[,7]))
     
   return(pltdata)
}

