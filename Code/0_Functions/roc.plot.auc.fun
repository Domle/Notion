roc.plot.auc<-function(x){      
	X<-data.frame(sens=x$se,spec=x$sp,eenminspec=1-x$sp,
                        diffsens=rep(0,length=length(x$se)),
         
               diffspec=rep(0,length=length(x$sp)),
                        area=rep(0,length=length(x$se)))
        X$diffsens[-(nrow(x))]<-X$sens[2:(nrow(X))]-X$sens[1:(nrow(X)-1)]
	X$diffspec[-(nrow(x))]<-X$eenminspec[2:(nrow(X))]-X$eenminspec[1:(nrow(X)-1)]
        X$area<-X$sens*X$diffspec+0.5*X$diffsens*X$diffspec
        return(sum(X$area))
}

