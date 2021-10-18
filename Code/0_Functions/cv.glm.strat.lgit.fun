cv.glm.strat.lgit<-function(data, glm.obj, K, lvo.cv = F) 
{
 
   # CONTROL IF THE NUMBER OF OBSERVATIONS IS SUFFICIENT
    
   n <- nrow(data)
   if ((K > n) || (K <= 1)) 
        stop("K outside allowable range")
    
   n.0 <- nrow(data[glm.obj$y==0,])
   n.1 <- nrow(data[glm.obj$y==1,])
   

   # IF LVC.CV IS DISABLED
   if (lvo.cv == F)
        {
        if (n.0 < K || n.1 < K)
                {
                K <- min(n.0,n.1)       
                
                warning(paste("K has been set to",K,"(0/1 not sufficient)"))
                }       
        }
    

    df.res <- data.frame()
    if (lvo.cv == T && n.1 < 30)
    {
        K <- nrow(data)
        warning(paste("K has been set to",K,"(leave-one-out CV is enabled!)"))  
        
        for (i in 1:K)
        {
                temp.glm <- update(glm.obj, data = data[-i,])
                df.res[i,"id"] <- i
                df.res[i,"predicted"] <- predict(temp.glm, data[i,],type="response")
                
        }
        return(df.res)
    }   else 
    {
    

    id.0 <- as.vector(row.names(data[glm.obj$y==0,]),mode = "numeric")
    K.0 <- K
    K.lst.0 <- round(K)
    kvals.0 <- unique(round(n.0/(1:floor(n.0/2))))
    temp.0 <- abs(kvals.0 - K.lst.0)
    if (!any(temp.0 == 0)) 
    {
        K.lst.0 <- kvals.0[temp.0 == min(temp.0)][1]
    }
    if (K.lst.0 != K.0) 
    { 
        warning(paste("K.0 has been set to", K.lst.0))
    }


    id.1 <- as.vector(row.names(data[glm.obj$y==1,]), mode = "numeric")
    K.1 <- K
    K.lst.1 <- round(K)
    kvals.1 <- unique(round(n.1/(1:floor(n.1/2))))
    temp.1 <- abs(kvals.1 - K.lst.1)
    if (!any(temp.1 == 0)) 
    { 
        K.lst.1 <- kvals.1[temp.1 == min(temp.1)][1]
    }
    if (K.lst.1 != K.1) 
    { 
        warning(paste("K.1 has been set to", K.lst.1))
     }   
     if (K.lst.0 != K.lst.1)
    {
         warning(paste("P/A stratifications have not the same values",sep=""))
         min.K <- min(K.lst.0,K.lst.1)
         K.lst.0 <- min.K
         K.lst.1 <- min.K
    }
    

    warning(paste("K has been finally set to",K.lst.0))
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
    
    # RESPONSE PREDICTION
    for (i in 1:K.lst) 
    {
        j.out <- id.0[(s.0 == i)]
        j.out <- append(j.out, id.1[(s.1 == i)], after=length(j.out)) 
        j.out <- sort(j.out)
        j.in <- id.0[(s.0 != i)]
        j.in <- append(j.in, id.1[(s.1 != i)], after=length(j.in)) 
        j.in <- sort(j.in)
        glm.cal <- update(glm.obj, data = data[j.in, , drop = FALSE])
        glm.val <- predict(glm.cal, data[j.out, , drop = FALSE], type = "response")
        
         if (i == 1)
        {
                vect.id <- j.out
                vect.predicted <- as.vector(glm.val)
        } else if (i > 1)
        {       
                vect.id <- append(vect.id, j.out, after=length(vect.id))
                vect.predicted <- append(vect.predicted, as.vector(glm.val), after=length(vect.predicted))
        }       
        
        
    }
    df.res[1:length(vect.id),"id"] <- vect.id
    df.res[1:length(vect.predicted),"predicted"] <- vect.predicted
    #df.res <- data.frame(sort.col(df.res,columns.to.sort="@ALL",columns.to.sort.by=1))
    df.res <-data.frame(df.res[order(df.res[,1]),])


    return(df.res)
    }
}

