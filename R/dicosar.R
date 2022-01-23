#' Testing the equality of Pearson correlation coefficients between two groups
#' 
#' These are further details.
#' 
#' @section The value of \emph{code}:
#' 
#' \describe{ 
#' \item{1}{No issue found.}
#' \item{-2}{Sample size is too small and NA returned.}
#' \item{-3}{Too many identical values and NA returned.}
#' \item{-4}{At least one of the PCCs is one and NA returned.}
#' \item{-5}{The algorithm in DICOSAR does not converge.}
#' \item{2 or 3}{The algorithm converges, but the simplified version of DICOSAR (i.e., he first-order approximation of the signed root of the likelihood ratio statistic) is used.}
#' }
#'
#' @param x1 an N1 by P matrix, in which N1 is the number of samples in the first group, and P is the number of variables. If P>2, tests will be run pairwisely.
#' @param x2 an N2 by P matrix, in which N2 is the number of samples in the second group.
#' @param delta a logical value. If TRUE, the p-value from the Delta method will also be returned when P=2.
#' @return r1: The Pearson correlation coefficient of x1.
#' @return r2: The Pearson correlation coefficient of x2. 
#' @return p: The p-value of the test for the equality of Pearson correlation coefficients using DICOSAR.
#' @return code: More information about the convergence, optimization, etc. A negative value indicates a problem. See details.
#' @return p_delta: The p-value of the test for the equality of Pearson correlation coefficients using the Delta method.
#' @return p_global: A p-value of the global test for the equality of two correlation matrices for each pair of the P variables in the input data. Only returned when P>2. 
#' @export
#' @examples
#' library(dicosar)
#' n = 100
#' x1 = matrix(rnorm(2*n),nrow=n)
#' x2 = matrix(rnorm(2*n),nrow=n)
#' dicosar(x1,x2)
#' 


dicosar = function(x1,x2,delta=FALSE)
{
  
  nd <- ncol(x1)
  if((nd!=ncol(x2))|(nd<2))
  {
    stop('The two groups must have the same number of variables.')
  }
  
  imp <- TRUE
  code <- 1
  
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  
  if(min(n1,n2)<7)
  {
    warnings('Too small sample size. The sample size in each group must be larger than 6.')
    code <- -2
  }
  
  x1 <- scale(x1)
  x2 <- scale(x2)
  xa <- rbind(x1,x2)
  
  if(nd==2)
  {
    nuniq <- c(length(unique(xa[,1])),length(unique(xa[,2])))
    uniqrow <- max(nuniq)
    if(uniqrow<7)
    {
      uniqrow <- nrow(unique(xa))
    }
    
    rho1 <- cor(x1)[1,2]
    r1 <- (log(1+rho1) - log(1-rho1))/2
    rho2 <- cor(x2)[1,2]
    r2 <- (log(1+rho2) - log(1-rho2))/2
    
    if((min(nuniq)<4) | (uniqrow<7) | is.na(rho1) | is.na(rho2))
    {
      warnings('Too many identical values in the data.')
      code <- -3
    }
    
    if(code %in% c(-2,-3))
    {
      return(list(r1 = rho1, r2 = rho2, p=NA,p_delta=NA,code=code))
    }else{
      
      if((abs(rho1)==1) | (abs(rho2)==1))
      {
        code <- -4
        return(list(r1 = rho1, r2 = rho2, p=NA,p_delta=NA,code=code))
      }
      
      p_delta <- NA
      if(delta==TRUE)
      {
        d1 <- ld_z_ci_u(x1,rho1)
        d2 <- ld_z_ci_u(x2,rho2)
        vtst <- d1 + d2
        p_delta <- pchisq((r1-r2)^2/vtst,1,lower.tail=FALSE)
      }
      
      info1 <- cbind(xa[,1],xa[,2],xa[,1]^2,xa[,2]^2,xa[,1]*xa[,2])
      x0 <- c(apply(info1,2,mean),apply(info1,2,mean))
      
      p <- rep(NA,3)
      cp <- cum_dist(r1-r2,info1,info1,x0,imp,n1,n2)
      if(is.na(cp[1]))
      {
        code <- -5
        return(list(r1 = rho1, r2 = rho2, p=NA,p_delta=p_delta,code=code))
      }else{
        p[1] <- 2*min(cp[1],1-cp[1])
        if((!is.na(cp[3]))&(cp[3]>0)&(cp[3]<1))
        {
          p[3] <- 2*min(cp[3],1-cp[3])
          if((p[1]>0.9)&(p[3]<0.1))
          {
            p[3] <- p[1]
            code <- 3
          }
        }else{
          p[3] <- p[1]
          code = 2
        }
        return(list(r1 = rho1, r2 = rho2, p=p[3],p_delta=p_delta,code=code))
      }
      
    }
    
  }else{
    
    pm <- matrix(NA,nd,nd)
    rhom1 <- matrix(NA,nd,nd)
    rhom2 <- matrix(NA,nd,nd)
    p_global <- NA
    
    nuniq <- sapply(1:nd,function(x)length(unique(xa[,x])))
    
    for(i in 1:(nd-1))
    {
      for(j in (i+1):nd)
      {
        info1 <- cbind(xa[,i],xa[,j],xa[,i]^2,xa[,j]^2,xa[,i]*xa[,j])
        x0 <- c(apply(info1,2,mean),apply(info1,2,mean))
        rho1 <- cor(x1[,i],x1[,j])
        r1 <- (log(1+rho1) - log(1-rho1))/2
        rho2 <- cor(x2[,i],x2[,j])
        r2 <- (log(1+rho2) - log(1-rho2))/2
        
        uniqrow <- max(nuniq[c(i,j)])
        if(uniqrow<7)
        {
          uniqrow <- nrow(unique(xa[,c(i,j)]))
        }
        
        p <- rep(NA,3)
        
        if((!is.na(rho1)) & (!is.na(rho2)))
        {
          if((abs(rho1)<1) & (abs(rho2)<1) & (uniqrow>6) & (min(nuniq[c(i,j)])>3))
          {
            cp <- cum_dist(r1-r2,info1,info1,x0,imp,n1,n2)
            
            if(!is.na(cp[1]))
            {
              p[1] <- 2*min(cp[1],1-cp[1])
            }
            
            if((is.na(cp[3])==FALSE)&(cp[3]>0)&(cp[3]<1))
            {
              p[3] <- 2*min(cp[3],1-cp[3])
              if((p[1]>0.9)&(p[3]<0.1))
              {
                p[3] <- p[1]
              }
            }else{
              p[3] <- p[1]
            }
          }
        }
        
        rhom1[i,j] <- rho1
        rhom2[i,j] <- rho2
        pm[i,j] <- p[3]
      }
    }
    
    ptemp <- c(pm)
    ptemp <- ptemp[!is.na(ptemp)]
    cct <- sum(tan(pi*(0.5-ptemp)))/length(ptemp)
    p_global <- pcauchy(cct,lower.tail=FALSE)
    
    list(r1 = rhom1, r2 = rhom2, p_single=pm,p_global=p_global)
  }
}


