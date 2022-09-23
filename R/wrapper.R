
CGF = function(t,info)
{
  CGF_c(info,t)
}

del_CGF = function(t,info)
{
  del_CGF_c(info,t)
}

seq = function(t,info1,info2,zeta)
{
  seq_c(info1,info2,t,zeta)
}

solve_seq = function(zeta,info1,info2)
{
  suppressWarnings(rootSolve::multiroot(seq,rep(0,10),info1=info1,info2=info2,zeta=zeta,rtol = 1e-10,atol = 1e-10, ctol = 1e-10)$root)
}

solve_seq2 = function(zeta,info1,info2)
{
  BBsolve(rep(0,10),seq,info1=info1,info2=info2,zeta=zeta, control = list(tol = 1e-10),quiet=TRUE)$par
}

ell = function(z,info1,info2,n1,n2)
{
  t <- solve_seq(z,info1,info2)
  (CGF(t[1:5],info1)-sum(t[1:5]*z[1:5]))*n1 + (CGF(t[6:10],info2)-sum(t[6:10]*z[6:10]))*n2
}

g = function(z)
{
  g_c(z)
}

del_g = function(z)
{
  del_g_c(z)
}


lagrange = function(z_lam,w,info1,info2,n1,n2,solver=0)
{
  npa <- length(z_lam) - 1
  if(solver==0)
  {
    t <- solve_seq(z_lam[1:npa],info1,info2)
  }else{
    t <- solve_seq2(z_lam[1:npa],info1,info2)
  }
  t <- c(n1*t[1:5],n2*t[6:10])
  c(z_lam[length(z_lam)]*del_g(z_lam[1:npa])+t,g(z_lam[1:npa])-w)
}
  
cum_dist = function(w,info1,info2,x0,imp=FALSE,n1,n2,solver)
{
  invisible(capture.output( 
   solv1 <- tryCatch(
    {
      nleqslv(c(x0,0),lagrange,solver=0,w=w,info1=info1,info2=info2,n1=n1,n2=n2,method='Newton',global='hook',xscalm='auto',jacobian=FALSE,control=list(allowSingular=TRUE,maxit=500))
    },
    error=function(cond){NA}
    )
  ))
  
  if(!is.na(solv1[[1]][1]))
  { 
    if(max(abs(solv1$fvec))>0.01)
    {solv1 <- NA}
  }
  
  if(is.na(solv1[[1]][1]))
  {
    if(solver==TRUE)
    {
      solv1 <- tryCatch({
        nleqslv(c(x0,0),lagrange,solver=1,w=w,info1=info1,info2=info2,n1=n1,n2=n2,method='Newton',global='hook',xscalm='auto',jacobian=FALSE,control=list(allowSingular=TRUE,maxit=500))
      },
      error=function(cond){NA}
      )
    
      if(is.na(solv1[[1]][1]))
      {
        return(rep(NA,3))
      }
      
      if(max(abs(solv1$fvec))>0.01)
      {
        return(rep(NA,3))
      }
    }else{
      return(rep(NA,3))
    }
  }
  
  zetas <- solv1$x[1:10]
  #t <- solve_seq(zetas,info1,info2)
  t = solv1$fvec[1:10] - solv1$x[11]*del_g(zetas)
  t[1:5] = t[1:5]/n1
  t[6:10] = t[6:10]/n2
  
  if(max(abs(zetas-x0))<1e-3)
  {
    imp <- FALSE
  }
  
  jtol = 1e-10
  
  if(imp==TRUE)
  {
    grad = c(-n1*t[1:5],-n2*t[6:10])
    # hes = numDeriv::jacobian(solve_seq,zetas,info1=info1,info2=info2,'simple',method.args=list(eps=jtol))
    hes = jacift_c(info1,info2,t,zetas)
    hes[1:5,1:5] = n1*hes[1:5,1:5]
    hes[6:10,6:10] = n2*hes[6:10,6:10]
    
    grad_g = del_g(zetas)
    hes_g = numDeriv::jacobian(del_g,zetas,'simple',method.args=list(eps=jtol))
    mj = hes + solv1$x[11]*hes_g
    qw = as.numeric(t(grad_g)%*%solve(mj,grad_g))
    
    # hes0 <- numDeriv::jacobian(solve_seq,x0,info1=info1,info2=info2,'simple',method.args=list(eps=jtol))
    hes0 <- jacift_c(info1,info2,rep(0,10),x0)
    hes0[1:5,1:5] <- n1*hes0[1:5,1:5]
    hes0[6:10,6:10] <- n2*hes0[6:10,6:10]
    
    detjac <- determinant(mj)
    qw = detjac$sign*qw
    dy = qw*exp(detjac$modulus-determinant(hes0)$modulus)
   
    if(qw<=0)
    {
      imp <- FALSE
    }else{
      dw <- 1/sqrt(dy)
      bd0 <- 1/sqrt(abs(det(numDeriv::jacobian(del_CGF,rep(0,5),info=info1,'simple')))*abs(det(numDeriv::jacobian(del_CGF,rep(0,5),info=info2,'simple'))))
    }
  }
  
  # ell_y_tilde <- ell(zetas,info1,info2,n1,n2)
  ell_y_tilde <- (CGF(t[1:5],info1)-sum(t[1:5]*zetas[1:5]))*n1 + (CGF(t[6:10],info2)-sum(t[6:10]*zetas[6:10]))*n2
  # ell_y <- ell(x0,info1,info2,n1,n2)
  ell_y <- 0
  
  if((ell_y-ell_y_tilde)<0)
  {
    return(rep(NA,3))
  }
  
  r_w <- if(w==g(x0)){0} else{sign(w-g(x0))*sqrt(2*(ell_y-ell_y_tilde))}
  
  res <- pnorm(r_w)
  if(imp==TRUE)
  {
    bd <- 1/sqrt(abs(det(numDeriv::jacobian(del_CGF,t[1:5],info=info1,'simple')))*abs(det(numDeriv::jacobian(del_CGF,t[6:10],info=info2,'simple'))))
    # c = dw*grad_g[nz]/grad[nz]*bd/bd0
    c <- dw*bd/bd0/solv1$x[11]
    res <- c(res, pnorm(r_w+0.5*(1/r_w+c)*(3+r_w*c)),pnorm(r_w)+dnorm(r_w)*(1/r_w+c))
  }else{
    res <- c(res, NA,NA)
  }
  res
}


ld_z_ci_u = function(ds,r)
{
  # data = na.omit(data)
  # data = as.matrix(data)
  # r = cor(data[,1],data[,2])
  # z = (log(1+r)-log(1-r))/2
  # ds = scale(data)
  m22 <- mean(ds[,1]^2*ds[,2]^2)
  m31 <- mean(ds[,1]^3*ds[,2])
  m13 <- mean(ds[,1]*ds[,2]^3)
  m40 <- mean(ds[,1]^4)
  m04 <- mean(ds[,2]^4)
  n <- nrow(ds)
  var_z <- ((m40+2*m22+m04)*r^2-4*r*(m13+m31)+4*m22)/(4*(1-r^2)^2)
  return(var_z/n)
}
  