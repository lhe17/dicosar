// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]


//
// via the exports attribute we tell Rcpp to make this function
// available from R
//
// [[Rcpp::export]]
double CGF_c(const Eigen::Map<Eigen::MatrixXd> & info,const Eigen::VectorXd & t) {
  
  int nr = info.rows();
  double res = (info*t).array().exp().sum();
  res = log(res/nr);
  return res;
  //log(sum(exp(as.matrix(info)%*%t[1:5]))/n)
}

//
// [[Rcpp::export]]
Eigen::VectorXd del_CGF_c(const Eigen::Map<Eigen::MatrixXd> & info,const Eigen::VectorXd & t) {
  
  Eigen::VectorXd res = (info*t).array().exp();
  double deno = res.sum();
  res = info.transpose()*res;
  return res/deno;
  //ext = exp(as.matrix(info)%*%t[1:5])
  //deno = sum(ext)
  //t(info)%*%ext/deno
}

//
// [[Rcpp::export]]
double g_c(const Eigen::VectorXd & z) {
  
  double r = (z[4]-z[0]*z[1])/sqrt(z[2]-z[0]*z[0])/sqrt(z[3]-z[1]*z[1]);
  double z1 = (log(1+r)-log(1-r))/2;
  r = (z[9]-z[5]*z[6])/sqrt(z[7]-z[5]*z[5])/sqrt(z[8]-z[6]*z[6]);
  double z2 = (log(1+r)-log(1-r))/2;
  double zd = z1 - z2;
  return zd;
}


// [[Rcpp::export]]
Eigen::VectorXd del_g_c(const Eigen::VectorXd & z) {
  
  Eigen::VectorXd dz(10);
  double r = (z[4]-z[0]*z[1])/sqrt(z[2]-z[0]*z[0])/sqrt(z[3]-z[1]*z[1]);
  r = 1-r*r;
  dz(0) = (z[4]*z[0]-z[2]*z[1])/sqrt(z[3]-z[1]*z[1])/pow(z[2]-z[0]*z[0],1.5)/r;
  dz(1) = -(z[3]*z[0]-z[4]*z[1])/sqrt(z[2]-z[0]*z[0])/pow(z[3]-z[1]*z[1],1.5)/r;
  dz(2) = -0.5*(z[4]-z[0]*z[1])/sqrt(z[3]-z[1]*z[1])/pow(z[2]-z[0]*z[0],1.5)/r;
  dz(3) = -0.5*(z[4]-z[0]*z[1])/sqrt(z[2]-z[0]*z[0])/pow(z[3]-z[1]*z[1],1.5)/r;
  dz(4) = 1/sqrt(z[3]-z[1]*z[1])/sqrt(z[2]-z[0]*z[0])/r;
  
  double r2 = (z[9]-z[5]*z[6])/sqrt(z[7]-z[5]*z[5])/sqrt(z[8]-z[6]*z[6]);
  r2 = 1-r2*r2;
  dz(5) = -(z[9]*z[5]-z[7]*z[6])/sqrt(z[8]-z[6]*z[6])/pow(z[7]-z[5]*z[5],1.5)/r2;
  dz(6) = (z[8]*z[5]-z[9]*z[6])/sqrt(z[7]-z[5]*z[5])/pow(z[8]-z[6]*z[6],1.5)/r2;
  dz(7) = 0.5*(z[9]-z[5]*z[6])/sqrt(z[8]-z[6]*z[6])/pow(z[7]-z[5]*z[5],1.5)/r2;
  dz(8) = 0.5*(z[9]-z[5]*z[6])/sqrt(z[7]-z[5]*z[5])/pow(z[8]-z[6]*z[6],1.5)/r2;
  dz(9) = -1/sqrt(z[8]-z[6]*z[6])/sqrt(z[7]-z[5]*z[5])/r2;
  
  return dz;
}


// [[Rcpp::export]]
Eigen::VectorXd seq_c(const Eigen::Map<Eigen::MatrixXd> & info1,const Eigen::Map<Eigen::MatrixXd> & info2,
             const Eigen::VectorXd & t,const Eigen::VectorXd & zeta) {
  
  Eigen::VectorXd s(10);
  Eigen::VectorXd temp = (info1*t.head(5)).array().exp();
  for(int i=0;i<5;i++)
  {
    s(i) = temp.dot((info1.col(i).array() - zeta(i)).matrix());
  }
  
  temp = (info2*t.tail(5)).array().exp();
  for(int i=5;i<10;i++)
  {
    s(i) = temp.dot((info2.col(i-5).array() - zeta(i)).matrix());
  }
    
  return s;
  
}
