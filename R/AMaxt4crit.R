#' Find test statistic value and corresponding critical value
#'
#' More detailed description
#'
#' @param g1 a real vector
#' @param g2 a real vector
#' @param g3 a real vector
#' @param g4 a real vector
#' @param alpha a real number between 0 and 1
#'
#' @return a numerical vector
#'
#' @examples
#' g1=c(1,2,4,1,1,3,4,0,2,4,1,3,5,2,1,4)
#' g2=c(0.4,2,-4,0.9,4.6,8,9.8,9,4,5,6,3)
#' g3=c(0,1,2,5,7,3,7.9,3,2)
#' g4=c(0,-3,-4,5,1,3,5,7,2,1,2,0)
#' AMaxt4crit(g1,g2,g3,g4,0.05)
#'
#' @export
AMaxt4crit<-function(g1,g2,g3,g4,alpha)
{
  fun1<-function(n1,n2,n3,n4,v1,v2,v3,v4)
  {
    N=n1+n2+n3+n4
    nu1=n1/N;nu2=n2/N;nu3=n3/N;nu4=n4/N
    d1=sqrt(nu1*v2+nu2*v1);d2=sqrt(nu2*v3+nu3*v2);d3=sqrt(nu3*v4+nu4*v3)
    D=matrix(c(d1,0,0,0,d2,0,0,0,d3),nrow=3,byrow=TRUE)
    s1=sqrt(nu1*nu3)*v2;s2=sqrt(nu2*nu4)*v3
    S=matrix(c(d1^2,-s1,0,-s1,d2^2,-s2,0,-s2,d3^2),nrow=3,byrow=TRUE)
    P=solve(D)%*%S%*%solve(D);mu=c(0,0,0)
    T<-mvrnorm(1,mu,P)
    A=max(T[1],T[2],T[3],na.rm = FALSE)
    return(A)
  }
  fun2<-function(n1,n2,n3,n4,alpha,v1,v2,v3,v4)
  {
    x<-replicate(5000,fun1(n1,n2,n3,n4,v1,v2,v3,v4))
    y<-sort(x,decreasing=FALSE)
    m=(1-alpha)*5000
    c<-y[m]
    return(c)
  }
  fun3<-function(n1,n2,n3,n4,alpha,v1,v2,v3,v4)
  {
    z=replicate(10,fun2(n1,n2,n3,n4,alpha,v1,v2,v3,v4))
    cri=mean(z)
    return(cri)
  }
  fun4<-function(g1,g2,g3,g4)
  {
    X1=mean(g1);X2=mean(g2);X3=mean(g3);X4=mean(g4)
    S1=var(g1);S2=var(g2);S3=var(g3);S4=var(g4)
    n1=length(g1);n2=length(g2);n3=length(g3);n4=length(g4)
    v1<-sqrt(S1/n1+S2/n2);v2<-sqrt(S2/n2+S3/n3);v3<-sqrt(S3/n3+S4/n4)
    T1=(X2-X1)/v1;T2=(X3-X2)/v2;T3=(X4-X3)/v3
    T=max(T1,T2,T3,na.rm = FALSE)
    return(T)
  }
  v1=var(g1);v2=var(g2);v3=var(g3);v4=var(g4)
  n1=length(g1);n2=length(g2);n3=length(g3);n4=length(g4)
  set.seed(5)
  statistic_value<-fun4(g1,g2,g3,g4)
  crit_value=fun3(n1,n2,n3,n4,alpha,v1,v2,v3,v4)
  result=c(statistic_value, crit_value)
  return(result)
}
#data<-lapply(data1, function(col)col[!is.na(col)])
