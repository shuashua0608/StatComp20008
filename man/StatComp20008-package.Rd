\name{StatComp20008-package}
\alias{StatComp20008-package}
\alias{StatComp20008}
\docType{package}
\title{
 Estimate Bart Simpson density by histogram and compare different bins.
}
\description{
  Estimate Bart Simpson density f(x)=1/2 phi(x;0,1)+1/10 sum_{j=0}^4 phi(x;j/2-1,1/10) by histogram and naive density estimator using different bandwidth.
}
\details{
  Function BS(x) draw ture density of Bart Simpson density function and function HBS(n,h) uses n random points generated from BS density with bandwidth of h. NB(x,y,h) uses naive estimator with x being estimated points and y being samples generated from BS density,h being bandwidth.
}
\author{
Yangshuhua.

Maintainer: Yangshuhua <shuashua0608@mail.ustc.edu.cn>
}
\references{
 https://www.bb.ustc.edu.cn/bbcswebdav/pid-66546-dt-content-rid-7067954_1/courses/STAT6125P.01.2020FA/Lec7.pdf
 ----Nonparametric Statistics
}
\keyword{ package }

