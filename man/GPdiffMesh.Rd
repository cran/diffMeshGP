\name{GPdiffMesh}
\alias{GPdiffMesh}
\docType{package}
\title{Multi-Fidelity Computer Experiments Using The Tuo-Wu-Yu Model}
\description{
This R function implements the nonstationary Kriging model proposed by Tuo, Wu and Yu (2014) for analyzing multi-fidelity computer outputs. This function computes the maximum likelihood estimates for the model parameters as well as the predictive means and variances of the exact solution (i.e., the conceptually highest fidelity).
}
\usage{
GPdiffMesh(x,X,meshT,Y,regFunX = function(x){  return(0*matrix(x[,1]))},
                       regFunT = function(x){  return(1*matrix(x[,1]))},
                       phi1 = 1,
                       sigma12 = 1,
                       sigma22 = 1,
                       phi2 = 1,
                       mybeta = FALSE,
                       l = 4)
}
\arguments{
\item{x}{An \eqn{(n_0 x d)} matrix of control points on which the values of the exact solution will be predicted. Each row represents one input point.}
        \item{X}{An \eqn{(n x d)} matrix of control variables.}
        \item{meshT}{An \eqn{(n x 1)} mesh density matrix. A same row of X and meshT represents the inputs from one trial of the computer experiment.}
        \item{Y}{An \eqn{(n x 1)} computer output matrix corresponding to (X, meshT).}
        \item{regFunX}{A scalar or a vector-valued regression function \eqn{f_1(x)}. The default value is \eqn{f_1(x)=1}. See Details below.}
        \item{regFunT}{A scalar or a vector-valued regression function \eqn{f_2(x,t)}. The default value is \eqn{f_2(x,t) = t}. See Details below.}
        \item{phi1,sigma12,sigma22,phi2,mybeta}{Initial values of the parameters for the maximum likelihood estimation. The default values are \eqn{phi1 = 1, sigma12 = 1, sigma22 = 1, phi2 = 1, mybeta = 1}. See Details below.}
        \item{l}{The value of parameter \eqn{l} in the covariance function of the nonstationary Kriging model. (We used fixed \eqn{l} in this R function.) The default value is \eqn{l=4}. See Details below.}
}
\value{
\item{outy}{An \eqn{(n_0 x 1)} matrix of predictive means corresponding to control points. A same row of x and outy represents a pair of one control point and corresponding predictive mean.}
        \item{sigy}{An \eqn{(n_0 x 1)} matrix of predictive variances corresponding to control points. A same row of x and sigy represents a pair of one control point and corresponding predictive variances.}
        \item{estipar}{A list of parameter estimates.}
}
\details{
This R function implements the nonstationary Kriging model proposed by Tuo, Wu and Yu (2014) for the modeling and analysis of multi-fidelity computer experiments. Denote an input-output pair from the computer simulation by \eqn{(x, t, y)}, where \eqn{x} is the vector of input variables, \eqn{t} is the mesh density, and \eqn{y} is the corresponding computer output. Tuo, Wu and Yu (2014) use the following Gaussian process model:

        \eqn{
        y = \eta(x, t) = \eta(x, 0) + \delta(x, t),
        }

        where \eqn{\eta(x, 0)} and \eqn{\delta(x, t)} are realizations of two mutually independent Gaussian processes \eqn{V(x)} and \eqn{Z(x, t)}. Assume

        \eqn{
        E(V(x))=f_{1}^{T}(x)\beta_1, E(Z(x, t))=f^{T}_{2}(x,t)\beta_2,
        }

        where \eqn{f_{1}(x)} can be set in regFunX, and \eqn{f_{2}(x,t)} can be set in regFunT. Both \eqn{f_{1}(x)} and \eqn{f_{2}(x,t)} can be vector-valued function. \eqn{f_{2}(x,t)} should satisfy for any \eqn{x}, \eqn{\lim_{t\rightarrow 0}f_{2}(x,t) = 0}. The default functions are \eqn{f_{1}(x) = 1} and \eqn{f_{2}(x,t) = t}.
        Assume the covariance of \eqn{V(x)} is

        \eqn{
            Cov(V(x_{1}),V(x_{2}))=\sigma_{1}^{2}e^{-\phi^2_1\|x_{1}-x_{2}\|^{2}_{2}},
        }

         and the covariance of \eqn{Z(x, t)} is

         \eqn{
            Cov(Z(x_{1}, t_{1}),Z(x_{2}, t_{2}))=\sigma_{2}^{2}\min(t_{1},t_{2})^{l}e^{-\phi_{2}^{2}\|x_{1}-x_{2}\|^{2}_{2}},
         }

         where \eqn{l} is a fixed parameter and is not estimated in this function.

          The goal of the Tuo-Wu-Yu model is to predict for the \eqn{exact solution} \eqn{\varphi(.):=\eta(.,0)}.
          For a set of control points \eqn{x=(x_{1},\ldots,x_{n_0})}, this R function predicts \eqn{\eta(x_{i},0)=\varphi(x_{i})} for \eqn{i=1,\ldots,n_{0}}.

    This R function uses maximum likelihood method to estimate the model parameters.
          Nelder-Mead method is used to maximize the likelihood function.
          The solution may depend on the choice of initial values.
          Users may specify the initial values or use the default values.
}
\examples{
hig02 <- function(s)
{
  #The test function is from [2].
  y <- s*sin(s) / 10
  return(y)
}
myX <- matrix(c(seq(from = 0,to = 10, by = 1),
    seq(from = 0,to = 10, by = 1)),ncol = 2)
myy <- hig02(matrix(myX[,1]))
myT <- matrix(c(0.01,0.5,0.01,0.02,0.02,0.01,0.01,0.02,0.002,
    0.003,0.03))
myregf <- function(x){
  return(x)
}
myregfn <- function(s){
  return(cbind((matrix(s[,1])^2*matrix(s[,2])),
    (matrix(s[,1])*matrix(s[,2]))))
}
#Here s=cbind(t,x), where x is a matrix of input variables and
#     t is the corresponding mesh density matrix.
x <- matrix(c(seq(from = 0,to = 10, by = 0.1),
    seq(from = 0,to = 10, by = 0.1)),ncol = 2)
myploty <- hig02(matrix(x[,1]))
y <- GPdiffMesh(x, myX, myT, myy, regFunX = myregf, regFunT = myregfn)
# The regression function is beta_0 + beta_1 x + (t^2x, tx)' beta_2
y$outy
y$sigy
y$estipar
plot(x[,1], myploty,"l")
lines(x[,1],y$outy, type="o", pch=22, lty=2, col="red")
}
