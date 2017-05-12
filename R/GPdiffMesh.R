GPdiffMesh <- function(x,X,meshT,Y,regFunX = function(x){  return(0*matrix(x[,1]))},
                       regFunT = function(x){  return(1*matrix(x[,1]))},
                       phi1 = 1,
                       sigma12 = 1,
                       sigma22 = 1,
                       phi2 = 1,
                       mybeta = FALSE,
                       l = 4) {
  ################
  #
  #X is the data points, T is the mesh size, y is the function value, x is the point want to predict
  #
  #For now it is linear regression
  #
  #MLE
  #
  ################
  #X <- matrix(X)
  nrowX <- dim(X)[1]
  d <- dim(X)[2]
  myreX <- regFunX(X)
  reDX <- dim(myreX)[2]
  myreT <- regFunT(cbind(meshT,X))
  reDT <- dim(myreT)[2]

  myphi1 = phi1
  mysigma22 = sigma22
  mysigma12 = sigma12
  myphi2 = phi2
  if(length(mybeta)== 1 && mybeta == FALSE){
    mybeta0 = 1
    mybeta1 = matrix(data = 1, nrow = reDX, ncol = 1)
    mybeta2 = matrix(data = 1, nrow = reDT, ncol = 1)
  }
  else{
    mybeta0 = 1
    mybeta1 = matrix(mybeta[1:reDX])
    mybeta2 = matrix(mybeta[(reDX+1):(reDX+reDT)])
  }


  myMLENSM <- function(pars){
    #function calculate MLE function for nonstationary Gaussian process
    phi1 <- pars[1]
    phi2 <- pars[2]
    sigma22 <- pars[3]
    sigma12 <- pars[4]
    beta0 <- pars[5]
    beta1 <- matrix(pars[6:(reDX+5)])
    beta2 <- matrix(pars[(reDX+6):(reDX+reDT+5)])
    t1 <- matrix(data = 1, nrow = nrowX, ncol = 1)
    mya <- X
    myh <- meshT
    mu <- beta0 * t1  + myreX %*% beta1 + myreT %*% beta2

    mySigmam <- matrix(data = 1, nrow = nrowX, ncol = nrowX)
    for(i in 1:nrowX){
      for(j in 1:nrowX){
        mySigmam[i,j] = sigma12^2 * exp(-phi1^2 * (sum((mya[i,]-mya[j,])^2))) + sigma22^2 * min(myh[i],myh[j])^l * exp(-phi2^2 * (sum((mya[i,]-mya[j,])^2)))
      }
    }
    myval <- 1/2*log(det(mySigmam))+1/2*t(Y-mu) %*% solve(mySigmam,Y-mu,tol=0)
    return(myval)
  }


  #b <- c(1,1/80,90,1,0,1,0)


  b <- c(myphi1, myphi2, mysigma22, mysigma12, mybeta0, mybeta1, mybeta2)
  myoptNSM <- stats::optim(b, myMLENSM, method = "Nelder-Mead")
  parss <- myoptNSM$par

  phi1 <- parss[1]
  phi2 <- parss[2]
  sigma22 <- parss[3]
  sigma12 <- parss[4]
  beta0 <- parss[5]
  beta1 <- parss[6:(reDX+5)]
  beta2 <- parss[(reDX+6):(reDX+reDT+5)]
  mySigmam <- matrix(data = 1, nrow = nrowX, ncol = nrowX)

  mya <- X
  h <- meshT

  for(i in 1:nrowX){
    for(j in 1:nrowX){
      mySigmam[i,j] = sigma12^2 * exp(-phi1^2 * (sum((mya[i,]-mya[j,])^2))) + sigma22^2 * min(h[i],h[j])^l * exp(-phi2^2 * (sum((mya[i,]-mya[j,])^2)))
    }
  }
  myouty <- matrix(data = 1, nrow = dim(x)[1], ncol = 1)
  mysigy <- matrix(data = 1, nrow = dim(x)[1], ncol = 1)

  for(i in 1:dim(x)[1]){
    myrr = matrix(data = 1, nrow = nrowX, ncol = 1)
    for(j in 1:nrowX){
      myrr[j] = sigma12^2 * exp(-phi1^2 * (sum((x[i,]-mya[j,])^2)))
    }
    temp = solve(mySigmam, myrr, tol=0)
    myCalm11 = Y - beta0 - myreX %*% beta1
    myouty[i] = beta0 +  regFunX(x[i,,drop=FALSE]) %*% beta1 + t(myCalm11) %*% temp
    mysigy[i] = sigma12^2 - t(myrr) %*% temp
  }

  estiP = list(phi1 = parss[1],
               phi2 = parss[2],
               sigma22 = parss[3],
               sigma12 = parss[4],
               beta0 = parss[5],
               beta1 = parss[6:(reDX+5)],
               beta2 = parss[(reDX+6):(reDX+reDT+5)])

  return(list(outy = myouty,sigy = mysigy,estipar = estiP))
}

#
# ####test functions ####
# hig02 <- function(s)
# {
#   y <- s*sin(s) / 10
#   return(y)
# }
# ####
# myX <- matrix(c(seq(from = 0,to = 10, by = 1),seq(from = 0,to = 10, by = 1)),ncol = 2)
# myy <- hig02(matrix(myX[,1]))
# myT <- matrix(c(0.01,0.5,0.01,0.02,0.02,0.01,0.01,0.02,0.002,0.003,0.03))
# ttttt <- function(x){
#   return(cbind((matrix(x[,1])^2*matrix(x[,2])),(matrix(x[,1])*matrix(x[,2]))) )
# }
# myregf <- function(x){
#   return(cbind(matrix(data = 1, nrow = dim(x)[1], ncol = 1),x))
# }
# x <- matrix(c(seq(from = 0,to = 10, by = 0.1),seq(from = 0,to = 10, by = 0.1)),ncol = 2)
# myploty <- hig02(matrix(x[,1]))
# y <- GPdiffMesh(x, myX, myT, myy, regFunT = ttttt, regFunX = myregf)
# #y$outy
# #y$sigy
# y$estipar
# plot(x[,1], myploty,"l")
# lines(x[,1],y$outy, type="o", pch=22, lty=2, col="red")
# #plot(x,y)
