library(L0TFinv)
library(AMIAS)
library(genlasso)
library(wbs)
library(not)
library(changepoint)
library(cpop)
library(ggpubr)
library(ggplot2)
library(foreach)
library(doParallel)

### Different trend filtering methods---------------------------

tau1 = c(0.1, 0.25, 0.3, 0.35, 0.4, 0.55, 0.65, 0.7, 0.85, 0.95)
h1 = c(-1, 5, 3, 0, -1, 2, 0, -1, -3, 2, 4)

tau2 = c(17, 25, 40, 50, 58, 69, 90, 100, 110, 115)/120
h2 = c(0, 2, -2, 1, -1, 4, 6, 2, -1, -4, -2)

tau3 = c(256, 512, 768, 1024, 1152, 1280, 1344)/1408
h3 = c(0, 14, -14, 28, -28, 42, -42, 56)

EvalMetrics <- function(beta, cpts=NULL, y0, tcpt = NULL){
  mse <- mean((beta-y0)^2)
  mad <- mean(abs(beta-y0))
  if(is.null(tcpt)){
    tab <- data.frame(MSE=mse, MAD=mad)
  }else{
    n <- length(beta)
    n.cpts <- length(cpts)
    segments.endpoints.true <- sort(unique(tcpt))
    segments.endpoints.est <- sort(unique(cpts))
    distm <- abs(matrix(rep(segments.endpoints.est, length(segments.endpoints.true)), nrow=length(segments.endpoints.est))
                 -matrix(rep(segments.endpoints.true, length(segments.endpoints.est)), nrow=length(segments.endpoints.est), byrow=TRUE))
    screening.dist <- max(apply(distm, 2, min)) * 100 / n # min distance for each true cpt
    precision.dist <- max(apply(distm, 1, min)) * 100 / n # min distance for each estimated cpt
    haus.dist <- max(screening.dist, precision.dist)
  }
  tab <- data.frame(MSE=mse, MAD=mad, dH=haus.dist, nknot=n.cpts,dP=precision.dist,dS=screening.dist)
  return(tab)
}

SimuL0Inv <- function(q = 0, n = 300, sigma = 0.1, seed = NA, tau = tau, h = h){
  if (q == 0){
    data = L0TFinv::SimuBlocksInv(n=n, sigma=sigma, seed=seed,tau = tau, h = h)
    kmax=20;
  }
  if (q == 1){
    data = L0TFinv::SimuWaveInv(n=n, sigma=sigma, seed=seed,tau = tau, h = h, a0 = 0)
    kmax=20;
  }
  tic = Sys.time()
  resL0 = L0TFinv::L0TFinv.opt(as.numeric(data$y),kmax=kmax,q=q,first = 0,last=1,penalty="bic")
  toc = Sys.time()
  metric = EvalMetrics(resL0$yopt, resL0$Aopt/n, data$y0, data$tau)
  metric$time = as.numeric(toc-tic,  units = "secs")
  return(metric)
}

SimuL0AMIAS <- function(q = 0, n = 300, sigma = 0.1, seed = NA, tau = tau, h = h){
  if (q == 0){
    data = L0TFinv::SimuBlocksInv(n=n, sigma=sigma, seed=seed,tau = tau, h = h)
    kmax=20;
    tic = Sys.time()
    resL0 = AMIAS::samias(as.numeric(data$y),kmax=kmax,q=q,D_type = "tf0")
    toc = Sys.time()
  }
  if (q == 1){
    data = L0TFinv::SimuWaveInv(n=n, sigma=sigma, seed=seed,tau = tau, h = h, a0 = 0)
    kmax=20;
    tic = Sys.time()
    resL0 = AMIAS::samias(as.numeric(data$y),kmax=kmax,q=q,D_type = "tfq")
    toc = Sys.time()
  }
  metric = EvalMetrics(resL0$alpha, resL0$A/n, data$y0, data$tau)
  metric$time = as.numeric(toc-tic,  units = "secs")
  return(metric)
}

SimuCpop <- function(q = 0, n = 300, sigma = 0.1, seed = NA, tau = tau, h = h){
  if (q == 1){
    data = L0TFinv::SimuWaveInv(n=n, sigma=sigma, seed=seed,tau = tau, h = h, a0 = 0)
    kmax=20;
  }
  tic = Sys.time()
  resCpop = cpop::cpop(data$y)
  toc = Sys.time()
  cpt <- resCpop@changepoints
  if(length(cpt)==2){
    cpt <- c()
  }else{
    cpt <- cpt[2:(length(cpt)-1)]
  }
  metric = EvalMetrics(cpop::estimate(resCpop)[,2], cpt/n, data$y0, data$tau)
  metric$time = as.numeric(toc-tic,  units = "secs")
  return(metric)
}

L1TF <- function(data, n, q, maxdf) {
  resL1  <- genlasso::trendfilter(pos=data$x, y=data$y, ord=q)
  idx <- which(resL1$df<=maxdf)
  bicL1 <- apply(resL1$beta[,idx], 2, function(beta) n*log(mean((data$y-beta)^2))) + 2*log(n)*resL1$df[idx]
  betaL1 <- resL1$beta[,which.min(bicL1)]
  knotL1 <- data$x[which(abs(diff(betaL1, diff=q+1))>1e-5)+1]
  return(list(beta=betaL1, knot=knotL1))
}

SimuL1 <- function(q = 0, n = 300, sigma = 0.1, seed = NA, tau = tau, h = h){
  if (q == 0){
    data = L0TFinv::SimuBlocksInv(n=n, sigma=sigma, seed=seed,tau = tau, h = h)
    kmax=20;
  }
  if (q == 1){
    data = L0TFinv::SimuWaveInv(n=n, sigma=sigma, seed=seed,tau = tau, h = h, a0 = 0)
    kmax=20;
  }
  tic = Sys.time()
  resL1 = L1TF(data, n, q, maxdf=kmax)
  toc = Sys.time()
  metric = EvalMetrics(resL1$beta, resL1$knot, data$y0, data$tau)
  metric$time = as.numeric(toc-tic,  units = "secs")
  return(metric)
}

Simuwbs <- function(q = 0, n = 300, sigma = 0.1, seed = NA, tau = tau, h = h){
  if (q == 0){
    data = L0TFinv::SimuBlocksInv(n=n, sigma=sigma, seed=seed,tau = tau, h = h)
    kmax=20;
  }
  if (q == 1){
    data = L0TFinv::SimuWaveInv(n=n, sigma=sigma, seed=seed,tau = tau, h = h, a0 = 0)
    kmax=20;
  }
  tic = Sys.time()
  resWbs = wbs::changepoints(wbs::wbs(data$y),penalty="ssic.penalty")
  toc = Sys.time()
  obj <- resWbs$cpt.th[[1]]
  metric = EvalMetrics(wbs::means.between.cpt(data$y, obj), obj/n, data$y0, data$tau)
  metric$time = as.numeric(toc-tic,  units = "secs")
  return(metric)
}

SimuPelt <- function(q = 0, n = 300, sigma = 0.1, seed = NA, tau = tau, h = h){
  if (q == 0){
    data = L0TFinv::SimuBlocksInv(n=n, sigma=sigma, seed=seed,tau = tau, h = h)
    kmax=20;
  }
  if (q == 1){
    data = L0TFinv::SimuWaveInv(n=n, sigma=sigma, seed=seed,tau = tau, h = h, a0 = 0)
    kmax=20;
  }
  tic = Sys.time()
  resPelt = changepoint::cpt.mean(data$y,method="PELT",class=TRUE,penalty="SIC",Q=kmax)
  toc = Sys.time()
  obj = resPelt@cpts
  lis = c()
  value = resPelt@param.est$mean
  pivot = obj[1]
  for(ii in 1:length(obj)){
    lis <- c(lis,rep(value[ii], pivot))
    pivot <- obj[ii+1] - obj[ii]
  }
  metric = EvalMetrics(lis, obj[1:(length(obj)-1)]/n, data$y0, data$tau)
  metric$time = as.numeric(toc-tic,  units = "secs")
  return(metric)
}

SimuNot <- function(q = 0, n = 300, sigma = 0.1, seed = NA, tau = tau, h = h){
  if (q == 0){
    data = L0TFinv::SimuBlocksInv(n=n, sigma=sigma, seed=seed,tau = tau, h = h)
    kmax=20;
    tic = Sys.time()
    resNot = not::not(data$y, contrast="pcwsConstMean")
    toc = Sys.time()
  }
  if (q == 1){
    data = L0TFinv::SimuWaveInv(n=n, sigma=sigma, seed=seed,tau = tau, h = h, a0 = 0)
    kmax=20;
    tic = Sys.time()
    resNot = not::not(data$y, contrast="pcwsLinMean")
    toc = Sys.time()
  }
  obj = not::features(resNot,penalty="sic",q.max=kmax)$cpt
  lis = predict(resNot)
  metric = EvalMetrics(lis, obj/n, data$y0, data$tau)
  metric$time = as.numeric(toc-tic,  units = "secs")
  return(metric)
}


## Table A1---------------------------

iter = 100
nn = c(100,1000,10000)
ResultL0Inv <- ResultL0AMIAS <- ResultCpop <- NULL
sigma = 0.2
cl <- makeCluster(20)
registerDoParallel(cl)

for (n in nn) {
  print(n)
  TabL1 = foreach(mcseed = 1:iter, .combine = rbind) %dopar%
    SimuL0Inv(q = 1, n = n, sigma = sigma, seed = mcseed, tau = tau1, h = h1)
  TabL2 = foreach(mcseed = 1:iter, .combine = rbind) %dopar%
    SimuL0AMIAS(q = 1, n = n, sigma = sigma, seed = mcseed, tau = tau1, h = h1)
  TabL3 = foreach(mcseed = 1:iter, .combine = rbind) %dopar%
    SimuCpop(q = 1, n = n, sigma = sigma, seed = mcseed, tau = tau1, h = h1)
  ResultL0Inv <- rbind(ResultL0Inv, c(colMeans(TabL1),apply(TabL1,2,sd)))
  ResultL0AMIAS <- rbind(ResultL0AMIAS, c(colMeans(TabL2),apply(TabL2,2,sd)))
  ResultCpop <- rbind(ResultCpop, c(colMeans(TabL3),apply(TabL3,2,sd)))
}
stopCluster(cl)


iter = 100
nn = c(500,5000,50000)
ResultL0Inv <- ResultL0AMIAS <- NULL
sigma = 0.2
cl <- makeCluster(20)
registerDoParallel(cl)

for (n in nn) {
  print(n)
  TabL1 = foreach(mcseed = 1:iter, .combine = rbind) %dopar%
    SimuL0Inv(q = 0, n = n, sigma = sigma, seed = mcseed, tau = tau1, h = h1)
  TabL2 = foreach(mcseed = 1:iter, .combine = rbind) %dopar%
    SimuL0AMIAS(q = 0, n = n, sigma = sigma, seed = mcseed, tau = tau1, h = h1)
  ResultL0Inv <- rbind(ResultL0Inv, c(colMeans(TabL1),apply(TabL1,2,sd)))
  ResultL0AMIAS <- rbind(ResultL0AMIAS, c(colMeans(TabL2),apply(TabL2,2,sd)))
}
stopCluster(cl)


## Table A2---------------------------

iter = 100
n = 500
sig = c(1,2)
ResultL0Inv <- ResultL0AMIAS <- ResultL1TF <- ResultNot <- ResultWbs <- ResultPelt <- NULL
cl <- makeCluster(20)
registerDoParallel(cl)

for (sigma in sig) {
  TabL1 = foreach(mcseed = 1:iter, .combine = rbind) %dopar%
    SimuL0Inv(q = 0, n = n, sigma = sigma, seed = mcseed, tau = tau2, h = h2)
  TabL2 = foreach(mcseed = 1:iter, .combine = rbind) %dopar%
    SimuL0AMIAS(q = 0, n = n, sigma = sigma, seed = mcseed, tau = tau2, h = h2)
  TabL3 = foreach(mcseed = 1:iter, .combine = rbind) %dopar%
    SimuL1(q = 0, n = n, sigma = sigma, seed = mcseed, tau = tau2, h = h2)
  TabL4 = foreach(mcseed = 1:iter, .combine = rbind) %dopar%
    Simuwbs(q = 0, n = n, sigma = sigma, seed = mcseed, tau = tau2, h = h2)
  TabL5 = foreach(mcseed = 1:iter, .combine = rbind) %dopar%
    SimuNot(q = 0, n = n, sigma = sigma, seed = mcseed, tau = tau2, h = h2)
  TabL6 = foreach(mcseed = 1:iter, .combine = rbind) %dopar%
    SimuPelt(q = 0, n = n, sigma = sigma, seed = mcseed, tau = tau2, h = h2)
  ResultL0Inv <- rbind(ResultL0Inv, c(colMeans(TabL1),apply(TabL1,2,sd)))
  ResultL0AMIAS <- rbind(ResultL0AMIAS, c(colMeans(TabL2),apply(TabL2,2,sd)))
  ResultL1TF <- rbind(ResultL1TF, c(colMeans(TabL3),apply(TabL3,2,sd)))
  ResultWbs <- rbind(ResultWbs, c(colMeans(TabL4),apply(TabL4,2,sd)))
  ResultNot <- rbind(ResultNot, c(colMeans(TabL5),apply(TabL5,2,sd)))
  ResultPelt <- rbind(ResultPelt, c(colMeans(TabL6),apply(TabL6,2,sd)))
}
stopCluster(cl)


## Figure 3---------------------------

iter = 100
n = 1000
sig = c(1,2)
ResultL0Inv <- ResultL0AMIAS <- ResultL1TF <- ResultNot <- ResultWbs <- ResultPelt <- ResultCpop <- NULL
cl <- makeCluster(20)
registerDoParallel(cl)

for (sigma in sig) {
  TabL1 = foreach(mcseed = 1:iter, .combine = rbind) %dopar%
    SimuL0Inv(q = 1, n = n, sigma = sigma, seed = mcseed, tau = tau3, h = h3)
  TabL2 = foreach(mcseed = 1:iter, .combine = rbind) %dopar%
    SimuL0AMIAS(q = 1, n = n, sigma = sigma, seed = mcseed, tau = tau3, h = h3)
  TabL3 = foreach(mcseed = 1:iter, .combine = rbind) %dopar%
    SimuL1(q = 1, n = n, sigma = sigma, seed = mcseed, tau = tau3, h = h3)
  TabL4 = foreach(mcseed = 1:iter, .combine = rbind) %dopar%
    Simuwbs(q = 1, n = n, sigma = sigma, seed = mcseed, tau = tau3, h = h3)
  TabL5 = foreach(mcseed = 1:iter, .combine = rbind) %dopar%
    SimuNot(q = 1, n = n, sigma = sigma, seed = mcseed, tau = tau3, h = h3)
  TabL6 = foreach(mcseed = 1:iter, .combine = rbind) %dopar%
    SimuPelt(q = 1, n = n, sigma = sigma, seed = mcseed, tau = tau3, h = h3)
  TabL7 = foreach(mcseed = 1:iter, .combine = rbind) %dopar%
    SimuCpop(q = 1, n = n, sigma = sigma, seed = mcseed, tau = tau3, h = h3)
  ResultL0Inv <- rbind(ResultL0Inv, c(colMeans(TabL1),apply(TabL1,2,sd)))
  ResultL0AMIAS <- rbind(ResultL0AMIAS, c(colMeans(TabL2),apply(TabL2,2,sd)))
  ResultL1TF <- rbind(ResultL1TF, c(colMeans(TabL3),apply(TabL3,2,sd)))
  ResultWbs <- rbind(ResultWbs, c(colMeans(TabL4),apply(TabL4,2,sd)))
  ResultNot <- rbind(ResultNot, c(colMeans(TabL5),apply(TabL5,2,sd)))
  ResultPelt <- rbind(ResultPelt, c(colMeans(TabL6),apply(TabL6,2,sd)))
  ResultCpop <- rbind(ResultCpop, c(colMeans(TabL7),apply(TabL7,2,sd)))
}
stopCluster(cl)


MSEL0Inv <- TabL1$MSE
MSEL0AMIAS <- TabL2$MSE
MSEL1 <- TabL3$MSE
MSEwbs <- TabL4$MSE
MSEnot <- TabL5$MSE
MSEpelt <- TabL6$MSE
MSEcpop <- TabL7$MSE
name <- c(rep("L0TFinv",iter),rep("AMIAS",iter),rep("genlasso",iter),rep("wbs",iter),
          rep("not",iter),rep("changepoint",iter),rep("cpop",iter))
mse <- c(MSEL0Inv,MSEL0AMIAS,MSEL1,MSEwbs,MSEnot,MSEpelt,MSEcpop)
data <- data.frame(group=name,MSE=mse)
ggboxplot(data, x = "group", y = "MSE", fill = "group",palette = "npg",add = "jitter")


dHL0Inv <- TabL1$dH
dHL0AMIAS <- TabL2$dH
dHL1 <- TabL3$dH
dHwbs <- TabL4$dH
dHnot <- TabL5$dH
dHpelt <- TabL6$dH
dHcpop <- TabL7$dH
name <- c(rep("L0TFinv",iter),rep("AMIAS",iter),rep("genlasso",iter),rep("wbs",iter),
          rep("not",iter),rep("changepoint",iter),rep("cpop",iter))
dH <- c(dHL0Inv,dHL0AMIAS,dHL1,dHwbs,dHnot,dHpelt,dHcpop)
data <- data.frame(group=name,dH=dH)
ggboxplot(data, x = "group", y = "dH", fill = "group",palette = "npg",add = "jitter")


nknotL0Inv <- TabL1$nknot
nknotL0AMIAS <- TabL2$nknot
nknotL1 <- TabL3$nknot
nknotwbs <- TabL4$nknot
nknotnot <- TabL5$nknot
nknotpelt <- TabL6$nknot
nknotcpop <- TabL7$nknot
name <- c(rep("L0TFinv",iter),rep("AMIAS",iter),rep("genlasso",iter),rep("wbs",iter),
          rep("not",iter),rep("changepoint",iter),rep("cpop",iter))
nknot <- c(nknotL0Inv,nknotL0AMIAS,nknotL1,nknotwbs,nknotnot,nknotpelt,nknotcpop)
data <- data.frame(group=name,nknot=nknot)
ggboxplot(data, x = "group", y = "nknot", fill = "group",palette = "npg",add = "jitter")


timeL0Inv <- TabL1$time
timeL0AMIAS <- TabL2$time
timeL1 <- TabL3$time
timewbs <- TabL4$time
timenot <- TabL5$time
timepelt <- TabL6$time
timecpop <- TabL7$time
name <- c(rep("L0TFinv",iter),rep("AMIAS",iter),rep("genlasso",iter),rep("wbs",iter),
          rep("not",iter),rep("changepoint",iter),rep("cpop",iter))
time <- c(timeL0Inv,timeL0AMIAS,timeL1,timewbs,timenot,timepelt,timecpop)
data <- data.frame(group=name,time=time)
ggboxplot(data, x = "group", y = "time", fill = "group",palette = "npg",add = "jitter")







