cluster_3D_BIC <- function(y,times,kmax=20,kint=10,maxit=2e2,max.order=7,plott=TRUE){
  
  
  m <- lapply(2:kmax,function(c) cluster_3D(y,c,times,kint,maxit,max.order))
  
  BIC <- sapply(m,function(c)c$BIC)
  print(BIC)
  if(plott==TRUE){
    plot(2:kmax,BIC,type="l",col="red",lwd=2,xlab="Number of Clusters")
  }
  
  ii <- c(2:kmax)[which.min(BIC)]
  
  ret <- m[[ii-1]]
  
  return(ret)
  
}

cluster_3D <- function(y,k,times,kint=10,maxit = 2e2,max.order=7){
  
  
  if(!is.list(times)){
    stop("Parameters times is not a list!")
  }
  n <- dim(y)[1]
  d <- dim(y)[2]
  #取的是ct有多少c个组成
  nt <- length(times)
  #一共有多少个时间点
  nts <- length(unlist(times))
  cs <- sapply(1:nt,function(i,x){
    length(x[[i]])
  },x=times)
  #cumsum参数传递的向量的累积和，12：1，3
  css <- c(0,cumsum(cs))
  
  init.va <- c()
  init.par1 <- list()
  for(i in 1:kint){
    #最优参数等结果
    init.t <- get.int(y, k,times)
    #模型预测值
    init.v <- get.mu.L(init.t$initial_mu_params,times,k)
    #最佳参数方差、系数
    init.cov <-init.t$initial_cov_params
    init.sigma <- init.cov[1:nt];init.rho <- init.cov[-c(1:nt)]
    #总协方差矩阵
    init.covv <- get.SAD1(init.sigma,init.rho,times)
    init.cond <- sapply(1:k, function(c) dmvnorm(y, init.v[c,], init.covv)* init.t$initial_probibality[c])
    #模型值
    init.va <- c(init.va,-sum(log(rowSums(init.cond))))
    init.par1[[i]] <- init.t
  }
  init.va 
  ii <- which(init.va==min(init.va))[1]
  
  init.par <- init.par1[[ii]]
  #init.par
  
  para.mu <- init.par$initial_mu_params
  para.cov <- init.par$initial_cov_params
  prob <-  init.par$initial_probibality
  npt <- init.par$np
  LL <- list()
  for(i in 1:length(times)){
    LL[[i]] <- Legendre.coef(times[[i]],normalize = T)
  }
  
  delta.ll <- 1
  iter <- 0
  LL.old <- 0
  
  while(delta.ll > 1e-5 && iter <= maxit ){
    #E step
    meanv <- get.mu.L(para.mu,times,k)
    sigma <- para.cov[1:nt];rho <- para.cov[-c(1:nt)]
    covv <- get.SAD1(sigma,rho,times)
    cond <- sapply(1:k, function(c) dmvnorm(y, meanv[c,], covv)* prob[c])
    if(any(cond==0)){
      cond[which(cond==0)] <- 1e-300
    }
    omega <- cond/rowSums(cond)
    
    #M step
    prob <- colSums(omega)/sum(omega)
    
    estv <- sapply(1:dim(omega)[2],function(c) colSums(y*matrix(rep(omega[,c],d),nrow=n,byrow=F)))/
      matrix(rep(colSums(omega),nts),nts,byrow=T)
    for(i in 1:nt){
      estv1 <- estv[(css[i]+1):css[i+1],]
      para.mu[[i]] <- sapply(1:k,function(iii,c,x,L){
        L.solve(c[iii],x[,iii],L)
      },c=npt[,i],x=estv1,L=LL[[i]])
    }
    prob1 <- matrix(rep(prob,n),n,byrow=T)
    
    #对得到参数进一步进行优化
    Q1 <- optim(par = para.cov, Q.M, 
                omega = omega, 
                v=estv,
                y = y,
                k = k,
                times=times,
                prob1=prob1,
                nt=nt,
                method = "Nelder-Mead",control=list(trace=F))
    para.cov <- Q1$par
    sigma <- para.cov[1:nt];rho <- para.cov[-c(1:nt)]
    covv <- get.SAD1(sigma,rho,times)
    cond <- sapply(1:k, function(c) dmvnorm(y, estv[,c], covv)* prob[c])
    if(any(cond==0)){
      cond[which(cond==0)] <- 1e-300
    }
    LL.new <- -sum(log(rowSums(cond)))
    #观察似然函数值之间的差异，判断是否收敛
    delta.ll <-  abs(LL.new - LL.old)
    LL.old <- LL.new
    
    iter <- iter + 1
    
    #cat("\n", "iter =", iter, "\n", "diff = ",delta.ll, "\n")
    
  }
  
  BIC <- 2*(LL.new) + log(n*nts)*(length(unlist(para.mu))+length(para.cov))
  
  
  clu.ii <- apply(omega,1,which.max)
  
  obj <- list(cluster_number = k,
              y <- y,
              BIC_value = BIC,
              cov_par = para.cov,
              mu_par = para.mu,
              probibality = prob,
              cluster = clu.ii)
  return(obj)
  
}

get.int <- function(y, k,times,max.order=7){
  
  n <- dim(y)[1]
  
  #set.seed((ss))
  #分为k个组
  init.cluster <- kmeans(y,centers = k)
  #统计每个组有多少个的概率，比如有两组时1：0.4，2：0.6
  prob <- table(init.cluster$cluster)/n 
  #得到均值
  mv1 <- sapply(1:k,function(x,y,ii){
    colMeans(y[which(ii==x),])
  },y=y,ii=init.cluster$cluster)
  
  #ts=c(0,9)
  ts <- cumsum(c(0,unlist(lapply(times,length))))
  #内容为0，k行ts-1列
  np.m <- matrix(0,k,length(ts)-1)
  muall <- list()
  for(i in 1:(length(ts)-1)){
    #取每组的数据
    mv1_s <- mv1[(ts[i]+1):ts[i+1],]
    L <- Legendre.coef(times[[i]],normalize = T)
    npt <- apply(mv1_s,2,get.order,max.order=max.order,L=L)
    np.m[,i] <- npt
    #进一步的筛选最佳阶数
    muall[[i]] <- sapply(1:k,function(i,c,x,L){
      L.solve(c[i],x[,i],L)
    },c=npt,x=mv1_s,L=L)
  }
  #初始化初始方差、系数
  init.sigma <- rep(1,length(times))
  if(length(init.sigma)==1){
    init.rho <- c(0.2)
  }else{
    nc <- combn(length(init.sigma),2)
    init.rho <- c(rep(0.2,length(init.sigma)),rep(0.5,dim(nc)[2]))
  }
  
  #迭代200次得到最优的参数
  init.cov <- SAD1.est(init.par=c(init.sigma,init.rho),mv=colMeans(y),ct=times,y=y,iter=200)
  return_obj <- list(initial_cov_params = init.cov$par,
                     initial_mu_params = muall,
                     initial_probibality = prob,
                     np=np.m)
  return(return_obj)
}

get.order <- function(x,max.order,L){
  
  xl <- length(x)
  llo <- sapply(2:max.order,L.solve,x=x,L=L)
  AIC1 <- c()
  for(i in 1:length(llo)){
    npar1 <- length(llo[[i]])
    if(npar1 >= max.order)
      break;
    #基于AIC准则和最小二乘法拟合不同阶数的线性模型
    RSS <- sum((x-L[,1:npar1]%*%as.matrix(llo[[i]]))^2)
    aic_cal <- 2*npar1+xl*log(RSS/xl)+2*npar1*(npar1+1)/(xl-npar1-1)
    AIC1 <- c(AIC1,aic_cal)
  }
  c(2:max.order)[which(AIC1== min(AIC1), arr.ind = TRUE)]
}

L.solve <- function(i,x,L){
  #调用包matlib，求解Ax=B形式的线性方程组，得到x
  gaussianElimination(L[,1:i], x, verbose = F, fractions = F)[1:i,i+1]
  
}

Legendre.coef <-function( t,normalize=F)
{
  
  if(normalize){
    u <- -1;
    v <- 1;
    tmin<-min(t);
    tmax<-max(t);
    #标准化的状态（归一化
    ti <- u + ((v-u)*(t-tmin))/(tmax - tmin);
  }else{
    ti <- t
  }
  
  ti <- as.matrix(ti)
  #勒让德多项式系数拟合
  L <- function(ti){
    tmp <- c(1,ti,0.5*(3*ti*ti-1),0.5*(5*ti^3-3*ti),0.125*(35*ti^4-30*ti^2+3),
             0.125*(63*ti^5-70*ti^3+15*ti),(1/16)*(231*ti^6-315*ti^4+105*ti^2-5),
             (1/16)*(429*ti^7-693*ti^5+315*ti^3-35*ti),(1/128)*(6435*ti^8-12012*ti^6+6930*ti^4-1260*ti^2+35),
             (1/128)*(12155*ti^9-25740*ti^7+18018*ti^5-4620*ti^3+315*ti),
             (1/256)*(46189*ti^10-109395*ti^8+90090*ti^6-30030*ti^4+3465*ti^2-63))
    as.matrix(tmp)
  } 
  
  lmm <- t(apply(ti,1, L))
  return(lmm);
}

SAD1.est <- function(init.par,mv,ct,y,iter=200){
  
  SAD1.mle <- function(spara,y,mv,ct){
    
    n <- length(ct)
    sigma <- spara[1:n];rho <- spara[-c(1:n)]
    #得到MM总协方差矩阵
    covm <- get.SAD1(sigma=sigma,rho=rho,ct=ct)
    LS <- -sum(dmvnorm(y,mv,covm,log=T))
    LS
  }
  
  #根据概率密度函数得到最优参数
  ls_f <- optim(init.par,SAD1.mle,mv=mv,y=y,ct=ct,method="Nelder-Mead",control = list(maxit=iter))
  #返回参数和模型值
  estp <- list(para=abs(ls_f$par),value=ls_f$value)
  estp 
}

get.SAD1 <- function(sigma,rho,ct){
  
  if(!is.list(ct)){
    stop("Parameters ct is not a list!")
  }
  
  n <- length(ct)
  ts <- unlist(lapply(ct,length))
  tsc <- c(0,cumsum(ts))
  nt <- length(unlist(ct))
  rho1 <- rho[c(1:n)]
  
  MM <- matrix(NA,nt,nt)
  
  for(i in 1:n){
    para <- c(sigma[i],rho1[i])
    rcii <- (tsc[i]+1):(tsc[i+1])
    
    #不同时间之间的协方差矩阵
    MM[rcii,rcii] <- get.SAD1_diag(para=para,times=ct[[i]])
  }
  
  if(n>1){
    nc <- combn(length(sigma),2)
    rho2 <- rho[-c(1:n)]
    for(i in 1:dim(nc)[2]){
      para <- c(sigma[nc[,i]],rho1[nc[,i]],rho2[i])
      ii1 <- (tsc[nc[1,i]]+1):(tsc[nc[1,i]+1])
      ii2 <- (tsc[nc[2,i]]+1):(tsc[nc[2,i]+1])
      #不同时间下的协方差矩阵
      sad_cor <- get.SAD1_corr(para=para,tim1=ct[[nc[1,i]]],tim2=ct[[nc[2,i]]])
      MM[ii1,ii2] <- sad_cor 
      MM[ii2,ii1] <- t(sad_cor) 
    }
  }
  
  MM
  
}

get.SAD1_diag <- function(para,times){
  
  nt1 <- length(times)
  nc1 <- combn(nt1,2)
  sig <- para[1]; phi <- para[2]
  
  MM1 <- matrix(NA,nt1,nt1)
  #保持恒定在时间点t的方差
  diag(MM1) <- (1-phi^(2*times))/(1-phi^2)
  #t2-t1
  ntc <- times[nc1[2,]]-times[nc1[1,]]
  #不同时间点之间的协方差
  MM1[lower.tri(MM1)] <- phi^ntc*sqrt((1-phi^(2*times[nc1[1,]]))/(1-phi^(2*times[nc1[2,]])))
  MM1t <- t(MM1)
  MM1[upper.tri(MM1)] <- MM1t[upper.tri(MM1t)]
  
  MM1*sig^2
  
}

get.SAD1_corr <- function(para,tim1,tim2){
  
  
  sig <- para[1]*para[2]; phi1 <- para[3];phi2 <- para[4];phi3 <- para[5];
  
  nt1 <- length(tim1)
  nt2 <- length(tim2)
  
  tc <- cbind(rep(1:nt1,each=nt2),rep(1:nt2,nt1))
  
  corf <- function(x,pi1,pi2,tim1,tim2){
    t1 <- tim1[x[1]];t2 <- tim2[x[2]];
    #不同组织与不同环境在不同时间点之间的协方差
    if(t2>=t1){
      tmp <- (pi2^(t2-t1)-pi1^t1*pi2^t2)/(1-pi1*pi2)*(sqrt(((1-pi1^(2*t1))*(1-pi2^(2*t2)))/((1-pi1^2)*(1-pi2^2))))^-1
    }else{
      tmp <- (pi1^(t1-t2)-pi1^t1*pi2^t2)/(1-pi1*pi2)*(sqrt(((1-pi1^(2*t1))*(1-pi2^(2*t2)))/((1-pi1^2)*(1-pi2^2))))^-1
    }
    tmp
  }
  
  MM2 <- matrix(apply(tc,1,corf,pi1=phi1,pi2=phi2,tim1=tim1,tim2=tim2),nrow=nt1,ncol=nt2,byrow=T)
  
  MM2*sig*phi3
  
}

Q.M <- function(para.cov,v,omega, y, k,times,prob1,nt){
  n <- dim(y)[1]
  sigma <- para.cov[1:nt];rho <- para.cov[-c(1:nt)]
  covv <- get.SAD1(sigma,rho,times)
  mvn.log <- sapply(1:k, function(c) dmvnorm(y, v[,c], covv, log = TRUE))+log(prob1)
  Q.LL <- -sum(omega * mvn.log)
  Q.LL
}

get.mu.L <- function(para.mu,times,k){
  
  muc <- c()
  for(i in 1:k) {
    fc <- c()
    for(j in 1:length(times)){
      mu1 <- para.mu[[j]]
      nt1 <- times[[j]]
      L <- Legendre.coef(times[[j]],normalize = T)
      
      if(is.matrix((mu1))){
        mu1.s <- mu1[,i]
      }
      if(is.list((mu1))){
        mu1.s <- mu1[[i]]
      }
      #计算预测值，根据多项式系数
      fc.t <- L[,1:length(mu1.s)]%*%as.matrix(mu1.s)
      fc <- c(fc,fc.t)
    }
    muc <- rbind(muc,fc)
  }
  muc
}

library(ggplot2)
library(mvtnorm)
library(plyr)
library(matlib)
library(readxl)

dat <- read.csv("D:/research/item2/论文写作/跑数据/filtered_datazg.csv", header = TRUE)
col_names <- dat[, 1]
dats <- dat[, -1]

# ct <- list(c(0.1,1.6,2.3,3.4,4.4,5.1),
#            c(1,3.9,4.6,5.7,6.8,7.4),
#            c(2,3.9,4.7,5.8,6.9,7.5),
#            c(3,4.8,5.9,7.0,7.6))
# ct <- list(c(1.1,1.6,2.3,3.4,4.4,5.1,1,3.9,4.6,5.7),
#            c(6.8,7.4,2,3.9,4.7,5.8,6.9,7.5,3,4.8,5.9,7.0,7.6)
#            )
ct=list(c(1,2,3,6,12,13,15,18,22),
        c(2,7,9,10,12,13,15),
        c(1,1.3,1.6,2,2.4,3.3,4.5)
       
        )
# ct=list(c(1,2,3,6,12,13),
#         c(2.1,7,9,10,12.1,13.1),
#         c(1.1,1.3,1.6,2.2,2.4,3.3),
#         c(4,7.1,8,10.1,11)
# )
# ct=list(c(1,2,3,6,12,13,14,15,16,17,18,19),
#         c(1.1,1.3,1.6,2.2,2.4,3.3,4,7.1,8,10.1,11)
# )
rr1 <- cluster_3D_BIC(y=dats,times=ct)

rr2 <- capture.output(rr1)
write.table(rr2, file = "cluster_data27102.txt", sep = "\t", row.names = FALSE)