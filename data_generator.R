library(MASS)

generate_basis = function(n){
  while (TRUE) {
    raw_matrix = matrix(rnorm(n^2), ncol=n)
    a = eigen(raw_matrix)
    if( all(a$values !=0) ){
      break
    }
  }
  eigen(raw_matrix %*% t(raw_matrix))$vector
}

generate_data = function(base, generator_type){
  if(generator_type == 'MDA-1') {
    return(MDA_1(base))
    } else if(generator_type == 'MDA-2') {
      return(MDA_2(base))
    } else if(generator_type == 'LDA-1-Ind') {
      return(LDA_1(base))
    } else if(generator_type == 'LDA-1-Cor') {
      return(LDA_1(base,sigma_type = "Cor"))
    } else if(generator_type == 'LDA-2-Ind') {
      return(LDA_1(base,sig_rate = 0.1))
    } else if(generator_type == 'LDA-2-Cor') {
      return(LDA_1(base,sigma_type = "Cor",sig_rate = 0.1))
      
      
    } else if(generator_type == 'SIM-1') {
      return(SIM_1(base))
    } else if(generator_type == 'SIM-2') {
      return(SIM_2(base))
    } else if(generator_type == 'LR-Ind') {
      return(LR(base))
    } else if(generator_type == 'LR-Cor') {
      return(LR(base,sigma_type = "Cor"))
      
      
    } else if(generator_type == 'MIM-1') {
      return(MIM_1(base))
    } else if(generator_type == 'MIM-2') {
      return(MIM_2(base))
    } else if(generator_type == 'MIM-3') {
      return(MIM_3(base))
    } else if(generator_type == 'LDA-Ind') {
      return(LDA(base))
    } else if(generator_type == 'LDA-Cor') {
      return(LDA(base,sigma_type = "Cor"))
    }
}

AR = function(n, rate){
  H = abs(outer(1:n,1:n,"-"))
  rate^H
}

logit = function(x){
  1/{1+exp(-x)}
}


MDA_1 = function(ortho_basis,d =1, N1 = 100, N2 =100){
  Data_YX = c()
  for (i in 1:N1){
    s = rbinom(1,1,0.5)
    betaX_y1 = rnorm(1,-1,sqrt(10))*s+rnorm(1,1,sqrt(0.1))*(1-s)
    B0X = mvrnorm(n = 1,mu= rep(0,dim(ortho_basis)[1]-d), Sigma = diag(rep(1,dim(ortho_basis)[1]-d)))
    X =   matrix(ortho_basis[,1],ncol = 1) %*% betaX_y1+  ortho_basis[,-1] %*% B0X
    Data_YX = rbind(Data_YX, c(1, t(X)))
  }
  for (i in 1:N2){
    s = rbinom(1,1,0.5)
    betaX_y1 = rnorm(1)*s+rnorm(1,2,1)*(1-s)
    B0X = mvrnorm(n = 1,mu= rep(0,dim(ortho_basis)[1]-d), Sigma = diag(rep(1,dim(ortho_basis)[1]-d)))
    X =   matrix(ortho_basis[,1],ncol = 1) %*% betaX_y1+  ortho_basis[,-1] %*% B0X
    Data_YX = rbind(Data_YX, c(2, t(X)))
  }
  Data_YX
}

MDA_2 = function(ortho_basis,d =1, N1 = 100, N2 =100){
  Data_YX = c()
  for (i in 1:N1){
    s = rbinom(1,1,0.5)
    betaX_y1 = rnorm(1,-2,sqrt(0.1))*s+rnorm(1,2,sqrt(0.1))*(1-s)
    B0X = mvrnorm(n = 1,mu= rep(0,dim(ortho_basis)[1]-d), Sigma = diag(rep(1,dim(ortho_basis)[1]-d)))
    X =   matrix(ortho_basis[,1],ncol = 1) %*% betaX_y1+  ortho_basis[,-1] %*% B0X
    Data_YX = rbind(Data_YX, c(1, t(X)))
  }
  for (i in 1:N2){
    s = rbinom(1,1,0.5)
    betaX_y1 = rnorm(1)*s+rnorm(1,5,1)*(1-s)
    B0X = mvrnorm(n = 1,mu= rep(0,dim(ortho_basis)[1]-d), Sigma = diag(rep(1,dim(ortho_basis)[1]-d)))
    X =   matrix(ortho_basis[,1],ncol = 1) %*% betaX_y1+  ortho_basis[,-1] %*% B0X
    Data_YX = rbind(Data_YX, c(2, t(X)))
  }
  Data_YX
}

LDA_1 = function(ortho_basis,d =1, N1 = 100, N2 =100,sigma_type = "Ind",sig_rate = 1){
  Data_YX = c()
  n = dim(ortho_basis)[1]
  sig = diag(rep(1,n))*sig_rate
  if(sigma_type == "Cor"){
    sig = AR(n,0.8)*sig_rate
  }
  for (i in 1:N1){
    X =   mvrnorm(n = 1,mu= rep(0,n), Sigma = sig)
    Data_YX = rbind(Data_YX, c(1, t(X)))
  }
  for (i in 1:N2){
    X =   mvrnorm(n = 1,mu= 2*sig %*% ortho_basis[,1:1]/sig_rate, Sigma = sig)
    Data_YX = rbind(Data_YX, c(2, t(X)))
  }
  Data_YX
}


SIM_1 = function(ortho_basis,d = 1,N = 200){
  Data_YX = c()
  n = dim(ortho_basis)[1]
  sig = AR(n,0.8)
  for(i in 1:N){
    X = mvrnorm(n = 1,mu= rep(0,n), Sigma = sig)
    value = 5*ortho_basis[,1] %*% X
    Data_YX = rbind(Data_YX, c(rbinom(1,1,logit(sin(value*pi/4)+0.1*value^2))+1, t(X)))
  }
  Data_YX
}


SIM_2 = function(ortho_basis, d =1, N =200){
  Data_YX = c()
  n = dim(ortho_basis)[1]
  sig = AR(n, 0.8)
  for(i in 1:N){
    X = mvrnorm(n = 1,mu= rep(0,n), Sigma = sig)
    X[3]= abs(X[1])+abs(X[2])+abs(X[1])*rnorm(1)
    X[4]=(X[1]+X[2])^2+abs(X[2])*rnorm(1)
    X[5]=rbinom(1,5,logit(X[2]))
    X[6]=rbinom(1,5,pnorm(X[2]))
    value = ortho_basis[,1] %*% X
    Data_YX = rbind(Data_YX, c(rbinom(1,1,logit(value+0.1*value^3))+1,t(X)))
  }
  Data_YX
}


LR = function(ortho_basis,d = 1,N = 200, sigma_type = "Ind"){
  Data_YX = c()
  n = dim(ortho_basis)[1]
  sig = diag(rep(1,n))
  if(sigma_type == "Cor"){
    sig = AR(n,0.8)
  }
  for(i in 1:N){
    X = mvrnorm(n = 1,mu= rep(0,n), Sigma = sig)
    value = ortho_basis[,1] %*% X
    Data_YX = rbind(Data_YX, c(rbinom(1,1,logit(2*value))+1, t(X)))
  }
  Data_YX
}


MIM_1 = function(ortho_basis,d=2,N=200){
  Data_YX = c()
  for (i in 1:N) {
    beta1X_y=runif(1,min=-4,max=4)
    beta2X_y=runif(1,min=-4,max=4)
    B0X = mvrnorm(n = 1,mu= rep(0,dim(ortho_basis)[1]-d), Sigma = diag(rep(1,dim(ortho_basis)[1]-d)))
    X =   matrix(ortho_basis[,1],ncol = 1) %*% beta1X_y+matrix(ortho_basis[,2],ncol = 1) %*%
      beta2X_y+  ortho_basis[,-c(1,2)] %*% B0X
    Y=2
    if (abs(beta1X_y)<=3 && abs(beta1X_y)>=1 && abs(beta2X_y)<=3 && abs(beta2X_y)>=1){
      Y = 1
    }
    Data_YX = rbind(Data_YX, c(Y, t(X)))
  }
  Data_YX
}

MIM_2 = function(ortho_basis,d=2,N=200){
  Data_YX = c()
  for (i in 1:N) {
    Y = rbinom(1,1,0.5)+1
    if (Y==1){
      betaX_y = mvrnorm(n=1, mu = rep(0,d),Sigma = diag(rep(1,2)))  
    } else {
      s = sample(1:4,1)
      betaX_y = mvrnorm(n=1, mu = c(0,3),Sigma = diag(rep(1,2)))*(s==1) +
        mvrnorm(n=1, mu = c(0,-3),Sigma = diag(rep(1,2)))*(s==2)+
        mvrnorm(n=1, mu = c(3,0),Sigma = 0.1*diag(rep(1,2)))*(s==3)+
        mvrnorm(n=1, mu = c(-3,0),Sigma = 0.1*diag(rep(1,2)))*(s==4)
    }
    B0X = mvrnorm(n = 1,mu= rep(0,dim(ortho_basis)[1]-d), Sigma = diag(rep(1,dim(ortho_basis)[1]-d)))
    X = ortho_basis[,1:2] %*% betaX_y + ortho_basis[,-(1:2)] %*% B0X
    Data_YX = rbind(Data_YX, c(Y, t(X)))
  }
  Data_YX
}


MIM_3 = function(ortho_basis,d=2,N=300){
  Data_YX = c()
  for (i in 1:N) {
    Y = sample(1:3,1)
    if (Y==1){
      betaX_y = mvrnorm(n=1, mu = rep(0,d),Sigma = diag(rep(1,2)))  
    } else if (Y==2){
      s = sample(0:1,1)
      betaX_y = mvrnorm(n=1, mu = c(3,0),Sigma = diag(rep(1,2)))*s +
        mvrnorm(n=1, mu = c(-3,0),Sigma = diag(rep(1,2)))*(1-s)
    } else {
      s = sample(0:1,1)
      betaX_y =mvrnorm(n=1, mu = c(0,3),Sigma = diag(c(5,1)))*s+
        mvrnorm(n=1, mu = c(0,-3),Sigma = diag(c(5,1)))*(1-s)
    }
    B0X = mvrnorm(n = 1,mu= rep(0,dim(ortho_basis)[1]-d), Sigma = diag(rep(1,dim(ortho_basis)[1]-d)))
    X = ortho_basis[,1:2] %*% betaX_y + ortho_basis[,-(1:2)] %*% B0X
    Data_YX = rbind(Data_YX, c(Y, t(X)))
  }
  Data_YX
}

LDA = function(ortho_basis,d=2,N=300,sigma_type = "Ind"){
  Data_YX = c()
  n = dim(ortho_basis)[1]
  sig = diag(rep(1,n))
  if(sigma_type == "Cor"){
    sig = AR(n,0.8)
  }
  for (i in 1:N) {
    Y = sample(1:3,1)
    if (Y==1){
      X = mvrnorm(n=1, mu = rep(0,n),Sigma = sig)  
    } else if (Y==2){
      X = mvrnorm(n=1, mu = 3*ortho_basis[,1],Sigma = sig)
    } else {
      X =mvrnorm(n=1, mu = 3*ortho_basis[,2],Sigma = sig)
    }
    Data_YX = rbind(Data_YX, c(Y, t(X)))
  }
  Data_YX
}
