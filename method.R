library("GrassmannOptim")
library("sSDR")

Simulation = function(truebasis,method,Data_generator,loop){
  setwd('.')
  C = 2
  if(Data_generator %in% c('MIM-3','LDA-Ind','LDA-Cor')) {C =3}
  d = 1
  if(Data_generator %in% c('MIM-1','MIM-2','MIM-3','LDA-Ind','LDA-Cor')) {d =2}
  for (i in 1:loop) {
    Data=generate_data(truebasis,Data_generator)
    Y = Data[,1]
    X = Data[,-1]
    for (usemethod in method) {
      print(usemethod)
      beta_hat = estimate(usemethod,X,Y,C,d)
      print(beta_hat)
      err = Cal_err(truebasis,beta_hat,d)
      write(paste(Data_generator,toString(err),sep = ","),paste("DATA//",usemethod,"data.csv",sep = "_"),append = TRUE,sep=",")
    }
  }
}
normalize = function(m,d){
  m = matrix(m,ncol = d)
  for (i in 1:d) {
    m[,i] = m[,i]/sqrt(sum((m[,i]^2)))
  }
  m
}

Cal_err = function(truebasis, beta_hat, d){
  err =NULL
  truebasis = normalize(truebasis[,1:d],d)
  beta_hat = normalize(beta_hat[,1:d],d)
  if(d == 1){
    err =acos(abs(t(truebasis)%*%beta_hat))*2/pi*90
  }else{
    Ptrue= truebasis[,1:d]%*% solve(t(truebasis[,1:d])%*% truebasis[,1:d])%*% t(truebasis[,1:d])
    Phat= beta_hat%*% solve(t(beta_hat)%*% beta_hat)%*% t(beta_hat)
    err = norm(Ptrue-Phat,"F")
  }
  return(err)
}

estimate = function(method, X, Y, C, d){
  if (method == "sir"){
    return(SIR(X,Y,C,d))
  }else if (method == "save"){
    return(SAVE(X,Y,C,d))
  }else if (method == "dr"){
    return(DR(X,Y,C,d))
  }else if (method == "mases"){
    return(MASES(X,Y,C,d))
  }else if (method == "energy"){
    return(Energy(X,Y,d))
  }
}

matpower = function(a,alpha){
  a = (a + t(a))/2
  tmp = eigen(a)
  return(tmp$vectors%*%diag((tmp$values)^alpha)%*%t(tmp$vectors))
}

DR = function(x,y,h,r,ytype="categorical"){
  p=ncol(x);n=nrow(x)
  signrt=matpower(var(x),-1/2)
  xc=t(t(x)-apply(x,2,mean))
  xst=xc%*%signrt
  if(ytype=="continuous") ydis=discretize(y,h)
  if(ytype=="categorical") ydis=y
  yless=ydis;ylabel=numeric()
  for(i in 1:n) {if(var(yless)!=0) {ylabel=c(ylabel,yless[1]);yless=yless[yless!=yless[1]]}}
  ylabel=c(ylabel,yless[1])
  prob=numeric() 
  for(i in 1:h) prob=c(prob,length(ydis[ydis==ylabel[i]])/n)
  vxy = array(0,c(p,p,h));exy=numeric()
  for(i in 1:h) {
    vxy[,,i]=var(xst[ydis==ylabel[i],])
    exy=rbind(exy,apply(xst[ydis==ylabel[i],],2,mean))}
  mat1 = matrix(0,p,p);mat2 = matrix(0,p,p)
  for(i in 1:h){
    mat1 = mat1+prob[i]*(vxy[,,i]+exy[i,]%*%t(exy[i,]))%*%
      (vxy[,,i]+exy[i,]%*%t(exy[i,]))
    mat2 = mat2+prob[i]*exy[i,]%*%t(exy[i,])}
  out = 2*mat1+2*mat2%*%mat2+2*sum(diag(mat2))*mat2-2*diag(p)
  return(signrt%*%eigen(out)$vectors[,1:r])
}

SAVE=function(x,y,h,r,ytype="categorical"){
  p=ncol(x)
  n=nrow(x)
  signrt=matpower(var(x),-1/2)
  xc=t(t(x)-apply(x,2,mean))
  xst=xc%*%signrt
  if(ytype=="continuous") ydis=discretize(y,h)
  if(ytype=="categorical") ydis=y
  yless=ydis
  ylabel=numeric()
  for(i in 1:n) {
    if(var(yless)!=0) {
      ylabel=c(ylabel,yless[1])
      yless=yless[yless!=yless[1]]
    }
  }
  ylabel=c(ylabel,yless[1])
  prob=numeric() 
  for(i in 1:h) prob=c(prob,length(ydis[ydis==ylabel[i]])/n)
  vxy = array(0,c(p,p,h))
  for(i in 1:h) vxy[,,i] = var(xst[ydis==ylabel[i],]) 
  savemat=0
  for(i in 1:h){
    savemat=savemat+prob[i]*(vxy[,,i]-diag(p))%*%(vxy[,,i]-diag(p))
  }
  return(signrt%*%eigen(savemat)$vectors[,1:r])
}

SIR=function(x,y,h,r,ytype="categorical"){
  p=ncol(x)
  n=nrow(x)
  signrt=matpower(var(x),-1/2)
  xc=t(t(x)-apply(x,2,mean))
  xst=xc%*%signrt
  if (ytype=="continuous") ydis=discretize(y,h)
  if(ytype=="categorical") ydis=y
  yless=ydis;ylabel=numeric()
  for(i in 1:n) {
    if(var(yless)!=0) {
      ylabel=c(ylabel,yless[1])
      yless=yless[yless!=yless[1]]
    }
  }
  ylabel=c(ylabel,yless[1])
  prob=numeric()
  exy=numeric()
  for(i in 1:h) 
    prob=c(prob,length(ydis[ydis==ylabel[i]])/n) 
  for(i in 1:h) exy=rbind(exy,apply(xst[ydis==ylabel[i],],2,mean))
  sirmat=t(exy)%*%diag(prob)%*%exy
  return(signrt%*%eigen(sirmat)$vectors[,1:r])
}

MASES = function(x,y,C,d){
  p = dim(x)[2]
  signrt=matpower(var(x),-1/2)
  xc=t(t(x)-apply(x,2,mean))
  z=xc%*%signrt
  Data = cbind(y,z)
  W = list(Qt = generate_basis(p), dim = c(d, p), sigmas = Data)
  init = initial(W,100)
  W = list(Qt = init$maxbasis, dim = c(init$d_hat, p), sigmas = Data)
  eps_err = 1e-3
  result = GrassmannOptim(objfun = objfun, W, eps_conv = eps_err,verbose = TRUE,eps_f = 1e-4)
  trans = signrt%*%result$Qt[,1:d]
  return(trans)
}

initial = function(W, n= 10){
  maxbasis = W$Qt
  d_hat = W$dim[1]
  maxv = f_gra(W, calculate_gradient = FALSE)$value
  if(W$dim[1] == 1){
    for(i in 1:n){
      print(i)
      currentbasis = generate_basis(p)
      W1 = list(Qt =currentbasis, dim = W$dim, sigmas = W$sigmas)
      currentvalue =f_gra(W1,calculate_gradient = FALSE)$value
      if(currentvalue<=maxv) {
        maxbasis=currentbasis
        maxv = currentvalue
      }
    }
  } else {
    seqmeth=sequential_method(W)
    maxbasis = seqmeth$B
    d_hat = seqmeth$d_hat
  }
  return(list(maxbasis=maxbasis,d_hat = d_hat))
}

sequential_method = function(W){
  p = W$dim[2]
  B = diag(rep(1,p))
  B_orthogal = diag(rep(1,p))
  value = rep(0,p)
  for (i in 1:(p-1)) {
    for (j in 1:2) {
      init_point = scale(rnorm(p+1-i))
      opt = optim(init_point,f_init,gr_init,bases=B_orthogal,W =W,method ="BFGS",
                    control = list(maxit = 20000, temp = 20,REPORT = 500, trace = TRUE))
      b1 = opt$par
      v = 1 - opt$value
      if(v>value[i]){
        value[i]=v
        b_hat = b1
      }
    }
    B[,i] = B_orthogal%*%b_hat
    B = orthnormal(B)
    B_orthogal = matrix(B[,(i+1):p],nrow = p)
  }
  d_hat = choose_d(value)
  B[,p] = B_orthogal
  return(list(B=B, d_hat =d_hat))
}

choose_d = function(v){
  a = c()
  for (i in 1:(length(v)-1)) {
    a = c(a, v[i+1]/v[i])
  }
  return(order(a)[1])
}

f_init = function(b,bases,W){
  d = dim(bases)[2]
  p = W$dim[2]
  dataset = W$sigmas
  n = dim(dataset)[1]
  U = bases%*%b
  X = dataset[,-1]
  Y = dataset[,1]
  Cata = unique(Y)
  C = length(Cata)
  yp = matrix(0,nrow = 1,ncol = C)
  for(i in Cata){
    yp[i] = sum(Y==i)/n
  }
  sigma_hat = 1
  epn =0
  hn = 1.06 * sigma_hat * n^(-0.2)
  Uij = matrix(0,ncol = n,nrow = n)
  for (i in 1:n){
    for (j in 1:n) {
      Uij[i,j]=exp(-(2*hn^(2))^(-1)*(norm(t(U)%*%matrix(X[i,]-X[j,],ncol = 1),type = "2")^2))
    }
  }
  f_hat = matrix(0,ncol = n,nrow = C)
  for (i in 1:C) {
    for (j in 1:n) {
      f_hat[i,j] = (2*pi)^(-1/2)*(sum(Uij[j,Y==i])-Uij[j,j]*(Y[j]==i))/(n-1)/hn
    }
  }
  
  value =sum(apply(f_hat, 2, function(x){allsqrt(x)/(C*(C-1)/2)/(yp%*%x+epn)}))
  return(value/n)
}

gr_init = function(b,bases,W){
  d = dim(bases)[2]
  p = W$dim[2]
  dataset = W$sigmas
  n = dim(dataset)[1]
  U = bases%*%b
  X = dataset[,-1]
  Y = dataset[,1]
  Cata = unique(Y)
  C = length(Cata)
  yp = matrix(0,nrow = 1,ncol = C)
  for(i in Cata){
    yp[i] = sum(Y==i)/n
  }
  sigma_hat = 1
  epn =0
  hn = 1.06 * sigma_hat * n^(-0.2)
  Uij = matrix(0,ncol = n,nrow = n)
  for (i in 1:n){
    for (j in 1:n) {
      Uij[i,j]=exp(-(2*hn^(2))^(-1)*(norm(t(U)%*%matrix(X[i,]-X[j,],ncol = 1),type = "2")^2))
    }
  }
  f_hat = matrix(0,ncol = n,nrow = C)
  for (i in 1:C) {
    for (j in 1:n) {
      f_hat[i,j] = (2*pi)^(-1/2)*(sum(Uij[j,Y==i])-Uij[j,j]*(Y[j]==i))/(n-1)/hn
    }
  }
   
  Df = matrix(0,ncol = 1,nrow = d)
  for (i in 1:n) {
    for (j in 1:C) {
      p1 = (sum(sqrt(f_hat[-j,i]))*(yp%*%f_hat[,i]+epn)/2/sqrt(f_hat[j,i])-allsqrt(f_hat[,i])*yp[j])/(yp%*%f_hat[,i]+epn)^2
      p2 = -1*(2*pi)^(-1/2)/(n-1)/hn^(3)
      p3 = matrix(0,ncol = 1,nrow = d)
      for (m in 1:n) {
        if(Y[m]==j&&m!=i){
          p3 = p3 + Uij[m,i]*t(bases)%*%(X[m,]-X[i,])%*%t(X[m,]-X[i,])%*%U
        }
      }
      Df = Df + p1[1,1]*p2*p3/(C*(C-1)/2)
      }
    }
  return(Df/n)
}



f_gra = function(W,calculate_gradient = TRUE){
  Qt = W$Qt
  d = W$dim[1]
  p = W$dim[2]
  dataset = W$sigmas
  n = dim(dataset)[1]
  U = matrix(Qt[,1:d],ncol = d)
  V = matrix(Qt[,(d+1):p], ncol = p-d)
  X = dataset[,-1]
  Y = dataset[,1]
  Cata = unique(Y)
  C = length(Cata)
  yp = matrix(0,nrow = 1,ncol = C)
  for(i in Cata){
    yp[i] = sum(Y==i)/n
  }
  sigma_hat = d
  epn =0
  hn = 1.06 * sigma_hat * n^(-0.2)
  Uij = matrix(0,ncol = n,nrow = n)
  for (i in 1:n){
    for (j in 1:n) {
      Uij[i,j]=exp(-(2*hn^(2*d))^(-1)*(norm(t(U)%*%matrix(X[i,]-X[j,],ncol = 1),type = "2")^2))
    }
  }
  f_hat = matrix(0,ncol = n,nrow = C)
  for (i in 1:C) {
    for (j in 1:n) {
      f_hat[i,j] = (2*pi)^(-d/2)*(sum(Uij[j,Y==i])-Uij[j,j]*(Y[j]==i))/(n-1)/hn^d
    }
  }
  
  value =sum(apply(f_hat, 2, function(x){allsqrt(x)/(C*(C-1)/2)/(yp%*%x+epn)}))
  
  Df = matrix(0,ncol = d,nrow = p)
  if(calculate_gradient){
    for (i in 1:n) {
      for (j in 1:C) {

          p1 = (sum(sqrt(f_hat[-j,i]))*(yp%*%f_hat[,i]+epn)/2/sqrt(f_hat[j,i])-allsqrt(f_hat[,i])*yp[j])/(yp%*%f_hat[,i]+epn)^2
          p2 = -1*(2*pi)^(-d/2)/(n-1)/hn^(3*d)
          p3 = matrix(0,ncol = d,nrow = p)
          for (m in 1:n) {
            if(Y[m]==j&&m!=i){
              p3 = p3 + Uij[m,i]*(X[m,]-X[i,])%*%t(X[m,]-X[i,])%*%U
            }
          }
          Df = Df + p1[1,1]*p2*p3/(C*(C-1)/2)
        
      }
    }
  }
  return(list(value = value/n,Gradient =t(Df)%*%V/n))
}

allsqrt = function(x){
  sqrtsum = 0 
  for (i in 1:(length(x)-1)) {
    for (j in (i+1):length(x)) {
      sqrtsum = sqrtsum + sqrt(x[i]*x[j])  
    }
  }
  sqrtsum
}

objfun = function(W){
  Value = f_gra(W)
  list(value =1- Value$value, gradient = -Value$Gradient)
}


energy_distance = function(W){
  Qt = W$Qt
  d = W$dim[1]
  p = W$dim[2]
  dataset = W$sigmas
  n = dim(dataset)[1]
  U = matrix(Qt[,1:d],ncol = d)
  V = matrix(Qt[,(d+1):p], ncol = p-d)
  
  Y =dataset[,1]
  X = dataset[,-1]
  energy = 0
  gr = matrix(0,nrow = p, ncol = d)
  Cata = unique(Y)
  C = length(Cata)
  for (i in 1:(C-1)){
    for (j in (i+1):C) {
      func = cal_energy(U, X[Y == Cata[i],],X[Y == Cata[j],],d)
      energy = energy + func$value
      gr = gr+func$Gradient
    }
  }
  return(list(value =energy,Gradient = t(gr)%*%V))
}


cal_energy = function(U,X,Y,d){
  n1 = dim(X)[1]
  n2 = dim(Y)[1]
  s = 0
  gr =matrix(0,nrow = dim(X)[2], ncol = d)
  for (i in 1:n1) {
    for (j in 1:n2) {
      h = sqrt(sum((t(U)%*%(X[i,]-Y[j,]))^2))
      s = s + h/n1/n2*2
      gr =  gr + 2/h/n1/n2*(X[i,]-Y[j,])%*%t(X[i,]-Y[j,])%*%U
    }
  }
  for (i in 1:n1) {
    for (j in 1:n1) {
      if(i!=j){
        h = sqrt(sum((t(U)%*%(X[i,]-X[j,]))^2))
        s = s - h/n1/n1
        gr =  gr - 1/h/n1/n1*(X[i,]-X[j,])%*%t(X[i,]-X[j,])%*%U
      }
      
    }
  }
  for (i in 1:n2) {
    for (j in 1:n2) {
      if(i!=j){
        h = sqrt(sum((t(U)%*%(Y[i,]-Y[j,]))^2))
        s = s - h/n2/n2
        gr =  gr - 1/h/n2/n2*(Y[i,]-Y[j,])%*%t(Y[i,]-Y[j,])%*%U
      }
      
    }
  }
  return(list(value = s*n1*n2/(n1+n2),Gradient =gr*n1*n2/(n1+n2)))
}

objfun_energy = function(W){
  Value = energy_distance(W)
  list(value = Value$value, gradient = Value$Gradient)
}


Energy = function(x,y,d){
  p = dim(x)[2]
  signrt=matpower(var(x),-1/2)
  xc=t(t(x)-apply(x,2,mean))
  z=xc%*%signrt
  Data = cbind(y,z)
  W = list(Qt = generate_basis(p), dim = c(d, p), sigmas = Data)
  W = list(Qt = inital_Energy(W,50), dim = c(d, p), sigmas = Data)
  eps_err = 1e-3
  sim_annel = NULL
  if (d == 1){
    sim_annel = FALSE
  }else{
    sim_annel = FALSE
  }
  result = GrassmannOptim(objfun = objfun_energy, W, sim_anneal = sim_annel, temp_init = 20,
                          cooling_rate = 4, eps_conv = eps_err,verbose = TRUE,eps_f = 1e-4)
  trans = signrt%*%result$Qt[,1:d]
  return(trans)
}

inital_Energy = function(W, n=100){
  maxbasis = W$Qt
  maxv = energy_distance(W)$value
  for (i in 1:n) {
    cbasis = generate_basis(p)
    W1 = list(Qt =cbasis, dim = c(W$dim[1], W$dim[2]), sigmas = W$sigma)
    cv = energy_distance(W1)$value
    if(cv>maxv){
      maxbasis=cbasis
      maxv = cv 
    }
  }
  return(maxbasis)
}
