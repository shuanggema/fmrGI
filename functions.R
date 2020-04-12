fnorm=function(sigma,gamma){
  f=(1/(sqrt(2*pi)*sigma))*exp(-1/(2*sigma^2)*gamma*gamma)
}

e=function(sigma,delta,x){
  
  
  
  sum(delta[,k]%*%x[,j]^2)/(n*sigma[k]^2)
}

S <- function(sigma,delta,gamma,x,y,eta,alpha) {
  
  s=1/n*sum(delta*(sigma*y-gamma)*x)+eta*alpha
  
}

#u <- function(alpha,tau,beta,c) {

# u=1/(alpha+tau)*sum(c[j,which(beta[,k]!=0)])

#}
###main###

sim <- function(x,z,prob,alpha,beta,ssd){
  ## Purpose: generates data from FMR model
  ## ----------------------------------------------------------------------
  ## Arguments: x (design-matrix);prob (component probabilities)
  ##            beta ((p+1)*k Matix of regression coef);ssd (standard dev.)
  ## ----------------------------------------------------------------------
  k <- length(prob)
  n <- dim(x)[1]
  y <- numeric(n)
  s <- sample(1:k,n,replace=TRUE,p=prob)
  for (i in 1:n){
    y[i] <- rnorm(1,mean=alpha[,s[i]]%*%x[i,]+beta[,s[i]]%*%z[i,],sd=ssd[s[i]])
  }
  # for (i in 1:n){
  #   y[i] <- rnorm(1,mean=alpha[,s[i]]%*%x[i,],sd=ssd[s[i]])
  # }
  output<-list(y=y,s=s)
}


frm_corr<-function(x,z,y,K,lambda1,lambda2,class_old,ksi=6,tau=1e-4,maxiter=200,iter_v=1e-4){
  
  n<-dim(x)[1]
  p<-dim(x)[2]
  q<-dim(z)[2]
  
  
  delta<-matrix(0,n,K)
  for (i in 1:n) {
    delta[i,class_old[i]]=1
    
  }
  
  
  

  miu.old <- rep(1/K,K)
  for (k in 1:K) {
    miu.old[k]=sum(delta[,k]/n)
  }
  miu.new <- miu.old
  
  alpha.old<-matrix(0,p,K)
  beta.old<-matrix(0,q,K)
  
  for (k in 1:K){
    fit1<-cv.glmnet(cbind(x[class_old==k,],z[class_old==k,]),y[class_old==k])
    temp<-fit1$glmnet.fit$beta[,fit1$lambda==fit1$lambda.min]
    alpha.old[,k]<-as.matrix(temp[1:p],p,1)
    beta.old[,k]<-as.matrix(temp[(p+1):(p+q)],q,1)
  }
  
  
 
  c=cor(x,z,method = 'pearson')
  c=abs(c)*(abs(c)>0.15)
  
  # c=matrix(0,p,q)
  # c[1:pp,1:pp]=1
  
  sigma<-matrix(0,K,1)
  for (k in 1:K) {
    sigma[k]= sqrt(sum(delta[,k]*(y-x%*%alpha.old[,k]-z%*%beta.old[,k])^2)/sum(delta[,k]))
  }
  
  
  sigma <-1/sigma
  for (k in 1:K) {
    alpha.old[,k]<-alpha.old[,k]*sigma[k]
    beta.old [,k]<- beta.old[,k]*sigma[k]
  }
  
  gamma<-matrix(0,n,K)
  for (k in 1:K){
    gamma[,k]=x%*%alpha.old[,k]+z%*%beta.old[,k]
  }
  f<-matrix(0,n,2)
  
  alpha<-alpha.old
  beta<-beta.old
  diff_v=1
  t<-1 ################################
  while ((diff_v>iter_v) && (t<maxiter)){
    
    alpha.old<-alpha
    beta.old<-beta
    miu.old<-miu.new 
    
    
    # for (k in 1:K) {
    #   f[,k]=fnorm(sigma = sigma[k],gamma[,k])
    # }
    # g=rowSums(f*(matrix(1,n,1)%*%miu.old))
    # delta1<-matrix(0,n,2)
    # 
    # 
    # for (k in 1:K) {
    #   delta1[,k]=miu.old[k]*f[,k]/g
    # }
    
    
    ss <- matrix(1,n,K)
    tt <- matrix(1,n,K)
    for(k in 1:K){
      ss[,k]<-rep(0,n)
      for(kk in 1:K){
        tt[,kk] <-exp(-1/2*((sigma[kk]*y-gamma[,kk])^2-(sigma[k]*y-gamma[,k])^2)) 
        ss[,k] <- ss[,k]+miu.old[kk]*(sigma[kk]/sigma[k])*tt[,kk]
      }
      if(miu.old[k]==0) {delta[,k]=0} else{
        delta[,k] <- miu.old[k]/ss[,k]
      }
    }
    a<-rep(1,2)
    b<-rep(1,2)
    c2<-rep(1,2)
    
    for (k in 1:K) {
      miu.new[k]=sum(delta[,k])/n
      a[k]=sum(delta[,k]*y^2)
      b[k]=sum(delta[,k]*gamma[,k]*y)
      c2[k]=sum(delta[,k])
      sigma[k]= (b[k]+sqrt(b[k]^2+4*a[k]*c2[k]))/(2*a[k])
    }
    
    
    
    #print(miu.new)
    
    eta <- matrix(0,p,2)
    for (k in 1:K) {
      u<-matrix(0,p,K)
      for (j in 1:p) {
        eta[j,k]= sum(delta[,k]*x[,j]^2)/n
        s=S(sigma = sigma[k],delta[,k],gamma = gamma[,k],x[,j],y=y,eta[j,k],alpha[j,k])
        
        
        u[j,k]=2/(tau)*exp(-alpha[j,k]^2/tau)*sum(c[j,which(beta[,k]!=0)])
        
        
        if(abs(s)>lambda1*ksi*(eta[j,k]+lambda2*u[j,k])){
          alpha[j,k]=(s)/(eta[j,k]+lambda2*u[j,k])
        } else {
          if (abs(s)>lambda1){
            alpha[j,k]=(s-sign(s)*lambda1)/(eta[j,k]+lambda2*u[j,k]-1/ksi)
          } else {
            alpha[j,k]=0
          }
        }
        
        alpha[abs(alpha)<1e-3]=0
        gamma[,k]=gamma[,k]+x[,j]*alpha[j,k]-x[,j]*alpha.old[j,k]
      }
      
      u<-matrix(0,q,K)
      for (l in 1:q) {
        eta[l,k]= sum(delta[,k]*z[,l]^2)/n
        s=S(sigma = sigma[k],delta[,k],gamma = gamma[,k],z[,l],y=y,eta[l,k],beta[l,k])
        
      
        
        u[l,k]=2/(tau)*exp(-beta[l,k]^2/tau)*sum(c[which(alpha[,k]!=0),l])
        
        if(abs(s)>lambda1*ksi*(eta[l,k]+lambda2*u[l,k])){
          beta[l,k]=(s)/(eta[l,k]+lambda2*u[l,k])
        } else {
          if(abs(s)>lambda1){
            beta[l,k]=(s-sign(s)*lambda1)/(eta[l,k]+lambda2*u[l,k]-1/ksi)
          } else {
            beta[l,k]=0
          }
        }
        
        beta[abs(beta)<1e-3]=0
        
        gamma[,k]=gamma[,k]+z[,l]*beta[l,k]-z[,l]*beta.old[l,k]
      }
    }
    
    diff_v=max(max(abs(miu.new-miu.old)),max(abs(alpha-alpha.old)),max(abs(beta-beta.old)))
    
    #print(max(abs(alpha-alpha.old)))
    
    t<-t+1
  }
  for (k in 1:K) {
    alpha[,k]<-alpha[,k]/sigma[k]
    beta[,k]<- beta[,k]/sigma[k]
  }
  
  alpha[abs(alpha)<1e-3]=0
  beta[abs(beta)<1e-3]=0
  
  sigma=1/sigma
  output<-list(alpha=alpha,beta=beta,sigma=sigma,delta=delta,mu=miu.new)
}


frm_corr_lasso<-function(x,z,y,K,lambda1,lambda2,class_old,tau=1e-4,maxiter=200,iter_v=1e-4){
  
  n<-dim(x)[1]
  p<-dim(x)[2]
  q<-dim(z)[2]
  
  
  delta<-matrix(0,n,K)
  for (i in 1:n) {
    delta[i,class_old[i]]=1
    #delta[i,output$s[i]]=1
  }
  
  
  
  #class_old<-kmresults$cluster
  #class_old<-output$s
  miu.old <- rep(1/K,K)
  for (k in 1:K) {
    miu.old[k]=sum(delta[,k]/n)
  }
  miu.new <- miu.old
  
  alpha.old<-matrix(0,p,K)
  beta.old<-matrix(0,q,K)
  
  for (k in 1:K){
    fit1<-cv.glmnet(cbind(x[class_old==k,],z[class_old==k,]),y[class_old==k])
    temp<-fit1$glmnet.fit$beta[,fit1$lambda==fit1$lambda.min]
    alpha.old[,k]<-as.matrix(temp[1:p],p,1)
    beta.old[,k]<-as.matrix(temp[(p+1):(p+q)],q,1)
  }
  
  
  
  c=cor(x,z,method = 'pearson')
  c=abs(c)*(abs(c)>0.2)
  
  # c=matrix(0,p,q)
  # c[1:pp,1:pp]=1
  
  sigma<-matrix(0,K,1)
  for (k in 1:K) {
    sigma[k]= sqrt(sum(delta[,k]*(y-x%*%alpha.old[,k]-z%*%beta.old[,k])^2)/sum(delta[,k]))
  }
  
  
  sigma <-1/sigma
  for (k in 1:K) {
    alpha.old[,k]<-alpha.old[,k]*sigma[k]
    beta.old [,k]<- beta.old[,k]*sigma[k]
  }
  
  gamma<-matrix(0,n,K)
  for (k in 1:K){
    gamma[,k]=x%*%alpha.old[,k]+z%*%beta.old[,k]
  }
  f<-matrix(0,n,2)
  
  alpha<-alpha.old
  beta<-beta.old
  diff_v=1
  t<-1 ################################
  while ((diff_v>iter_v) && (t<maxiter)){
    
    alpha.old<-alpha
    beta.old<-beta
    miu.old<-miu.new 
    
    
    # for (k in 1:K) {
    #   f[,k]=fnorm(sigma = sigma[k],gamma[,k])
    # }
    # g=rowSums(f*(matrix(1,n,1)%*%miu.old))
    # delta1<-matrix(0,n,2)
    # 
    # 
    # for (k in 1:K) {
    #   delta1[,k]=miu.old[k]*f[,k]/g
    # }
    
    
    ss <- matrix(1,n,K)
    tt <- matrix(1,n,K)
    for(k in 1:K){
      ss[,k]<-rep(0,n)
      for(kk in 1:K){
        tt[,kk] <-exp(-1/2*((sigma[kk]*y-gamma[,kk])^2-(sigma[k]*y-gamma[,k])^2)) 
        ss[,k] <- ss[,k]+miu.old[kk]*(sigma[kk]/sigma[k])*tt[,kk]
      }
      if(miu.old[k]==0) {delta[,k]=0} else{
        delta[,k] <- miu.old[k]/ss[,k]
      }
    }
    a<-rep(1,2)
    b<-rep(1,2)
    c2<-rep(1,2)
    
    for (k in 1:K) {
      miu.new[k]=sum(delta[,k])/n
      a[k]=sum(delta[,k]*y^2)
      b[k]=sum(delta[,k]*gamma[,k]*y)
      c2[k]=sum(delta[,k])
      sigma[k]= (b[k]+sqrt(b[k]^2+4*a[k]*c2[k]))/(2*a[k])
    }
    
    
    
    #print(miu.new)
    
    eta <- matrix(0,p,2)
    for (k in 1:K) {
      u<-matrix(0,p,K)
      for (j in 1:p) {
        eta[j,k]= sum(delta[,k]*x[,j]^2)/n
        s=S(sigma = sigma[k],delta[,k],gamma = gamma[,k],x[,j],y=y,eta[j,k],alpha[j,k])
        
        
        # u[j,k]=1/(alpha[j,k]+tau)*sum(c[j,which(beta[,k]!=0)])
        # 
        # if(abs(s-lambda2*u[j,k])>lambda1*eta[j,k]*ksi){
        #   alpha[j,k]=(s-lambda2*u[j,k])/eta[j,k]
        # } else {
        #   if (abs(s-lambda2*u[j,k])>lambda1){
        #     alpha[j,k]=(s-lambda2*u[j,k]-sign(s-lambda2*u[j,k])*lambda1)/(eta[j,k]-1/ksi)
        #   } else {
        #     alpha[j,k]=0
        #   }
        # }
        u[j,k]=2/(tau)*exp(-alpha[j,k]^2/tau)*sum(c[j,which(beta[,k]!=0)]/(abs(beta[beta[,k]!=0,k])))/(abs(alpha[j,k])+1e-10)#1/(alpha[j,k]+tau)*sum(c[j,which(beta[,k]!=0)])
        
        #u[j,k]=2/(tau)*exp(-alpha[j,k]^2/tau)*sum(c[j,which(beta[,k]!=0)])#1/(alpha[j,k]+tau)*sum(c[j,which(beta[,k]!=0)])
        
        
        if (abs(s)>lambda1){
            alpha[j,k]=(s-sign(s)*lambda1)/(eta[j,k]+lambda2*u[j,k])
          } else {
            alpha[j,k]=0
          }
        
        alpha[abs(alpha)<1e-4]=0
        gamma[,k]=gamma[,k]+x[,j]*alpha[j,k]-x[,j]*alpha.old[j,k]
      }
      
      u<-matrix(0,q,K)
      for (l in 1:q) {
        eta[l,k]= sum(delta[,k]*z[,l]^2)/n
        s=S(sigma = sigma[k],delta[,k],gamma = gamma[,k],z[,l],y=y,eta[l,k],beta[l,k])
        
        # u[l,k]=1/(beta[l,k]+tau)*sum(c[which(alpha[,k]!=0),l])
        # 
        # if(abs(s-lambda2*u[l,k])>lambda1*eta[l,k]*ksi){
        #   beta[l,k]=(s-lambda2*u[l,k])/eta[l,k]
        # } else {
        #   if(abs(s-lambda2*u[l,k])>lambda1){
        #     beta[l,k]=(s-lambda2*u[l,k]-sign(s-lambda2*u[l,k])*lambda1)/(eta[l,k]-1/ksi)
        #   } else {
        #     beta[l,k]=0
        #   }
        # }
        
        u[l,k]=2/(tau)*exp(-beta[l,k]^2/tau)*sum(c[which(alpha[,k]!=0),l]/(abs(alpha[alpha[,k]!=0,k])))/(abs(beta[l,k])+1e-10)
        
        #u[l,k]=2/(tau)*exp(-beta[l,k]^2/tau)*sum(c[which(alpha[,k]!=0),l])
        
        if (abs(s)>lambda1){
            beta[l,k]=(s-sign(s)*lambda1)/(eta[l,k]+lambda2*u[l,k])
          } else {
            beta[l,k]=0
          }
        
        
        beta[abs(beta)<1e-4]=0
        
        gamma[,k]=gamma[,k]+z[,l]*beta[l,k]-z[,l]*beta.old[l,k]
      }
    }
    
    diff_v=max(max(abs(miu.new-miu.old)),max(abs(alpha-alpha.old)),max(abs(beta-beta.old)))
    
    #print(max(abs(alpha-alpha.old)))
    
    t<-t+1
  }
  for (k in 1:K) {
    alpha[,k]<-alpha[,k]/sigma[k]
    beta[,k]<- beta[,k]/sigma[k]
  }
  
  sigma=1/sigma
  output<-list(alpha=alpha,beta=beta,sigma=sigma,delta=delta,mu=miu.new)
}

matrix_vec = function(P, cov.str = "ar", rho = 0.3){
  tmp = matrix(rep(0, P * P), nrow = P)
  for (i in 1:P) 
    for (j in 1:P)
      if (i == j) tmp[i, j] = 1 else 
        if (cov.str == "ar") tmp[i, j] = (rho) ^ abs(i - j) else 
          if (cov.str == "flat") tmp[i, j] = rho
          return(tmp)
}
GetFPTP<-function(theta,theta_hat){
  # to get TNR (True Negative Rate ) and TPR (True Positive Rate) 
  thea = abs(theta) > 0   # transform coefficients to binary values
  thea_hat = abs(theta_hat) > 1e-8  # convert estimated coefficients to binary values
  A = sum((!thea)*(!thea_hat))  # A: TN
  B = sum((!thea)*thea_hat)   # B: FP
  C = sum(thea*(!thea_hat))   # C: FN
  D = sum(thea*thea_hat)    # D: TP
  TPR = D/(D+C)    # TPR=TP/(TP+FN)  true positive rate (TPR) sensitivity
  FPR = B/(B+A)    # FPR=FP/(TN+FP)  false positive rate     
  result=list(TPR= TPR, FPR = FPR,TP=D,FP=B)
  return(result)
}






mcp_corr<-function(x,z,y,K,lambda1,lambda2,ksi=6,tau=1e-4,maxiter=200,iter_v=1e-4){
  
  n<-dim(x)[1]
  p<-dim(x)[2]
  q<-dim(z)[2]
  
  
  c=cor(x,z,method = 'pearson')
  c=abs(c)*(abs(c)>0.15)
  
  # c=matrix(0,p,q)
  # c[1:pp,1:pp]=1
  
  
  
  alpha<-matrix(0,p,1)
  beta<-matrix(0,q,1)
  
  gamma=y-x%*%alpha-z%*%beta
  
  diff_v=1
  t<-1 ################################
  while ((diff_v>iter_v) && (t<maxiter)){
    
    alpha.old<-alpha
    beta.old<-beta
    
      eta <- matrix(0,p,1)
      u<-matrix(0,p,1)
      for (j in 1:p) {
        eta[j]= sum(x[,j]^2)/n
        s=sum(gamma*x[,j])/n+eta[j]*alpha[j]
        u[j]=2/(tau)*exp(-alpha[j]^2/tau)*sum(c[j,which(beta!=0)])
        
        
        if(abs(s)>lambda1*ksi*(eta[j]+lambda2*u[j])){
          alpha[j]=(s)/(eta[j]+lambda2*u[j])
        } else {
          if (abs(s)>lambda1){
            alpha[j]=(s-sign(s)*lambda1)/(eta[j]+lambda2*u[j]-1/ksi)
          } else {
            alpha[j]=0
          }
        }
        
        alpha[abs(alpha)<1e-3]=0
        gamma=gamma+x[,j]*alpha.old[j]-x[,j]*alpha[j]
      }
      
      u<-matrix(0,q,1)
      eta <- matrix(0,p,1)
      for (l in 1:q) {
        eta[l]= sum(z[,l]^2)/n
        s=sum(gamma*z[,l])/n+eta[l]*beta[l]
        u[l]=2/(tau)*exp(-beta[l]^2/tau)*sum(c[which(alpha!=0),l])
        
        if(abs(s)>lambda1*ksi*(eta[l]+lambda2*u[l])){
          beta[l]=(s)/(eta[l]+lambda2*u[l])
        } else {
          if(abs(s)>lambda1){
            beta[l]=(s-sign(s)*lambda1)/(eta[l]+lambda2*u[l]-1/ksi)
          } else {
            beta[l]=0
          }
        }
        beta[abs(beta)<1e-3]=0
        gamma=gamma+z[,l]*beta.old[l]-z[,l]*beta[l]
      }
    
    diff_v=max(max(abs(alpha-alpha.old)),max(abs(beta-beta.old)))
    
    
    t<-t+1
  }
  
  df1=sum(alpha!=0)
  df2=sum(beta!=0)
  df=df1+df2
  
  
  xx=cbind(x[,alpha!=0],z[,beta!=0])
  aa=solve(t(xx)%*%xx+0.001*diag(1,df))%*%t(xx)%*%y
  
  alpha[alpha!=0]=aa[1:df1]
  beta[beta!=0]=aa[(df1+1):df]
  
  RSS=sum((y-x%*%alpha-z%*%beta)^2)
  
  # if (df<n){
  #   fit=lm.fit(cbind(x[,alpha!=0],z[,beta!=0]),y)
  #   alpha[alpha!=0]=fit$coefficients[1:df1]
  #   beta[beta!=0]=fit$coefficients[(df1+1):df]
  #   RSS=sum((fit$residuals)^2)
  # }
  print(RSS)
  BIC=log(RSS)+df*2/n
  output<-list(alpha=alpha,beta=beta,BIC=BIC,RSS=RSS,df=df)
}


L<-function(y,x,z,alpha,beta,miu,mysigma){
  b=0
  a=0
  for (i in 1:n) {
    
    for (k in 1:K) {
      a=a+miu[k]*dnorm(y[i],mean = x[i,]%*%alpha[,k]+z[i,]%*%beta[,k],sd=mysigma[k])
      
    }
    b=b+log(a)
  }
  return(b/n)
}

# check if the group indicator
shat=function(output,s.hat){
  if(var(s.hat)!=0){
    if((table(output$s)[1]>table(output$s)[2])&(table(s.hat)[1]>table(s.hat)[2])|((table(output$s)[2]>table(output$s)[1])&(table(s.hat)[2]>table(s.hat)[1])))
    { s.hat=s.hat}else
    {s.hat[which(s.hat==1)]=0
    s.hat[which(s.hat==2)]=1
    s.hat[which(s.hat==0)]=2
    }
    return(s.hat)
  }
  else{return(s.hat)}
}

mse=function(ini,x,z,y,kmrandom){
  mse=0
  for (i in 1:n) {
    if(kmrandom[[ini]][[1]]$delta[i,1]>kmrandom[[ini]][[1]]$delta[i,2]){s.hat[i]=1}else{s.hat[i]=2}
    
    mse=mse+y[i]-kmrandom[[ini]][[1]]$alpha[,s.hat[i]]%*%x[i,]+kmrandom[[ini]][[1]]$beta[,s.hat[i]]%*%z[i,]
  }
  return(list(mse=mse,s.hat=s.hat))
}