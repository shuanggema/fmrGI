
#source("E:/code2017/R/2019FRM/main.R") 

rm(list=ls())
old<-Sys.time()
library('glmnet')
library('MASS')
library('ncvreg')

source("FMR_FUN190813.R") 


for (case in 1){
  
  repli<-100
  case1<-matrix(0,repli*4+4,13)
  colnames(case1)<-c('tp','alpha.mse','alpha.tp','alpha.fp','beta.mse','beta.tp','beta.fp','df','alpha.tpr','alpha.fpr','beta.tpr','beta.fpr','pmse')
  repli.out<-matrix(0,4,14)
  colnames(repli.out)<-c('method',colnames(case1))
  #model specifications
  
  n <- 300
  n_test<-100
  p <- 1000
  q <- 500
  init_time<-20
  
  
  
  if (case%%12==1){
    
    cov.str='ar'
    rho=0.7
    theta=matrix(0,200,200)
    for (i in 1:20){
      theta[((i-1)*10+1):(i*10),((i-1)*10+1):(i*10)]=1
    }
    
  } else if (case%%12==2){
    
    ### same as case 1 except for that rho=0.5
    
    cov.str='ar'
    rho=0.5
    theta=matrix(0,200,200)
    for (i in 1:20){
      theta[((i-1)*10+1):(i*10),((i-1)*10+1):(i*10)]=1
    }
    
    
  } else if (case%%12==3){
    
    
    cov.str='ar'
    rho=0.3
    theta=matrix(0,200,200)
    for (i in 1:20){
      theta[((i-1)*10+1):(i*10),((i-1)*10+1):(i*10)]=1
    }
    
    
  } else if (case%%12==4){
    ### same as case 1 except for that each block of theta is AR(0.7).
    cov.str='ar'
    rho=0.7
    theta=matrix(0,200,200)
    for (i in 1:20){
      theta[((i-1)*10+1):(i*10),((i-1)*10+1):(i*10)]=matrix_vec(10,cov.str = cov.str,rho = 0.7)
    }
    
  } else if (case%%12==5){
    
    cov.str='ar'
    rho=0.5
    theta=matrix(0,200,200)
    for (i in 1:20){
      theta[((i-1)*10+1):(i*10),((i-1)*10+1):(i*10)]=matrix_vec(10,cov.str = cov.str,rho = 0.7)
    }
    
  } else if (case%%12==6){
    
    cov.str='ar'
    rho=0.3
    theta=matrix(0,200,200)
    for (i in 1:20){
      theta[((i-1)*10+1):(i*10),((i-1)*10+1):(i*10)]=matrix_vec(10,cov.str = cov.str,rho = 0.7)
    }
    
  } else if (case%%12==7){
    
    ### same as case 1 except for that there are only ten correlated x-z blocks and the size of each block is 20. 
    cov.str='ar'
    rho=0.7
    theta=matrix(0,200,200)
    for (i in 1:10){
      theta[((i-1)*20+1):(i*20),((i-1)*20+1):(i*20)]=1
    }
    
  } else if (case%%12==8){
    
    cov.str='ar'
    rho=0.5
    theta=matrix(0,200,200)
    for (i in 1:10){
      theta[((i-1)*20+1):(i*20),((i-1)*20+1):(i*20)]=1
    }
    
  } else if (case%%12==9){
    cov.str='ar'
    rho=0.3
    theta=matrix(0,200,200)
    for (i in 1:10){
      theta[((i-1)*20+1):(i*20),((i-1)*20+1):(i*20)]=1
    }
    
  } else if (case%%12==10){
    ### same as case 7 except for that each block of theta is AR(0.7).
    cov.str='ar'
    rho=0.7
    theta=matrix(0,200,200)
    for (i in 1:10){
      theta[((i-1)*20+1):(i*20),((i-1)*20+1):(i*20)]=matrix_vec(20,cov.str = cov.str,rho = 0.7)
    }
    
    
  } else if (case%%12==11){
    
    cov.str='ar'
    rho=0.5
    theta=matrix(0,200,200)
    for (i in 1:10){
      theta[((i-1)*20+1):(i*20),((i-1)*20+1):(i*20)]=matrix_vec(20,cov.str = cov.str,rho = 0.7)
    }
    
  } else if (case%%12==0){
    
    cov.str='ar'
    rho=0.3
    theta=matrix(0,200,200)
    for (i in 1:10){
      theta[((i-1)*20+1):(i*20),((i-1)*20+1):(i*20)]=matrix_vec(20,cov.str = cov.str,rho = 0.7)
    }
  } 
  
  if (case<=12 || (case<=48 && case>36)){
    # two clusters have the same important markers but with different magnitudes
    cc<-11:210
    p0=10
    alpha.set <- cbind(c(rep(1.5,5),rep(0.5,5),rep(0,p-10)),c(rep(-0.5,10),rep(0,p-10)))
    beta.set<-cbind(c(rep(1.5,2),rep(0.5,8),rep(0,q-10)),c(rep(-0.5,10),rep(0,q-10)))
  } else if ((case<=24 && case>12) || (case<=60 && case>48)) {
    cc<-21:220
    p0=13
    # two clusters have the different important markers 
    alpha.set <- cbind(c(rep(1.5,5),rep(0.5,5),rep(0,p-10)),c(rep(0,3),rep(-0.5,10),rep(0,p-13)))
    beta.set<-cbind(c(rep(1.5,2),rep(0.5,8),rep(0,q-10)),c(rep(0,3),rep(-0.5,10),rep(0,q-13)))
  } else if ((case<=36 && case>24)|| (case<=72 && case>60)) {
    cc<-11:210
    p0=10
    # two clusters have some important markers with the same magnitudes
    alpha.set <- cbind(c(rep(1.5,5),rep(0.5,5),rep(0,p-10)),c(rep(-0.5,5),rep(0.5,2),rep(-0.5,3),rep(0,p-10)))
    beta.set<-cbind(c(rep(1.5,2),rep(0.5,8),rep(0,q-10)),c(rep(-0.5,5),rep(0.5,2),rep(-0.5,3),rep(0,q-10)))
  }
  
  if (case<=36){
    prob<-c(0.5,0.5)
  } else {
    prob<-c(0.4,0.6)
  }
  
  for (re in 1:repli) {
    print(paste0('The ',re,'th replicate'))
    set.seed(re)
    
    Sigma1=matrix(matrix_vec(p0,cov.str = cov.str,
                             rho = rho), ncol = p0)
    Sigma2 = matrix(matrix_vec(p-p0,
                               cov.str = cov.str,
                               rho = rho), ncol = p-p0)
    Sigma=matrix(0,p,p)
    Sigma[1:p0,1:p0]=Sigma1
    Sigma[(p0+1):p,(p0+1):p]=Sigma2
    
    x <-  matrix(mvrnorm(n, rep(0, p),Sigma),n,p)
    x_test<-matrix(mvrnorm(n_test, rep(0, p),Sigma),n_test,p)
    
    
    Sigma1=matrix(matrix_vec(p0,cov.str = cov.str,
                             rho = rho), ncol = p0)
    Sigma2 = matrix(matrix_vec(q-p0,
                               cov.str = cov.str,
                               rho = rho), ncol = q-p0)
    
    Sigma=matrix(0,q,q)
    Sigma[1:p0,1:p0]=Sigma1
    Sigma[(p0+1):q,(p0+1):q]=Sigma2
    
    z <-  matrix(mvrnorm(n, rep(0, q),Sigma),n,q)
    z_test <-  matrix(mvrnorm(n_test, rep(0, q),Sigma),n_test,q)
    
    z[,cc]=x[,cc]%*%theta
    z_test[,cc]=x_test[,cc]%*%theta
    ssd <- c(0.5,0.5)
    
    #generate data
    output <- sim(x=x,z=z,ssd=ssd,alpha=alpha.set,beta=beta.set,prob=prob)
    y<-output$y
    
    output_test <- sim(x=x_test,z=z_test,ssd=ssd,alpha=alpha.set,beta=beta.set,prob=prob)
    y_test<-output_test$y
    
    K <-length(prob)
    
    lambda1_set <- c(seq(from=0.05,to=0.1,by=0.01))#0.09
    lambda2_set <-c(seq(from=0.001,to=0.01,by=0.001))
    lll=1
    l1=length(lambda1_set)
    l2=length(lambda2_set)
    para=matrix(0,l1*l2,3)
    for (m1 in 1:l1){
      for (m2 in 1:l2){
        para[lll,1]=lll
        para[lll,2]=lambda1_set[m1]
        para[lll,3]=lambda2_set[m2]
        lll=lll+1
      }
    }
    maxiter <- 200#500
    ksi <- 6
    tau <-1e-2
    iter_v<-1e-3
    
    
    method_result<-list()
    init_set<-list()
    result<-list()
    
    method=c('proprosed','FMR_MCP','kmeans_MCP','MCP')
    
    for (m in 1:init_time){
      print(paste0('The ',m,'th initialization'))
      if (m==1){
        class_old<-sample(1:2,n,replace=T)
        init_set[[1]]<-class_old
      } else {
        aa=0
        while (aa<0.5*n){
          for (ii in 1:(m-1)){
            class_old<-sample(1:2,n,replace=T)
            aa=sum(abs(class_old-init_set[[ii]]))
            if (aa>0.5*n){
              break
            }
          }
        }
      }
      
      init_set[[m]]<-class_old
    }
    
    
    for (tt in 1:4){
      
      if (tt==1){
        l1=length(lambda1_set)
        l2=length(lambda2_set)
        lambda1=lambda1_set[ceiling(l1/2)]
        lambda2=lambda2_set[ceiling(l2/2)]
        bic_init<-matrix(1e+10,init_time,1)
        for (m in 1:init_time){
          # print(paste0('The ',m,'th initialization'))
          class_old=init_set[[m]]
          temp<-frm_corr(x,z,y,K,lambda1,lambda2,class_old,ksi=6,tau=tau,maxiter=maxiter,iter_v=iter_v)
          delta=temp$delta
          alpha=temp$alpha
          beta=temp$beta
          mu=temp$mu
          sigma=temp$sigma
          lik<-L(y,x,z,alpha,beta,mu,sigma)
          bic_init[m] <- -2*lik+(sum(beta!=0)+sum(alpha!=0))*log(n)/n
        }
        id=which.min(bic_init)
        class_old=init_set[[id]]
        bic=matrix(1e+20,l1*l2,1)
        result_temp<-list()
        lll=1
        for (m1 in 1:l1){
          print(paste0('The ',m1,'th lambda'))
          lambda1=lambda1_set[m1]
          for (m2 in 1:l2){
            lambda2=lambda2_set[m2]
            temp<-frm_corr(x,z,y,K,lambda1,lambda2,class_old,ksi=6,tau=tau,maxiter=maxiter,iter_v=iter_v)
            result_temp[[lll]]<-temp
            delta=temp$delta
            alpha=temp$alpha
            beta=temp$beta
            mu=temp$mu
            sigma=temp$sigma
            lik<-L(y,x,z,alpha,beta,mu,sigma)
            bic[lll] <- -2*lik+(sum(beta!=0)+sum(alpha!=0))*log(n)/n
            lll=lll+1
          }
        }
        id=which.min(bic)
        result[[tt]]<-result_temp[[id]]
        
        result_temp1<-result_temp
        
        delta=result[[tt]]$delta
        s.hat<-rep(0,n)
        for (i in 1:n) {
          if(delta[i,1]>delta[i,2]){s.hat[i]=1}else{s.hat[i]=2}
        }
        
        result[[tt]]$s.hat=s.hat
      } else if (tt==2) {
        
        l1=length(lambda1_set)
        lambda1=lambda1_set[ceiling(l1/2)]
        bic_init<-matrix(1e+10,init_time,1)
        for (m in 1:init_time){
          #print(paste0('The ',m,'th initialization'))
          class_old=init_set[[m]]
          temp<-frm_corr(x,z,y,K,lambda1,0,class_old,ksi=6,tau=tau,maxiter=maxiter,iter_v=iter_v)
          delta=temp$delta
          alpha=temp$alpha
          beta=temp$beta
          mu=temp$mu
          sigma=temp$sigma
          lik<-L(y,x,z,alpha,beta,mu,sigma)
          bic[m] <- -2*lik+(sum(beta!=0)+sum(alpha!=0))*log(n)/n
        }
        id=which.min(bic_init)
        class_old=init_set[[id]]
        bic=matrix(1e+20,l1,1)
        result_temp<-list()
        for (m in 1:l1){
          #print(paste0('The ',m,'th lambda'))
          lambda1=lambda1_set[m]
          temp<-frm_corr(x,z,y,K,lambda1,0,class_old,ksi=6,tau=tau,maxiter=maxiter,iter_v=iter_v)
          result_temp[[m]]<-temp
          delta=temp$delta
          alpha=temp$alpha
          beta=temp$beta
          mu=temp$mu
          sigma=temp$sigma
          lik<-L(y,x,z,alpha,beta,mu,sigma)
          bic[m] <- -2*lik+(sum(beta!=0)+sum(alpha!=0))*log(n)/n
        }
        id=which.min(bic)
        result[[tt]]<-result_temp[[id]]
        
        result_temp2<-result_temp
        
        delta=result[[tt]]$delta
        s.hat<-rep(0,n)
        for (i in 1:n) {
          if(delta[i,1]>delta[i,2]){s.hat[i]=1}else{s.hat[i]=2}
        }
        
        result[[tt]]$s.hat=s.hat
        
        
      } else if (tt==3) {
        
        y3<-kmeans(y,2)
        class_3<-y3$cluster
        y3<-cbind(y,class_3)
        yid1=which(y3[,2]==1)
        y3_1<-y3[yid1,1]
        
        yid2=which(y3[,2]==2)
        y3_2<-y3[yid2,1]
        X <-cbind(x,z)
        X3_1<-X[yid1,]
        X3_2<-X[yid2,]
        cvfita3_1 <- ncvreg(X3_1,y3_1)
        cvfita3_2 <- ncvreg(X3_2,y3_2)
        df1<-colSums(cvfita3_1$beta!=0)
        n1<-dim(X3_1)[1]
        BIC1<-log(cvfita3_1$loss)+df1*log(n1)/n1
        n2<-dim(X3_2)[1]
        df2<-colSums(cvfita3_2$beta!=0)
        BIC2<-log(cvfita3_2$loss)+df2*log(n2)/n2
        
        lambda3_1=which.min(BIC1)
        lambda3_2=which.min(BIC2)
        temp1=cvfita3_1$beta[,lambda3_1]
        temp2=cvfita3_2$beta[,lambda3_2]
        alpha=matrix(0,p,K)
        beta=matrix(0,q,K)
        
        alpha[,1]=temp1[2:(p+1)]
        alpha[,2]=temp2[2:(p+1)]
        
        beta[,1]=temp1[(p+2):(p+q+1)]
        beta[,2]=temp2[(p+2):(p+q+1)]
        
        s.hat<-class_3
        
        result[[tt]]<-list(alpha=alpha,beta=beta,s.hat=s.hat)
      } else if (tt==4) {
        ###A4
        X <-cbind(x,z)
        fita4 <- ncvreg(X,y)
        df1<-colSums( fita4 $beta!=0)
        BIC1<-log(fita4$loss)+df1*log(n)/n
        lambda3_1=which.min(BIC1)
        temp1=cvfita3_1$beta[,lambda3_1]
        
        alpha=matrix(0,p,K)
        beta=matrix(0,q,K)
        
        alpha[,1]=temp1[2:(p+1)]
        alpha[,2]=temp1[2:(p+1)]
        
        beta[,1]=temp1[(p+2):(p+q+1)]
        beta[,2]=temp1[(p+2):(p+q+1)]
        
        s.hat<-rep(1,n)
        
        result[[tt]]<-list(alpha=alpha,beta=beta,s.hat=s.hat)
      }
      
      
      alpha=result[[tt]]$alpha
      beta=result[[tt]]$beta
      s.hat=result[[tt]]$s.hat
      ##True positive rate of delta
      
      indicator<-rep(0,n)
      for (i in 1:n) {
        if(s.hat[i]==output$s[i]){indicator[i]=1}else{indicator[i]=0}
      }
      
      if (sum(indicator)/n<1-sum(indicator)/n){
        alpha=alpha[,c(2,1)]
        beta=beta[,c(2,1)]
      }
      
      if (tt==1 || tt==2){
        y_hat1=x_test%*%alpha[,1]+z_test%*%beta[,1]
        y_hat2=x_test%*%alpha[,2]+z_test%*%beta[,2]
        mu=result[[tt]]$mu
        class.hat=rbinom(n_test,1,mu[1])
        #class.hat=output_test$s
        y_hat=y_hat1
        y_hat[class.hat==0]=y_hat2[class.hat==0]
        
      } else if (tt==3){
        y_hat1=x_test%*%alpha[,1]+z_test%*%beta[,1]
        y_hat2=x_test%*%alpha[,2]+z_test%*%beta[,2]
        
        y_mean1=mean(y[s.hat==1])
        y_mean2=mean(y[s.hat==2])
        y_hat=y_hat1
        y_hat[abs(y_hat1-y_mean1)>abs(y_hat2-y_mean2)]=y_hat2[abs(y_hat1-y_mean1)>abs(y_hat2-y_mean2)]
      } else if (tt==4){
        y_hat=x_test%*%alpha[,1]+z_test%*%beta[,1]
      }
      
      pmse=median((y_test-y_hat)^2)
      
      
      df<-sum(beta!=0)+sum(alpha!=0)
      tp<-max(sum(indicator)/n,1-sum(indicator)/n)
      alpha.mse<-sqrt(sum((alpha-alpha.set)^2)) ####################
      beta.mse<-sqrt(sum((beta-beta.set)^2)) #############################
      alpha.tp<-GetFPTP(as.vector(alpha.set),as.vector(alpha))
      beta.tp<-GetFPTP(as.vector(beta.set),as.vector(beta))
      
      case1[(tt-1)*repli+re,]<-c(tp,alpha.mse,alpha.tp$TP,alpha.tp$FP,beta.mse,beta.tp$TP,beta.tp$FP,df,alpha.tp$TPR,alpha.tp$FPR,beta.tp$TPR,beta.tp$FPR,pmse)
      
    }
    print('---------------------------------')
    for (tt in 1:4){
      print(c(method[tt],round(case1[(tt-1)*repli+re,],2)),quote = FALSE)
      repli.out[tt,]<-c(method[tt],round(case1[(tt-1)*repli+re,],2))
    }
    write.csv(repli.out,file =  paste0("repli",re,'case',case,".csv") )
    print('---------------------------------')
  }
  
  
  case1[4*repli+1,]=colMeans(case1[1:repli,])
  case1[4*repli+2,]=colMeans(case1[(repli+1):(2*repli),])
  case1[4*repli+3,]=colMeans(case1[(2*repli+1):(3*repli),])
  case1[4*repli+4,]=colMeans(case1[(3*repli+1):(4*repli),])
  
  
  print('---------------------------------')
  print('FINAL RESULTS')
  for (tt in 1:4){
    print(c(method[tt],round(case1[4*repli+tt,],2)),quote = FALSE)
    
  }
  
  write.csv(case1,file=paste0('result_case',case,'.csv'))
}
cost<-Sys.time()-old
print(cost)
