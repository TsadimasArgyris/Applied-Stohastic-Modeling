#Assignment 2
#NR for MLE of λ
whales <- scan(r"(C:\Users\rg_ts\OneDrive\Υπολογιστής\MsC\Applied Stohastic Modeling\1st lect\whales.txt)")
whales

n <- length(whales)
sum(whales)

#Suppose x~exp(λ)

logLikelihood<-function(lambda){
  n*log(lambda)-lambda*sum(whales)
}

gradll <- function(x,lambda){
  #first derivative of log likelihood
  (length(x)/lambda)-sum(x)
}

hessianll<-function(x,lambda){
  #second derivative of log likelihood
  -length(x)/(lambda^2)
}

# Newton - Raphson for MLE of lambda, when X~exp(λ)

newton_raphson<-function(x,lambda){
  iter<-0
  t<-0
  while(t==0){
    iter<-iter+1
    lambdanew<-lambda - gradll(x,lambda)/hessianll(x,lambda)
    print(c(lambdanew,iter))
    if(abs(lambdanew-lambda)<=1e-6){
      t<-1
      lambda <- lambdanew
    }
    lambda <- lambdanew
  }
  return(lambda)
}

newton_raphson(whales,0.1)

#FPI
###we set the second derivative equal to -1 #######
###
fixed_point<-function(x,lambda){
  iter<-0
  t<-0
  while(t==0){
    iter<-iter+1
    lambdanew<-lambda - gradll(x,lambda)/(-1)
    print(c(lambdanew,iter))
    if(abs(lambdanew-lambda)<=1e-6){
      t<-1
      lambda <- lambdanew
    }
    lambda <- lambdanew
    if (iter > 5000) break
  }
  return(lambda)
}

fixed_point(whales,0.1)
#does not converge

#scaled fixed point iteration

a<-0.01
scaled_fixed_point<-function(x,lambda){
  iter<-0
  t<-0
  while(t==0){
    iter<-iter+1
    lambdanew<-lambda - gradll(x,lambda)/(-1/a)
    print(c(lambdanew,iter))
    if(abs(lambdanew-lambda)<=1e-6){
      t<-1
      lambda <- lambdanew
    }
    lambda <- lambdanew
    if (iter > 5000) break
  }
  return(lambda)
}
scaled_fixed_point(whales,0.1)

#now the lazy methods
library(numDeriv)

?grad
?hessian

lazy_newton_raphson<-function(x,lambda){
  iter<-0
  t<-0
  while(t==0){
    iter<-iter+1
    lambdanew<-lambda - grad(logLikelihood,lambda)/hessian(logLikelihood,lambda)[1,1]
    print(c(lambdanew,iter))
    if(abs(lambdanew-lambda)<=1e-6){
      t<-1
      lambda <- lambdanew
    }
    lambda <- lambdanew
  }
  return(lambda)
}

lazy_newton_raphson(whales,0.1)

a<-0.01
lazy_scaled_fixed_point<-function(x,lambda){
  iter<-0
  t<-0
  while(t==0){
    iter<-iter+1
    lambdanew<-lambda - grad(logLikelihood,lambda)/(-1/a)
    print(c(lambdanew,iter))
    if(abs(lambdanew-lambda)<=1e-6){
      t<-1
      lambda <- lambdanew
    }
    lambda <- lambdanew
    if (iter > 5000) break
  }
  return(lambda)
}
lazy_scaled_fixed_point(whales,0.1)


##Now the same for x~Gamma(a,b)

GammalogLikelihood<-function(x,theta){
  n <- length(x)
  ll<- n*theta[1]*log(theta[2])-n*lgamma(theta[1])+(theta[1]-1)*sum(log(x))-theta[2]*sum(x)
  return(ll)
}

Gamma_gradll <- function(x,theta){
  n = length(x)
  #first derivative of log likelihood
  gdash<- c(n*log(theta[2])-n*digamma(theta[1])+sum(log(x)),
            n*theta[1]/theta[2] - sum(x))
  return(gdash)
}


Gamma_hessianll<-function(x,theta){
  n = length(x)
  #second derivative of log likelihood
  gddash<- matrix(c(-n*trigamma(theta[1]),n/theta[2],
                    n/theta[2],-n*theta[1]/theta[2]^2), byrow=T,nrow = 2)
  return(gddash)
}

# Newton - Raphson for MLE of alpha,beta when X~Gamma(a,b)

gamma_newton_raphson<-function(x,theta){
  iter<-0
  t<-0
  while(t==0){
    iter<-iter+1
    thetanew<-theta -  solve(Gamma_hessianll(x,theta))%*%Gamma_gradll(x,theta)
    print(c(thetanew,iter))
    if(sqrt(sum((thetanew - theta)^2))<=1e-6){
      #my convergence criterion here is the euclidean distance of the vectors
      t<-1
      theta <- thetanew
    }
    theta <- thetanew
  }
  return(theta)
}

gamma_newton_raphson(whales,c(0.1,0.2))

#now fpi

#scaled FPI
#we set the second derivative equal to -1
a<-0.003
gamma_fixed_point<-function(x,theta){
  iter<-0
  t<-0
  while(t==0){
    iter<-iter+1
    thetanew<-theta -  Gamma_gradll(x,theta)/(-1/a)
    print(c(thetanew,iter))
    if(sqrt(sum((thetanew - theta)^2))<=1e-6){
      t<-1
      theta <- thetanew
    }
    theta <- thetanew
    if (iter > 5000) break
  }
  return(theta)
}

gamma_fixed_point(whales,c(0.1,0.2))

#now the lazy methods

library(numDeriv)
#lazy NR
lazy_gamma_newton_raphson<-function(x,theta){
  iter<-0
  t<-0
  while(t==0){
    iter<-iter+1
    logLikelihoodWrapper <- function(theta) {
      GammalogLikelihood(x, theta)
    }
    thetanew <- theta - solve(hessian(logLikelihoodWrapper, theta))%*%
      grad(logLikelihoodWrapper,theta)
    print(c(thetanew,iter))
    if(sqrt(sum((thetanew - theta)^2))<=1e-6){
      #my convergence criterion here is the euclidean distance of the vectors
      t<-1
      theta <- thetanew
    }
    theta <- thetanew
  }
  return(theta)
}

lazy_gamma_newton_raphson(whales,c(0.1,0.2))

# lazy scaled fpi
a<-0.003
lazy_gamma_fixed_point<-function(x,theta){
  iter<-0
  t<-0
  while(t==0){
    iter<-iter+1
    logLikelihoodWrapper <- function(theta) {
      GammalogLikelihood(x, theta)
    }
    thetanew <- theta - grad(logLikelihoodWrapper,theta)/(-1/a)
    print(c(thetanew,iter))
    if(sqrt(sum((thetanew - theta)^2))<=1e-6){
      t<-1
      theta <- thetanew
    }
    theta <- thetanew
    if (iter > 5000) break
  }
  return(theta)
}

lazy_gamma_fixed_point(whales,c(0.1,0.2))
