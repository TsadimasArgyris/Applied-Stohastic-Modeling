#constructing ci LRT-based , profile likelihood

whales <- scan(r"(C:\Users\rg_ts\OneDrive\Υπολογιστής\MsC\Applied Stohastic Modeling\1st lect\whales.txt)")
whales

#x~Gamma(a,b)
gamma_loglikelihood = function(theta,k){
  n = length(k)
  a = theta[1]
  b = theta[2]
  y = n*a*log(b) - n*lgamma(a) + (a-1)*sum(log(k)) - b*sum(k)
  return(y)
}

#gamma log-likelihood - using dgamma 
gamma_ll_lazy <- function(theta,k){
  sum(dgamma(k,shape=theta[1],rate=theta[2],log=TRUE))
}

?optim

outgam <- optim(c(1,2),gamma_ll_lazy,method='BFGS',k=whales,hessian=T,control=list(fnscale=-1))
outgam

sqrt(diag(solve(-outgam$hessian)))   
# this the se of a 0.1423873 and this is the se of b 0.2754907

#or

-solve(outgam$hessian)
#[,1]       [,2]
#[1,] 0.02027415 0.03345575
#[2,] 0.03345575 0.07589512
var<- .Last.value
var
sqrt(var[1,1]) #se of a 0.1423873
sqrt(var[2,2])  #se of b 0.2754907 

#the variance of a is 0.02027415
#the variance of b is 0.07589512

#-------CONFIDENCE INTERVAL for a ------
outgam$par[1] + c(-1,1)*1.96*sqrt(diag(solve(-outgam$hessian)))[1]   
##-------------- CONFIDENCE INTERVAL for b -----------
outgam$par[2] + c(-1,1)*1.96*sqrt(diag(solve(-outgam$hessian)))[2]  


#now the contour plot
alpha <- seq(0.5,2,0.01)
beta <- seq(1.5,3.5,0.01)
#intialize the log likelihood
#it is a matrix length(a)*length(b) of zeros initially 
logL <- matrix(0,length(alpha),length(beta))


#calculate for each pair of a and b the values of the log likelihood
i <- 0
for (a in alpha){
  i<-i+1; j<-0
  for (b in beta){
    j<-j+1
    logL[i,j] <- gamma_loglikelihood(c(a,b),k=whales)
  }
}

#plot the loglikelihood function
contour(alpha,beta,logL)

#the MLEs of a and b
outgam$par
##[1] 1.595409 2.632693

# + where maximum loglikelihood is
points(1.595409, 2.632693,pch=3)

#this is the maximum log likelihood of gamma and it is -92.72379
outgam$value
##[1] -92.72379

#qchisq(0.95,2) returns the 95th percentile of the chi-squared distribution with 2 df
# outgam$value-0.5*qchisq(0.95,2) in order to see the interval in which log(l(θ,δ)) exceeds X^2 1-a
#add=T:This argument indicates that the contour line should be added to the existing plot
contour(alpha,beta,logL,levels=outgam$value-0.5*qchisq(0.95,2),col='red',add=T)



#--- calculating the 95% confidence intervals  ------
#------diag(solve(-outgam$hessian)) computes the variance-covariance matrix for the MLE ----
#by taking the inverse of the negative Hessian matrix from outgam.
#The diag function then extracts the variances (diagonal elements) for each parameter.
#bazei sqrt gia na parei tin tipiki apoklisi σ
outgam$par[1] + c(-1,1)*1.96*sqrt(diag(solve(-outgam$hessian)))[1]
##[1] 1.316330 1.874488
outgam$par[2] + c(-1,1)*1.96*sqrt(diag(solve(-outgam$hessian)))[2]
##[1] 2.092731 3.172655

#The lty=3 argument specifies the line type to be dashed (διακεκομμενη)
#adds vertical lines to the plot at the values representing the bounds of the CI for a
abline(v=c(1.316330, 1.874488),lty=3)
#adds vertical lines to the plot at the values representing the bounds of the CI for b
abline(h=c(2.092731, 3.172655),lty=3)


#finally calculating the LRT -based CI s


#log likelihood for beta given alpha
alpha_gamma_loglikelihood = function(b,a){
  gamma_loglikelihood(c(a,b),whales)
}
#log likelihood for alpha  given beta
beta_gamma_loglikelihood = function(a,b){
  gamma_loglikelihood(c(a,b),whales)
}

a = seq(1,3,by = 0.01)
b = seq(2,4,by = 0.01)

profile_likelihood_alpha = NULL
for (i in 1:length(a)){
  maxb = optim(1,fn = alpha_gamma_loglikelihood,method='BFGS', a = a[i],
               hessian=T,control=list(fnscale=-1))
  profile_likelihood_alpha = c(profile_likelihood_alpha,maxb$value)
}

plot(a,profile_likelihood_alpha,type='l',xlab = "Alpha Values", 
     ylab = "loglikelihood of alpha", 
     main = "Profile Likelihood of Alpha")
abline(h = max(profile_likelihood_alpha) - 1.92, col="red")


log_diff = function(a, ...){
  mle_beta = optim(1,alpha_gamma_loglikelihood,a = a,method='BFGS',
                   hessian=T,control=list(fnscale=-1))$par
  loglikelihood = alpha_gamma_loglikelihood(mle_beta,a)
  difference = max(profile_likelihood_alpha) - 1.92
  loglikelihood - difference
}

lower = uniroot(log_diff,interval = c(1,1.59))$root
upper = uniroot(log_diff,interval = c(1.59,3))$root

lower
upper

cat("Confidence Interval for mle of a is:", lower,",", upper)

plot(a,profile_likelihood_alpha,type='l',xlab = "Alpha Values", 
     ylab = "loglikelihood of alpha", 
     main = "Profile Likelihood of Alpha")
abline(h = max(profile_likelihood_alpha) - 1.92, col="red")
abline(v = lower,col = "blue")
abline(v = upper,col = "blue")


profile_likelihood_beta = NULL
for (i in 1:length(b)){
  maxa = optim(1,fn = beta_gamma_loglikelihood,method='BFGS',b = b[i],
               hessian=T,control=list(fnscale=-1))
  profile_likelihood_beta = c(profile_likelihood_beta,maxa$value)
}

plot(b,profile_likelihood_beta,type='l',xlab = "Beta Values", 
     ylab = "loglikelihood of beta", 
     main = "Profile Likelihood of Beta")
abline(h = max(profile_likelihood_beta) - 1.92, col="red")

log_diff_beta = function(b, ...){
  mle_alpha = optim(1,beta_gamma_loglikelihood,b = b,method='BFGS',
                    hessian=T,control=list(fnscale = -1))$par
  loglikelihood_b = beta_gamma_loglikelihood(mle_alpha,b)
  difference_b = max(profile_likelihood_beta) - 1.92
  loglikelihood_b - difference_b
}


lower_b = uniroot(log_diff_beta,interval = c(2,2.5))$root
upper_b = uniroot(log_diff_beta,interval = c(2.5,4))$root

lower_b
upper_b

cat("Confidence Interval for mle of b is:",lower_b,",",upper_b)

plot(b,profile_likelihood_beta,type='l',xlab = "Beta Values", 
     ylab = "loglikelihood of beta", 
     main = "Profile Likelihood of Beta")
abline(h = max(profile_likelihood_beta) - 1.92, col="red")
abline(v = lower_b,col = "blue")
abline(v = upper_b,col = "blue")


##############################################
##############################################
###            for exp(lambda)        ########
##############################################



logl <- function(lambda) {
  #exponential log-likelihood - hard-wired
  length(whales)*log(lambda)-lambda*sum(whales)
}
#run quasi - Newton to find the max od logl
outexp <- optim(0.5, logl, method='BFGS', control=list(fnscale=-1), hessian=TRUE)
outexp

outexp$par #this is my MLE
outexp$hessian #my hessian at the MLE

#hence var(λ) is -solve(outexp$hessian)
-solve(outexp$hessian)


lmax <- outexp$value

sqrt(diag(solve(-outexp$hessian)))                     
outexp$par + c(-1,1)*1.96*sqrt(diag(solve(-outexp$hessian)))
#this is the CI based on the asymptotic normal [1.427000 1.873387]



# LRT-based  CI

log_likelihood_difference<-function(lambda){
  logl(lambda)-lmax+1.92
}

plot(seq(1,3,0.01),logl(seq(1,3,0.01)),type='l')
abline(h=outexp$value-1.92,lty=2,col="red")

outexp$value-1.92 # that is the cutoff point

lower_exp = uniroot(log_likelihood_difference,interval = c(1,1.65))$root
upper_exp = uniroot(log_likelihood_difference,interval = c(1.65,2))$root

lower_exp
upper_exp

abline(v = lower_exp ,col = "blue")
abline(v = upper_exp,col = "blue")

