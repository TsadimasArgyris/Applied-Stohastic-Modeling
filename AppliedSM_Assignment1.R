#Assignment 1
#1
#reading the data , my x 's
data = c(1.77, -0.23, 2.76, 3.80, 3.47, 56.75, -1.34, 4.24, -2.44, 3.29
                  , 3.71, -2.40, 4.53, -0.07,
                  -1.05, -13.87, -2.53, -1.75, 0.27, 43.21)

pi
mean_data <- mean(data) #finding the mean of my dataset
mean_data

#my pdf cauchy(θ,1)
cauchy_pdf <- function(x,theta){
  1/(pi*(1+(x-theta)^2))
}

#and its logLikelihood
logLikelihood <- function(x,theta){
  n <- length(x)
  -n*log(pi)-sum(log((1+(x-theta)^2)))
}

#creating values to use to draw the line/loglikelihood on the x axis
theta_values = seq(-20,60,by = 0.1)

#now taking values for the y axis for its x (theta) in the graph
log_likelihoods <- sapply(theta_values, logLikelihood, x = data)#(argument,function,extra arg)

#now plotting the graph
plot(theta_values, log_likelihoods, type = 'l', col = 'skyblue1',
     main = 'Log Likelihood of Cauchy Distribution Cauchy(θ,1)',
     xlab = 'Theta', ylab = 'Log Likelihood', lwd = 2)

#now the Newton Raphson for the MLE

#first derivative
grad = function(theta){
  sum(2*(data-theta)/(1+(data - theta)^2))
}

#second derivative
hessian = function(theta){
  sum(2*((data-theta)^2 - 1)/(1+(data-theta)^2)^2)
}

# NR Algorithm

newton_raphson = function(x0){
  iter = 0  #counter
  t=0      #checks when NR stops
  while(t==0){
    iter = iter + 1
    xnew = x0 - grad(x0)/hessian(x0) #NR formula
    print(c(xnew,iter))
    if(abs(x0-xnew)<=1e-6){ #convergence criterion
      t=1
      x0=xnew
    }
    x0=xnew
  }
  return(x0)
}

newton_raphson(-11)
newton_raphson(8)
newton_raphson(mean(data))
x0 = c(-1,0,1.5,4,4.7,7,38)

results = sapply(x0,newton_raphson)

results_df = data.frame(Initial_Theta = x0, Final_Theta = results)
results_df

ismle = sapply(x0,hessian)
ismle

outcome = NULL
for(i in 1:length(results)){
  if (ismle[i]<0){
    outcome = c(outcome,"yes")
  }else{
    outcome =  c(outcome,"no")
  }
}
outcome

mle_verif = data.frame(Initial_Theta = x0,solution = results, hessian_value = ismle, is.mle = outcome)
mle_verif

#2

newton_raphson_mod = function(x0){
  iter = 0  #counter
  t = 0     #checks when NR stops
  while(t == 0){
    iter = iter + 1
    hessianValue = hessian(x0)  # Store the hessian value in a variable
    if (hessianValue > 0){
      hessianValue <- -abs(hessianValue)  # Modify the variable if needed
    }
    delta <- -grad(x0) / hessianValue  # Use the variable here
    xnew = x0 + delta #NR formula
    while (logLikelihood(data, xnew) < logLikelihood(data, x0)){#checks if i jump over the maximum and i end up in a lower position,because i only expect to go up
      delta <- delta / 2  #so i cut the step in half
      xnew = x0 + delta
    }
    print(c(xnew, iter))
    if(abs(x0 - xnew) <= 1e-6){ #convergence criterion
      t = 1
    }
    x0 = xnew
  }
  return(x0)
}


modified_results = sapply(x0,newton_raphson_mod)

modified_results_df = data.frame(Initial_Theta = x0, Final_Theta = modified_results)
modified_results_df

