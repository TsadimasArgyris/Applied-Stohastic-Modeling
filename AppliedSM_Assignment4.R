#read in Fecundability data
fec <- read.table("C:/Users/rg_ts/OneDrive/Υπολογιστής/MsC/Applied Stohastic Modeling/data/fecundability.txt", header=F, sep = "", col.names = c("cycle", "fsmoke", "fnsmoke"))

xns <- rep(fec$cycle,fec$fnsmoke)              #focus on non-smokers, for illustration

hist(xns)
plot(xns)

xs<-rep(fec$cycle,fec$fsmoke)
table(xs)

hist(xs)
plot(xs)

#the mass function of a geometric distribution
geo_pmf = function(x,theta){
  p = plogis(theta)
  dgeom(x-1,prob=p)
}

#its LogLikelihood function
lgeo <- function(thet,x){ 
  p <- plogis(thet)  
  sum(dgeom(x-1,prob=p,log=TRUE)) 
}
#the mixture geometric function
mixturegeo1<-function(thet,x){
  p1 <- plogis(thet[1])
  p2 <- plogis(thet[2])
  w <- thet[3]
  w*dgeom(x-1,prob=p1)+(1-w)*dgeom(x-1,prob = p2)
}

#the LogLikelihood of the mixture
llmixturegeo <- function(thet, x){
  p1 <- plogis(thet[1])
  p2 <- plogis(thet[2])
  w <- plogis(thet[3])
  
  
  logLikelihood <- sum(log(w * dgeom(x-1, prob=p1) + (1-w) * dgeom(x-1, prob=p2)))
}

llmixturegeo(c(7,13, 0.5),xns) #just checkin

#fitting my model
out <- optim(c(0.2,0.3, 0.5),llmixturegeo,method='BFGS',control=list(fnscale=-1),x=xns)     
out
#$value
#[1] -907.7368

#and now i transform back my parameters
out$par[1:3] <- plogis(out$par[1:3])
out$par
#[1] 0.2433558 0.6096634 0.5365670


# this is the saturated model
lsatgeo <- function(x){
  #geometric saturated log-likelihood 
  lgeo(qlogis(1/x),x=x) 
} 

#'fit' model 2
out2value <- lsatgeo(x=xns)
out2value
#[1] -599.0804


#i will use deviance to compare the models ,since they are nested.
deviance <- 2 * (out2value - out$value)
deviance
#[1] 617.3129

glm(cbind(1,xns-1)~1,family=binomial())$deviance

#if i had another model i would compare their deviances
#although my deviance seems pretty large ,i need context
p1<-out$par[1]
p2<-out$par[2]
prob_w<-out$par[3]

mixgeom_sample = function(){
  mysample = NULL
  for (i in 1:length(xns)){
    b = rbinom(1,1,prob_w)
    f = b*(rgeom(1,p1)+1) + (1-b)*(rgeom(1,p2)+1)
    mysample = c(mysample,f)
  }
  mysample
}

Deviance_sim<-numeric(250)
for (i in 1:250){
  newsample = mixgeom_sample()
  out_value_sim<-optim(c(0.1,0.2,0.35),fn = llmixturegeo,method='L-BFGS-B',x = round(newsample),
                       hessian=T,control=list(fnscale=-1))$value
  out_2_value_sim<-lsatgeo(x=newsample)
  Deviance_sim[i]<-2 * (out_2_value_sim - out_value_sim)
}


pvalue = sum(Deviance_sim>=deviance)/length(Deviance_sim)
pvalue
#[1] 0.588 
#the model fits the data well

hist(Deviance_sim , col = "skyblue", xlim=c(min(Deviance_sim),max(max(Deviance_sim),deviance)))
abline(v=deviance,col='red')
