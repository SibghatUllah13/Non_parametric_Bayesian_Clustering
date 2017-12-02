
# Question 2 --------------------------------------------------------------




require(Compositional)
require(MCMCpack)

#Stick Breaking

#G0 = normal ( 0,1)



#draw from G0

G0 = function(n){return(rnorm(n,0,1))} #random points on the real line


draw_betas = function(n,M){return(rbeta(n,1,M))}




stick_breaking = function(nn = 10,M = 2,thetas = G0(nn)){
  
   
  betas = draw_betas(nn,M)
  
  weights = c()
  weights[1] = betas[1]
  
  for(i in 2:nn){
    weights[i] = betas[i] * prod(1 - betas[1:(i-1)])
  }
  
  theta <- sample(thetas, prob = weights, replace = TRUE)
  df = data.frame(weights,thetas)
  return(df)
}



#plots the cdfs and generates them
gen_cdf_mean = function(runs = 10, nn = 100, M = 50,theta = thetas ){
  
  weights_df = data.matrix(matrix(ncol = nn, nrow = runs))
  samples_df = data.matrix(matrix(ncol = nn, nrow = runs))
  for(i in 1:runs)
  {
    out = stick_breaking(nn , M, theta)
    samples = out[,2]
    samples_df[i,] = samples
    weights = out[,1]
    weights_df[i,] = cumsum(weights)
    #plot((samples),(weights), type = 'h')
  }
  
  matplot(sort(samples),t(weights_df),col = 1:5,type = 'b')
  curve(pnorm(x),add = T, title(paste('M = ',M)))
  
  means = apply(weights_df,2,mean)
  variance = apply(weights_df,2,var)
  return(list(weights_df,means,variance))
}



#Realizations of DP for different values of M , along with prior mean and prior variances
# M = 0.2 -----------------------------------------------------------------
thetas = G0(100)
runs = 10
lst = gen_cdf_mean(runs = 10, nn = 100,M=0.2, theta = thetas)
wt = lst[[1]]
means = lst[[2]]
variance = lst[[3]]
par(mfrow=c(1,2))
plot(means,col='orchid')
plot(variance,col='blue')

# M = 5 -------------------------------------------------------------------
lst = gen_cdf_mean(runs = 10, nn = 100,M=5, theta = thetas)
wt = lst[[1]]
means = lst[[2]]     #this mean vector will hold the mean of generated realizations from DP i.e. E[G(B)]
variance = lst[[3]] 
par(mfrow=c(1,2))
plot(means,col='orchid')
plot(variance,col='blue')


# M = 25 ------------------------------------------------------------------
lst = gen_cdf_mean(runs = 10, nn = 100,M=25, theta = thetas)
wt = lst[[1]]
means = lst[[2]]
variance = lst[[3]]
par(mfrow=c(1,2))
plot(means,col='orchid')
plot(variance,col='blue')

# M = 50 ------------------------------------------------------------------

lst = gen_cdf_mean(runs = 10, nn = 100,M=50, theta = thetas)
wt = lst[[1]]
means = lst[[2]]
variance = lst[[3]]
par(mfrow=c(1,2))
plot(means,col='orchid')
plot(variance,col='blue')



# M = 1000 ----------------------------------------------------------------
lst = gen_cdf_mean(runs = 10, nn = 100,M=1000, theta = thetas)
wt = lst[[1]]
means = lst[[2]]
variance = lst[[3]]
par(mfrow=c(1,2))
plot(means,col='orchid')
plot(variance,col='blue')






# Ferguson's definition ---------------------------------------------------
require(Compositional)
require(MCMCpack)


sample.dir<-function(nn=10, M=2, ln = 11){
  x <- seq(-4,4, length=ln)
  y<-c()
  y[1]<-pnorm(x[1])
  for(i in 2:ln) y[i]<-pnorm(x[i])-pnorm(x[(i-1)])
  y<-c(y, 1-pnorm(x[ln]))
  param<-M*y
  #return(param)
  sample.dir<-rdirichlet(nn,param)
  draw<-apply(t(sample.dir), 2, cumsum)
  return(draw)}



# M = 0.2 -----------------------------------------------------------------

draws<-sample.dir(nn = 10, M = 2 , ln = 10)
xx<-c(seq(-4,4, length=10),5)
matplot(xx, draws, col=1:10, type="b")
curve(pnorm(x), add=T)
mean_vec = apply(draws,2,mean) #vector to hold values of mean for all the sets B1,B2,.....Bk, first element is the mean of observations
var_vec = apply(draws,2,var)
par(mfrow=c(1,2))
plot(mean_vec,col='orchid')
plot(var_vec,col='blue')



# M = 5 -------------------------------------------------------------------

draws<-sample.dir(nn = 100, M = 5 , ln = 10)
xx<-c(seq(-4,4, length=10),5)
matplot(xx, draws, col=1:10, type="b")
curve(pnorm(x), add=T)
mean_vec = apply(draws,2,mean) #vector to hold values of mean for all the sets B1,B2,.....Bk, first element is the mean of observations
var_vec = apply(draws,2,var)
par(mfrow=c(1,2))
plot(mean_vec,col='orchid')
plot(var_vec,col='blue')





# M = 25 ------------------------------------------------------------------
draws<-sample.dir(nn = 10, M = 25 , ln = 10)
xx<-c(seq(-4,4, length=10),5)
matplot(xx, draws, col=1:10, type="b")
curve(pnorm(x), add=T)
mean_vec = apply(draws,2,mean) #vector to hold values of mean for all the sets B1,B2,.....Bk, first element is the mean of observations
var_vec = apply(draws,2,var)
par(mfrow=c(1,2))
plot(mean_vec,col='orchid')
plot(var_vec,col='blue')


# M = 50 ------------------------------------------------------------------

draws<-sample.dir(nn = 10, M = 50 , ln = 10)
xx<-c(seq(-4,4, length=10),5)
matplot(xx, draws, col=1:10, type="b")
curve(pnorm(x), add=T)
mean_vec = apply(draws,2,mean) #vector to hold values of mean for all the sets B1,B2,.....Bk, first element is the mean of observations
var_vec = apply(draws,2,var)
par(mfrow=c(1,2))
plot(mean_vec,col='orchid')
plot(var_vec,col='blue')




# adding a prior on M  ----------------------------------------------------


M_vec = rgamma(5,3,3)
par(mfrow = c(5,3))
for (i in M_vec){
  lst = gen_cdf_mean(M=i)   #stick breaking  used here
  wt = lst[[1]]
  means = lst[[2]]
  variance = lst[[3]]
  plot(means)
  plot(variance)
}



# Question 3 --------------------------------------------------------------
library(coda)
library(distr)

#generate mixture of normals
myMix <- UnivarMixingDistribution(Norm(mean=2.5, sd=0.5), 
                                  Norm(mean=0.5, sd=0.7),
                                  Norm(mean=1.5, sd=2),
                                  mixCoeff=c(0.5, 0.3, 0.2))




plot(myMix, to.draw.arg="p") 





#posterior function to calculate dp posterior
dp_post = function(data,M,nn= 1000,thetas){
  
  nsize = length(data)
  Ghat = nsize/(nsize+M) * data + M/(nsize+M) * thetas   #this is the weighted distribution of empirical data and prior distro
  return(gen_cdf_mean(runs = 10, nn = nsize,M = nsize + M, theta = Ghat))
}

# n = 20, M = 5 -----------------------------------------------------------
par(mfrow = c(3,2))

nn=20
data = rnorm(nn)
theta_points <- r(myMix)(nn)



lst_dp = dp_post(data,M = 5, nn ,theta = theta_points)
wt = lst_dp[[1]]
post_mean = lst_dp[[2]]
post_var = lst_dp[[3]]
plot(post_mean)
plot(post_var)

#HPDinterval(mcmc(wt), 0.95)
#This will compute credible intervals for each of the runs of the stick breaking process
credible_int <- apply(wt, 2, function(col) HPDinterval(as.mcmc(col), prob=0.95))
credible_int
# n = 200, M = 5 ----------------------------------------------------------

nn=200
data = rnorm(nn)
theta_points <- r(myMix)(nn)



lst_dp = dp_post(data,M = 5, nn ,theta = theta_points)
wt = lst_dp[[1]]
post_mean = lst_dp[[2]]
post_var = lst_dp[[3]]
plot(post_mean)
plot(post_var)

#HPDinterval(mcmc(wt), 0.95)
#This will compute credible intervals for each of the runs of the stick breaking process
credible_int <- apply(wt, 2, function(col) HPDinterval(as.mcmc(col), prob=0.95))
credible_int

# n = 2000, M = 5 ---------------------------------------------------------

nn=2000
data = rnorm(nn)
theta_points <- r(myMix)(nn)



lst_dp = dp_post(data,M = 5, nn = 10,theta = theta_points)
wt = lst_dp[[1]]
post_mean = lst_dp[[2]]
post_var = lst_dp[[3]]
plot(post_mean)
plot(post_var)

#HPDinterval(mcmc(wt), 0.95)
#This will compute credible intervals for each of the runs of the stick breaking process
credible_int <- apply(wt, 2, function(col) HPDinterval(as.mcmc(col), prob=0.95))
credible_int
# n = 20, M = 100 ---------------------------------------------------------
nn=20
data = rnorm(nn)
theta_points <- r(myMix)(nn)



lst_dp = dp_post(data,M = 100, nn ,theta = theta_points)
wt = lst_dp[[1]]
post_mean = lst_dp[[2]]
post_var = lst_dp[[3]]
plot(post_mean)
plot(post_var)

#HPDinterval(mcmc(wt), 0.95)
#This will compute credible intervals for each of the runs of the stick breaking process
credible_int <- apply(wt, 2, function(col) HPDinterval(as.mcmc(col), prob=0.95))
credible_int
# n = 200, M = 100 --------------------------------------------------------
nn=200
data = rnorm(nn)
theta_points <- r(myMix)(nn)



lst_dp = dp_post(data,M = 100, nn = 10,theta = theta_points)
wt = lst_dp[[1]]
post_mean = lst_dp[[2]]
post_var = lst_dp[[3]]
plot(post_mean)
plot(post_var)

#HPDinterval(mcmc(wt), 0.95)
#This will compute credible intervals for each of the runs of the stick breaking process
credible_int <- apply(wt, 2, function(col) HPDinterval(as.mcmc(col), prob=0.95))
credible_int

# n = 2000, M = 100 -------------------------------------------------------

nn=2000
data = rnorm(nn)
theta_points <- r(myMix)(nn)



lst_dp = dp_post(data,M = 100, nn ,theta = theta_points)
wt = lst_dp[[1]]
post_mean = lst_dp[[2]]
post_var = lst_dp[[3]]
plot(post_mean)
plot(post_var)



#HPDinterval(mcmc(wt), 0.95)
#This will compute credible intervals for each of the runs of the stick breaking process
credible_int <- apply(wt, 2, function(col) HPDinterval(as.mcmc(col), prob=0.95))
credible_int


# For higher values of M, the distribution will shift towards the prior that is the mixture of normals
# and for higher values of n i.e. the sample size, the distributions will shift towards the empirical distribution, 
# which is generated from pnorm(x)


