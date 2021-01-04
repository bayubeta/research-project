library(mvtnorm)

mu = matrix(c(4,4), ncol = 1)
Sigma = matrix(c(1,0,0,1), ncol=2)
N = 1000
Fmat = matrix(c(1,0,0,1, -1, 1, 1.1, -1), byrow = TRUE, ncol = 2)

m = nrow(Fmat)
d = ncol(Fmat)
g = matrix(rep(0, m), ncol = 1)

initial_X = matrix(c(2.1, 2.2), ncol = 1)


######### transform to whitened frame ########
initial_X = initial_X - mu
g = g + Fmat%*%mu # shift g

# whitening
R = chol(Sigma)
W = solve(t(R)) # solve(Sigma) = t(W)%*%W
initial_X = W%*%initial_X
Fmat = Fmat%*%t(R) # t(R) = solve(W)


llik = function(x, Fmat, g){
  # assuming F and g are already transformed
  # to the whitened frame
  x = as.vector(x)
  L = dmvnorm(x, log = TRUE)
  constr = prod(Fmat%*%x + g > 0)
  
  if (constr){
    return(L)
  }
  else{
    return(-1000) # returns a really low value
  }
}

last_X = initial_X

# placeholder matrix for values:
Xs = last_X

iter_times = c()

set.seed(12345)

# begin loop:
while(ncol(Xs)<N){
  ######### start time #########
  start.time = Sys.time()
  ##############################
  
  # prior ellipse:
  v = t(rmvnorm(1, mean = rep(0, d)))
  
  # log-likelihood threshold:
  U = runif(1, 0, 1)
  log_y = llik(last_X, Fmat = Fmat, g = g) + log(U)
  
  # initial proposal:
  t = runif(1, 0, 2*pi)
  t_min = t - 2*pi
  t_max = t
  
  accept = 0
  
  while (accept == 0){
    # check new proposal:
    x_prime = last_X*cos(t) + v*sin(t)
    log_x_prime = llik(x_prime, Fmat = Fmat, g = g)
    
    if (log_x_prime > log_y){
      # accept
      new_X = x_prime
      accept = 1
    }
    else{
      # shrink the bracket bracket:
      if (t<0){
        t_min = t
      }
      else{
        t_max = t
      }
      # sample new angle
      t = runif(1, t_min, t_max)
    }
  }
  
  # check if satisfy constraints
  if (all(Fmat %*% new_X + g > 0)){
    
    # append new X value
    Xs = cbind(Xs, new_X)
    # set new_X as last_X
    last_X = new_X
    if (ncol(Xs)%%100 == 0){
      cat("sample no:", ncol(Xs),"\n")
    }
  }
  
  ########## end time ##########
  end.time = Sys.time()
  iter_times = cbind(iter_times, as.vector(end.time - start.time))
  ##############################
  
}


# transform back
Xs = t(R)%*%Xs + matrix(mu, nrow = length(mu), ncol = N)

plot(Xs[1,], Xs[2,], cex = 0.5,
     xlab = "X1", ylab = "X2",
     xlim = c(1.5,6.5), ylim = c(1.5,6.5),
     main = "Elliptical Slice Sampling")


# moments of predicted distributions
# predicted means
# rowMeans(Xs)
# 
# # predicted cov
# cov(t(Xs))
# 
# # theoretical moments
# tmvtnorm::mtmvnorm(mean = as.vector(mu),
#                    lower = c(0,0), upper = c(Inf, Inf),
#                    sigma = Sigma)

saveRDS(Xs, file = "x_ESS.rds")
saveRDS(iter_times, file = "t_ESS.rds")