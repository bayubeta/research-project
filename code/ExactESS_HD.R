library(mvtnorm)
library(sets)
library(rootSolve)

rm(list=ls())

lb = rep(-1, 10)
ub = rep(1, 10)
d = length(lb)

Fmat = rbind(diag(1,d), diag(-1,d))
g = matrix(c(-lb,ub), ncol = 1)

m = nrow(Fmat)
M = diag(d)

mu = matrix(rep(0,d), ncol = 1)

Sigma = matrix(rep(0.8, d*d), ncol = d, nrow = d)
diag(Sigma) = rep(1, d)

N = 10000
J = 10

initial_X = matrix(rep(0.5, d), ncol = 1)

######### transform to whitened frame ########
initial_X = initial_X - mu
g = g + Fmat%*%mu # shift g

# whitening
R = chol(Sigma)
W = solve(t(R)) # W satisfies solve(Sigma) = t(W)%*%W
initial_X = W%*%initial_X
Fmat = Fmat%*%t(R) # t(R) = solve(W)
##############################################


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


multRunif = function(Intv){
  intv = as.list(Intv)
  intv_vec = lapply(intv, FUN = function(x){r = as.vector(range(x)); c(r[1],r[2])})
  thr = cumsum(sapply(intv_vec, function(x){x[2]-x[1]}))
  max_thr = max(thr)
  thr = c(0, thr[1:(length(thr)-1)])
  u = runif(1, 0, max_thr)
  
  pos = sum(u>thr)
  
  offset = u - thr[pos]
  
  return(intv_vec[[pos]][1] + offset)
}


angles_to_intervals = function(t, x, Fmat, g){
  t_new = t %% (2*pi)
  # sort
  t_new = sort(t_new)
  eps = 10^(-3)
  t1 = t_new[1] - eps
  t2 = mean(t_new)
  t3 = t_new[2] + eps
  
  ts = c(t1, t2, t3)
  
  ps = sapply(ts, FUN = function(t){Fmat %*% (x*cos(t) + v*sin(t)) + g >0})
  ps_region = apply(ps, 2, all)
  
  bounds = matrix(c(0,t_new[1],
                    t_new[1], t_new[2],
                    t_new[2],2*pi),
                  ncol = 3)
  bounds = bounds[,ps_region]
  
  if (!is.null(dim(bounds))){
    theta_js = apply(bounds, 2, FUN = function(x){interval(x[1], x[2])})
    Theta_js = Reduce('|', theta_js)
    return(Theta_js)
  }
  else{
    Theta_js = interval(bounds[1], bounds[2])
    return(Theta_js)
  }
}


characterize_slice = function(x, v){
  # x: last_X, x[i-1]
  # find the slices for constraints
  fx = Fmat %*% x
  fv = Fmat %*% v
  
  U = sign(fx)*sqrt(fx^2 + fv^2)
  phi = atan(-fv/fx)
  phi_rev = atan(fv/fx)
  
  # find wall(s) that may be hit
  pn = abs(g/U) <= 1
  
  if (any(pn)){
    t0 = -phi[pn] + acos(-g[pn]/U[pn])
    
    t0_rev = -phi_rev[pn] + acos(-g[pn]/U[pn])
    t0_rev = 2*pi-t0_rev # because going the other way
    
    # create sets of intervals
    if (sum(pn)>1){
      theta_j1 = angles_to_intervals(t0, x = last_X, Fmat = Fmat, g = g)
      theta_j2 = angles_to_intervals(t0_rev, x = last_X, Fmat = Fmat, g = g)
      Theta_js = theta_j1 | theta_j2
    }
    else{
      Theta_js = angles_to_intervals(c(t0, t0_rev), x = last_X, Fmat = Fmat, g = g)
    }
  }
  else{
    # if no trajectories are hitting the walls,
    # the interval is the whole angle
    Theta_js = interval(0, 2*pi)
  }
  
  Slice = function(log_y){
    # find the slice for likelihood given y
    # and given Theta_js
    
    A = as.vector(t(x) %*% x)
    B = as.vector(t(x) %*% v)
    C = as.vector(t(v) %*% v)
    
    D = -(d/2)*log(2*pi) - log_y
    
    f = function(t){
      D - 0.5*(A*cos(t)^2 + 2*B*cos(t)*sin(t) + C*sin(t)^2)
    }
    
    roots = uniroot.all(f, c(0,2*pi))
    
    Theta_y = interval(0, 2*pi)
    
    if (!identical(roots, numeric(0))){
      # trim intervals if roots found
      l_0 = roots[c(TRUE, FALSE)]
      r_0 = roots[c(FALSE, TRUE)]
      for (i in 1:length(l_0)){
        Theta_y = interval_symdiff(Theta_y, interval(l_0[i], r_0[i]))
      }
    }
    
    Intv = Theta_js & Theta_y
    
    return(Intv)
  }
  
  return(Slice)
}


last_X = initial_X

# placeholder matrix for values:
Xs = last_X

iter_times = c()

set.seed(12345)

# begin loop:
while(ncol(Xs)<=N){
  ######### start time #########
  start.time = Sys.time()
  ##############################
  
  # sample ellipse from prior:
  v = t(rmvnorm(1, mean = rep(0, d)))
  
  # characterize slice:   
  # using last_X and v
  Slice = characterize_slice(last_X, v)
  
  
  # sample U to be recycled J times:
  U = runif(1, 0, 1)
  
  new_Xs = matrix(nrow = d, ncol = J)
  for (j in 1:J){
    # log-likelihood threshold:
    log_y = log(j-U) - log(J) + llik(last_X, Fmat = Fmat, g = g)
    
    # input log_y to Slice:
    Intv = Slice(log_y)
    
    if (interval_is_empty(Intv)){
      break()
    }
    
    # uniformly random sample from the interval:
    t_new = multRunif(Intv)
    x_prime = last_X*cos(t_new) + v*sin(t_new)
    
    new_Xs[,j] = x_prime
    
  }
  
  
  
  # check if satisfy constraints
  res = all(Fmat %*% new_Xs + matrix(g, nrow = length(g), ncol = J) > 0)
  if (!is.na(res) && res){
    if (J>1){
      # random shuffle new_Xs
      ind = sample(1:J, J)
      newXs_randperm = new_Xs[,ind]
      
      # append new X value
      Xs = cbind(Xs, newXs_randperm)
      # set the last element as last_X
      last_X = as.matrix(newXs_randperm[,J])
    }
    else{
      # append new X value
      Xs = cbind(Xs, x_prime)
      # set the last element as last_X
      last_X = as.matrix(x_prime)
    }
    
    if ((ncol(Xs)-1)%%100 == 0){
      cat("sample no:", ncol(Xs)-1,"\n")
    }
    
    ########## end time ##########
    end.time = Sys.time()
    iter_times = cbind(iter_times, as.vector(end.time - start.time))
    ##############################
  }
}

# remove first element
Xs = Xs[,-1]

# transform back
Xs = t(R)%*%Xs + matrix(mu, nrow = length(mu), ncol = N)

{
  plot(Xs[1,], Xs[3,], cex = 0.5,
       xlab = "X1", ylab = "X2",
       xlim = c(-1.5,1.5), ylim = c(-1.5,1.5),
       main = "Exact HMC")
  abline(v=-1)
  abline(v=1)
  abline(h=-1)
  abline(h=1)
}

plot(Xs[1,1:1000], type = "l", ylim = c(-2,2))

saveRDS(Xs, file = "../report/x_EESS_HD_J10.rds")
saveRDS(iter_times, file = "../report/t_EESS_HD_J10.rds")