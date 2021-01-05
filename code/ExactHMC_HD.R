# X: d x 1 vector
# F: m x d matrix
# g: m x 1 vector
rm(list=ls())

Fmat = diag(2)
m = nrow(Fmat)
d = ncol(Fmat)
g = matrix(rep(0, m), ncol = 1)
M = diag(d)
L = 1000

mu = matrix(rep(1.5,d), ncol = 1)

Sigma = matrix(rep(0.75, d*d), ncol = d, nrow = d)
diag(Sigma) = rep(1, d)

# M = solve(Sigma)
# r = M %*% mu

initial_X = matrix(rep(2, d), ncol = 1)

# set max travel time
max_T = pi/2


######## transform to whitened frame ########
# translate X and g by mu
initial_X = initial_X - mu
g = g + Fmat%*%mu
# transform by W
R = chol(Sigma) # R = upper triangle matrix = t(L), t(R) = L
Fmat = Fmat%*%t(R)
initial_X = solve(t(R), initial_X)


# verify that X is feasible
# c = (Fmat %*% initial_X) + g
# any(c<0)

nearzero = 1e-10


# assign last_X as the last X value (initial_X)
last_X = initial_X

# create placeholder matrix for Xs
Xs = matrix(initial_X, ncol = 1)

iter_times = c()

set.seed(12345)
# loop to sample X
while(ncol(Xs) < L){
  ######### start time #########
  start.time = Sys.time()
  ##############################
  
  # initialize initial velocity
  V0 = rnorm(d, 0, 1)
  
  # use the last X as initial X
  X = last_X
  
  # restart traveling time
  
  tt = 0
  
  # indicates which wall the particle hits
  # h = 0 -> no walls are hit
  h = 0
  
  stop_t = 0 # stopping condition for particle
  
  # loop to sample 1 value
  nbounce = 0
  while(1){
    # set values of initial conditions a and b
    # a and b can be results from a bounce
    a = V0
    b = X
    
    # necessary terms
    fa = Fmat%*%a
    fb = Fmat%*%b
    
    sign_fb = sign(fb)
    sign_fb[which(sign_fb == 0)] = 1
    
    U = sign_fb*sqrt(fa^2 + fb^2)
    phi = atan(-fa/fb)
    
    # find wall(s) that may be hit
    pn = abs(g/U) <= 1
    
    # find h, the first constraint that becomes zero
    if (any(pn)){
      inds = which(pn) # indices of constraint that wil be hit
      # find t which the particle hits the wall
      t_hit = -phi[pn] + acos(-g[pn]/U[pn])
      
      if (h>0){
        if (pn[h] == 1){
          cs = cumsum(pn)
          indj = cs[h]
          tt1 = t_hit[indj]
          if (is.nan(tt1)){
            t_hit[indj] = Inf
          }
        }
      }
      
      mt = min(t_hit) # (smallest) moving time
      h = inds[which.min(t_hit)] # which constraint is mt
    }
    else{
      mt = max_T
    }
    
    # update traveling time with moving time
    tt = tt + mt
    
    # check if total traveling time exceeds max_T
    if (tt >= max_T){
      mt = mt - (tt - max_T) # substract excess time from mt
      stop_t = 1 # change stop_t so no more iterations
    }
    
    # move particle an mt time
    X = a*sin(mt) + b*cos(mt)
    V = a*cos(mt) - b*sin(mt)
    
    if (stop_t){
      break() # break the loop if stop_t = 1
    }
    
    # if stop_t != 1, continue the trajectory
    # compute reflected velocity as new V0
    alpha_h = as.vector((Fmat[h,] %*% V)/sum(Fmat[h,]^2))
    F_h = matrix(Fmat[h,], ncol = 1)
    V0 = V - 2*alpha_h*F_h
    
    nbounce = nbounce + 1
  }
  
  # check if all constraints are satisfied
  if (all(Fmat%*%X + g > 0)){
    # if yes, keep current X and change the value of last_X
    Xs = cbind(Xs, X)
    last_X = X
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
Xs = t(R) %*% Xs + matrix(mu, nrow = length(mu), ncol = L)


plot(Xs[1,], Xs[2,], cex = 0.5,
     xlab = "X1", ylab = "X2",
     xlim = c(-1.5,5.5), ylim = c(-1.5,5.5),
     main = "Exact HMC")
abline(h=0)
abline(v=0)

plot(Xs[1,], type = "l", ylim = c(0,6))

# saveRDS(Xs, file = "../report/x_EHMC_HD.rds")
# saveRDS(iter_times, file = "../report/t_EHMC_HD.rds")