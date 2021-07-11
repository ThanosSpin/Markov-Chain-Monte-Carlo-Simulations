#### Import Libraries ###### 

library(svMisc)# progress for loops
library(doParallel)# parallel processing
library(zoo)#Rolling means (moving average)
library(magicfor)

########## Set parallel computations ################
getwd()

print("----- Parallel processing -----")
source("./parallel_processing.R")


options(digits = 10) # control number of digits in numeric format

########## Monte Carlo (k dimensions) ##########

###MC function k - dimensions (without plots)###
mc.d_wp <-function(n, d, k)
{  
  accept <- 0
  
  x <- (2*runif(k, 0, d) - d)
  data <- c(x)
  accepts <- matrix(data, 1, k)
  for (i in 1:n) {
    x <- t(matrix(2*runif(k, 0, d) - d)) 
    if(sum(x^2) <= 1) { 
      accepts <- rbind(accepts, c(x))
      accept <- accept + 1
    }
    else {
      x <- x
    }
  }
  alpha <- accept/n
  return(alpha)
}

###Running MC for d (dimensions)####
t0 <- Sys.time()
set.seed(32)
matrix_mc <- c()  
for (i in 1:50) {
 matrix_mc[i] <- matrix(i, ncol = 1)
 matrix_mc[i] <- mc.d_wp(500000, 1, i)
}

l <- seq(from = 75, to = 200, by = 25)

set.seed(33)
for (i in 1:length(l)) {
  x <- mc.d_wp(500000, 1, l[i])
  matrix_mc <- c(matrix_mc, x)
}

t1 <- Sys.time()
t1 - t0

matrix_mc <- data.frame(matrix_mc)
colnames(matrix_mc)[1] <- c("Prob")

###B sphere area (true volume)####
B_true <- c()
for (i in 1:50) {
  B_true[i] <- matrix(i, ncol = 1)
  (B_true[i] <- (pi^(i/2))/gamma(i/2 + 1))
}

for (i in 1:length(l)) {
  x <- (pi^(l[i]/2))/gamma(l[i]/2 + 1)
  B_true <- c(B_true, x)
}



###Join the two columns (B_true, prob)####
matrix_mc$B_True <- B_true
matrix_mc$Dimensions <- c(seq(1, 50, by = 1), l)
benchmark_matrix <- matrix_mc[ , c(3, 1:2)]
# names(benchmark_matrix)[which(names(benchmark_matrix) %in% c("B_MC"))] <- c("B_True")

####Estimate the integrate of B (d dimensions)####
benchmark_matrix$B_MC <- benchmark_matrix$Prob*(2^benchmark_matrix$Dimensions)


###The relative error of estimation of area B (d dimensions)####
(benchmark_matrix$Rel_Err_MC <- round(abs(1 - (benchmark_matrix$B_MC/benchmark_matrix$B_True)), 4)*100)

###Plot acceptance ratio of MC####
plot(benchmark_matrix$Dimensions[1:20], matrix_mc$Prob[1:20], type = "l", col = "red", ylab = "Probability", xlab = "Dimensions", main = "Monte Carlo Estimator", xaxt = "n")
axis(1, at = c(seq(0, 20, by = 1)), labels = c(seq(0, 20, by = 1)))


#### Monte Carlo Markov Chain (d dimensions) ####

###Initial point function (run the chain with k - 1 dimensions for 5.000 steps)
initial_point <- function(k) {
  y <- mh.d(c(rep(0, k - 1)), 5000, k)
  return(y)}



####MCMC function####
mh.d<-function(starting_point, n, k) 
{
  accept <- 0
  x <- starting_point
  for (i in 2:n) {
      can <- (runif(k - 1, -1/(k - 1), 1/(k - 1))) + x
    sphere <- sum(can^2)
    if (sphere <= 1) { 
        x <- can 
        # vec <- x
    }
  }
  last_step <- c((runif(1, -1, 1)))
  if (sum(last_step^2, x^2) < 1) {
    accept <- accept + 1
  }
  alpha <- accept
  return(list(x, alpha, mean(x)))
}



####MCMC function with storing the simulated values of Uniform Distribution####
mh.d_dim <-function(starting_point, n, k) 
{
  accept <- 0

  x <- starting_point
  vec <- c(x)
  vec <- data.frame(t(vec))
  for (i in 2:n) {
    can <- (runif(k - 1, -1/(k - 1), 1/(k - 1))) + x
    sphere <- sum(can^2)
    if (sphere <= 1) { 
      x <- can
      vec[i, ] <- x   
    }
    else {
      vec[i, ] <- x
    }
    }
    last_step <- c(runif(1, -1,  1))
    if ((sum(last_step^2, vec$X1[n]^2)) < 1) {
    accept <- accept + 1
    
    alpha <- accept
  }
  return(list(vec, alpha, mean(vec$X1)))
}



### Set the initial point ofthe markov chain ####
set.seed(1234)

  
### Iterating MCMC process - Initial points of the chain the ending points of the previous (markov chain) #### 

start.time <-  Sys.time()
alpha_mcmc <- c()
alpha_mcmc <- list(mh.d(c(initial_point(10)[[1]]), 10000, 10))
for (i in 2:1000) {
  alpha_mcmc[i] <- list(mh.d(c(alpha_mcmc[[i - 1]][[1]]), 10000, 10))
}
end.time <- Sys.time()


alpha_prob <- c()
for (i in 1:length(alpha_mcmc)) {alpha_prob[i] <- ((alpha_mcmc[[i]][[2]]))}
sum(alpha_prob)/length(alpha_mcmc)


###Running MCMC for k (dimensions)####

  matrix_mcmc <- c()  
  probability <- c()
  alpha <- c()
  
  
t2 <- Sys.time()#current time before running


####dimesions - k = 2,...,10#####
set.seed(152)
for (k in 2:10) {
    alpha[1] <- list(mh.d(c(initial_point(k)[[1]]), 1000, k))
    probability[1] <- alpha[[1]][[2]]
    for (i in 2:10000) {
      progress(i)
      alpha[i] <- list(mh.d(c(alpha[[i - 1]][[1]]), 10000, k))
      probability[i] <- alpha[[i - 1]][[2]]
      }
      matrix_mcmc[k - 1] <- sum(as.numeric(probability))/i
      if (k == k[length(k)]) cat("Done!\n")  
  }
  
t3 <- Sys.time()
time.taken <- t3 - t2
  
  
mean_vl <- c()


t2 <- Sys.time()#current time before running


####dimesions - k = 11,...,20#####
set.seed(155)
  for (k in 11:20) {
    alpha[1] <- list(mh.d(c(initial_point(k)[[1]]), 1000, k))
    probability[1] <- alpha[[1]][[2]]
    for (i in 2:10000) {
      progress(i)
      alpha[i] <- list(mh.d(c(alpha[[i - 1]][[1]]), 10000, k))
      probability[i] <- alpha[[i - 1]][[2]]
      mean_vl[i] <- alpha[[i - 1]][[3]]
    }
    matrix_mcmc[k - 1] <- sum(as.numeric(probability))/i
    if (k == k[length(k)]) cat("Done!\n")  
  }

t3 <- Sys.time()
time.taken <- t3 - t2


#####Volume of shpere (dimensions- k = 2,..,20)####
for (i in 2:20) {
print(2*matrix_mcmc[i - 1]*(benchmark_matrix$B_True[i - 1]))}#for start we calculate the volume of next dimension(k)
#depending on the real volume of sphere of previous dimension (k - 1) in order to evalueate the error of estimation
#of expected mean value


t4 <- Sys.time()

####dimesions - k = 21,...,30#####
set.seed(165)
for (k in 21:30) {
  alpha[1] <- list(mh.d(c(initial_point(k)[[1]]), 1000, k))
  probability[1] <- alpha[[1]][[2]]
  for (i in 2:10000) {
    progress(i)
    alpha[i] <- list(mh.d(c(alpha[[i - 1]][[1]]), 10000, k))
    probability[i] <- alpha[[i - 1]][[2]]
  }
  matrix_mcmc[k - 1] <- sum(as.numeric(probability))/i
  if (k == k[length(k)]) cat("Done!\n")  
}

t5 <- Sys.time()
time.taken <- t5 - t4


t6 <- Sys.time()

####dimesions - k = 31,...,40#####
set.seed(166)

for (k in 31:40) {
  alpha[1] <- list(mh.d(c(initial_point(k)[[1]]), 1000, k))
  probability[1] <- alpha[[1]][[2]]
  for (i in 2:10000) {
    progress(i)
    alpha[i] <- list(mh.d(c(alpha[[i - 1]][[1]]), 10000, k))
    probability[i] <- alpha[[i - 1]][[2]]
  }
  matrix_mcmc[k - 1] <- sum(as.numeric(probability))/i
  if (k == k[length(k)]) cat("Done!\n")  
}


t7 <- Sys.time()
time.taken <- t7 -t6


t8 <- Sys.time()

####dimesions - k = 41,...,50#####
set.seed(168)

for (k in 41:50) {
  alpha[1] <- list(mh.d(c(initial_point(k)[[1]]), 1000, k))
  probability[1] <- alpha[[1]][[2]]
  for (i in 2:100000) {
    progress(i)
    alpha[i] <- list(mh.d(c(alpha[[i - 1]][[1]]), 1000, k))
    probability[i] <- alpha[[i - 1]][[2]]
  }
  matrix_mcmc[k - 1] <- sum(as.numeric(probability))/i
  if (k == k[length(k)]) cat("Done!\n")  
}

t9 <- Sys.time()
time.taken <- t9 - t8



for (i in 2:50) {
print(2*matrix_mcmc[i - 1]*(benchmark_matrix$B_True[i-1]))}#for start we calculate the volume of next dimension(k)
#depending on the real volume of sphere of previous dimension (k - 1) in order to evalueate the error of estimation
#of expected mean value
                                                    

# matrix_mcmc <-  c(0.78650, 0.65860, 0.59080, 0.52730, 0.48340, 0.45580, 0.43770, 0.41150, 0.39110, 0.37030,
# 0.35040, 0.33520, 0.32180, 0.31760, 0.30430, 0.30860, 0.29780, 0.29420, 0.28080, 0.26893,
# 0.26769, 0.25870, 0.25309, 0.24922, 0.24738, 0.23904, 0.23692, 0.23040, 0.22846, 0.22443,
# 0.21797, 0.21855, 0.21214, 0.21067, 0.20690, 0.20537, 0.20419, 0.20030, 0.19875, 0.19740,
# 0.19220, 0.18380, 0.19410, 0.18750, 0.17660, 0.19200, 0.18550, 0.18370, 0.17920) 
# 
# matrix_mcmc <- c(matrix_mcmc, 0.12600, 0.12800, 0.10300, 0.09500, 0.09400, 0.07400)


alpha <- c()
probability <- c()


t12 <- Sys.time()

####dimesions - k = 75,...,200#####
set.seed(244)

for (k in l) {
  alpha[1] <- list(mh.d(c(initial_point(k)[[1]]), 1000, k))
  probability[1] <- alpha[[1]][[2]]
  for (i in 2:10000) {
    progress(i)
    alpha[i] <- list(mh.d(c(alpha[[i - 1]][[1]]), 10000, k))
    probability[i] <- alpha[[i - 1]][[2]]
  }
  matrix_mcmc[k - 1] <- sum(as.numeric(probability))/i
  if (k == k[length(k)]) cat("Done!\n")  
}

t13 <- Sys.time()
time.taken <- t13 - t12


# benchmark_matrix <- benchmark_matrix[-c(51), ] ## removing a redundant row of the data frame

###mcmc function:: If we divide pi with number of dimensions (k) for the generation of the values (uniform interval), 
###as we increase dimensions, the probability remains stable despite the increase of the dimensions

matrix_mcmc <- matrix_mcmc[complete.cases(matrix_mcmc)]
matrix_mcmc <- data.frame(matrix_mcmc)
colnames(matrix_mcmc)[1] <- c("Prob_MCMC")

B_estimated_MCMC <- c()
B_estimated_MCMC[1] <- 2
# B_estimated_MCMC <- c(B_estimated_MCMC, 2*matrix_mcmc$Prob_MCMC[1]*(B_estimated_MCMC))

# names(matrix_mcmc)

for (i in 2:50) {
B_estimated_MCMC[i] <- 2*matrix_mcmc$Prob_MCMC[i - 1]*(B_estimated_MCMC[i - 1])
}


####Volume of k - 1 hyper - sphere (for MCMC estimations over 50 dimensions####

B_MCMC <- c()

for (i in 1:length(l)) {
 B_MCMC[i] <- ((pi^((l - 1)[i]/2))/gamma(((l - 1)[i]/2) + 1))
}
 
 for (j in 51:56) {
    B_estimated_MCMC[j] <- 2*matrix_mcmc$Prob_MCMC[j - 1]*(B_MCMC[j - 50])
 }

benchmark_matrix$B_MCMC <- B_estimated_MCMC


###Join the two methods (MC, MCMC)####
# benchmark_matrix$Prob_MCMC <- matrix_mcmc$Prob_MCMC
# names(benchmark_matrix)[which(names(benchmark_matrix) %in% c("Prob"))] <- c("Prob_MC")



####Estimate the integrate of B (d dimensions, MCMC)####
# B_estimated_MCMC <- 2*benchmark_matrix$Prob_MCMC
# for (i in 2:nrow(benchmark_matrix)) {
#   B_estimated_MCMC[i] <- 2*benchmark_matrix$Prob_MCMC[i]*(benchmark_matrix$B_estimated_MCMC[i-1])
# }
# 


###The relative error of estimation of area B (d dimensions, MCMC)####
(benchmark_matrix$Rel_Err_MCMC <- round(abs(1 - (benchmark_matrix$B_MCMC/benchmark_matrix$B_True)), 4)*100)
for (i in 1:nrow(benchmark_matrix)) {
if (benchmark_matrix$Rel_Err_MCMC[i] > 100) {
  benchmark_matrix$Rel_Err_MCMC[i] <- 100
}
  else {benchmark_matrix$Rel_Err_MCMC[i]}
}

#### Add MCMC estimator in benchmark matrix####
benchmark_matrix$Prob_MCMC <- c(1, matrix_mcmc$Prob_MCMC)

###Sumamry statistics for benchmark matrix####
summary(benchmark_matrix)
benchmark_matrix <- benchmark_matrix[c(1:5, 8, 6:7)]

  
###Plot probability of MCMC####
plot(benchmark_matrix$Dimensions, benchmark_matrix$Prob_MCMC, type = "l", col = "red", ylab = "Probability", xlab = "Dimensions", main = "MCMC estimator - MH Algorithm", xaxt = "n")
axis(1, at = c(seq(0, 200, by = 10)), labels = c(seq(0, 200, by = 10)))

benchmark_matrix$B_True[20]##true value of integral calculation in 20 dimensions 
benchmark_matrix$B_MCMC[20]##estimated value of integral calculation in 20 dmensions through MCMC - MH method
round(abs(benchmark_matrix$B_True[20]/benchmark_matrix$B_MCMC[20] - 1)*100, 2)

##We observe that the error is below 1.5%. 
##We can estimate with high accuracy the integrals of a given function in high dimensional space (until 200 dimensions - relative error MCMC below 10% in benchmark matrix)
##by simulating values from uniform random variables Ui, i = 1,..,n with the appropriate delta 
##MC method fails in estimating integrals after 17 dimensions (relative error MC equals to 100% in benchmark matrix). 


###Optimal value of delta (step of each simulation)#### 
d <- seq(0.1, 5, length.out = 50)

source("mh_d_opt_delta.R")

mh.d_opt_delta(c(initial_point(10)[[1]]), 1000, 10, d[3])


delta <- c(); accept <- c(); mh <- c();

set.seed(157)
for (k in 11:20) {
  for (j in d) {
  progress(j)
  mh <- mh.d_opt_delta(c(initial_point(k)[[1]]), 20000, k, j)
  delta <- c(delta, mh[[1]])
  accept <- c(accept, mh[[2]])
   }
  if (k == k[length(k)]) cat("Done!\n")  
}


delta[seq(3, 50, by = 10)]
delta[seq(2, 500, by = 50)]
delta[seq(5, 500, by = 50)]
1/c(10:19)
benchmark_matrix$Prob_MCMC[11:20]


plot(delta, type = "l", col = "red", xlab = "", xaxt = "n")
axis(1, at = c(seq(0, 500, by = 50)), labels = c(seq(0, 500, by = 50)))

plot(accept, type = "l",  xlab = "", xaxt = "n", ylab = "MCMC Estimator")
axis(1, at = c(seq(0, 500, by = 50)), labels = c(seq(0, 500, by = 50)))

plot(seq(0.1, 5, by = 0.1), delta[1:50], type = "l", col = "red", ylab = "(%)", xlab = "delta", ylim = c(0,0.8), xaxt = "n")
axis(1, at = c(seq(0, 5, by = 0.1)), labels = c(seq(0, 5, by = 0.1)))

plot(seq(0.1, 5, by = 0.1), accept[1:50], type = "l", col = "red", ylab = "MCMC estimator", xlab = "", xaxt = "n")
axis(1, at = c(seq(0, 5, by = 0.1)), labels = c(seq(0, 5, by = 0.1)))


d[1]
delta[seq(1, 500, by = 50)]
accept[seq(1, 500, by = 50)]
benchmark_matrix$Prob_MCMC[11:20]




##### Comparison of the running means of Uniform distribution used for sphere volume ######

#### Assign the values of previous run of the Markov Chain (k = 4) 
#### (Burn In period included)


####dimension -> k = 4#####
set.seed(50)
for (k in 4) {
  alpha <- (mh.d_dim(c(initial_point(k)[[1]]), 1000, k))
  # alpha <- list(alpha[[1]][1000, ], alpha[[2]], alpha[[3]])
  probability[1] <- (alpha[[2]])
  for (i in 2:1000) {
    progress(i)
    alpha <- (mh.d_dim(c(t(alpha[[1]][1000, ])), 1000, k))
    probability <- c(probability, (t(alpha[[2]])))
    mean_vl <- c(probability, (t(alpha[[3]])))
  }
  matrix_mcmc[k - 1] <- sum(as.numeric(probability))/i
  if (k == k[length(k)]) cat("Done!\n")  
}

run_mean <- c()#create an empty vector

for (i in 1:10000) {run_mean[i] <- as.data.frame(t(alpha[[i]][[3]]))}#assign the values
#from the repeated process to the run_mean vector (dimension -> k = 4)

run_mean <- t(run_mean)#transpose vector run_mean
run_mean <- as.data.frame(apply(run_mean, 2, as.numeric))#transform run_mean to a data frame and tranform the values to numeric 
names(run_mean) <- "mean" #set the name mean to run_mean data frame 

plot(ts(run_mean$mean), ylab = 'Running Mean')#plot variable mean as a timeseries  

hist(run_mean$mean, probability = TRUE, main = "Mean of 10.000 iterations", ylim = c(0, 1.5), xlab = "")#histogram of variable mean

x2 <- seq(min(c(-1, run_mean$mean)), max(run_mean$mean), length.out=100)#generate 100 values with min and max value of run_mean variable
lines(x2, dnorm(x2, mean(run_mean$mean), sd(run_mean$mean)), lty = 2, col = 2)#adding a normality line to plot of variable mean

alpha <- c()# empty alpha vector  

##### Running the Markov Chain without Burn - In period #####
set.seed(65)
for (k in 4) {
  alpha[1] <- list(mh.d_dim(runif(k - 1, -1/(k - 1), 1/(k - 1)), 100, k))
  probability[1] <- alpha[[1]][[2]]
  for (i in 2:10000) {
    progress(i)
    alpha[i] <- list(mh.d(c(alpha[[i - 1]][[1]]), 1000, k))
    probability[i] <- alpha[[i - 1]][[2]]
    mean_vl[i] <- alpha[[i - 1]][[3]]
  }
  matrix_mcmc[k - 1] <- sum(as.numeric(probability))/i
  if (k == k[length(k)]) cat("Done!\n") 
}


#### Running Average for MCMC (dimension -> k = 3)####
run_mean_wburn <- c()

for (i in 1:10000) {run_mean_wburn[i] <- data.frame(t(alpha[[i]][[3]]))}#assign the values
#from the repeated process to the test vector (dimension -> k = 3)

run_mean_wburn <- t(run_mean_wburn)#transpose vector run_mean_wburn
run_mean_wburn <- as.data.frame(apply(run_mean_wburn, 2, as.numeric))#transform run_mean_wburn to a data frame and tranform the values to numeric 
names(run_mean_wburn) <- "mean" #set the name mean to run_mean_wburn data frame 

plot(ts(run_mean_wburn)) #plot variable run_mean_wburn as a timeseries 

hist(run_mean_wburn[ , 1], probability = TRUE, ylim = c(0, 1.5))#histogram of variable mean

x2 <- seq(min(c(-1, run_mean_wburn[ , 1])), max(run_mean_wburn[, 1]), length.out=100)#generate 100 values with min and max value of test variable
lines(x2, dnorm(x2, mean(run_mean_wburn[ , 1]), sd(run_mean_wburn[ , 1])), lty = 2, col =2)#adding a normality line in plot of variable mean

par(mfrow = c(1, 2)) #two graphs per window
run_mean <-  cumsum(run_mean)/seq(along = run_mean)#running mean of variable run_mean
run_mean_wburn <- cumsum(run_mean_wburn)/seq(along = run_mean_wburn)#running mean of variable run_mean_wburn

####Plot the graphs of running means with and without burn - in period ####
plot(ts(run_mean), xlim=c(0, 10000), ylim = c(-8, 6), xaxt='n', yaxt = 'n', ylab = '', main = "Running mean of MCMC with burn - in period of 5.000 steps")#plot the running mean of run_mean variable
axis(side = 1, at = seq(0, 10000, by = 1000))#set in axis - x the desired scale
axis(side = 2, at = seq(-8, 6, by = 2))

plot(ts(run_mean_wburn), xlim=c(0, 10000), xaxt='n', yaxt = 'n',ylab = '', main = "Running mean of MCMC without burn - in period")#plot the running mean of run_mean_wburn variable
axis(side = 1, at = seq(0, 10000, by = 1000))#set in axis - x the desired scale
axis(side = 2, at = seq(-8, 6, by = 2))

rm(list = ls(pattern = "^run_")) # remove test objects from environment



####MCMC function of hypercylinder (k - dimensions) with height {-1, 1}####

###MCMC function - hypercylinder (map one vector from [-1, 1] to create the hypersphere)###
mh.cylinder <- function(starting_point, n, k) 
{
    accept <- 0
    x <- starting_point
    for (i in 2:n) {
      can <- (runif(k - 1, -1/(k - 1), 1/(k - 1))) + x
      sphere <- sum(can^2)
      if (sphere <= 1) { 
        x <- can 
      }
    last_step <- c((runif(1, -1, 1)))
    if (sum(last_step^2, x^2) < 1) {
      accept <- accept + 1
      }
    }
    alpha <- accept
    return(list(x, alpha, mean(x)))
}


alpha_cyl <- c()
probability_cyl <- c()

###Running MCMC cyl for k (dimensions)####
t4 <- Sys.time()
matrix_mcmc_cyl <- c()  

set.seed(444)
for (k in 2:50) {
  alpha_cyl[1] <- list(mh.cylinder(c(initial_point(k)[[1]]), 1000, k))
  probability_cyl[1] <- alpha_cyl[[1]][[2]]/1000
  for (i in 2:10000) {
    progress(i)
    alpha_cyl[i] <- list(mh.cylinder(c(alpha_cyl[[i - 1]][[1]]), 10000, k))
    probability_cyl[i] <- alpha_cyl[[i - 1]][[2]]/10000
  }
  matrix_mcmc_cyl[k - 1] <- sum(as.numeric(probability_cyl))/i
  if (k == k[length(k)]) cat("Done!\n")  
}

t5 <- Sys.time()
t5 - t4


alpha_cyl <- c()
probability_cyl <- c()


t6 <- Sys.time()

set.seed(445)
for (k in l) {
  alpha_cyl[1] <- list(mh.cylinder(c(initial_point(k)[[1]]), 1000, k))
  probability_cyl[1] <- alpha_cyl[[1]][[2]]/1000
  for (i in 2:10000) {
    progress(i)
    alpha_cyl[i] <- list(mh.cylinder(c(alpha_cyl[[i - 1]][[1]]), 10000, k))
    probability_cyl[i] <- alpha_cyl[[i - 1]][[2]]/10000
  }
  matrix_mcmc_cyl[k - 1] <- sum(as.numeric(probability_cyl))/i
  if (k == k[length(k)]) cat("Done!\n")  
}


t7 <- Sys.time()
t7 - t6

matrix_mcmc_cyl <- c(1, matrix_mcmc_cyl)
matrix_mcmc_cyl_df <- data.frame(matrix_mcmc_cyl)
colnames(matrix_mcmc_cyl_df) <- c("Prob_MCMC_cyl")

matrix_mcmc_cyl_df <- matrix_mcmc_cyl_df$Prob_MCMC_cyl[complete.cases(matrix_mcmc_cyl_df$Prob_MCMC_cyl)]


###Join the two methods (MC, MCMC)####
benchmark_matrix$Prob_MCMC_cyl <- matrix_mcmc_cyl_df


####Estimate the integrate of B (d dimensions, MCMC)####
B_estimated_MCMC_cyl <- c()
B_estimated_MCMC_cyl[1] <- 2
# B_estimated_MCMC <- c(B_estimated_MCMC, 2*matrix_mcmc$Prob_MCMC[1]*(B_estimated_MCMC))

# names(matrix_mcmc)

for (i in 2:50) {
  B_estimated_MCMC_cyl[i] <- 2*matrix_mcmc_cyl_df[i]*(B_estimated_MCMC_cyl[i - 1])
}


####Volume of k - dimensions hyper - sphere (for MCMC estimations over 50 dimensions####

B_MCMC_cyl <- c()

for (i in 1:length(l)) {
  B_MCMC_cyl[i] <- ((pi^((l - 1)[i]/2))/gamma(((l - 1)[i]/2) + 1))
}

for (j in 51:56) {
  B_estimated_MCMC_cyl[j] <- 2*matrix_mcmc_cyl_df[j]*(B_MCMC[j - 50])
}

benchmark_matrix$B_MCMC_et <- B_estimated_MCMC_cyl


###The relative error of estimation of area B (k dimensions, MCMC) - 3rd method####
(benchmark_matrix$Rel_Err_MCMC_et <- round(abs(1 - (benchmark_matrix$B_MCMC_et/benchmark_matrix$B_True)), 4)*100)
for (i in 1:nrow(benchmark_matrix)) {
  if (benchmark_matrix$Rel_Err_MCMC_et[i] > 100) {
    benchmark_matrix$Rel_Err_MCMC_et[i] <- 100
  }
  else {benchmark_matrix$Rel_Err_MCMC_et[i]}
}


names(benchmark_matrix)[which(names(benchmark_matrix) %in% c("Prob_MCMC_cyl"))] <- "Prob_MCMC_et"

###Plot of MCMC - 3rd method####
plot(benchmark_matrix$Dimensions, matrix_mcmc_cyl_df, type = "l", col = "red", ylab = "Probability (alpha)", xlab = "Dimensions", main = "Monte Carlo Markov Chain (cylinder)", xaxt = "n")
axis(1, at = c(seq(0, 200, by = 10)), labels = c(seq(0, 200, by = 10)))

benchmark_matrix$B_True[20]##true value of integral calculation in 20 dmensions through MCMC method
benchmark_matrix$B_MC[20]##MC estimated value of integral calculation in 20 dimensions 
benchmark_matrix$B_MCMC[20]##MCMC estimated value of integral calculation in 20 dimensions 
benchmark_matrix$B_MCMC_et[20]##MCMC (ET) value of integral calculation in 20 dimensions 
### The better volume estimation in k=20 dimensions  is with 3rd method, stationary distribution (Ergodic Theorem) ###

###Plot relative error (%) of volumes estimation with three methods####
plot(benchmark_matrix$Dimensions, benchmark_matrix$Rel_Err_MC, type = "l", col = "red", ylab = "%", xlab = "Dimensions", main = "Relative Error of Volumes Estimation, MC vs MCMC (CLT vs Ergodicity)", xaxt = "n", yaxt = "n")
axis(1, at = c(seq(0, 200, by = 10)), labels = c(seq(0, 200, by = 10)))
axis(2, at = c(seq(0, 100, by = 10)), labels = c(seq(0, 100, by = 10)))
lines(benchmark_matrix$Dimensions, benchmark_matrix$Rel_Err_MCMC, col = "green")
lines(benchmark_matrix$Dimensions, benchmark_matrix$Rel_Err_MCMC_et, col = "blue")
legend("topright",  c("Rel_Err_MC_CLT", "Rel_Err_MCMC_CLT", "Rel_Err_MCMC_ET"), cex = 0.7, lty = c(1, 1), lwd = c(1.5, 1.5), col = c("red", "green", "blue"))

