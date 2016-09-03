#==============================================
# Simulation data functions 
#==============================================
# univariate
makeData1 <- function(n){
  x1 <- runif(n,-4,4)
  y <- x1/20 - 0.42*x1^2 + rnorm(n)
  return(data.frame(x1=x1,y=y))
}

# jumps
makeData2 <- function(n){
  x1 <- runif(n,-4,4)
  y <- -2.7*as.numeric(x1< -3) + 2.5*as.numeric(x1 > -2) - 
    2*as.numeric(x1>0) + 4*as.numeric(x1>2) - 3*as.numeric(x1>3) + rnorm(n)
  return(data.frame(x1=x1,y=y))
}

# sins
makeData3 <- function(n){
  x1 <- runif(n,-4,4)
  y <- 2*sin(pi/2*abs(x1)) + 2*cos(pi/2*abs(x1)) + rnorm(n)
  
  return(data.frame(x1=x1,y=y))
}

# trivariate 
# smooth 
makeData4 <- function(n){
  x1 <- runif(n,-4,4)
  x2 <- runif(n,-4,4)
  x3 <- rbinom(n, 1, 0.5)
  y <- x1/15 - 0.28*x1^2 + 0.5*x2 + 0.25*x3*x2 + rnorm(n)
  return(data.frame(x1=x1,x2=x2,x3=x3,y=y))
}

# jumps + interactions
makeData5 <- function(n){
  x1 <- runif(n,-4,4)
  x2 <- runif(n,-4,4)
  x3 <- rbinom(n, 1, 0.5)
  y <- -2*as.numeric(x1< -3)*x3 + 2.5*as.numeric(x1 > -2) - 
    2*as.numeric(x1>0) + 2.5*as.numeric(x1>2)*x3 - 2.5*as.numeric(x1>3) + 
    1*as.numeric(x2 > -1) - 
    4*as.numeric(x2>1)*x3 + 2*as.numeric(x2>3) + rnorm(n)
  return(data.frame(x1=x1,x2=x2,x3=x3,y=y))
}

# really hard
makeData6 <- function(n){
  x1 <- runif(n,-4,4)
  x2 <- runif(n,-4,4) 
  x3 <- rbinom(n, 1, 0.5)
  y <- x3*(4*sin(pi/2*abs(x1))*(x2 < 0) + 4.1*cos(pi/2*abs(x1))*(x2 > 0)) + rnorm(n)
  return(data.frame(x1=x1,x2=x2,x3=x3,y=y))
}

### five-variate
# smooth
makeData7 <- function(n){
  x1 <- runif(n,-4,4)
  x2 <- runif(n,-4,4)
  x3 <- rbinom(n, 1, 0.5)
  x4 <- rnorm(n)
  x5 <- rgamma(n, 2, 1)
  
  y <- x1/10 - 0.3*x1^2 + 0.25*x2 + 0.5*x3*x2 - 0.5*x4 + 0.2*x5^2/5 - 0.1*x5 + rnorm(n)
  return(data.frame(x1=x1,x2=x2,x3=x3,x4=x4,x5=x5,y=y))
}

# jumps + interactions
makeData8 <- function(n){
  x1 <- runif(n,-4,4)
  x2 <- runif(n,-4,4)
  x3 <- rbinom(n, 1, 0.5)
  x4 <- rnorm(n)
  x5 <- rgamma(n, 2, 1)
  
  y <- -1*as.numeric(x1< -3)*x3 + 0.5*as.numeric(x1 > -2) - 
    1*as.numeric(x1>0) + 2*as.numeric(x1>2)*x3 - 3*as.numeric(x1>3) + 
    1.5*as.numeric(x2 > -1) - 
    5*as.numeric(x2>1)*x3 + 2*as.numeric(x2>3) + as.numeric(x4 < 0)*2 - 
    as.numeric(x5 > 5)*1 - as.numeric(x4<0)*as.numeric(x1<0) + 2*x3 + 
    rnorm(n)
  return(data.frame(x1=x1,x2=x2,x3=x3,x4=x4,x5=x5,y=y))
}

# really hard
makeData9 <- function(n){
  x1 <- runif(n,-4,4)
  x2 <- runif(n,-4,4) 
  x3 <- rbinom(n, 1, 0.5)
  x4 <- rnorm(n)
  x5 <- rgamma(n, 2, 1)
  
  y <- x3*(3.75*sin(pi/2*abs(x1))*(x2 < 0) + 4*cos(pi/2*abs(x1))*(x2 > 0)) + 
    sin(x4*pi)*x5/10 + cos(abs(x4-x5))*x3 + rnorm(n)
  return(data.frame(x1=x1,x2=x2,x3=x3,x4=x4,x5=x5,y=y))
}

### seven variate
# smooth
makeData10 <- function(n){
  x1 <- runif(n,-4,4)
  x2 <- runif(n,-4,4)
  x3 <- rbinom(n, 1, 0.5)
  x4 <- rnorm(n)
  x5 <- rgamma(n, 2, 1)
  x6 <- rbinom(n, 1, 0.25)
  x7 <- rpois(n, 4)
  
  y <- 0.5*x1 - 0.24*x1^2 + 0.25*x2 + 0.25*x3*x2 - 0.25*x4 + 0.05*x5^2 - 0.25*x5 + 0.05*x6*x3*2 + 0.05*x7 + 0.02*x7*x1 + rnorm(n)
  return(data.frame(x1=x1,x2=x2,x3=x3,x4=x4,x5=x5,x6=x6, x7=x7,y=y))
}

# jumps + interactions
makeData11 <- function(n){
  x1 <- runif(n,-4,4)
  x2 <- runif(n,-4,4)
  x3 <- rbinom(n, 1, 0.5)
  x4 <- rnorm(n)
  x5 <- rgamma(n, 2, 1)
  x6 <- rbinom(n, 1, 0.25)
  x7 <- rpois(n, 4)
  
  y <- -1*as.numeric(x1< -3)*x3 + 0.5*as.numeric(x1 > -2) - 
    0.5*as.numeric(x1>0) + 2*as.numeric(x1>2)*x3 - 0.5*as.numeric(x1>3) + 
    1.5*as.numeric(x2 > -1) - 
    1*as.numeric(x2>1)*x3 + 1*as.numeric(x2>3)  + 
    as.numeric(x4 < 0)*2 - as.numeric(x5 > 5)*1 - 
    as.numeric(x4<0)*as.numeric(x1<0) + 1.75*x3 + 
    2*as.numeric(x7 > 1) - 0.5*x6 + as.numeric(x7 > x5)*2 + rnorm(n)
  return(data.frame(x1=x1,x2=x2,x3=x3,x4=x4,x5=x5,x6=x6, x7=x7,y=y))
}

# really hard
makeData12 <- function(n){
  x1 <- runif(n,-4,4)
  x2 <- runif(n,-4,4) 
  x3 <- rbinom(n, 1, 0.5)
  x4 <- rnorm(n)
  x5 <- rgamma(n, 2, 1)
  x6 <- rbinom(n, 1, 0.25)
  x7 <- rpois(n, 4)
  
  y <- x3*(2.75*sin(pi/2*abs(x1))*(x2<0) + 4*cos(pi/2*abs(x1))*(x2 > 0)) + 
    sin(x4*pi)*x5/10 + cos(abs(x4-x5))*x3 + 
    x7*sin(abs(x5)*x6)/5 + x5*cos(abs(x7)*x6)*x3/1.5 + rnorm(n)
  return(data.frame(x1=x1,x2=x2,x3=x3,x4=x4,x5=x5,x6=x6, x7=x7,y=y))
}
