
makeData <- function(n = n){
    L0 <- data.frame(x.0 = runif(n,-1,1))
    A0 <- rbinom(n, 1, plogis(L0$x.0^2))
    L1 <- data.frame(x.1 = L0$x.0^2*A0 + runif(n))
    A1 <- rbinom(n, 1, plogis(L0$x.0*L1$x.1))
    L2 <- rnorm(n, L0$x.0^2*A0*A1 + L1$x.1)
    return(list(L0 = L0, L1 = L1, L2 = L2, A0 = A0, A1 = A1))
}

getTruth <- function(n = 1e6, abar = c(1,1)){
    L0 <- data.frame(x.0 = runif(n,-1,1))
    L1 <- data.frame(x.1 = L0$x.0^2*abar[1] + runif(n))
    L2 <- rnorm(n, L0$x.0^2*abar[1]*abar[2] + L1$x.1)
    return(mean(L2))
}