SL.glm3 <- function(Y, X, newX, family, obsWeights, ...){
    n <- length(Y)
    Xmat <- cbind(rep(1,n), X[,1], X[,1]^2, X[,1]^3)
    newXmat <- cbind(rep(1,n), newX[,1], newX[,1]^2, newX[,1]^3)
    beta_hat <- crossprod(tcrossprod(Xmat,solve(crossprod(Xmat))), matrix(Y))
    pred <- newXmat%*%beta_hat
    out <- list(fit = list(object = beta_hat), pred = pred)
    class(out$fit) <- "SL.glm3"
    return(out)
}

predict.SL.glm3 <- function(object, newdata){
    newXmat <- cbind(rep(1,nrow(newdata)), newdata[,1], newdata[,1]^2, newdata[,1]^3)
    return(newXmat%*%object$object)
}


# SL.glm3_slow <- function (Y, X, newX, family, obsWeights, model = TRUE, ...) 
# {
#     if (is.matrix(X)) {
#         X = as.data.frame(X)
#     }
#     fit.glm <- glm(Y ~ x + I(x^2) + I(x^3), data = X, family = family, weights = obsWeights, 
#         model = model)
#     if (is.matrix(newX)) {
#         newX = as.data.frame(newX)
#     }
#     pred <- predict(fit.glm, newdata = newX, type = "response")
#     fit <- list(object = fit.glm)
#     class(fit) <- "SL.glm"
#     out <- list(pred = pred, fit = fit)
#     return(out)
# }

# n <- 5000
# X <- data.frame(x = runif(n, -1, 1))
# Y <- rnorm(n, X$x, 1)
# system.time(SL.glm3_slow(Y = Y, X = X, newX = X, family = "gaussian", obsWeights = rep(1, n)))
# system.time(SL.glm3(Y = Y, X = X, newX = X, family = "gaussian", obsWeights = rep(1, n)))


SL.hal9001 <- function(Y, X, newX, family, obsWeights, ...){
    # get rid of non-unique X columns
    if(dim(X)[2] > 1){
        dropCol <- apply(X, 2, function(x){ length(unique(X)) == 1})
    }else{
        dropCol <- length(unique(X[,1])) == 1
    }
    if(all(dropCol)){
        hal_out <- "all columns of X the same"
        pred <- rep(mean(Y), nrow(newX))
        out <- list(fit = list(object = hal_out, m = mean(Y)), 
                    pred = pred)
    }else{
        hal_out <- fit_hal(Y = Y, X = as.matrix(X)[,!dropCol], yolo = FALSE,
                           fit_type = "glmnet")
        pred <- predict(hal_out, new_data = as.matrix(newX)[,!dropCol])
        out <- list(fit = list(object = hal_out, dropCol = dropCol), pred = pred)
    }
    class(out$fit) <- "SL.hal9001"
    return(out)
}

predict.SL.hal9001 <- function(object, newdata, ...){
    if(!class(object$object) == "character"){
        return(predict(object$object, new_data = as.matrix(newdata)[,!object$dropCol]))        
    }else{
        return(rep(object$m, nrow(newdata)))
    }
}

# makeData <- function(n = n){
#     L0 <- data.frame(x.0 = runif(n,-1,1))
#     A0 <- rbinom(n, 1, plogis(L0$x.0^2))
#     L1 <- data.frame(x.1 = L0$x.0^2*A0 + runif(n))
#     A1 <- rbinom(n, 1, plogis(L0$x.0*L1$x.1))
#     L2 <- rnorm(n, L0$x.0^2*A0*A1 + L1$x.1)
#     return(list(L0 = L0, L1 = L1, L2 = L2, A0 = A0, A1 = A1))
# }

# getTruth <- function(n = 1e6, abar = c(1,1)){
#     L0 <- data.frame(x.0 = runif(n,-1,1))
#     L1 <- data.frame(x.1 = L0$x.0^2*abar[1] + runif(n))
#     L2 <- rnorm(n, L0$x.0^2*abar[1]*abar[2] + L1$x.1)
#     return(mean(L2))
# }

makeData <- function(n = n){
    L0 <- data.frame(x.0 = runif(n,-2,2), x.1 = rbinom(n, 1, 1/2))
    A0 <- rbinom(n, 1, plogis(L0$x.0 - 2*L0$x.0*L0$x.1))
    L1 <- data.frame(x.2 = L0$x.0 - 2 * L0$x.0*L0$x.1 + 0.2 * A0 + runif(n))
    A1 <- rbinom(n, 1, plogis(-L1$x.2 + 2*L1$x.2*L0$x.1 + 0.2*L0$x.0 - 0.1*A0))
    L2 <- rnorm(n, - L0$x.0 + 2*L0$x.0*L0$x.1 - L1$x.2 + A1 + A0)
    return(list(L0 = L0, L1 = L1, L2 = L2, A0 = A0, A1 = A1))
}

getTruth <- function(n = 1e6, abar = c(1,1)){
    L0 <- data.frame(x.0 = runif(n,-2,2), x.1 = rbinom(n, 1, 1/2))
    L1 <- data.frame(x.2 = L0$x.0 - 2 * L0$x.0*L0$x.1 + 0.2 * abar[1] + runif(n))
    L2 <- rnorm(n, - L0$x.0 + 2*L0$x.0*L0$x.1 - L1$x.2 + abar[1] + abar[2])
    return(mean(L2))
}

compare_stuff <- function(n = 1e6){
    grbg <- makeData(n = n)
    truth <- getTruth(n = n)
    c(mean(grbg$L2[grbg$A0==1 & grbg$A1==1]), truth)
}




