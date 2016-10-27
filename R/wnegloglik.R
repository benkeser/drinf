#' wnegloglik
#' 
#' A function that computes the weighted negative log-likelihood loss of the 
#' intercept only fluctuation model. 
#' 
#' @param epsilon The scalar parameter of the fluctuation submodel. 
#' @param weight The \code{vector} of weights, i.e., the clever covariates.
#' @param Y The \code{vector} of regression outcomes. 
#' @param offset The \code{vector} of offsets. 
#' 
#' @return A \code{numeric} value of the negative log-likelihood loss

wnegloglik <- function(epsilon, weight, Y, offset){
    X <- as.matrix(cbind(offset, rep(1, length(offset))))
    mu <- plogis(X%*%as.matrix(c(1,epsilon)))
    mu[mu == 0] <- .Machine$double.neg.eps
    mu[mu == 1] <- 1 - .Machine$double.neg.eps
    wloglik <- sum(weight * (Y * log(mu) + (1-Y)*log(1-mu)))
    return(-wloglik)
}

# 
# optim(
#     par = 0.25, fn = wnegloglik, gr = gradient.wnegloglik,
#     method = "L-BFGS-B", lower = -100, upper = 100,
#     control = list(maxit = 10000),
#     Y = L2s, offset = flucOff, weight = flucCov2
# )

offnegloglik <- function(epsilon, weight, Y, offset){
    X <- as.matrix(cbind(offset, weight))
    mu <- plogis(X%*%as.matrix(c(1,epsilon)))
    mu[mu == 0] <- .Machine$double.neg.eps
    mu[mu == 1] <- 1 - .Machine$double.neg.eps
    ologlik <- sum(weight * (Y * log(mu) + (1-Y)*log(1-mu)))
    return(-ologlik)
}

gradient.offnegloglik <- function(epsilon, weight, Y, offset){
    X <- as.matrix(cbind(offset, weight))
    mu <- plogis(X%*%matrix(c(1,epsilon)))
    mu[mu == 0] <- .Machine$double.neg.eps
    mu[mu == 1] <- 1 - .Machine$double.neg.eps
    grad <- crossprod(weight, Y - mu)
    return(-grad)
}

# optim(
#     par = 0.25, fn = offnegloglik, gr = gradient.offnegloglik,
#     method = "L-BFGS-B", lower = -100, upper = 100,
#     control = list(maxit = 10000),
#     Y = L2s, offset = flucOff, weight = flucCov2
# )


#' gradient.wnegloglik
#' 
#' @param epsilon The scalar parameter of the fluctuation submodel. 
#' @param weight The \code{vector} of weights, i.e., the clever covariates.
#' @param Y The \code{vector} of regression outcomes. 
#' @param offset The \code{vector} of offsets. 
#' 
#' @return A \code{vector} of the gradient of the loss evaluated at the 
#' data inputs. 
gradient.wnegloglik <- function(epsilon, weight, Y, offset){
    X <- cbind(offset, rep(1,length(offset)))
    mu <- plogis(X%*%matrix(c(1,epsilon)))
    mu[mu == 0] <- .Machine$double.neg.eps
    mu[mu == 1] <- 1 - .Machine$double.neg.eps
    grad <- crossprod(weight, Y - mu)
    return(-grad)
}
