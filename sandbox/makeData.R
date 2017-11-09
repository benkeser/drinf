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
        hal_out <- fit_hal(Y = Y, X = as.matrix(X)[,!dropCol], yolo = FALSE)
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