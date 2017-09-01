SL.hal9001 <- function(Y, X, newX, family, obsWeights, ...){
    hal_out <- fit_hal(Y = Y, X = as.matrix(X))
    pred <- predict(hal_out, newdata = as.matrix(newX))
    out <- list(fit = list(object = hal_out), pred = pred)
    class(out$fit) <- "SL.hal9001"
    return(out)
}

predict.SL.hal9001 <- function(object, newdata, ...){
	return(predict(object$object, newdata = as.matrix(newdata)))
}