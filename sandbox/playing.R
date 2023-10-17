# playing around locally
source("~/Dropbox/R/drinf/sandbox/makeData.R")
library(drinf)
library(SuperLearner)
library(hal9001)
set.seed(1234124)
# dat <- makeData(n = 10000)

# do.one <- function(){
	dat <- makeData(n = 2000)

	# faster to call mean.tmle
	# debug(drinf.tmle)
	object <- drinf.tmle(
	L0 = dat$L0, L1 = dat$L1, L2 = dat$L2, A0 = dat$A0, A1 = dat$A1, 
	abar = c(1,1), 
	cvFolds = 1, 
	SL.Q = "SL.glm",
	SL.g = "SL.glm.interaction", 
	SL.Qr = "SL.gam",
	SL.gr = "SL.gam",
	universal = TRUE, 
	universalStepSize = 1e-8,
	tolIF = 1/length(dat$L0[,1]),
	flucOrd = c("targetg0","targetg1","redReg",
	           "targetQ2","targetQ1", "redReg"),
	return.models = TRUE,
	verbose = FALSE,
	maxIter = 20,
	return.ltmle = TRUE,
	allatonce = FALSE,
	tolg = 1e-2,
	tolQ = 1e-2, stratify = TRUE
	)
	# object$est_trace
	c(object$est, object$est_trace[1], object$est.ltmle)
}

rslt <- replicate(50, do.one())



# visualizing g misspecification
debug(drinf:::estimateG)
g00 <- plogis(2*L0$x.0*L0$x.1)
plot(g00 ~ g0n, col = L0$x.1 + 1)
g10 <- plogis(2*L0$x.0*L0$x.1 + 0.2*L1$x.2 - 0.1 * 1)
plot(g10 ~ g1n, col = L0$x.1 + 1)

# visualizing Q2 hal fit
Q20 <- L0$x.0*L0$x.1 + L1$x.2 + 2
plot(Q20 ~ Q2n, col = L0$x.1 + 1)
abline(a = 0, b = 1, col = 4)
# vary x.0
pred_seq <- seq(-2,2,length=1000)
pred_data_1 <- data.frame(x.0 = pred_seq,
                          x.1 = 1, x.2 = 0)
pred_1 <- predict(Q2mod$fit, newdata = pred_data_1)
pred_data_0 <- data.frame(x.0 = pred_seq,
                          x.1 = 0, x.2 = 0)
pred_0 <- predict(Q2mod$fit, newdata = pred_data_0)
plot(pred_1 ~ pred_seq, col = 1, ylim = c(-2,4))
points(pred_0 ~ pred_seq, col = 2)
# add truth
lines(y = pred_seq*0 + 0 + 2, x = pred_seq, col = 2, lwd = 2)
lines(y = pred_seq*1 + 0 + 2, x = pred_seq, col = 1, lwd = 2)


# visualizing Q2r fit
g_seq <- seq(min(g1n), max(g1n), length=1000)
pred <- predict(Q2rmod$fit, newdata = data.frame(g1n = g_seq))
plot(rQ2 ~ g1n, col = "gray90")
lines(pred ~ g_seq, type = "l")

# visualizing g0r fit
Q_seq <- seq(min(Q1n), max(Q1n), length = 1000)
pred <- predict(g0rmod$fit, newdata = data.frame(Q1n = Q_seq))
plot(A0 ~ Q1n, col = "gray90")
lines(pred ~ Q_seq, type = "l")

# visualizing g1r fit
Q_seq <- seq(min(Q2n), max(Q2n), length = 1000)
pred <- predict(g0rmod$fit, newdata = data.frame(Q2n = Q_seq))
plot(I(A0*A1) ~ Q2n, col = "gray90")
lines(pred ~ Q_seq, type = "l")

# visualizing h0rfit
Q_seq <- seq(min(Q1n), max(Q1n), length = 1000)
pred <- predict(h0rmod$fit, newdata = data.frame(Q1n = Q_seq))
plot(rg0 ~ Q1n, col = "gray90")
lines(pred ~ Q_seq, type = "l")

# try re-fit with something else
otherAlg <- "SL.gam"
h0rmod_mod <- do.call(ifelse(multiAlgos, getFromNamespace("SuperLearner", 
    "SuperLearner"), otherAlg), args = list(Y = rg0, X = data.frame(Q1n = Q1n), 
    newX = data.frame(Q1n = Q1n), SL.library = SL.gr, obsWeights = rep(1, 
        length(Q1n)), family = gaussian(), verbose = verbose))

Q_seq <- seq(min(Q1n), max(Q1n), length = 1000)
pred <- predict(h0rmod_mod$fit, newdata = data.frame(Q1n = Q_seq))
plot(rg0 ~ Q1n, col = "gray90")
lines(pred ~ Q_seq, type = "l")

# visualizing h1rfit
Q_seq <- seq(min(Q2n), max(Q2n), length = 1000)
pred <- predict(h1rmod$fit, newdata = data.frame(Q2n = Q_seq))
plot(rg1 ~ Q2n, col = "gray90")
lines(pred ~ Q_seq, type = "l")

# visualizing hbar fit
Q_seq <- seq(min(Q2n), max(Q2n), length = 1000)
pred <- predict(hbarrmod$fit, newdata = data.frame(Q2n = Q_seq))
plot(y = A0/g0nr * h0nr, x = Q2n, col = "gray90")
lines(pred ~ Q_seq, type = "l")


# checking g0n vs. g0nstar fit
plot(gnstar$g0n ~ gn$g0n)