L0 <- data.frame(x.0 = rnorm(100))
A0 <- rbinom(100, 1, plogis(L0$x.0))
L1 <- data.frame(x.1 = rnorm(100, L0$x.0 + A0))
A1 <- rbinom(100, 1, plogis(L0$x.0 - A0 + L1$x.1))
L2 <- rnorm(100, L0$x.0 + A0 + A1 + L1$x.1)

abar <- c(1,1)

SL.Q.options = list(family = gaussian())
SL.g.options = list(family = binomial())
glm.Q.options = list(family = gaussian())
verbose <- TRUE

SL.g <- c("SL.glm","SL.mean")
SL.Q <- c("SL.glm","SL.mean")
glm.g <- glm.Q <- "."

verbose <- TRUE
tolg <- 1e-8


