library(INLA)

set.seed(100)
n = 3000
x <- c()
beta_x <- c()
nivel <- c()
y <- c()

x <- rnorm(1)
beta_x[1] <- rnorm(1)
nivel[1] <- rnorm(1)
y[1] <- rnorm(1)

for(i in 2:n){
  beta_x[i] <- beta_x[i-1] + rnorm(1, sd = 0.1)
  nivel[i] <- nivel[i-1] + rnorm(1, sd = 0.1)
  x[i] <- x[i-1] + rnorm(1, sd = 0.05)
}

x <- x/10

y <- rpois(n, lambda = exp(nivel + beta_x*x))
ts.plot(y)

id <- id1 <- 1:n

formula <- y ~ f(id, model="rw1",param=c(1,0.001), initial= 0, 
                 constr = F) +
  x + 
  +f(id1, x, model="rw1",param=c(0.001,0.001), initial = 0) #+
  #x2 + f(id2, x2, model="rw1",param=c(1,0.001),initial=5)

r = inla(formula, family="poisson", data = data.frame(id,id1,y),
         control.predictor=list(compute=TRUE))
summary(r)

ts.plot(r$summary.random$id1$`0.5quant`)
lines(r$summary.random$id1$`0.975quant`)
lines(r$summary.random$id1$`0.025quant`)
abline(h=0)
lines(beta_x - mean(beta_x), col = "green")


ts.plot(r$summary.random$id$`0.5quant`, ylim=c(-6,6))
lines(r$summary.random$id$`0.975quant`)
lines(r$summary.random$id$`0.025quant`)
abline(h=0)
lines(nivel - mean(nivel), col = "green")







