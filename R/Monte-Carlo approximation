The following codes show the Monte Carlo approximations of the mean and variance of gene expression data from pooled experiments ($Y_k$). 
The results presented in  Figure S1 of the supplementary file.

```{r, echo=FALSE, "Function required the chunk below"}

EA <- function(q, n){
  #q*choose(n-q, q-1)/choose(n,q)
  q/n
} 
VA <- function(q, n){
  (1-EA(q, n))*EA(q, n)
  
}
EW <- function(q){
  1/q
} 
VW <- function(q, alp){
  (q-1)/((q^2)*(q*alp+1))
}


cv <- function(x){sd(x)/mean(x)}


creatAmatV2 <- function(nr, nc, q, shuffle=TRUE){
  A <- matrix(0, nrow=nr, ncol=nc)
  i=1
  for(j in 1:nc){ 
    if(i>nr){
      i=1
    }
    if(sum(A[i,])<q){
      A[i, j] <- 1
      i <- i+1
    } 
  } 
  if(shuffle){
    A <- A[, sample(1:ncol(A), ncol(A), replace = FALSE)]
  }
  A
}
```


```{r, echo=FALSE, eval=FALSE, "MC approximation"} 
library(MCMCpack) # for Drichilet distribution


# ---------- dispersion = 0.5 ------------------------------------------------------------

n    <- 60
muj  <- rgamma(n, 1, 0.1)
disp <- 0.5
alp=1  #Drichilet parameter
 
Q <- c(2:n)[sapply(2:n, function(i) n%%i==0 & n!=i)]

tht <- as.data.frame(t(sapply(Q, function(q){
  #var.y = ((alp+1)/(n*(q*alp+1)))*sum(muj+(disp+1)*(muj^2)) - sum((muj^2))/(n^2) 
  var.y = sum((VA(q, n)+EA(q, n)^2)*(VW(q, alp)+EW(q)^2)*(muj+(disp+1)*(muj^2)) 
                  - (EA(q, n)^2)*(EW(q)^2)*(muj^2))
  c(mean.y=mean(muj), var.y=var.y, cv.y=sqrt(var.y)/mean(muj), q=q, m=n/q)
})))
  
B <- 2000
mcs <- as.data.frame(do.call(rbind, lapply(Q, function(q){
  m <- n/q
   
 est=t(sapply(1:B, function(i){  
   U <- sapply(muj, function(m) rnbinom(1, mu=m, size=1/disp))
   A <- creatAmatV2(nr = m, nc = n, q = q, shuffle = TRUE)
   y <- numeric(m)  
   for(k in 1:m){ 
     w <- rdirichlet(1, alpha = rep(alp, q))
     y[k] <- sum(U[A[k,]==1]*w)
   }
   #y
   c(mean.y=mean(y), var.y=var(y), cv.y=sd(y)/mean(y), q=q, m=m)
 }))
 est
 #data.frame(mean.y=colMeans(est), var.y=apply(est, 2, sd), cv.y=apply(est, 2, cv), q=q, m=m)
})))
 
saveRDS(list(mcs=mcs, tht=tht), ".../MonteCarloSim1.rds")


# ---------- dispersion = 2 ------------------------------------------------------------

n    <- 60
muj  <- rgamma(n, 1, 0.1)
disp <- 2
alp=1  #Drichilet parameter
 
Q <- c(2:n)[sapply(2:n, function(i) n%%i==0 & n!=i)]

tht <- as.data.frame(t(sapply(Q, function(q){
  #var.y = ((alp+1)/(n*(q*alp+1)))*sum(muj+(disp+1)*(muj^2)) - sum((muj^2))/(n^2) 
  var.y = sum((VA(q, n)+EA(q, n)^2)*(VW(q, alp)+EW(q)^2)*(muj+(disp+1)*(muj^2)) 
                  - (EA(q, n)^2)*(EW(q)^2)*(muj^2))
  c(mean.y=mean(muj), var.y=var.y, cv.y=sqrt(var.y)/mean(muj), q=q, m=n/q)
})))


B <- 2000
mcs <- as.data.frame(do.call(rbind, lapply(Q, function(q){
  m <- n/q
   
 est=t(sapply(1:B, function(i){ 
   U <- sapply(muj, function(m) rnbinom(1, mu=m, size=1/disp))
   A <- creatAmatV2(nr = m, nc = n, q = q)
   y <- numeric(m)
   for(k in 1:m){ 
     w <- rdirichlet(1, alpha = rep(alp, q))
     y[k] <- sum(U[A[k,]==1]*w)
   }
  #y
   c(mean.y=mean(y), var.y=var(y), cv.y=sd(y)/mean(y), q=q, m=m)
 }))
 est
 #data.frame(mean.y=colMeans(est), var.y=apply(est, 2, var), cv.y=apply(est, 2, cv), q=q, m=m)
})))

saveRDS(list(mcs=mcs, tht=tht), ".../MonteCarloSim2.rds")


# ---------- dispersion = 5 ------------------------------------------------------------

n    <- 60
muj  <- rgamma(n, 1, 0.1)
disp <- 5
alp=1  #Drichilet parameter
 
Q <- c(2:n)[sapply(2:n, function(i) n%%i==0 & n!=i)]

tht <- as.data.frame(t(sapply(Q, function(q){
  #var.y = ((alp+1)/(n*(q*alp+1)))*sum(muj+(disp+1)*(muj^2)) - sum((muj^2))/(n^2) 
  var.y = sum((VA(q, n)+EA(q, n)^2)*(VW(q, alp)+EW(q)^2)*(muj+(disp+1)*(muj^2)) 
                  - (EA(q, n)^2)*(EW(q)^2)*(muj^2))
  c(mean.y=mean(muj), var.y=var.y, cv.y=sqrt(var.y)/mean(muj), q=q, m=n/q)
})))


B <- 2000
mcs <- as.data.frame(do.call(rbind, lapply(Q, function(q){
  m <- n/q
   
 est=t(sapply(1:B, function(i){ 
   U <- sapply(muj, function(m) rnbinom(1, mu=m, size=1/disp))
   A <- creatAmatV2(nr = m, nc = n, q = q)
   y <- numeric(m)
   for(k in 1:m){ 
     w <- rdirichlet(1, alpha = rep(alp, q))
     y[k] <- sum(U[A[k,]==1]*w)
   }
   #y
   c(mean.y=mean(y), var.y=var(y), cv.y=sd(y)/mean(y), q=q, m=m)
 }))
 est
 #data.frame(mean.y=colMeans(est), var.y=apply(est, 2, var), cv.y=apply(est, 2, cv), q=q, m=m)
})))

saveRDS(list(mcs=mcs, tht=tht), ".../MonteCarloSim3.rds")



#-------------------------------------------------------------------------------
# Generte plots
mcs0p5phi <- readRDS(".../MonteCarloSim1.rds")
mcs2phi   <- readRDS(".../MonteCarloSim2.rds")
mcs5phi   <- readRDS(".../MonteCarloSim3.rds")
 
Q <- mcs0p5phi$tht$q

par(mfrow=c(3,2), cex.lab=1.5, mar=c(4, 5, 4, 1))
boxplot(mcs0p5phi$mcs$mean.y~mcs0p5phi$mcs$q, xlab="q", ylab=expression(mu[Y[k]]), 
        main=expression(mu[j]~"~"~Gamma(1, 0.1)~~~~phi==0.5)) 
points(1:length(Q), tapply(mcs0p5phi$mcs$mean.y, mcs0p5phi$mcs$q, mean), col=4, pch=16, cex=1.5)
lines(1:length(Q), mcs0p5phi$tht$mean.y, col=2, type = 'b', pch=16, cex=1.5)
legend("topleft", c("MC approximation", "Theoretical"), col=c(1,2), lty=1, lwd=2, cex=1.25)

boxplot(mcs0p5phi$mcs$var.y~mcs0p5phi$mcs$q, xlab="q", ylab=expression(sigma[Y[k]]^2),
        main=expression(mu[j]~"~"~Gamma(1, 0.1)~~~~phi==0.5), ylim=c(0, 1000)) 
points(1:length(Q), tapply(mcs0p5phi$mcs$var.y, mcs0p5phi$mcs$q, mean), col=4, pch=16, cex=1.5)
lines(1:length(Q), mcs0p5phi$tht$var.y, col=2, type = 'b', pch=16, cex=1.5)
 


boxplot(mcs2phi$mcs$mean.y~mcs2phi$mcs$q, xlab="q", ylab=expression(mu[Y[k]]), 
        main=expression(mu[j]~"~"~Gamma(1, 0.1)~~~~phi==2)) 
points(1:length(Q), tapply(mcs2phi$mcs$mean.y, mcs2phi$mcs$q, mean), col=4, pch=16, cex=1.5)
lines(1:length(Q), mcs2phi$tht$mean.y, col=2, type = 'b', pch=16, cex=1.5)

boxplot(mcs2phi$mcs$var.y~mcs2phi$mcs$q, xlab="q", ylab=expression(sigma[Y[k]]^2),
        main=expression(mu[j]~"~"~Gamma(1, 0.1)~~~~phi==2), ylim=c(0, 1000)) 
points(1:length(Q), tapply(mcs2phi$mcs$var.y, mcs2phi$mcs$q, mean), col=4, pch=16, cex=1.5)
lines(1:length(Q), mcs2phi$tht$var.y, col=2, type = 'b', pch=16, cex=1.5)



boxplot(mcs5phi$mcs$mean.y~mcs5phi$mcs$q, xlab="q", ylab=expression(mu[Y[k]]), 
        main=expression(mu[j]~"~"~Gamma(1, 0.1)~~~~phi==5)) 
points(1:length(Q), tapply(mcs5phi$mcs$mean.y, mcs5phi$mcs$q, mean), col=4, pch=16, cex=1.5)
lines(1:length(Q), mcs5phi$tht$mean.y, col=2, type = 'b', pch=16, cex=1.5)

boxplot(mcs5phi$mcs$var.y~mcs5phi$mcs$q, xlab="q", ylab=expression(sigma[Y[k]]^2),
        main=expression(mu[j]~"~"~Gamma(1, 0.1)~~~~phi==5), ylim=c(0, 2000)) 
points(1:length(Q), tapply(mcs5phi$mcs$var.y, mcs5phi$mcs$q, mean), col=4, pch=16, cex=1.5)
lines(1:length(Q), mcs5phi$tht$var.y, col=2, type = 'b', pch=16, cex=1.5)

box("outer")
