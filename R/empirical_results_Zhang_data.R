# Unequal library sizes in the source data

# -------------------- prepare data -----------------------------
Zhang.data  <- readRDS("Data/Zhang_Full_Data.RData")
str(Zhang.data)
counts <- Zhang.data$counts
group  <- Zhang.data$group
table(group)

# Keep 80 and 80 samples in each group (keep only samples with large library size)
counts2 <- counts[, c(which(group==0), which(group==1))]
group2  <- group[c(which(group==0), which(group==1))]

LS.grp0 <- colSums(counts2[, group2==0])
sel.sample_g0 <- which(order(LS.grp0)>=(length(LS.grp0)-80+1))
LS.grp1 <- colSums(counts2[, group2==1])
sel.sample_g1 <- sum(group==0)+which(order(LS.grp1)>=(length(LS.grp1)-80+1))

counts2 <- counts2[, c(sel.sample_g0, sel.sample_g1)]
dim(counts2)
group2 <- group2[c(sel.sample_g0, sel.sample_g1)]
table(group2)

saveRDS(list(counts=counts2, group=group2), ".../Zhang.data2.rds")

 



# -------------------- generate  data for each scenario -----------------------------
rm(list = ls())
source("Global files/subSampling.R")
library(MCMCpack)
library(RNAsub)


Zhang.data2 <- readRDS(".../Zhang.data2.rds")
counts0 <- Zhang.data2$counts
group0  <- Zhang.data2$group
table(group0)

round(sum(counts0)/1e6)  # we have now a count of 3642 M reads
summary(colSums(counts0)/1e6)



# ---------------- REF Scenario (max budget): ---------------- 
# One has 2 groups each with 40 samples and a budget of 1600 M sequencing reads. The classic scenario would be to make 80 libraries and sequence them so each sample would get 20 M reads (total cost is 80S+80L+1600R).


set.seed(6431)
sel.g0 <- sample(which(group0==0), 40)
sel.g1 <- sample(which(group0==1), 40)

counts_REF <- counts0[, c(sel.g0, sel.g1)]
group_REF  <- group0[c(sel.g0, sel.g1)]
sum(counts_REF)/1e6
summary(colSums(counts_REF)/1e6)

saveRDS(list(counts=counts_REF, group=group_REF), ".../data_REF.rds")

counts0 <- counts_REF
group0  <- group_REF

# ---------------- Scenario A: ---------------- 
# One has 2 groups each with 20 samples and a budget of 800 M sequencing reads. The classic scenario would be to make 40 libraries and sequence them so each sample would get 20 M reads (total cost is 40S+40L+800R).


set.seed(6431)
sel.g0 <- sample(which(group0==0), 20)
sel.g1 <- sample(which(group0==1), 20)

counts_A <- counts0[, c(sel.g0, sel.g1)]
group_A  <- group0[c(sel.g0, sel.g1)]
sum(counts_A)/1e6

saveRDS(list(counts=counts_A, group=group_A), ".../data_A.rds")




# ---------------- Scenario B: ---------------- 
# same as A but half sequencing depth (but no RNA sample pooling yet): One has 2 groups each with 20 samples and a budget of only 400 M sequencing reads. One makes 40 libraries and sequence them so each sample would get 10 M reads (total cost is 40S+40L+400R).

set.seed(62431)
sel.g0 <- sample(which(group0==0), 20)
sel.g1 <- sample(which(group0==1), 20)

counts_B <- counts0[, c(sel.g0, sel.g1)]
group_B  <- group0[c(sel.g0, sel.g1)]

counts_B2 <- RNAsub(counts_B, target.total.count = 400e6, 
                       seed = 62431, BPPARAM = BiocParallel::SnowParam())#round(0.5*counts_B)
sum(counts_B2)/1e6

summary(colSums(counts_B))/1e6 ;summary(colSums(counts_B2))/1e6
plot(colSums(counts_B)/1e6, colSums(counts_B2)/1e6, xlab="800M reads", ylab="400M reads",
     ylim=c(min(colSums(counts_B2)/1e6), max(colSums(counts_B)/1e6)))
abline(0, 1)

saveRDS(list(counts=counts_B2, group=group_B), ".../data_B.rds")


# ---------------- Scenario C: ---------------- 
# same as A but one was able to prepare twice as many RNA samples and libraries, hence sequencing deph is half. One has 2 groups each with 40 samples and a budget of 800 M sequencing reads. One makes 80 libraries and sequence them so each sample would get 10 M reads (total cost is 80S+80L+800R).
set.seed(6431)
sel.g0 <- sample(which(group0==0), 40)
sel.g1 <- sample(which(group0==1), 40)

counts_C <- counts0[, c(sel.g0, sel.g1)]
group_C  <- group0[c(sel.g0, sel.g1)]

counts_C2 <- RNAsub(counts_C, target.total.count = 800e6, 
                       seed = 62431, BPPARAM = BiocParallel::SnowParam())#round(0.5*counts_C)
sum(counts_C2)/1e6

summary(colSums(counts_C))/1e6 ;summary(colSums(counts_C2))/1e6
plot(colSums(counts_C)/1e6, colSums(counts_C2)/1e6, xlab="800M reads", ylab="400M reads",
     ylim=c(min(colSums(counts_C2)/1e6), max(colSums(counts_C)/1e6)), main = "scenario C")
abline(0, 1)


saveRDS(list(counts=counts_C2, group=group_C), ".../data_C.rds")



# ---------------- Scenario C2: ---------------- 
# same as A but one was able to prepare twice as many RNA samples and libraries, hence sequencing deph is half. One has 2 groups each with 40 samples and a budget of 800 M sequencing reads. One makes 80 libraries and sequence them so each sample would get 10 M reads (total cost is 80S+80L+800R).
set.seed(6431)
sel.g0 <- sample(which(group0==0), 40)
sel.g1 <- sample(which(group0==1), 40)

counts_C2 <- counts0[, c(sel.g0, sel.g1)]
group_C2  <- group0[c(sel.g0, sel.g1)]

counts_C22 <- RNAsub(counts_C2, target.total.count = 400e6, 
                       seed = 62431, BPPARAM = BiocParallel::SnowParam())#round(0.5*counts_C)
sum(counts_C22)/1e6

summary(colSums(counts_C2))/1e6 ;summary(colSums(counts_C22))/1e6 
saveRDS(list(counts=counts_C22, group=group_C2), ".../data_C2.rds")



# ---------------- Scenario D: ---------------- 
# same as A, but 40 pools of 2 RNA samples are made. One has 2 groups each with 40 samples and a budget of 800 M sequencing reads. Pools are made per 2 RNA samples, hence 40 libraries and sequenced so each sample would get 20 M reads (total cost is 80S+40L+800R).

set.seed(6431)
sel.g0 <- sample(which(group0==0), 40)
sel.g1 <- sample(which(group0==1), 40)

counts_D <- counts0[, c(sel.g0, sel.g1)]
group_D  <- group0[c(sel.g0, sel.g1)]
sum(counts_D)/1e6

set.seed(2486) 
n <- 40  # biological replicates per group
m <- 20  # number of pools per group to be made
q <- n/m  # pool size per group

pool2_D <- do.call(cbind, lapply(unique(group_D), function(g){
  counts_D.g  <- counts_D[, sample(which(group_D==g), sum(group_D==g))] 
  counts_D.gg <- counts_D.g[, 1:m]   # initialize a matrix
  #w  <- runif(m, 0.2, 0.8)  # random mixing fraction
  w <- rdirichlet(m, alpha = rep(1, q))
  for(i in 1:m){ 
    counts_D.gg[, i] <- round(w[i,1]*counts_D.g[, i]+w[i,2]*counts_D.g[, m+i])
  }
  counts_D.gg
}))
dim(pool2_D)
sum(pool2_D)/1e6

group_D2 <- rep(unique(group_D), each=m)
table(group_D2)

summary(colSums(pool2_D))/1e6 ; summary(colSums(counts_D))/1e6

saveRDS(list(counts=pool2_D, group=group_D2), ".../data_D.rds")




# ---------------- Scenario E: ---------------- 
# same as D, but half the sequencing depth. One has 2 groups each with 40 samples and a budget of 400 M sequencing reads. Pools are made per 2 RNA samples, hence 40 libraries and sequenced so each sample would get 10 M reads (total cost is 80S+40L+400R).

set.seed(6431)
sel.g0 <- sample(which(group0==0), 40)
sel.g1 <- sample(which(group0==1), 40)

counts_E <- counts0[, c(sel.g0, sel.g1)]
group_E  <- group0[c(sel.g0, sel.g1)]

counts_E2 <- RNAsub(counts_E, target.total.count = 800e6, seed = 6431, BPPARAM = BiocParallel::SnowParam()) #round(0.5*counts_E)
sum(counts_E2)/1e6

set.seed(2486) 
n <- 40  # biological replicates per group
m <- 20  # number of pools per group to be made
q <- n/m  # pool size per group

pool2_E <- do.call(cbind, lapply(unique(group_E), function(g){
  counts_E.g  <- counts_E2[, sample(which(group_E==g), sum(group_E==g))] 
  counts_E.gg <- counts_E.g[, 1:m]   # initialize a matrix
  w <- rdirichlet(m, alpha = rep(1, q))
  for(i in 1:m){ 
    counts_E.gg[, i] <- round(w[i,1]*counts_E.g[, i]+w[i,2]*counts_E.g[, m+i])
  }
  counts_E.gg
}))
dim(pool2_E)
sum(pool2_E)/1e6

group_E2 <- rep(unique(group_E), each=m)
table(group_E2)

summary(colSums(pool2_E))/1e6 ; summary(colSums(counts_E))/1e6 ; summary(colSums(counts_E2))/1e6

saveRDS(list(counts=pool2_E, group=group_E2), ".../data_E.rds")


 
# ---------------- Scenario F:  NEW ---------------- 
# same as A, but now with RNA sample pooling. One has 2 groups each with 20 samples and a budget of 800 M sequencing reads. Pools are made per 2 RNA samples, hence 20 libraries and sequenced so each sample would get 40 M reads (total cost is 40S+20L+800R).

set.seed(6431)
sel.g0 <- sample(which(group0==0), 20)
sel.g1 <- sample(which(group0==1), 20)

counts_F <- counts0[, c(sel.g0, sel.g1)]
group_F  <- group0[c(sel.g0, sel.g1)]

counts_F2 <- RNAsub(counts_F, target.total.count = 400e6, seed = 4315, BPPARAM = BiocParallel::SnowParam())
sum(counts_F2)/1e6

set.seed(2486) 
n <- 20  # biological replicates per group
m <- 10  # number of pools per group to be made
q <- n/m  # pool size per group

pool2_F <- do.call(cbind, lapply(unique(group_F), function(g){
  counts_F.g  <- counts_F2[, sample(which(group_F==g), sum(group_F==g))] 
  counts_F.gg <- counts_F.g[, 1:m]   # initialize a matrix
  w <- rdirichlet(m, alpha = rep(1, q)) #w   <- runif(m, 0.2, 0.8)  # random mixing fraction
  for(i in 1:m){ 
    counts_F.gg[, i] <- round(w[i, 1]*counts_F.g[, i]+ w[i, 2]*counts_F.g[, m+i])
  }
  counts_F.gg
}))
dim(pool2_F)
sum(pool2_F)/1e6

group_F2 <- rep(unique(group_F), each=m)
table(group_F2)

summary(colSums(pool2_F))/1e6 ; summary(colSums(counts_F))/1e6 ; summary(colSums(counts_F2))/1e6
saveRDS(list(counts=pool2_F, group=group_F2), ".../data_F.rds")



# ---------------- Scenario G: ---------------- 
# same as F, but half the sequencing deph. One has 2 groups each with 20 samples and a budget of 400 M sequencing reads. Pools are made per 2 RNA samples, hence 20 libraries and sequenced so each sample would get 20 M reads (total cost is 40S+20L+400R).

set.seed(6431)
sel.g0 <- sample(which(group0==0), 20)
sel.g1 <- sample(which(group0==1), 20)

counts_G <- counts0[, c(sel.g0, sel.g1)]
group_G  <- group0[c(sel.g0, sel.g1)]

sum(counts_G)/1e6

set.seed(2486) 
n <- 20  # biological replicates per group
m <- 10  # number of pools per group to be made
q <- n/m  # pool size per group

pool2_G <- do.call(cbind, lapply(unique(group_G), function(g){
  counts_G.g  <- counts_G[, sample(which(group_G==g), sum(group_G==g))] 
  counts_G.gg <- counts_G.g[, 1:m]   # initialize a matrix
  w <- rdirichlet(m, alpha = rep(1, q)) #w  <- runif(m, 0.2, 0.8)  # random mixing fraction
  for(i in 1:m){ 
    counts_G.gg[, i] <- round(w[i,1]*counts_G.g[, i]+w[i,2]*counts_G.g[, m+i])
  }
  counts_G.gg
}))
dim(pool2_G)
sum(pool2_G)/1e6

group_G2 <- rep(unique(group_G), each=m)
table(group_G2)

summary(colSums(pool2_G))/1e6 ; summary(colSums(counts_G))/1e6 
saveRDS(list(counts=pool2_G, group=group_G2), ".../data_G.rds")




 
# ---------------- Scenario H2: ---------------- 
set.seed(6431)
sel.g0 <- sample(which(group0==0), 40)
sel.g1 <- sample(which(group0==1), 40)

counts_H2 <- counts0[, c(sel.g0, sel.g1)]
group_H2  <- group0[c(sel.g0, sel.g1)]
 
sum(counts_H2)/1e6

set.seed(2486) 
n <- 40  # biological replicates per group
m <- 10  # number of pools per group to be made
q <- n/m  # pool size per group

pool2_H2 <- do.call(cbind, lapply(unique(group_H2), function(g){
  counts_H2.g  <- counts_H2[, sample(which(group_H2==g), sum(group_H2==g))] 
  counts_H2.gg <- counts_H2.g[, 1:m]   # initialize a matrix
  # w           <- matrix(rnorm(m*(q-1), 1/q, 0.05), nrow=q-1, ncol=m)  # random mixing fraction
  # w <- rbind(w, 1-colSums(w))
  w <- rdirichlet(m, alpha = rep(1, q))
  for(i in 1:m){ 
    counts_H2.gg[, i] <-   round(w[i,1]*counts_H2.g[, (0*m)+i] +
                                 w[i,2]*counts_H2.g[, (1*m)+i] +
                                 w[i,3]*counts_H2.g[, (2*m)+i] + 
                                 w[i,4]*counts_H2.g[, (3*m)+i])
  }
  counts_H2.gg
}))
dim(pool2_H2)
sum(pool2_H2)/1e6

group_H22 <- rep(unique(group_H2), each=m)
table(group_H22)

summary(colSums(pool2_H2))/1e6 ; summary(colSums(counts_H2))/1e6 
saveRDS(list(counts=pool2_H2, group=group_H22), ".../data_H2.rds")



# ---------------- Scenario H3: ---------------- 
set.seed(6431)
sel.g0 <- sample(which(group0==0), 40)
sel.g1 <- sample(which(group0==1), 40)

counts_H3 <- counts0[, c(sel.g0, sel.g1)]
group_H3  <- group0[c(sel.g0, sel.g1)]

counts_H3 <- RNAsub(counts_H3, target.total.count = 800e6, seed = 6431, BPPARAM = BiocParallel::SnowParam()) #round(0.5*counts_H3) 
sum(counts_H3)/1e6

set.seed(2486) 
n <- 40  # biological replicates per group
m <- 10  # number of pools per group to be made
q <- n/m  # pool size per group

pool2_H3 <- do.call(cbind, lapply(unique(group_H3), function(g){
  counts_H3.g  <- counts_H3[, sample(which(group_H3==g), sum(group_H3==g))] 
  counts_H3.gg <- counts_H3.g[, 1:m]   # initialize a matrix
  # w            <- matrix(rnorm(m*(q-1), 1/q, 0.05), nrow=q-1, ncol=m)  # random mixing fraction
  # w <- rbind(w, 1-colSums(w))
  w <- rdirichlet(m, alpha = rep(1, q))
  for(i in 1:m){ 
    counts_H3.gg[, i] <-   round(w[i,1]*counts_H3.g[, (0*m)+i] +
                                 w[i,2]*counts_H3.g[, (1*m)+i] +
                                 w[i,3]*counts_H3.g[, (2*m)+i] + 
                                 w[i,4]*counts_H3.g[, (3*m)+i])
  }
  counts_H3.gg
}))
dim(pool2_H3)
sum(pool2_H3)/1e6

group_H32 <- rep(unique(group_H3), each=m)
table(group_H32)

summary(colSums(pool2_H3))/1e6 ; summary(colSums(counts_H3))/1e6 
saveRDS(list(counts=pool2_H3, group=group_H32), ".../data_H3.rds")


 
# ---------------- Scenario I3: ---------------- 
set.seed(6431)
sel.g0 <- sample(which(group0==0), 20)
sel.g1 <- sample(which(group0==1), 20)

counts_I3 <- counts0[, c(sel.g0, sel.g1)]
group_I3  <- group0[c(sel.g0, sel.g1)]

#counts_I3 <- round(2*counts_I3) 
sum(counts_I3)/1e6

set.seed(2486) 
n <- 20  # biological replicates per group
m <- 5  # number of pools per group to be made
q <- n/m  # pool size per group

pool2_I3 <- do.call(cbind, lapply(unique(group_I3), function(g){
  counts_I3.g  <- counts_I3[, sample(which(group_I3==g), sum(group_I3==g))] 
  counts_I3.gg <- counts_I3.g[, 1:m]   # initialize a matrix
  # w           <- matrix(rnorm(m*(q-1), 1/q, 0.05), nrow=q-1, ncol=m)  # random mixing fraction
  # w <- rbind(w, 1-colSums(w))
  w <- rdirichlet(m, alpha = rep(1, q))
  for(i in 1:m){ 
    counts_I3.gg[, i] <- round(  w[i,1]*counts_I3.g[, (0*m)+i] +
                                 w[i,2]*counts_I3.g[, (1*m)+i] +
                                 w[i,3]*counts_I3.g[, (2*m)+i] + 
                                 w[i,4]*counts_I3.g[, (3*m)+i])
  }
  counts_I3.gg
}))
dim(pool2_I3)
sum(pool2_I3)/1e6

group_I32 <- rep(unique(group_I3), each=m)
table(group_I32)

summary(colSums(pool2_I3))/1e6 ; summary(colSums(counts_I3))/1e6 
saveRDS(list(counts=pool2_I3, group=group_I32), ".../data_I3.rds")




# ---------------- Scenario I4: ---------------- 
set.seed(6431)
sel.g0 <- sample(which(group0==0), 20)
sel.g1 <- sample(which(group0==1), 20)

counts_I4 <- counts0[, c(sel.g0, sel.g1)]
group_I4  <- group0[c(sel.g0, sel.g1)]

counts_I4 <- RNAsub(counts_I4, target.total.count = 400e6, seed=6431, BPPARAM = BiocParallel::SnowParam()) #round(0.5*counts_I4) 
sum(counts_I4)/1e6

set.seed(2486) 
n <- 20  # biological replicates per group
m <- 5  # number of pools per group to be made
q <- n/m  # pool size per group

pool2_I4 <- do.call(cbind, lapply(unique(group_I4), function(g){
  counts_I4.g  <- counts_I4[, sample(which(group_I4==g), sum(group_I4==g))] 
  counts_I4.gg <- counts_I4.g[, 1:m]   # initialize a matrix
  # w           <- matrix(rnorm(m*(q-1), 1/q, 0.05), nrow=q-1, ncol=m)  # random mixing fraction
  # w <- rbind(w, 1-colSums(w))
  w <- rdirichlet(m, alpha = rep(1, q))
  for(i in 1:m){ 
    counts_I4.gg[, i] <- round(w[i,1]*counts_I4.g[, (0*m)+i] +
                                 w[i,2]*counts_I4.g[, (1*m)+i] +
                                 w[i,3]*counts_I4.g[, (2*m)+i] + 
                                 w[i,4]*counts_I4.g[, (3*m)+i])
  }
  counts_I4.gg
}))
dim(pool2_I4)
sum(pool2_I4)/1e6

group_I42 <- rep(unique(group_I4), each=m)
table(group_I42)

summary(colSums(pool2_I4))/1e6 ; summary(colSums(counts_I4))/1e6 
saveRDS(list(counts=pool2_I4, group=group_I42), ".../data_I4.rds")

 
 
