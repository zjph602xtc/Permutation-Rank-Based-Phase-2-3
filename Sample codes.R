
# simulation power  ------------------------------------------------------
setwd("F:/BMS/Rong paper")
dunnett <- function(p, m){
  f <- AdjustPvalues(c(p, rep(0.999, m-1)),proc = "DunnettAdj",par = parameters(n = Inf))[1]
  return(f)
}

library(MASS)
library(mvtnorm)
library(Mediana)
mu_e <- c(0.3,0.3,0.3,0) # last is control
mu_b <- c(0.25,0.3,0.35,0)

# mu_e <- c(0,0,0,0) # last is control
n1 <- 50  # sample size for stage 1
n2 <- 185  # sample size for stage 2

rho_list <- seq(0,0.3,0.1)
rho_list <- seq(0.4,0.6,0.1)
rho_list <- seq(0.7,1,0.1)
for (rho in rho_list){
  
  # rho <- 0.3  # cor between Ye, Yb
  w1 <- sqrt(n1/(n1+n2))
  w2 <- sqrt(n2/(n1+n2)) # stage 1 and 2 weight
  
  sig_e <- matrix(nrow=length(mu_e), ncol=length(mu_e), 0)
  diag(sig_e) <- 1
  sig_b <- matrix(nrow=length(mu_b), ncol=length(mu_b), 0)
  diag(sig_b) <- 1
  sig_be <- matrix(nrow=length(mu_b), ncol=length(mu_b), 0)
  diag(sig_be) <- rho
  sig <- rbind(cbind(sig_e, sig_be), cbind(t(sig_be),sig_b)) # Ye first, then Yb
  
  
  iter <- 50000
  z_s1_e_all <- matrix(nrow=iter, ncol=3)
  z_s1_b_all <-  matrix(nrow=iter, ncol=3)
  z_s1_unadjsted <- matrix(nrow=iter, ncol=3) # rank 1, 2, 3
  p_s1_unadjusted <- p_s1_dunnett <- p_s1_perm <- p_s1_exact <- matrix(nrow=iter, ncol=3)
  z_s2_unadjusted_all <- z_s2_dunnett_all <- z_s2_perm_all <- z_s2_exact_all <- matrix(nrow=iter, ncol=3)
  for (i in 1:iter){
    
    dat_s1 <- mvrnorm(n1, c(mu_e, mu_b), Sigma = sig) # Ye first, then Yb
    dat_s2 <- mvrnorm(n2, mu_e, Sigma = sig_e)
    
    z_s1_e <- sapply(1:3, function(i)t.test(dat_s1[,i], dat_s1[, 4], var.equal = T)$statistic)
    z_s1_b <- sapply(5:7, function(i)t.test(dat_s1[,i], dat_s1[, 8], var.equal = T)$statistic)
    z_s1_e_all[i,] <- z_s1_e
    z_s1_b_all[i,] <- z_s1_b 
    
    m <- 3
    for (k in 1:m){
      selected <- which(rank(z_s1_b)==k)
      z_s1_unadjsted[i, k] <- z_s1_e[selected]
      raw_p <- pnorm(z_s1_unadjsted[i, k], lower.tail = F)
      p_s1_unadjusted[i, k] <- raw_p
      
      # Dunnett -----------------------------------------------------------------
      p_dunnett <- dunnett(raw_p, m=k)
      p_s1_dunnett[i, k] <- p_dunnett
      
      # permutation -------------------------------------------------------------
      z_rand <- rep(NA, 20000)
      for (rrr in 1:20000){
        dat_rand <- t(sapply(1:n1, function(x){tp <- sample.int(4,4); dat_s1[x,c(tp,tp+4)]}))
        tmp_e <- sapply(1:3, function(i)t.test(dat_rand[,i], dat_rand[, 4], var.equal = T)$statistic)
        tmp_b <- sapply(5:7, function(i)t.test(dat_rand[,i], dat_rand[, 8], var.equal = T)$statistic)
        z_rand[rrr] <- tmp_e[rank(tmp_b)==k]
      }
      p_perm <- mean(z_s1_unadjsted[i, k] < z_rand)
      p_s1_perm[i, k] <- p_perm

      # exact -------------------------------------------------------------------
      sig_z_b <- matrix(byrow = T, ncol=3, nrow=3, 0.5)
      diag(sig_z_b) <- 1
      rho_est <- cor(c(dat_s1[,selected],dat_s1[,4]), c(dat_s1[,selected+4], dat_s1[,8]))
      
      a <- matrix(c(1,0,0,0, 0,1,-1,0, 0,0,1,-1), byrow = T, nrow=3)
      if (k==2){
        sig_z_eb <- rbind(c(1, 0.5*rho_est,rho_est,0.5*rho_est), cbind(c(0.5*rho_est,rho_est,0.5*rho_est), sig_z_b))
      }else if (k==1){
        sig_z_eb <- rbind(c(1, rho_est,0.5*rho_est,0.5*rho_est), cbind(c(rho_est,0.5*rho_est,0.5*rho_est), sig_z_b))
      }else if (k==3){
        sig_z_eb <- rbind(c(1, 0.5*rho_est,0.5*rho_est,rho_est), cbind(c(0.5*rho_est,0.5*rho_est,rho_est), sig_z_b))
      }
      p_exact <- 1 - pmvnorm(upper=c(z_s1_unadjsted[i, k], 0, 0), sigma = a %*% sig_z_eb %*% t(a))*6
      p_s1_exact[i, k] <- p_exact
      
      
      # stage 2 -----------------------------------------------------------------
      z_s2_e <- t.test(dat_s2[,selected], dat_s2[, 4], var.equal = T)$statistic
      z_s2_unadjusted_all[i,k] <- w1*z_s1_unadjsted[i, k]+w2*z_s2_e
      z_s2_dunnett_all[i,k] <- w1*qnorm(p_dunnett, lower.tail = F)+w2*z_s2_e
      z_s2_perm_all[i,k] <- w1*qnorm(p_perm, lower.tail = F)+w2*z_s2_e
      z_s2_exact_all[i,k] <- w1*qnorm(min(1,p_exact), lower.tail = F)+w2*z_s2_e
      
    }
    
    
    if (i %% 1000==0)print(i)
  }
  colMeans(z_s2_unadjusted_all > 1.96)
  colMeans(z_s2_dunnett_all > 1.96)
  colMeans(z_s2_perm_all > 1.96)
  colMeans(z_s2_exact_all > 1.96, na.rm = T)
  print(rho)
  save.image(paste0('power_rho=',rho,'.RData'))
}


z_s2_unadjusted <- z_s2_dunnett <- z_s2_perm <- z_s2_exact <- NULL
for (rho in seq(0,1,0.1)){
  load(paste0("F:/BMS/Rong paper/power_rho=", rho, ".RData"))
  z_s2_unadjusted <- rbind(z_s2_unadjusted, colMeans(z_s2_unadjusted_all > 1.96, na.rm = T))
  z_s2_dunnett <- rbind(z_s2_dunnett, colMeans(z_s2_dunnett_all > 1.96, na.rm = T))
  z_s2_perm <- rbind(z_s2_perm, colMeans(z_s2_perm_all > 1.96, na.rm = T))
  z_s2_exact <- rbind(z_s2_exact, colMeans(z_s2_exact_all > 1.96, na.rm = T))
}



# simulation for appendix ---------------------------------------------
setwd("C:/Users/xut13/OneDrive - Bristol Myers Squibb/work/Rong paper/with perm/simu_appendix")


dunnett <- function(p, m){
  f <- AdjustPvalues(c(p, rep(0.999, m-1)),proc = "DunnettAdj",par = parameters(n = Inf))[1]
  return(f)
}
g <- function(N,mu,rate,r){
  corr <- matrix(c(1,r,r,1),nrow=2)
  Z = rmvnorm(N, sigma = corr)
  cdf <- pnorm(Z)
  res <- data.frame(Bio=qnorm(cdf[, 2], mean = mu, sd = 1),CMR=qbinom(cdf[, 1], 1, prob=rate))
  return(res)
}


library(MASS)
library(mvtnorm)
library(Mediana)
library(nph)
# mu_e <- c(0.78,0.78,0.78,0.78,0.5) # last is control
mu_b <- c(0.25,0.3,0.35,0.4,0)
mu_e <- c(0.5,0.5,0.5,0.5,0.5) # last is control

n1 <- 20  # sample size for stage 1
n2 <- 40  # sample size for stage 2

rho_list <- seq(0,0.3,0.1)
rho_list <- seq(0.4,0.6,0.1)
rho_list <- seq(0.7,1,0.1)
rho_list=0
for (rho in rho_list){
  
  # rho <- 0.3  # cor between Ye, Yb
  w1 <- sqrt(n1/(n1+n2))
  w2 <- sqrt(n2/(n1+n2)) # stage 1 and 2 weight
  
  sig_e <- matrix(nrow=length(mu_e), ncol=length(mu_e), 0)
  diag(sig_e) <- 1
  sig_b <- matrix(nrow=length(mu_b), ncol=length(mu_b), 0)
  diag(sig_b) <- 1
  sig_be <- matrix(nrow=length(mu_b), ncol=length(mu_b), 0)
  diag(sig_be) <- rho
  sig <- rbind(cbind(sig_e, sig_be), cbind(t(sig_be),sig_b)) # Ye first, then Yb
  
  
  iter <- 50000
  z_s1_e_all <- matrix(nrow=iter, ncol=4)
  z_s1_b_all <-  matrix(nrow=iter, ncol=4)
  z_s1_unadjsted <- matrix(nrow=iter, ncol=4) # rank 1, 2, 3, 4
  p_s1_unadjusted <- p_s1_dunnett <- p_s1_perm <- matrix(nrow=iter, ncol=4)
  z_s2_unadjusted_all <- z_s2_dunnett_all <- z_s2_perm_all <- matrix(nrow=iter, ncol=4)
  for (i in 1:iter){
    tp1 <- g(n1, mu_b[1], mu_e[1], rho)
    tp2 <- g(n1, mu_b[2], mu_e[2], rho)
    tp3 <- g(n1, mu_b[3], mu_e[3], rho)
    tp4 <- g(n1, mu_b[4], mu_e[4], rho)
    tp5 <- g(n1, mu_b[5], mu_e[5], rho)
    dat_s1 <- cbind(tp1[,2],tp2[,2],tp3[,2],tp4[,2],tp5[,2],tp1[,1],tp2[,1],tp3[,1],tp4[,1],tp5[,1])# Ye first, then Yb
    dat_s2 <- sapply(mu_e, function(x)rbinom(n2, 1, x))
    
    
    z_s1_e <- sapply(1:4, function(i){dd <- prop.test(c(sum(dat_s1[,i]),sum(dat_s1[, 5])), c(n1,n1)); sqrt(dd$statistic) * ifelse(diff(dd$estimate)<0, -1, 1)})
    z_s1_b <- sapply(6:9, function(i)t.test(dat_s1[,i], dat_s1[, 10], var.equal = T)$statistic)
    z_s1_e_all[i,] <- z_s1_e
    z_s1_b_all[i,] <- z_s1_b 
    
    m <- 4
    for (k in 1:m){
      selected <- which(rank(z_s1_b, ties.method = 'random')==k)
      z_s1_unadjsted[i, k] <- z_s1_e[selected]
      raw_p <- pnorm(z_s1_unadjsted[i, k], lower.tail = T)
      p_s1_unadjusted[i, k] <- raw_p
      
      ## Dunnett -----------------------------------------------------------------
      p_dunnett <- dunnett(raw_p, m=k)
      p_s1_dunnett[i, k] <- p_dunnett
      
      ## permutation -------------------------------------------------------------
      z_rand <- rep(NA, 20000)
      for (rrr in 1:20000){
        dat_rand <- t(sapply(1:n1, function(x){tp <- sample.int(5,5); dat_s1[x,c(tp,tp+5)]}))
        
        tmp_e <- sapply(1:4, function(i){dd <- prop.test(c(sum(dat_rand[,i]),sum(dat_rand[, 5])), c(n1,n1)); sqrt(dd$statistic) * ifelse(diff(dd$estimate)<0, -1, 1)})
        tmp_b <- sapply(6:9, function(i)t.test(dat_rand[,i], dat_rand[, 10], var.equal = T)$statistic)
        z_rand[rrr] <- tmp_e[rank(tmp_b, ties.method = 'random')==k]
      }
      p_perm <- mean(z_s1_unadjsted[i, k] > z_rand)
      p_s1_perm[i, k] <- p_perm
      
      
      ## stage 2 -----------------------------------------------------------------
      z_s2_e <- sapply(selected, function(i){dd <- prop.test(c(sum(dat_s2[,i]),sum(dat_s2[, 5])), c(n2,n2)); sqrt(dd$statistic) * ifelse(diff(dd$estimate)<0, -1, 1)})
      
      z_s2_unadjusted_all[i,k] <- w1*z_s1_unadjsted[i, k]+w2*z_s2_e
      z_s2_dunnett_all[i,k] <- w1*qnorm(p_dunnett, lower.tail = T)+w2*z_s2_e
      z_s2_perm_all[i,k] <- w1*qnorm(p_perm, lower.tail = T)+w2*z_s2_e
    }
    
    
    if (i %% 100==0)print(i)
  }
  colMeans(z_s2_unadjusted_all < -1.96, na.rm = T)
  colMeans(z_s2_dunnett_all < -1.96, na.rm = T)
  colMeans(z_s2_perm_all < -1.96, na.rm = T)
  print(rho)
  save.image(paste0('appendix_s2_null_rho=',rho,'.RData'))
  # save.image(paste0('appendix_s1_null_rho=',rho,'.RData'))
}


z_s2_unadjusted <- z_s2_dunnett <- z_s2_perm <- z_s2_exact <- NULL
for (rho in seq(0,1,0.1)){
  load(paste0("C:/Users/xut13/OneDrive - Bristol Myers Squibb/work/Rong paper/with perm/simu_appendix/appendix_s2_null_rho=", rho, ".RData"))
  z_s2_unadjusted <- rbind(z_s2_unadjusted, colMeans(z_s2_unadjusted_all < -1.96, na.rm = T))
  z_s2_dunnett <- rbind(z_s2_dunnett, colMeans(z_s2_dunnett_all < -1.96, na.rm = T))
  z_s2_perm <- rbind(z_s2_perm, colMeans(z_s2_perm_all < -1.96, na.rm = T))
}


