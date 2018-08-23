library('Rcpp')
library("mclust")
sourceCpp("/storage/home/g/gzx103/scratch/sc/funs_rcpp5_noeta.cpp") 

get_nloglik <- function(y,matx,parc_all,parg_all){
  nloglik <- 0
  for (i in 1:dim(y)[1]){
    nloglik <- nloglik + get_nloglik_c(parc_all[i],parg_all, y[i,], matx[i,])
  }
  return(nloglik)
}

mle <- function(y,x,parg0,parc0,iter){
  lik_trace <- get_nloglik(y,x,parc0,parg0)
  for (i in 1:iter){
    nc <- dim(y)[1]
    ng <- dim(y)[2]
    eparc <- matrix(NA,nc)
    eparg <- matrix(NA,nrow=dim(parg0)[1], ncol=dim(parg0)[2])
    
    for (ic in 1:nc){
      temp <- nlminb(start = parc0[ic], objective = get_nloglik_c, parg_all= parg0, y_c=y[ic,], matx_c=x[ic,], lower = 0.0000001, upper = 1000000)
      eparc[ic] <- temp$par
    }
    parc0 <- eparc
    
    for (ig in 1:ng){
      temp <- nlminb(start = parg0[ig,], objective = get_nloglik_g, parc_all= parc0, y_g=y[,ig], matx_g=x[,ig], lower = c(0.00000001,0.00000001,0.00000001,0.00000001), upper = c(1000,1000,10000,10000))
      eparg[ig,] <- temp$par
    }
    parg0 <- eparg
    lik_trace <- c(lik_trace, get_nloglik(y,x,parc0,parg0))
  }
  rm <- NULL
  rm$a0 <- parg0[,1]
  rm$a1 <- parg0[,2]
  rm$b0 <- parg0[,3]
  rm$b1 <- parg0[,4]
  rm$lam <- parc0
  rm$lik_trace <- lik_trace 
  return(rm)
}

Samx <- function(x,y,ac,bg,a0,a1,b0,b1,lam,gam_fix,phi_up){
  C = dim(x)[1]
  G = dim(x)[2]
  gam = gam_fix
  phi = phi_up
  K <- max(ac)
  L <- max(bg)
  nc = tabulate(ac, nbins = K)
  ng = tabulate(bg, nbins = L)
  nx = matrix(NA, K, L)
  for(k in 1:K)
  {for(l in 1:L){	
    nx[k,l] = sum(x[which(ac==k), which(bg==l)]) 
    }
  }
  newx <- matrix(NA, C, G)
  p <- matrix(NA, C, G)
  for (c in 1:C) {
    for (g in 1:G) {
      k = ac[c]
      l = bg[g]
      sx = nx[k,l] - x[c,g]
      sizex = nc[k] * ng[l] 
      if(y[c,g]!=0) {
        logp_cg <- lgamma(a0[g]+y[c,g]) + lgamma(a0[g]+a1[g]) + y[c,g]*log(b0[g]) + log(phi+sizex-1-sx) - (a0[g]+y[c,g])*log(b0[g]*1+1) + log(1 - ((b0[g]*1+1)/(b0[g]*(lam[c]+1)+1))^(a0[g]+y[c,g])) - lgamma(a0[g]) - lgamma(a0[g]+a1[g]+y[c,g]) - y[c,g]*log(b0[g]+b1[g]) - log(gam+sx) + (a0[g]+a1[g]+y[c,g])*log((b0[g]+b1[g])*1+1) - log(1 - (((b0[g]+b1[g])*1+1)/((b0[g]+b1[g])*(1+lam[c])+1))^(a0[g]+a1[g]+y[c,g])) 
      }else{
        if(lam[c]>1){
          logp_cg <- log(phi+sizex-1-sx) - log(gam+sx) - a0[g]*log(b0[g]*1+1) + log(1 - ((b0[g]*1+1)/(b0[g]*(lam[c]+1)+1))^a0[g] + ((b0[g]*1+1)/(b0[g]*lam[c]+1))^a0[g]) + (a0[g]+a1[g])*log((b0[g]+b1[g])*1+1) - log(1 - (((b0[g]+b1[g])*1+1)/((b0[g]+b1[g])*(1+lam[c])+1))^(a0[g]+a1[g]) + (((b0[g]+b1[g])*1+1)/((b0[g]+b1[g])*lam[c]+1))^(a0[g]+a1[g]))
        }else{
          logp_cg <- log(phi+sizex-1-sx) - log(gam+sx) - a0[g]*log(b0[g]*lam[c]+1) + log(((b0[g]*lam[c]+1)/(b0[g]*1+1))^a0[g] - ((b0[g]*lam[c]+1)/(b0[g]*(lam[c]+1)+1))^a0[g] + 1) + (a0[g]+a1[g])*log((b0[g]+b1[g])*lam[c]+1) - log((((b0[g]+b1[g])*lam[c]+1)/((b0[g]+b1[g])*1+1))^(a0[g]+a1[g]) -  (((b0[g]+b1[g])*lam[c]+1)/((b0[g]+b1[g])*(1+lam[c])+1))^(a0[g]+a1[g]) + 1)
        }
      }             
      p[c,g] <- exp(logp_cg)
      newx[c,g] <- rbinom(n=1, size=1, prob=1/(1+p[c,g]))
    }
  }
  sampos=NULL
  sampos$p = p
  sampos$newx = newx
  return(sampos)
}

scDir<-function(x, B, gam_fix, phi_up, ac_ini, bg_ini){	
  lc = dim(x)[1]
  lg = dim(x)[2]
  alpha1 <- 1
  alpha2 <- 1
  gam = gam_fix
  phi = phi_up
  ac <- ac_ini
  bg <- bg_ini
  Kc <- max(ac)
  Kg <- max(bg)
  nc = tabulate(ac, nbins = Kc)
  ng = tabulate(bg, nbins = Kg)
  nx = array(0, dim=c(Kc,Kg))
  for(k in unique(ac)){	
    for(l in unique(bg)){	
      nx[k, l] = sum(x[which(ac==k), which(bg==l)])
    }
  }
  
  for (b in 1:B) {
    for(c in 1:lc){	
      nc[ac[c]] = nc[ac[c]] - 1
      tn = rep(0, Kg)
      for(l in unique(bg)){	
        tn[l] = sum(x[c, which(bg == l)])
      }
      nx[ac[c],] = nx[ac[c],] - tn			
      
      lp = rep(-1e6, Kc)
      for(k in unique(ac[-c])){			
        lp[k] = sum(lgamma(tn + nx[k,] + gam) - lgamma(nx[k,] + gam)) + sum(lgamma(ng * (nc[k] + 1) - tn - nx[k,] + phi) - lgamma(ng * nc[k] - nx[k,] + phi)) - sum(lgamma(ng*(nc[k]+1) + gam + phi) - lgamma(ng*nc[k] + gam + phi)) + log(nc[k]) - log(lc+alpha1-1)
      }
      if(nc[ac[c]] == 0){	
        lp[ac[c]] = sum(lgamma(tn + gam) - lgamma(gam)) + sum(lgamma(ng - tn + phi) - lgamma(phi)) - sum(lgamma(ng + gam + phi) - lgamma(gam + phi)) + log(alpha1) - log(lc+alpha1-1)
      }else{ 
        Kc <- Kc+1
        lp[Kc] = sum(lgamma(tn + gam) - lgamma(gam)) + sum(lgamma(ng - tn + phi) - lgamma(phi)) - sum(lgamma(ng + gam + phi) - lgamma(gam + phi)) + log(alpha1) - log(lc+alpha1-1)
      }
      
      lp = exp(lp - max(lp))
      ac[c] = sample(Kc, size=1, prob = lp)

      ac <- as.factor(ac)
      levels(ac) <- c(1:length(unique(ac)))
      
      ac <- as.numeric(ac)
      Kc <- max(ac)
      nc = tabulate(ac, nbins = Kc)
      nx = array(0, dim=c(Kc,Kg));
      
      for(k in unique(ac)){	
        for(l in unique(bg)){	
          nx[k, l] = sum(x[which(ac==k), which(bg==l)])
        }
      }
    }
    
    
    for(i in 1:lg){	
      ng[bg[i]] = ng[bg[i]] - 1
      tn = rep(0, Kc)
      for(j in unique(ac)){	
        tn[j] = sum(x[which(ac == j), i])
      }
      nx[,bg[i]] = nx[,bg[i]] - tn			
      
      lp = rep(-1e6, Kg)
      for(j in unique(bg[-i])){	
        lp[j] = sum(lgamma(tn + nx[,j] + gam) - lgamma(nx[,j] + gam)) + sum(lgamma(nc * (ng[j] + 1) - tn - nx[,j] + phi) - lgamma(nc * ng[j] - nx[,j] + phi)) - sum(lgamma(nc*(ng[j]+1) + gam + phi) - lgamma(nc*ng[j] + gam + phi)) + log(ng[j]) - log(lg+alpha2-1)
      }
      if(ng[bg[i]] == 0){	
        lp[bg[i]] = sum(lgamma(tn + gam) - lgamma(gam)) + sum(lgamma(nc - tn + phi) - lgamma(phi)) - sum(lgamma(nc + gam + phi) - lgamma(gam + phi)) + log(alpha2) - log(lg+alpha2-1)
      }else {
        Kg <- Kg + 1
        lp[Kg] = sum(lgamma(tn + gam) - lgamma(gam)) + sum(lgamma(nc - tn + phi) - lgamma(phi)) - sum(lgamma(nc + gam + phi) - lgamma(gam + phi)) + log(alpha2) - log(lg+alpha2-1)
      }
      
      lp = exp(lp - max(lp))
      bg[i] = sample(Kg, size=1, prob = lp)

      bg <- as.factor(bg)
      levels(bg) <- c(1:length(unique(bg)))
      bg <- as.numeric(bg)
      Kg <- max(bg)
      ng = tabulate(bg, nbins = Kg)
      nx = array(0, dim=c(Kc,Kg))
      for(k in unique(ac)){	
        for(l in unique(bg)){	
          nx[k, l] = sum(x[which(ac==k), which(bg==l)])
        }
      }
    }
  }
  
  rt=NULL
  rt$ac = ac
  rt$bg = bg
  return(rt)
}

BHP_MCMC <- function(y, threshold_ini, I, iter_mle, B_Dir, parg0, parc0, gam, phi, numcellclu, numgeneclu, cell_cluster, logname){
  C <- dim(y)[1]
  G <- dim(y)[2]
  x0 <- matrix(NA, C, G)
  for (c in 1:C) {
    for (g in 1:G) {
      if(y[c,g] < threshold_ini){x0[c,g] <- 0}
      else{x0[c,g] <- 1}
    }
  }
  
  ac_iter <- matrix(NA, I, C)
  bg_iter <- matrix(NA, I, G)
  ARI_c <- matrix(NA, I)
  mle_lik_trace <- matrix(NA, I, iter_mle+1) 
  ac = sample(numcellclu, size=C, replace=T)
  bg = sample(numgeneclu, size=G, replace=T)

  
  for (i in 1:I) {
    re_clu <- scDir(x0, B_Dir, gam, phi, ac, bg)
    ac <- re_clu$ac
    bg <- re_clu$bg
    ac_iter[i,] <- ac
    bg_iter[i,] <- bg
    re_mle <- mle(y,x0,parg0,parc0,iter_mle)
    a0 <- re_mle$a0
    a1 <- re_mle$a1
    b0 <- re_mle$b0
    b1 <- re_mle$b1
    lam <- re_mle$lam
    parg0 <- cbind(a0,a1,b0,b1)
    parc0 <- lam
    resam <- Samx(x0,y,ac,bg,a0,a1,b0,b1,lam,gam,phi)
    x0 <- resam$newx
    mle_lik_trace[i,] <- re_mle$lik_trace
    ARI_c[i] <- adjustedRandIndex(cell_cluster, ac) 
    print(c(i, ARI_c[i]))
    write(c(i, ARI_c[i]),file=logname,append=TRUE)
  }
  
  re <- NULL
  re$ac_iter <- ac_iter
  re$bg_iter <- bg_iter
  re$mle_lik_trace <- mle_lik_trace
  re$ARI_c <- ARI_c
  re$ac <- ac
  re$bg <- bg
  re$a0 <- a0
  re$a1 <- a1
  re$b0 <- b0
  re$b1 <- b1
  re$lam <- lam
  re$x <- x0
  return(re)
}

ini_par <- function(a0,a1,b0,b1,lam,C,G){
  tpar0 <- NULL
  tpar0$a0 <- rep(a0, G)
  tpar0$a1 <- rep(a1, G)
  tpar0$b0 <- rep(b0, G)
  tpar0$b1 <- rep(b1, G)
  tpar0$lam <- rep(lam, C) 
  tpar0$parg0 <- cbind(tpar0$a0,tpar0$a1,tpar0$b0,tpar0$b1)
  return(tpar0)
}

y0_rate <- function(y){
  length(which(y==0))/(dim(y)[1]*dim(y)[2])
}
