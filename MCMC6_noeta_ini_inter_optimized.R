library('Rcpp')
library("mclust")
sourceCpp("/storage/home/gzx103/group/software/SCM/funs_rcpp5_noeta.cpp") 


############################################################################################################

get_nloglik <- function(y,matx,parc_all,parg_all){
  nloglik <- 0
  for (i in 1:dim(y)[1]){
    nloglik <- nloglik + get_nloglik_c(parc_all[i],parg_all, y[i,], matx[i,])
  }
  return(nloglik)
}

mle_g = function(x, parc0){
  nlminb(start = x[1:4], objective = get_nloglik_g, parc_all= parc0, y_g=x[5:1924], matx_g=x[1925:3844], lower = c(0.00000001,0.00000001,0.00000001,0.00000001), upper = c(1000,1000,10000,10000))
}

############################################################################################################

mle <- function(y,x,parg0,parc0,iter){
  lik_trace <- get_nloglik(y,x,parc0,parg0)
  for (i in 1:iter){
    nc <- dim(y)[1]
    ng <- dim(y)[2]
    eparc <- matrix(NA,nc)
    eparg <- matrix(NA,nrow=dim(parg0)[1], ncol=dim(parg0)[2])
    
    for (ic in 1:nc){
      #temp <- nlminb(start = parc0[ic], objective = get_nloglik_c, parg_all= parg0, y_c=y[ic,], matx_c=x[ic,], lower = 0.0000001, upper = 1000000)
      temp <- optim(par = parc0[ic], fn = get_nloglik_c, parg_all= parg0, y_c=y[ic,], matx_c=x[ic,], lower = 0.0000001, upper = 1000000)
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

############################################################################################################

scDir<-function(x, B, gam_fix, phi_up, ac_ini, bg_ini){ 
  lc = dim(x)[1]
  lg = dim(x)[2]
  alpha1 = 1 ### fix 1
  alpha2 = 1 ### fix 1
  gam = 1 ### fix 1
  phi = 1 ### fix 1
  ac = ac_ini
  bg = bg_ini
  Kc = max(ac)
  Kg = max(bg)
  ### get the number of each cell / gene cluster
  nc = tabulate(ac, nbins = Kc)
  ng = tabulate(bg, nbins = Kg)
  ### (1) get the number of binary label of 1 of each cell-gene cluster-pair
  nx = array(0, dim=c(Kc,Kg))
  #print(ac)
  for(k in unique(ac)){ 
    for(l in unique(bg)){
      nx[k, l] = sum(x[which(ac==k), which(bg==l)])
    }
  }
  ################## 
  ###### DP for cell (loop each cell)
  for(c in 1:lc){   
    ### The number of cells in the leave-out-cell's cluster subtract DP removed cell
    rm_c = 0
    leave_out_cell_c = ac[c]
    ### subtract 1 for the leave-one-out cell's cluster
    nc_leave_out_cell_c = nc[leave_out_cell_c] - 1
    nc[leave_out_cell_c] = nc_leave_out_cell_c
    ### get the gene number in each gene cluster of the leave-out-cell
    tn = rep(0, Kg)
    for(l in unique(bg)){   
      tn[l] = sum(x[c, which(bg == l)])
    }
    ### The number of genes in the leave-out-cell's genes' clusters subtract DP removed genes
    nx[leave_out_cell_c,] = nx[leave_out_cell_c,] - tn    
    ### initialize likelihood-p
    lp = rep(-1e6, Kc)
    ### (1) get lp of the leave-out-cell of existing clusters (loop cluster after remove the DP removed cell)
    for(k in unique(ac[-c])){           
      lp[k] = sum(lgamma(tn + nx[k,] + gam) - lgamma(nx[k,] + gam)) + sum(lgamma(ng * (nc[k] + 1) - tn - nx[k,] + phi) - lgamma(ng * nc[k] - nx[k,] + phi)) - sum(lgamma(ng*(nc[k]+1) + gam + phi) - lgamma(ng*nc[k] + gam + phi)) + log(nc[k]) - log(lc+alpha1-1)
    }
    ### (2) get lp of the leave-out-cell of empty clusters
    if(nc[leave_out_cell_c] == 0){ 
      rm_c=1
      ### if the cell cluster has 0 member
      lp[leave_out_cell_c] = sum(lgamma(tn + gam) - lgamma(gam)) + sum(lgamma(ng - tn + phi) - lgamma(phi)) - sum(lgamma(ng + gam + phi) - lgamma(gam + phi)) + log(alpha1) - log(lc+alpha1-1)
    }else{ 
      ### add a new cluster 
      Kc <- Kc+1
      lp[Kc] = sum(lgamma(tn + gam) - lgamma(gam)) + sum(lgamma(ng - tn + phi) - lgamma(phi)) - sum(lgamma(ng + gam + phi) - lgamma(gam + phi)) + log(alpha1) - log(lc+alpha1-1)
    }
    ### (3) get the new lp of the leave-out-cell of each cluster
    #lp = exp(lp - max(lp))### get the new cluster of the leave-out-cell by sampling #new_c = sample(Kc, size=1, prob = lp)
    ### replacing sampling by maximization
    new_c = which.max(lp)
    ### update new cluster id
    ac[c] = new_c
    ### if the leave-one-out cell go back to the original cell cluster then not removing
    if ((rm_c == 1) & (new_c==leave_out_cell_c)){
        rm_c=0
    }
    ### (4) get new cluster list of all cells and convert to factor
    ### check if add new cluster
    if (new_c>dim(nx)[1]){
        nx = rbind(nx, tn)
        nc = c(nc, 1)
    } else if (nc_leave_out_cell_c!=0) {
        Kc = Kc-1
        nx[new_c, ] = nx[new_c, ] + tn
        nc[new_c] = nc[new_c]+1
    } else {
        nx[new_c, ] = nx[new_c, ] + tn
        nc[new_c] = nc[new_c]+1      
    }

    ### check if remove cell
    if (rm_c==1){
        Kc = Kc - 1
        nx = nx[-leave_out_cell_c,]
        nc = nc[-leave_out_cell_c]
        ac[ac>leave_out_cell_c] = ac[ac>leave_out_cell_c]-1
    }
    ### update DONE
  ac_new = ac
  }
  ##################
  ###### DP for gene (loop each gene)
  for(i in 1:lg){ 
    rm_g_c = 0
    leave_out_gene_c = bg[i]
    ### The number of genes in the leave-out-gene's cluster subtract DP removed gene
    ng_leave_out_gene_c = ng[leave_out_gene_c] - 1
    ng[leave_out_gene_c] = ng_leave_out_gene_c
    ### get the cell number in each cell cluster of the leave-out-gene
    tn = rep(0, Kc)
    for(j in unique(ac)){ 
      tn[j] = sum(x[which(ac == j), i])
    }
    ### The number of cell in the leave-out-gene's cells' clusters subtract DP removed cells 
    nx[,leave_out_gene_c] = nx[,leave_out_gene_c] - tn         
    ### initialize likelihood-p
    lp = rep(-1e6, Kg)
    ### get lp of the leave-out-gene of existing clusters (loop cluster after remove the DP removed gene)
    for(j in unique(bg[-i])){ 
      lp[j] = sum(lgamma(tn + nx[,j] + gam) - lgamma(nx[,j] + gam)) + sum(lgamma(nc * (ng[j] + 1) - tn - nx[,j] + phi) - lgamma(nc * ng[j] - nx[,j] + phi)) - sum(lgamma(nc*(ng[j]+1) + gam + phi) - lgamma(nc*ng[j] + gam + phi)) + log(ng[j]) - log(lg+alpha2-1)
    }
    ### get lp of the leave-out-gene of empty clusters
    if(ng[leave_out_gene_c] == 0){ 
      rm_g_c = 1
      ### if the gene cluster has 0 member
      lp[leave_out_gene_c] = sum(lgamma(tn + gam) - lgamma(gam)) + sum(lgamma(nc - tn + phi) - lgamma(phi)) - sum(lgamma(nc + gam + phi) - lgamma(gam + phi)) + log(alpha2) - log(lg+alpha2-1)
    }else {
      ### add a new cluster 
      Kg <- Kg + 1
      lp[Kg] = sum(lgamma(tn + gam) - lgamma(gam)) + sum(lgamma(nc - tn + phi) - lgamma(phi)) - sum(lgamma(nc + gam + phi) - lgamma(gam + phi)) + log(alpha2) - log(lg+alpha2-1)
    }   
    ### get the new lp of the leave-out-gene of each cluster
    new_g_c = which.max(lp)
    bg[i] = new_g_c
    ### if the leave-one-out gene go back to the original gene cluster then not removing
    if ((rm_g_c == 1) & (new_g_c==leave_out_gene_c)){
        rm_g_c=0
    }
    ### check if add new cluster
    if (new_g_c>dim(nx)[2]){
        nx = cbind(nx, tn)
        ng = c(ng, 1)
    } else if (ng_leave_out_gene_c != 0){
        Kg = Kg-1
        nx[,new_g_c] = nx[,new_g_c] + tn
        ng[new_g_c] = ng[new_g_c]+1
    } else {
        nx[,new_g_c] = nx[,new_g_c] + tn
        ng[new_g_c] = ng[new_g_c]+1
    }
    ### check if remove cell
    if (rm_g_c==1){
        Kg = Kg - 1
        nx = nx[,-leave_out_gene_c]
        ng = ng[-leave_out_gene_c]
        bg[bg>leave_out_gene_c] = bg[bg>leave_out_gene_c]-1
    }
    ###
  }
  ### return new cell/gene clusters
  rt=NULL
  rt$ac = ac
  rt$bg = bg
  return(rt)
}

############################################################################################################

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
      a0_g = a0[g]
      a1_g = a1[g]
      b0_g = b0[g]
      b1_g = b1[g]
      lam_c = lam[c]
      y_c_g = y[c,g]
      if(y_c_g!=0) {
        logp_cg <- lgamma(a0_g+y_c_g) + lgamma(a0_g+a1_g) + y_c_g*log(b0_g) + log(phi+sizex-1-sx) - (a0_g+y_c_g)*log(b0_g*1+1) + log(1 - ((b0_g*1+1)/(b0_g*(lam_c+1)+1))^(a0_g+y_c_g)) - lgamma(a0_g) - lgamma(a0_g+a1_g+y_c_g) - y_c_g*log(b0_g+b1_g) - log(gam+sx) + (a0_g+a1_g+y_c_g)*log((b0_g+b1_g)*1+1) - log(1 - (((b0_g+b1_g)*1+1)/((b0_g+b1_g)*(1+lam_c)+1))^(a0_g+a1_g+y_c_g)) 
      }else{
        if(lam_c>1){
          logp_cg <- log(phi+sizex-1-sx) - log(gam+sx) - a0_g*log(b0_g*1+1) + log(1 - ((b0_g*1+1)/(b0_g*(lam_c+1)+1))^a0_g + ((b0_g*1+1)/(b0_g*lam_c+1))^a0_g) + (a0_g+a1_g)*log((b0_g+b1_g)*1+1) - log(1 - (((b0_g+b1_g)*1+1)/((b0_g+b1_g)*(1+lam_c)+1))^(a0_g+a1_g) + (((b0_g+b1_g)*1+1)/((b0_g+b1_g)*lam_c+1))^(a0_g+a1_g))
        }else{
          logp_cg <- log(phi+sizex-1-sx) - log(gam+sx) - a0_g*log(b0_g*lam_c+1) + log(((b0_g*lam_c+1)/(b0_g*1+1))^a0_g - ((b0_g*lam_c+1)/(b0_g*(lam_c+1)+1))^a0_g + 1) + (a0_g+a1_g)*log((b0_g+b1_g)*lam_c+1) - log((((b0_g+b1_g)*lam_c+1)/((b0_g+b1_g)*1+1))^(a0_g+a1_g) -  (((b0_g+b1_g)*lam_c+1)/((b0_g+b1_g)*(1+lam_c)+1))^(a0_g+a1_g) + 1)
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

############################################################################################################

BHP_MCMC <- function(y, threshold_ini, I, iter_mle, B_Dir, parg0, parc0, gam, phi, numcellclu, numgeneclu, cell_cluster, logname){
  ### get cell /gene number 
  C <- dim(y)[1]
  G <- dim(y)[2]
  
  ### initialize binary matrix based on a given threshold
  x0 <- matrix(NA, C, G)
  x0[y < threshold_ini] = 0
  x0[y >= threshold_ini] = 1

  ### initialize output matrices
  ac_iter <- matrix(NA, I, C)
  bg_iter <- matrix(NA, I, G)
  ARI_c <- matrix(NA, I)
  mle_lik_trace <- matrix(NA, I, iter_mle+1) 

  ### initialized cell cluster
  ac = sample(numcellclu, size=C, replace=T)
  ### initialized gene cluster
  bg = sample(numgeneclu, size=G, replace=T)

  ### training
  for (i in 1:I) {
    ###### Dirichlet Process (DP)
    ### DP
    re_clu <- scDir(x0, B_Dir, gam, phi, ac, bg)
    ### update cell / gene clusters after DP
    ac <- re_clu$ac
    bg <- re_clu$bg
    ### add to output matrix
    ac_iter[i,] <- ac
    bg_iter[i,] <- bg

    ###### MLE
    ### MLE learning parameters
    re_mle <- mle(y,x0,parg0,parc0,iter_mle)
    ### update parameters after MLE
    a0 <- re_mle$a0
    a1 <- re_mle$a1
    b0 <- re_mle$b0
    b1 <- re_mle$b1
    lam <- re_mle$lam
    parg0 <- cbind(a0,a1,b0,b1)
    parc0 <- lam

    ###### update Binary labels
    ### Samx
    resam <- Samx(x0,y,ac,bg,a0,a1,b0,b1,lam,gam,phi)
    x0 <- resam$newx

    ###### calculate Adjusted Random Index
    ARI_c[i] <- adjustedRandIndex(cell_cluster, ac) 
    print(c(i, ARI_c[i]))
    write(c(i, ARI_c[i]),file=logname,append=TRUE)

    ###### Save trace
    mle_lik_trace[i,] <- re_mle$lik_trace
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

############################################################################################################

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

############################################################################################################

y0_rate <- function(y){
  length(which(y==0))/(dim(y)[1]*dim(y)[2])
}

############################################################################################################





