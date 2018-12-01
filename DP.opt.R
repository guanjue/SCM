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
    nc[leave_out_cell_c] = nc[leave_out_cell_c] - 1
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
    #lp = exp(lp - max(lp))
    ### get the new cluster of the leave-out-cell by sampling
    #new_c = sample(Kc, size=1, prob = lp)
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
    } else {
        Kc = Kc-1
        nx[new_c, ] = nx[new_c, ] + tn
        nc[new_c] = nc[new_c]+1
    }
    ### check if remove cell
    if (rm_c==1){
        nx = nx[-leave_out_cell_c,]
        nc = nc[-leave_out_cell_c]
        ac[ac>leave_out_cell_c] = ac[ac>leave_out_cell_c]-1
    }
    #if (c %% 500 == 0){
    #  print('iteration:')
    #  print(c)
    #  print(unique(ac))
    #}
    ### update DONE
  ac_new = ac
  }
  ##################
  ###### DP for gene (loop each gene)
  for(i in 1:lg){ 
    rm_g_c = 0
    leave_out_gene_c = bg[i]
    ### The number of genes in the leave-out-gene's cluster subtract DP removed gene
    ng[leave_out_gene_c] = ng[leave_out_gene_c] - 1
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
    #lp = exp(lp - max(lp))
    ### get the new cluster of the leave-out-gene by sampling
    #bg[i] = sample(Kg, size=1, prob = lp)
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
    } else {
        Kg = Kg-1
        nx[,new_g_c] = nx[,new_g_c] + tn
        ng[new_g_c] = ng[new_g_c]+1
    }
    ### check if remove cell
    if (rm_g_c==1){
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


set.seed(2018)
x = sample(c(0,1), replace=TRUE, size=2000)
for (i in 1:50){
  x = cbind(x, sample(c(0,1), replace=TRUE, size=2000))
}
ac_ini = sample(c(1:10), replace=TRUE, size=dim(x)[1])
bg_ini = sample(c(1:10), replace=TRUE, size=dim(x)[2])


s0 = Sys.time()
#re_clu = scDir(x, 1, 1, 1, ac_ini, bg_ini)
e0 = Sys.time()
print('All time:')
print(e0-s0)




set.seed(2018)
c_num = 2000
x = sample(c(0,1), replace=TRUE, size=c_num)
for (i in 1:19){
  x = cbind(x, sample(c(0,1), replace=TRUE, size=2000))
}

x0 = read.table('/Users/universe/Downloads/test.bimat.txt')
x = cbind(x, x0[c(1:c_num),-c(1:4)])

ac_ini = sample(c(1:10), replace=TRUE, size=dim(x)[1])
bg_ini = sample(c(1:10), replace=TRUE, size=dim(x)[2])


s0 = Sys.time()
re_clu = scDir(x, 1, 1, 1, ac_ini, bg_ini)
e0 = Sys.time()
print('All time:')
print(e0-s0)




