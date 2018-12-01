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
  s1 = Sys.time()
  nx = array(0, dim=c(Kc,Kg))
  for(k in unique(ac)){ 
    for(l in unique(bg)){
      nx[k, l] = sum(x[which(ac==k), which(bg==l)])
    }
  }
  e1 = Sys.time()
  #print('DP (1) time:')
  #print(e1-s1)
  ### 
 
 
  ###################
  ###################
  ###################
  s_cell = Sys.time()
  e_cell = Sys.time()
  t_cell1 = e_cell-s_cell
  t_cell2 = t_cell1
  t_cell3 = t_cell1
  t_cell4 = t_cell1
  t_cell5 = t_cell1
  t_cell6 = t_cell1
  ###### DP for cell (loop each cell)
  s_cell = Sys.time()
  for(c in 1:lc){	
    ### The number of cells in the leave-out-cell's cluster subtract DP removed cell
    nc[ac[c]] = nc[ac[c]] - 1

    ### get the gene number in each gene cluster of the leave-out-cell
    tn = rep(0, Kg)
    for(l in unique(bg)){	
      tn[l] = sum(x[c, which(bg == l)])
    }

    ### The number of genes in the leave-out-cell's genes' clusters subtract DP removed genes
    nx[ac[c],] = nx[ac[c],] - tn
    
    ### initialize likelihood-p
    lp = rep(-1e6, Kc)
    s = Sys.time()
    ### (1) get lp of the leave-out-cell of existing clusters (loop cluster after remove the DP removed cell)
    for(k in unique(ac[-c])){			
      lp[k] = sum(lgamma(tn + nx[k,] + gam) - lgamma(nx[k,] + gam)) + sum(lgamma(ng * (nc[k] + 1) - tn - nx[k,] + phi) - lgamma(ng * nc[k] - nx[k,] + phi)) - sum(lgamma(ng*(nc[k]+1) + gam + phi) - lgamma(ng*nc[k] + gam + phi)) + log(nc[k]) - log(lc+alpha1-1)
    }
    e = Sys.time()
    #print('(1) time:')
    t_cell1 = t_cell1 + (e-s)

    ### (2) get lp of the leave-out-cell of empty clusters
    s = Sys.time()
    if(nc[ac[c]] == 0){	
      ### if the cell cluster has 0 member
      lp[ac[c]] = sum(lgamma(tn + gam) - lgamma(gam)) + sum(lgamma(ng - tn + phi) - lgamma(phi)) - sum(lgamma(ng + gam + phi) - lgamma(gam + phi)) + log(alpha1) - log(lc+alpha1-1)
    }else{ 
      ### add a new cluster 
      Kc <- Kc+1
      lp[Kc] = sum(lgamma(tn + gam) - lgamma(gam)) + sum(lgamma(ng - tn + phi) - lgamma(phi)) - sum(lgamma(ng + gam + phi) - lgamma(gam + phi)) + log(alpha1) - log(lc+alpha1-1)
    }
    e = Sys.time()
    t_cell2 = t_cell2 + (e-s)
    
    ### (3) get the new lp of the leave-out-cell of each cluster
    s = Sys.time()
    lp = exp(lp - max(lp))
    ### get the new cluster of the leave-out-cell by sampling
    new_c = sample(Kc, size=1, prob = lp)
    ac[c] = new_c
    e = Sys.time()
    t_cell3 = t_cell3 + (e-s)

    ### (4) get new cluster list of all cells and convert to factor
    s = Sys.time()
    ac <- as.factor(ac)
    levels(ac) <- c(1:length(unique(ac)))
    e = Sys.time()
    t_cell4 = t_cell4 + (e-s)

    ### (5) update cell cluster numbers
    s = Sys.time()
    ac <- as.numeric(ac)
    Kc <- max(ac)
    nc = tabulate(ac, nbins = Kc)


    nx = array(0, dim=c(Kc,Kg));
    e = Sys.time()
    t_cell5 = t_cell5 + (e-s)
    
    ### (6) update the number of binary label of 1 of each cell-gene cluster-pair
    s = Sys.time()
    acu = unique(ac)
    bgu = unique(bg)
    #print(acu)
    #print(bgu)
    for(k in acu){	
      for(l in bgu){	
        #print(sum(x[which(ac==k), which(bg==l)]))
        nx[k, l] = sum(x[which(ac==k), which(bg==l)])
      }
    }
    e = Sys.time()
    t_cell6 = t_cell6 + (e-s)

    if (c %% 500 == 0){
      print('iteration:')
      print(c)
      print(unique(ac))
    }
    ### update DONE
  }

  e_cell = Sys.time()
  print('DP cell all time:')
  print(e_cell-s_cell)
  print('DP cell (1) time:')
  print(t_cell1)
  print('DP cell (2) time:')
  print(t_cell2)
  print('DP cell (3) time:')
  print(t_cell3)
  print('DP cell (4) time:')
  print(t_cell4)
  print('DP cell (5) time:')
  print(t_cell5)
  print('DP cell (6) time:')
  print(t_cell6)

  ###################
  ###################
  ###################


  ###### DP for gene (loop each gene)
  s_gene = Sys.time()
  for(i in 1:lg){ 
    ### The number of genes in the leave-out-gene's cluster subtract DP removed gene
    ng[bg[i]] = ng[bg[i]] - 1

    ### get the cell number in each cell cluster of the leave-out-gene
    tn = rep(0, Kc)
    for(j in unique(ac)){ 
      tn[j] = sum(x[which(ac == j), i])
    }

    ### The number of cell in the leave-out-gene's cells' clusters subtract DP removed cells 
    nx[,bg[i]] = nx[,bg[i]] - tn      
    
    ### initialize likelihood-p
    lp = rep(-1e6, Kg)
    ### get lp of the leave-out-gene of existing clusters (loop cluster after remove the DP removed gene)
    for(j in unique(bg[-i])){ 
      lp[j] = sum(lgamma(tn + nx[,j] + gam) - lgamma(nx[,j] + gam)) + sum(lgamma(nc * (ng[j] + 1) - tn - nx[,j] + phi) - lgamma(nc * ng[j] - nx[,j] + phi)) - sum(lgamma(nc*(ng[j]+1) + gam + phi) - lgamma(nc*ng[j] + gam + phi)) + log(ng[j]) - log(lg+alpha2-1)
    }

    ### get lp of the leave-out-gene of empty clusters
    if(ng[bg[i]] == 0){ 
      ### if the gene cluster has 0 member
      lp[bg[i]] = sum(lgamma(tn + gam) - lgamma(gam)) + sum(lgamma(nc - tn + phi) - lgamma(phi)) - sum(lgamma(nc + gam + phi) - lgamma(gam + phi)) + log(alpha2) - log(lg+alpha2-1)
    }else {
      ### add a new cluster 
      Kg <- Kg + 1
      lp[Kg] = sum(lgamma(tn + gam) - lgamma(gam)) + sum(lgamma(nc - tn + phi) - lgamma(phi)) - sum(lgamma(nc + gam + phi) - lgamma(gam + phi)) + log(alpha2) - log(lg+alpha2-1)
    }
    
    ### get the new lp of the leave-out-gene of each cluster
    lp = exp(lp - max(lp))
    ### get the new cluster of the leave-out-gene by sampling
    bg[i] = sample(Kg, size=1, prob = lp)

    ### get new cluster list of all genes and convert to factor
    bg <- as.factor(bg)
    levels(bg) <- c(1:length(unique(bg)))

    ### update gene cluster numbers
    bg <- as.numeric(bg)
    Kg <- max(bg)
    ng = tabulate(bg, nbins = Kg)
    nx = array(0, dim=c(Kc,Kg))

    ### update the number of binary label of 1 of each cell-gene cluster-pair
    acu = unique(ac) ### get uniq cell clusters
    bgu = unique(bg) ### get uniq gene clusters
    ### count the number of cell-gene-clusters
    for(k in acu){ 
      for(l in bgu){ 
        nx[k, l] = sum(x[which(ac==k), which(bg==l)])
      }
    }
    ###
  }
  e_gene = Sys.time()
  print('DP gene time:')
  print(e_gene-s_gene)  


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
re_clu = scDir(x, 1, 1, 1, ac_ini, bg_ini)
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






