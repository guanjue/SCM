rm(list=ls(all=TRUE))
args = commandArgs(trailingOnly=TRUE)
input_mat_rowCT_colGENE = args[1]
True_label_list = args[2]
source_script = args[3]
output = args[4]
random_seed = as.numeric(args[5])
iteration_num = as.numeric(args[6])
initial_thresh = as.numeric(args[7])

set.seed(as.numeric(random_seed))

#source("/storage/home/g/gzx103/scratch/sc/MCMC6_noeta_ini_inter.R")
source(source_script)

#y <- read.table(file = "/storage/home/g/gzx103/scratch/sc/Usoskin_filter10.txt")
y <- read.table(file = input_mat_rowCT_colGENE)
C <- dim(y)[1]
G <- dim(y)[2]
y <- as.numeric(unlist(y))
y <- matrix(y, C, G)
#cell_cluster <- read.table(file = "/storage/home/g/gzx103/scratch/sc/Usoskin_clu_t.txt")
cell_cluster <- read.table(file = True_label_list)
cell_cluster <- as.numeric(unlist(cell_cluster))
cell_cluster <- as.vector(cell_cluster)

par0 <- ini_par(1,1,10,10,1,C,G)

y0_rate(y)
iteration = as.numeric(iteration_num)
#logname = paste('Uso_filter10_seed', random_seed, '_iter', toString(iteration), '.log.txt', sep='')
logname = paste(output, random_seed, '_iter', toString(iteration), '.log.txt', sep='')
resim <- BHP_MCMC(y,initial_thresh,iteration,1,1,par0$parg0,par0$lam,1,1,10,10,cell_cluster, logname)

#save(file = paste('Uso_filter10_seed', random_seed, '_iter', toString(iteration), '.Rdata', sep=''), resim)
save(file = paste(output, random_seed, '_iter', toString(iteration), '.Rdata', sep=''), resim)
