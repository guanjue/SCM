### get parameters
args = commandArgs(trailingOnly=TRUE)

library(SingleCellExperiment)
library(SC3)
library(mclust)

counts_mat = args[1]
known_labels = args[2]
output_name = args[3]

counts_mat = 'GSE57872_GBM_data_matrix.sc.cvfiltered.txt'
counts_mat = 'Usoskin_filter10.txt'
counts_mat = 'GSE51372_readCounts.sc.cvfiltered.txt'
counts_mat = 'deng2014.sc.cvfiltered.txt'

counts_mat = read.table(counts_mat)
counts_mat = t(counts_mat)
#counts_mat = log2(counts_mat+1)
### run SingleCellExperiment
sce = SingleCellExperiment(assays = list(counts = counts_mat, logcounts = log2(counts_mat+1)))

### read real labels
known_labels = "sc_label.txt"
known_labels = "Usoskin_clu_t.txt"
known_labels = "ct_labels.txt"
known_labels = "sc_label.txt"

cell_cluster = read.table(file = known_labels)
cell_cluster = t(cell_cluster)

rowData(sce)$feature_symbol = rownames(counts_mat)
sce$feature_symbol = colnames(counts_mat)

k = sc3_estimate_k(sce)

### sc3
set.seed(2018)

sc3_re = sc3(sce, ks = 8, gene_filter = FALSE)
sc3_re = sc3(sce, ks = 11, gene_filter = FALSE)
sc3_re = sc3(sce, ks = 7, gene_filter = FALSE)
sc3_re = sc3(sce, ks = 10, gene_filter = FALSE)

col_data = colData(sc3_re)

cell_sc3 = as.vector(col_data$sc3_8_clusters)
cell_sc3 = as.vector(col_data$sc3_11_clusters)
cell_sc3 = as.vector(col_data$sc3_7_clusters)
cell_sc3 = as.vector(col_data$sc3_10_clusters)


sc3_dri = adjustedRandIndex(cell_cluster, cell_sc3) 

sc3_dri

write.table(sc3_dri, 'sc3_dri.txt', quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')
