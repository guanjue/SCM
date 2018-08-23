#PBS -l nodes=2:ppn=8
#PBS -l walltime=500:00:00
#PBS -j oe
#PBS -A yzz2_e_g_sc_default

module load gcc
module load python/2.7.14-anaconda5.0.1

cd /storage/home/gzx103/scratch/sc
#rm pbs_check.txt
#time Rscript R_code.R

#rm Uso_filter10_seed1_iter10000.log.txt
#time Rscript R_code_inter.R 100 100000

time Rscript R_code_inter.R Usoskin_filter10.txt Usoskin_clu_t.txt /storage/home/g/gzx103/scratch/sc/MCMC6_noeta_ini_inter.R Uso_filter10_seed 100 10000 128

