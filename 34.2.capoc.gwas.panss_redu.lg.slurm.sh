#!/bin/bash

#SBATCH -J gwas_panssredu_lg
#SBATCH -p cnl
#SBATCH -o out.%J.txt
#SBATCH -e err.%J.txt
#SBATCH -N 2
#SBATCH --ntasks-per-node=5

# 定义表型变量列表
phenos=("PANSS_reduce_rate_res" "PANSS_G_reduce_rate_res" "PANSS_N_reduce_rate_res" "PANSS_P_reduce_rate_res" "LDL_res" "HDL_res" "TC_res" "TG_res" "glu_res")

# 输入输出路径
plink2_path="/home/zhuyunqing/software/plink20/plink2"
bfile="/home/zhuyunqing/863/imputed_clean_qc"
pheno_file="/home/zhuyunqing/drugtarget_lipid/gwas/CAPOC.cis.capoc_cases.grs.panss_redu.lipid.glu.to_gwas.txt"
out_dir="/home/zhuyunqing/drugtarget_lipid/gwas"

# 循环执行每个表型的 GWAS
for pheno in "${phenos[@]}"; do
  $plink2_path \
    --bfile $bfile \
    --pheno $pheno_file \
    --pheno-name $pheno \
    --covar $pheno_file \
    --covar-name PC1 PC2 PC3 PC4 PC5 \
    --glm hide-covar \
    --out "${out_dir}/CAPOC.${pheno}.gwas"
done


