vi extract.lipid.target.glgc.slurm

#!/bin/bash

#SBATCH -J extract_lipid_glgc
#SBATCH -p cnl
#SBATCH -o out.%J.txt
#SBATCH -e err.%J.txt
#SBATCH -N 2
#SBATCH --ntasks-per-node=5


cd /home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/
Rscript extract.lipid.target.glgc.R 2>&1 | tee extract_output.log





vi extract.lipid.target.glgc.R

library(data.table)
setwd("/home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/")

# 读取bim文件并构建marker列
capoc.snp <- fread("/home/share1/check/863_qc/imputed_clean_qc.bim", header = FALSE)
capoc.snp <- capoc.snp[, .(chr = V1, rsid = V2, bp = V4)]
capoc.snp[, marker := paste0(chr, ":", bp)]
capoc.snp <- capoc.snp[, .(marker)]

# 基因信息表
gene_info <- data.table(
  gene  = c("ITGAL", "AHR", "PPARA", "PPARD", "PPARG", "NR1I2", "RXRA", "RXRB", "RXRG", "SLCO1B1", "SLCO2B1", "MMP25", "HCAR3", "HCAR2", "THRA", "ABCA1", "MTTP", "PCSK9", "ACLY", "ANGPTL3", "APOC3", "CETP", "APOB", "ANGPTL4", "LDLR", "LPL", "LPA", "HMGCR", "HDAC2", "DPP4", "SLC22A8", "SLCO1B3", "THRB", "NPC1L1", "SOAT1", "THRA"),
  chr   = c(16, 7, 22, 6, 3, 3, 9, 6, 1, 12, 11, 16, 12, 12, 17, 9, 4, 1, 17, 1, 11, 16, 2, 19, 19, 8, 6, 5, 6, 2, 11, 12, 3, 7, 1, 17),
  start = c(30484063, 17338276, 46546429, 35310335, 12328867, 119500948, 137218301, 33161365, 165370159, 21284128, 74862152, 3096562, 123199303, 123185840, 38218446, 107543287, 100485287, 55505221, 40023170, 63063191, 116700623, 56995862, 21224301, 8429039, 11200139, 19796764, 160952514, 74632993, 114254192, 162848755, 62760296, 20963639, 24158644, 44552134, 179262932, 38218446),
  end   = c(30534506, 17385771, 46639652, 35395955, 12475843, 119537334, 137332431, 33168630, 165414363, 21392730, 74917594, 3110727, 123201358, 123187904, 38250120, 107690436, 100545154, 55530525, 40075275, 63071984, 116703788, 57017757, 21266945, 8439254, 11244496, 19824770, 161085307, 74657941, 114292312, 162930725, 62783313, 21069845, 24537199, 44580929, 179327815, 38250120)
)

# 表型数据路径及一次性读取
traits <- list(
  LDL   = "deduplicated_LDL_GLGC_EAS_2021_hg19.txt",
  HDL   = "deduplicated_HDL_GLGC_EAS_2021_hg19.txt",
  logTG = "deduplicated_logTG_GLGC_EAS_2021_hg19.txt",
  TC    = "deduplicated_TC_GLGC_EAS_2021_hg19.txt"
)

base_path <- "/home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/"

# 一次性读取所有表型数据，按名称存入list
trait_data <- lapply(traits, function(file) fread(file.path(base_path, file)))

# 循环处理每个基因和每个表型
for (i in 1:nrow(gene_info)) {
  gene     <- gene_info[i, gene]
  chr_i     <- gene_info[i, chr]        # <-- 改名
  region_start <- gene_info[i, start] - 1e5
  region_end   <- gene_info[i, end]   + 1e5

  for (trait in names(trait_data)) {
    dat <- trait_data[[trait]]
    dat[, chr := as.numeric(chr)]
    dat[, bp19 := as.numeric(bp19)]
    dat[, p    := as.numeric(p)]

    # 现在用 chr == chr_i 去比对
    dat1 <- dat[
      chr   == chr_i &
      bp19 >= region_start &
      bp19 <= region_end &
      p    <= 5e-8
    ]


    dat1[, marker := paste0(chr, ":", bp19)]
    m <- merge(dat1, capoc.snp, by = "marker", all = FALSE)

    outfile <- file.path(base_path, paste0(gene, ".", trait, ".SNP.txt"))
    write.table(m, outfile, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  }
}










vi extract.glu.target.taiwan.slurm

#!/bin/bash

#SBATCH -J extract_glu_taiwan
#SBATCH -p cnl
#SBATCH -o out.%J.txt
#SBATCH -e err.%J.txt
#SBATCH -N 2
#SBATCH --ntasks-per-node=5


cd /home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/
Rscript extract.glu.target.taiwan.R 2>&1 | tee extract.glu.target.taiwan.log





vi extract.glu.target.taiwan.R

library(data.table)
setwd("/home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/")

# 读取bim文件并构建marker列
capoc.snp <- fread("/home/share1/check/863_qc/imputed_clean_qc.bim", header = FALSE)
capoc.snp <- capoc.snp[, .(chr = V1, rsid = V2, bp = V4)]
capoc.snp[, marker := paste0(chr, ":", bp)]
capoc.snp <- capoc.snp[, .(marker)]

# 基因信息表
gene_info <- data.table(
  gene  = c("ABCA1", "ABCB11", "ABCC8", "ACACB", "ACSL4", "AKR1B1", "AMY2A", "CALCR", "CFTR", "CPT1A", "DPP4", "ESRRA", "ETFDH", "GAA", "GANAB", "GANC", "GCK", "GIPR", "GLP1R", "GPD1", "HRH1", "IGF1R", "IGFBP7", "INSR", "KCNJ1", "KCNJ11", "KCNJ8", "MAOB", "MGAM", "PPARA", "PPARG", "PRKAA1", "PRKAB1", "PRKAG1", "RAMP1", "RAMP2", "RAMP3", "SERPINE1", "SI", "SLC29A1", "SLC5A1", "SLC5A2", "TRPM4"),
  chr   = c(9, 2, 11, 12, 23, 7, 1, 7, 7, 11, 2, 11, 4, 17, 11, 15, 7, 19, 6, 12, 3, 15, 4, 19, 11, 11, 12, 23, 7, 22, 3, 5, 12, 12, 2, 17, 7, 7, 3, 6, 22, 16, 19),
  start = c(107543287, 169777291, 17414045, 109554392, 108884564, 134127102, 104160049, 93053798, 117120079, 68522088, 162848755, 64072996, 159593448, 78075380, 62392301, 42565685, 44182812, 46171479, 39016557, 50497791, 11178924, 99191768, 57896939, 7112276, 128707915, 17406795, 21917889, 43625857, 141695679, 46546429, 12328867, 40759481, 120105938, 49396057, 238768266, 40913245, 45197390, 100770385, 164696686, 44187352, 32439248, 31494444, 49661049),
  end   = c(107690436, 169887834, 17498392, 109706031, 108976486, 134143991, 104168402, 93204042, 117308719, 68609384, 162930725, 64084215, 159630775, 78093680, 62414085, 42645864, 44229038, 46186980, 39059079, 50505096, 11305243, 99507759, 57976551, 7294425, 128737191, 17410893, 21927640, 43741696, 141806547, 46639652, 12475843, 40798297, 120119424, 49412559, 238820748, 40915059, 45223849, 100782528, 164796284, 44201879, 32509016, 31502090, 49715093)
)

# 表型数据路径及一次性读取
traits <- list(
  glu = "deduplicated_glucose_taiwan_summ_hg19.txt"
)

base_path <- "/home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/"

# 一次性读取所有表型数据，按名称存入list
trait_data <- lapply(traits, function(file) fread(file.path(base_path, file)))

# 循环处理每个基因和每个表型
for (i in 1:nrow(gene_info)) {
  gene <- gene_info[i, gene]
  chr <- gene_info[i, chr]
  region_start <- gene_info[i, start] - 100000
  region_end   <- gene_info[i, end] + 100000

  for (trait in names(trait_data)) {
    dat <- trait_data[[trait]]  # 从缓存中取出该表型数据
    dat[, chr := as.numeric(chr)]
    dat[, bp19 := as.numeric(bp19)]
    dat[, P := as.numeric(P)]
    dat1 <- dat[chr == chr & bp19 >= region_start & bp19 <= region_end & P <= 5e-8]
    dat1[, marker := paste0(chr, ":", bp19)]
    m <- merge(dat1, capoc.snp, by = "marker", all = FALSE)

    outfile <- file.path(base_path, paste0(gene, ".", trait, ".SNP.txt"))
    write.table(m, outfile, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  }
}











##########
#clumping

#!/bin/bash



cd /home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/
vi clump.lipid.target.glgc.slurm

#!/bin/bash

#SBATCH -J clump_lipid_glgc
#SBATCH -p cnl
#SBATCH -o out.%J.txt
#SBATCH -e err.%J.txt
#SBATCH -N 2
#SBATCH --ntasks-per-node=5


# 定义基因名列表
genes=("ITGAL" "AHR" "PPARA" "PPARD" "PPARG" "NR1I2" "RXRA" "RXRB" "RXRG" "SLCO1B1" "SLCO2B1" "MMP25" "HCAR3" "HCAR2" "THRA" "ABCA1" "MTTP" "PCSK9" "ACLY" "ANGPTL3" "APOC3" "CETP" "APOB" "ANGPTL4" "LDLR" "LPL" "LPA" "HMGCR" "HDAC2" "DPP4" "SLC22A8" "SLCO1B3" "THRB" "NPC1L1" "SOAT1" "THRA")  # 这里添加你想处理的所有基因

# 定义 trait 名列表
traits=("LDL" "HDL" "logTG" "TC")

# 循环处理每个基因和每个 trait
for gene in "${genes[@]}"; do
    for trait in "${traits[@]}"; do

        input_file="${gene}.${trait}"
        clump_file="/home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/${input_file}.SNP.txt"
        output_file="/home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/${input_file}.snp"

        echo "Running PLINK for ${input_file}..."

        /home/zhuyunqing/software/plink19/plink \
            --bfile /home/share1/Public_database/1000G_plink/1000G_EAS/1000G_merge_allchr_EAS/ALL \
            --clump "$clump_file" \
            --clump-p1 5e-8 \
            --clump-r2 0.1 \
            --clump-field p \
            --clump-snp-field name \
            --clump-range /home/share1/Public_database/glist-hg19.txt \
            --out "$output_file"

        echo "Processed ${input_file}"
    done
done


#有基因位点的药靶
PPARG.HDL
HCAR3.HDL
HCAR2.HDL
ABCA1.LDL
ABCA1.HDL
ABCA1.TC
MTTP.HDL
MTTP.TC
PCSK9.LDL
PCSK9.TC
ANGPTL3.logTG
ANGPTL3.TC
APOC3.HDL
APOC3.logTG
APOC3.TC
CETP.HDL
CETP.logTG
CETP.TC
APOB.LDL
APOB.HDL
APOB.logTG
APOB.TC
LDLR.LDL
LDLR.HDL
LDLR.TC
LPL.HDL
LPL.logTG
LPA.LDL
LPA.HDL
LPA.TC
HMGCR.LDL
HMGCR.TC






cd /home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/
vi clump.glu.target.taiwan.slurm

#!/bin/bash

#SBATCH -J clump_glu_taiwan
#SBATCH -p cnl
#SBATCH -o out.%J.txt
#SBATCH -e err.%J.txt
#SBATCH -N 2
#SBATCH --ntasks-per-node=5


# 定义基因名列表
genes=("ABCA1" "ABCB11" "ABCC8" "ACACB" "ACSL4" "AKR1B1" "AMY2A" "CALCR" "CFTR" "CPT1A" "DPP4" "ESRRA" "ETFDH" "GAA" "GANAB" "GANC" "GCK" "GIPR" "GLP1R" "GPD1" "HRH1" "IGF1R" "IGFBP7" "INSR" "KCNJ1" "KCNJ11" "KCNJ8" "MAOB" "MGAM" "PPARA" "PPARG" "PRKAA1" "PRKAB1" "PRKAG1" "RAMP1" "RAMP2" "RAMP3" "SERPINE1" "SI" "SLC29A1" "SLC5A1" "SLC5A2" "TRPM4")  # 这里添加你想处理的所有基因

# 定义 trait 名列表
traits=("glu")

# 循环处理每个基因和每个 trait
for gene in "${genes[@]}"; do
    for trait in "${traits[@]}"; do

        input_file="${gene}.${trait}"
        clump_file="/home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/${input_file}.SNP.txt"
        output_file="/home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/${input_file}.snp"

        echo "Running PLINK for ${input_file}..."

        /home/zhuyunqing/software/plink19/plink \
            --bfile /home/share1/Public_database/1000G_plink/1000G_EAS/1000G_merge_allchr_EAS/ALL \
            --clump "$clump_file" \
            --clump-p1 5e-8 \
            --clump-r2 0.1 \
            --clump-field p_pheno \
            --clump-snp-field name \
            --clump-range /home/share1/Public_database/glist-hg19.txt \
            --out "$output_file"

        echo "Processed ${input_file}"
    done
done


#######
ABCB11.glu
GCK.glu
GLP1R.glu








################
#整理clumping的结果

cd /home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/
vi clump.lipid.target.glgc.res.slurm

#!/bin/bash

#SBATCH -J clumpres_lipid_glgc
#SBATCH -p cnl
#SBATCH -o out.%J.txt
#SBATCH -e err.%J.txt
#SBATCH -N 2
#SBATCH --ntasks-per-node=5

cd /home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/
rm -f cis.allgenes.lipid.clumped.sum

genes=("PPARG" "HCAR3" "HCAR2" "ABCA1" "ABCA1" "ABCA1" "MTTP" "MTTP" "PCSK9" "PCSK9" "ANGPTL3" "ANGPTL3" "APOC3" "APOC3" "APOC3" "CETP" "CETP" "CETP" "APOB" "APOB" "APOB" "APOB" "LDLR" "LDLR" "LDLR" "LPL" "LPL" "LPA" "LPA" "LPA" "HMGCR" "HMGCR")   
traits=("HDL" "HDL" "HDL" "LDL" "HDL" "TC" "HDL" "TC" "LDL" "TC" "logTG" "TC" "HDL" "logTG" "TC" "HDL" "logTG" "TC" "LDL" "HDL" "logTG" "TC" "LDL" "HDL" "TC" "HDL" "logTG" "LDL" "HDL" "TC" "LDL" "TC")    

len=${#genes[@]}

for ((i=0; i<$len; i++)); do
    gene=${genes[$i]}
    trait=${traits[$i]}
    prefix="${gene}.${trait}"
    clumped_file="${prefix}.snp.clumped"

    if [[ -f "$clumped_file" ]]; then
        awk '($3 != "SNP" && $3 != "") {print $3}' "$clumped_file" >> cis.allgenes.lipid.clumped.sum
        awk '($3 != "SNP" && $3 != "") {print $3}' "$clumped_file" > "${prefix}.snp.clumped.txt"
        rm -f "${prefix}.snp.nosex" "${prefix}.snp.clumped.ranges" "${prefix}.snp.log"
        echo "Processed: $prefix"
    else
        echo "Warning: $clumped_file not found"
    fi
done







################
#整理clumping的结果

cd /home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/
vi clump.glu.target.taiwan.res.slurm

#!/bin/bash

#SBATCH -J clumpres_glu_taiwan
#SBATCH -p cnl
#SBATCH -o out.%J.txt
#SBATCH -e err.%J.txt
#SBATCH -N 2
#SBATCH --ntasks-per-node=5


cd /home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/
rm -f cis.allgenes.glu.clumped.sum

genes=("ABCB11" "GCK" "GLP1R")
trait="glu"

for gene in "${genes[@]}"; do
    prefix="${gene}.${trait}"
    clumped_file="${prefix}.snp.clumped"

    if [[ -f "$clumped_file" ]]; then
        awk '($3 != "SNP" && $3 != "") {print $3}' "$clumped_file" >> cis.allgenes.glu.clumped.sum
        awk '($3 != "SNP" && $3 != "") {print $3}' "$clumped_file" > "${prefix}.snp.clumped.txt"
        rm -f "${prefix}.snp.nosex" "${prefix}.snp.clumped.ranges" "${prefix}.snp.log"
        echo "Processed: $prefix"
    else
        echo "Warning: $clumped_file not found"
    fi
done











R
library(data.table)
lipid <- fread("/home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/cis.allgenes.lipid.clumped.sum", h=F)
glu <- fread("/home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/cis.allgenes.glu.clumped.sum", h=F)
combine <- rbind(lipid, glu)

duplicated_rows <- combine[duplicated(combine$V1), ]
data_unique <- combine[!duplicated(combine$V1), ]
write.table(data_unique, "/home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/cis.apoc3_insr.lipid_glu_sbp.clumped.sum.uniq", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")






#capoc中提取位点
/home/zhuyunqing/software/plink20/plink2 \
--bfile /home/share1/check/863_qc/imputed_clean_qc \
--extract /home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/cis.apoc3_insr.lipid_glu_sbp.clumped.sum.uniq \
--recode A \
--out /home/zhuyunqing/drugtarget_lipid/onesamplemr/capoc_snp/CAPOC.cis.lipid_glu_sbp.snp


R
library(data.table)
dosage <- fread("/home/zhuyunqing/drugtarget_lipid/onesamplemr/capoc_snp/CAPOC.cis.lipid_glu_sbp.snp.raw")
dosage <- dosage[,-c(2:6)]
write.table(dosage, "/home/zhuyunqing/drugtarget_lipid/onesamplemr/capoc_snp/CAPOC.cis.lipid_glu_sbp.snp.txt", row.names=FALSE, col.names=TRUE, quote=FALSE)






##############
#在summary数据中提取effect
vi effect.lipid.target.glgc.slurm

#!/bin/bash

#SBATCH -J effect_lipid_glgc
#SBATCH -p cnl
#SBATCH -o out.%J.txt
#SBATCH -e err.%J.txt
#SBATCH -N 2
#SBATCH --ntasks-per-node=5


cd /home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/
Rscript effect.lipid.target.glgc.R 2>&1 | tee effect.lipid.target.glgc.log





vi effect.lipid.target.glgc.R

R
library(data.table)

# 路径设置
base_path <- "/home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/"
effect_path <- "/home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/"

# 定义基因名和血脂 trait
genes <- c("PPARG", "HCAR3", "HCAR2", "ABCA1", "ABCA1", "ABCA1", "MTTP", "MTTP", "PCSK9", "PCSK9", "ANGPTL3", "ANGPTL3", "APOC3", "APOC3", "APOC3", "CETP", "CETP", "CETP", "APOB", "APOB", "APOB", "APOB", "LDLR", "LDLR", "LDLR", "LPL", "LPL", "LPA", "LPA", "LPA", "HMGCR", "HMGCR")  
traits <- c("HDL", "HDL", "HDL", "LDL", "HDL", "TC", "HDL", "TC", "LDL", "TC", "logTG", "TC", "HDL", "logTG", "TC", "HDL", "logTG", "TC", "LDL", "HDL", "logTG", "TC", "LDL", "HDL", "TC", "HDL", "logTG", "LDL", "HDL", "TC", "LDL", "TC") 

# 血脂trait对应的summary文件
effect_files <- list(
  LDL   = "deduplicated_LDL_GLGC_EAS_2021_hg19.txt",
  HDL   = "deduplicated_HDL_GLGC_EAS_2021_hg19.txt",
  logTG = "deduplicated_logTG_GLGC_EAS_2021_hg19.txt",
  TC    = "deduplicated_TC_GLGC_EAS_2021_hg19.txt"
)


effect_data <- list()
for (traitname in unique(traits)) {
  file_path <- file.path(effect_path, effect_files[[traitname]])
  if (file.exists(file_path)) {
    effect_data[[traitname]] <- fread(file_path)
    cat("Loaded effect file for", traitname, "\n")
  } else {
    warning("Missing summary file for trait:", traitname)
  }
}

# 收集所有合并结果
all_merged <- list()

# 逐对处理
if (length(genes) != length(traits)) stop("genes和traits长度需一致！")
for (i in seq_along(genes)) {
  gene <- genes[i]
  trait <- traits[i]
  clump_file <- file.path(base_path, paste0(gene, ".", trait, ".snp.clumped"))
  output_file <- file.path(base_path, paste0(gene, ".", trait, ".snp.clumped.csv"))
  
  if (file.exists(clump_file) && !is.null(effect_data[[trait]])) {
    dat <- fread(clump_file, header = TRUE)
    snp <- dat[, 3, with = FALSE]
    setnames(snp, "SNP")

    merged <- merge(snp, effect_data[[trait]], by.x = "SNP", by.y = "name", all = FALSE)
    merged$allele1 <- paste0("dat$", merged$SNP, "_", merged$EA, " = 2 - dat$", merged$SNP, "_", merged$OA)
    merged$allele2 <- paste0("dat$", merged$SNP, "_", merged$OA, " = 2 - dat$", merged$SNP, "_", merged$EA)
    merged$allele1_cal <- paste0("dat$", merged$SNP, "_", merged$EA, " * (", merged$b, ") + ")
    merged$allele2_cal <- paste0("dat$", merged$SNP, "_", merged$OA, " * (-1) * (", merged$b, ") + ")
    merged$gene <- gene
    merged$trait <- trait

    write.csv(merged, output_file, row.names = FALSE)
    all_merged[[paste(gene, trait, sep = "_")]] <- merged
    cat("Processed:", gene, trait, "- SNPs matched:", nrow(merged), "\n")
  } else {
    cat("Skipped:", gene, trait, "- missing input or effect data\n")
  }
}

# 合并所有结果，并保存总表
if (length(all_merged) > 0) {
  merged_all_df <- rbindlist(all_merged, fill = TRUE)
  fwrite(merged_all_df, file = file.path(base_path, "all_genes_traits_merged_lipid.csv"))
  cat("Saved merged result of all genes and traits\n")
}




cd /home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/

##############
#在summary数据中提取effect
vi effect.glu.target.taiwan.slurm

#!/bin/bash

#SBATCH -J effect_glu_taiwan
#SBATCH -p cnl
#SBATCH -o out.%J.txt
#SBATCH -e err.%J.txt
#SBATCH -N 2
#SBATCH --ntasks-per-node=5


cd /home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/
Rscript effect.glu.target.taiwan.R 2>&1 | tee effect.glu.target.taiwan.log





vi effect.glu.target.taiwan.R

R
library(data.table)

# 路径设置
base_path <- "/home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/"
effect_path <- "/home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/"

# 定义基因名和血脂 trait
genes <- c("ABCB11", "GCK", "GLP1R")  # 可自由扩展
traits <- c("glu")

# 血脂trait对应的summary文件
effect_files <- list(
  glu   = "deduplicated_glucose_taiwan_summ_hg19.txt"
)

# 读取每个trait的summary文件一次并缓存
effect_data <- list()
for (trait in traits) {
  file_path <- file.path(effect_path, effect_files[[trait]])
  if (file.exists(file_path)) {
    effect_data[[trait]] <- fread(file_path)
    cat("Loaded effect file for", trait, "\n")
  } else {
    warning("Missing summary file for trait:", trait)
  }
}

# 用于收集所有结果
all_merged <- list()

# 双重循环处理每个 基因 × trait 的组合
for (gene in genes) {
  for (trait in traits) {
    clump_file <- file.path(base_path, paste0(gene, ".", trait, ".snp.clumped"))
    output_file <- file.path(base_path, paste0(gene, ".", trait, ".snp.clumped.csv"))

    if (file.exists(clump_file) && !is.null(effect_data[[trait]])) {
      dat <- fread(clump_file, header = TRUE)
      snp <- dat[, 3, with = FALSE]
      setnames(snp, "SNP")

      merged <- merge(snp, effect_data[[trait]], by.x = "SNP", by.y = "name", all = FALSE)

      merged$allele1 <- paste0("dat$", merged$SNP, "_", merged$EA, " = 2 - dat$", merged$SNP, "_", merged$OA)
      merged$allele2 <- paste0("dat$", merged$SNP, "_", merged$OA, " = 2 - dat$", merged$SNP, "_", merged$EA)
      merged$allele1_cal <- paste0("dat$", merged$SNP, "_", merged$EA, " * (", merged$BETA, ") + ")
      merged$allele2_cal <- paste0("dat$", merged$SNP, "_", merged$OA, " * (-1) * (", merged$BETA, ") + ")

      # 添加基因名和trait名列
      merged$gene <- gene
      merged$trait <- trait

      write.csv(merged, output_file, row.names = FALSE, col.names = TRUE)

      # 收集合并结果
      all_merged[[paste(gene, trait, sep = "_")]] <- merged

      cat("Processed:", gene, trait, "- SNPs matched:", nrow(merged), "\n")
    } else {
      cat("Skipped:", gene, trait, "- missing input or effect data\n")
    }
  }
}

# 合并所有结果，并保存总表
if (length(all_merged) > 0) {
  merged_all_df <- rbindlist(all_merged, fill = TRUE)
  fwrite(merged_all_df, file = file.path(base_path, "all_genes_traits_merged_glu.csv"))
  cat("Saved merged result of all genes and traits\n")
}












######
#grs of glu

R
library(data.table)
glu.summ <- fread("/home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/deduplicated_glucose_taiwan_summ_hg19.txt")
glu.summ <- subset(glu.summ, P<5E-08)
capoc.snp <- fread("/home/share1/check/863_qc/imputed_clean_qc.bim", header = FALSE)
capoc.snp <- capoc.snp[,2]
colnames(capoc.snp) <- "SNP"
m <- merge(capoc.snp, glu.summ, by.x="SNP", by.y="name", all=F)

fwrite(m, "/home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/deduplicated_glucose_taiwan_summ_hg19.sig.txt",
       sep = "\t")



/home/zhuyunqing/software/plink19/plink \
--bfile /home/share1/Public_database/1000G_plink/1000G_EAS/1000G_merge_allchr_EAS/ALL \
--clump /home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/deduplicated_glucose_taiwan_summ_hg19.sig.txt \
--clump-field P \
--clump-snp-field SNP \
--clump-p1 5e-8 \
--clump-kb 10000 \
--clump-r2 0.001 \
--clump-range /home/share1/Public_database/glist-hg19.txt \
--out /home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/deduplicated_glucose_taiwan_summ_hg19


######
#iv merge effect size
R
library(data.table)
snp <- fread("/home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/deduplicated_glucose_taiwan_summ_hg19.clumped")
snp <- snp[,3]

gwas <- fread("/home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/deduplicated_glucose_taiwan_summ_hg19.sig.txt")
m <- merge(snp, gwas, by.x="SNP", by.y="SNP", all=F)

m$allele1 <- paste0("dat$", m$SNP, "_", m$EA, " = 2 - dat$", m$SNP, "_", m$OA)
m$allele2 <- paste0("dat$", m$SNP, "_", m$OA, " = 2 - dat$", m$SNP, "_", m$EA)
m$allele1_cal <- paste0("dat$", m$SNP, "_", m$EA, " * (", m$BETA, ") + ")
m$allele2_cal <- paste0("dat$", m$SNP, "_", m$OA, " * (-1) * (", m$BETA, ") + ")
m.snp <- m[,1]

fwrite(m, "/home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/deduplicated_glucose_taiwan_summ_hg19.sig.beta.txt",
       sep = "\t")

write.table(m.snp, "/home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/deduplicated_glucose_taiwan_summ_hg19.sig.snp", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


#capoc中提取位点
/home/zhuyunqing/software/plink20/plink2 \
--bfile /home/share1/check/863_qc/imputed_clean_qc \
--extract /home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/deduplicated_glucose_taiwan_summ_hg19.sig.snp \
--recode A \
--out /home/zhuyunqing/drugtarget_lipid/onesamplemr/capoc_snp/CAPOC.cis.glu_grs.snp


R
library(data.table)
dosage <- fread("/home/zhuyunqing/drugtarget_lipid/onesamplemr/capoc_snp/CAPOC.cis.glu_grs.snp.raw")
dosage <- dosage[,-c(2:6)]
write.table(dosage, "/home/zhuyunqing/drugtarget_lipid/onesamplemr/capoc_snp/CAPOC.cis.glu_grs.snp.txt", row.names=FALSE, col.names=TRUE, quote=FALSE)











######
#grs of lipid
R
library(data.table)

tg.summ <- fread("/home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/deduplicated_logTG_GLGC_EAS_2021_hg19.txt")
tg.summ$p <- as.numeric(tg.summ$p)
tg.summ <- subset(tg.summ, p<5E-08)
capoc.snp <- fread("/home/share1/check/863_qc/imputed_clean_qc.bim", header = FALSE)
capoc.snp <- capoc.snp[,2]
colnames(capoc.snp) <- "SNP"
m <- merge(capoc.snp, tg.summ, by.x="SNP", by.y="name", all=F)

fwrite(m, "/home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/deduplicated_logTG_GLGC_EAS_2021_hg19.sig.txt",
       sep = "\t")


tc.summ <- fread("/home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/deduplicated_TC_GLGC_EAS_2021_hg19.txt")
tc.summ$p <- as.numeric(tc.summ$p)
tc.summ <- subset(tc.summ, p<5E-08)
capoc.snp <- fread("/home/share1/check/863_qc/imputed_clean_qc.bim", header = FALSE)
capoc.snp <- capoc.snp[,2]
colnames(capoc.snp) <- "SNP"
m <- merge(capoc.snp, tc.summ, by.x="SNP", by.y="name", all=F)

fwrite(m, "/home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/deduplicated_TC_GLGC_EAS_2021_hg19.sig.txt",
       sep = "\t")


/home/zhuyunqing/software/plink19/plink \
--bfile /home/share1/Public_database/1000G_plink/1000G_EAS/1000G_merge_allchr_EAS/ALL \
--clump /home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/deduplicated_logTG_GLGC_EAS_2021_hg19.sig.txt \
--clump-field p \
--clump-snp-field SNP \
--clump-p1 5e-8 \
--clump-kb 10000 \
--clump-r2 0.001 \
--clump-range /home/share1/Public_database/glist-hg19.txt \
--out /home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/deduplicated_logTG_GLGC_EAS_2021_hg19


/home/zhuyunqing/software/plink19/plink \
--bfile /home/share1/Public_database/1000G_plink/1000G_EAS/1000G_merge_allchr_EAS/ALL \
--clump /home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/deduplicated_TC_GLGC_EAS_2021_hg19.sig.txt \
--clump-field p \
--clump-snp-field SNP \
--clump-p1 5e-8 \
--clump-kb 10000 \
--clump-r2 0.001 \
--clump-range /home/share1/Public_database/glist-hg19.txt \
--out /home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/deduplicated_TC_GLGC_EAS_2021_hg19


######
#iv merge effect size
R
library(data.table)
#LOGTG
snp <- fread("/home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/deduplicated_logTG_GLGC_EAS_2021_hg19.clumped")
snp <- snp[,3]

gwas <- fread("/home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/deduplicated_logTG_GLGC_EAS_2021_hg19.sig.txt")
m <- merge(snp, gwas, by.x="SNP", by.y="SNP", all=F)

m$allele1 <- paste0("dat$", m$SNP, "_", m$EA, " = 2 - dat$", m$SNP, "_", m$OA)
m$allele2 <- paste0("dat$", m$SNP, "_", m$OA, " = 2 - dat$", m$SNP, "_", m$EA)
m$allele1_cal <- paste0("dat$", m$SNP, "_", m$EA, " * (", m$b, ") + ")
m$allele2_cal <- paste0("dat$", m$SNP, "_", m$OA, " * (-1) * (", m$b, ") + ")
m.snp1 <- m[,1]

fwrite(m, "/home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/deduplicated_logTG_GLGC_EAS_2021_hg19.sig.beta.txt",
       sep = "\t")

write.table(m.snp1, "/home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/deduplicated_logTG_GLGC_EAS_2021_hg19.sig.snp", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


#TC
snp <- fread("/home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/deduplicated_TC_GLGC_EAS_2021_hg19.clumped")
snp <- snp[,3]

gwas <- fread("/home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/deduplicated_TC_GLGC_EAS_2021_hg19.sig.txt")
m <- merge(snp, gwas, by.x="SNP", by.y="SNP", all=F)

m$allele1 <- paste0("dat$", m$SNP, "_", m$EA, " = 2 - dat$", m$SNP, "_", m$OA)
m$allele2 <- paste0("dat$", m$SNP, "_", m$OA, " = 2 - dat$", m$SNP, "_", m$EA)
m$allele1_cal <- paste0("dat$", m$SNP, "_", m$EA, " * (", m$b, ") + ")
m$allele2_cal <- paste0("dat$", m$SNP, "_", m$OA, " * (-1) * (", m$b, ") + ")
m.snp2 <- m[,1]

fwrite(m, "/home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/deduplicated_TC_GLGC_EAS_2021_hg19.sig.beta.txt",
       sep = "\t")

write.table(m.snp2, "/home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/deduplicated_TC_GLGC_EAS_2021_hg19.sig.snp", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


m.snp <- rbind(m.snp1, m.snp2)
m.snp <- m.snp[!duplicated(SNP)]
write.table(m.snp, "/home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/deduplicated_TC_TG_GLGC_EAS_2021_hg19.sig.snp", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


#capoc中提取位点
/home/zhuyunqing/software/plink20/plink2 \
--bfile /home/share1/check/863_qc/imputed_clean_qc \
--extract /home/zhuyunqing/drugtarget_lipid/onesamplemr/GLGC/deduplicated_TC_TG_GLGC_EAS_2021_hg19.sig.snp \
--recode A \
--out /home/zhuyunqing/drugtarget_lipid/onesamplemr/capoc_snp/CAPOC.cis.tc_tg_grs.snp

R
library(data.table)
dosage <- fread("/home/zhuyunqing/drugtarget_lipid/onesamplemr/capoc_snp/CAPOC.cis.tc_tg_grs.snp.raw")
dosage <- dosage[,-c(2:6)]
write.table(dosage, "/home/zhuyunqing/drugtarget_lipid/onesamplemr/capoc_snp/CAPOC.cis.tc_tg_grs.snp.txt", row.names=FALSE, col.names=TRUE, quote=FALSE)





