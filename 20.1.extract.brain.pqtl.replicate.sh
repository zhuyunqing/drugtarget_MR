#在washington pqtl中提取与insr关联显著的snp
#数据来源：https://www.researchsquare.com/article/rs-2814616/v1
#alt-effect allele
#hg38

library(data.table)
library(tidyr)

setwd("/home/zhuyunqing/drugtarget_lipid/wustl_edu_pQTL/")
link <- fread("/home/share1/Public_database/snp_hg19_nodup.gz",h=T)

bim <- fread("/home/share1/check/863_qc/imputed_clean_qc.bim", h=F)
bim <- bim[,2]
colnames(bim) <- c("snp")


######
#apoc3
dat <- fread("phenocode-ApoC-III-CSF.tsv.gz", h=T)
#保留显著的
dat1 <- subset(dat, pval < 0.05)
#38版转37版
dat2 <- merge(link, dat1, by.x="name", by.y="rsids", all=F)
dat2 <- dat2 %>%
  separate(`chromosome:start`, into = c("chr", "bp19"), sep = ":", convert = TRUE)
dat2 <- dat2[,-c(4,5)]
#alt ref改名名
colnames(dat2)[4] <- "OA"
colnames(dat2)[5] <- "EA"
#截取apoc3区域的
apoc3.w <- subset(dat2, chr==11 & bp19 >= 116700623-100000 & bp19 <= 116703788+100000 & pval < 0.05000)
apoc3.w <- apoc3.w[,c(1,7)]
colnames(apoc3.w) <- c("SNP", "P")
#与863的bim文件merge 提出重合的
m <- merge(apoc3.w, bim, by.x = "SNP", by.y="snp", all=F)
write.table(m, "APOC3.wustl.cis.snp.capoc.txt", row.names=FALSE, col.names=TRUE, quote=FALSE)


/home/zhuyunqing/software/plink19/plink \
	--bfile /home/share1/Public_database/1000G_plink/1000G_EUR/1000G_merge_allchr_EUR/ALL \
	--clump /home/zhuyunqing/drugtarget_lipid/wustl_edu_pQTL/APOC3.wustl.cis.snp.capoc.txt \
	--clump-p1 1 \
	--clump-kb 500 \
	--clump-r2 0.1 \
	--clump-field P \
	--clump-range /home/share1/Public_database/glist-hg19.txt \
	--out /home/zhuyunqing/drugtarget_lipid/wustl_edu_pQTL/APOC3.wustl.cis.snp



#capoc中提取位点
/home/zhuyunqing/software/plink20/plink2 \
--bfile /home/share1/check/863_qc/imputed_clean_qc \
--extract /home/zhuyunqing/drugtarget_lipid/onesamplemr/snplist/brain.pqtl.extract.snp.list.txt \
--recode A \
--out /home/zhuyunqing/drugtarget_lipid/onesamplemr/capoc_snp/CAPOC.cis.brain.pqtl1

R
library(data.table)
dosage <- fread("/home/zhuyunqing/drugtarget_lipid/onesamplemr/capoc_snp/CAPOC.cis.brain.pqtl1.raw")
dosage <- dosage[,-c(2:6)]
write.table(dosage, "/home/zhuyunqing/drugtarget_lipid/onesamplemr/capoc_snp/CAPOC.cis.brain.pqtl1.txt", row.names=FALSE, col.names=TRUE, quote=FALSE)


