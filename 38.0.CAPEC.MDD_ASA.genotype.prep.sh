#######
#CAPEC

#genotype prepare
#CAPEC中提取位点
/home/zhuyunqing/software/plink20/plink2 \
--bfile /home/share1/check/drugtargetmr_check/data/CAPEC/imputedSnpsDataset \
--extract /home/share1/check/drugtargetmr_check/data/SNP.extract.APOC3.GCK.txt \
--recode A \
--out /home/zhuyunqing/drugtarget_lipid/CAPEC_replication/CAPEC.cis.rep.APOC3.GCK

R
library(data.table)
dosage <- fread("/home/zhuyunqing/drugtarget_lipid/CAPEC_replication/CAPEC.cis.rep.APOC3.GCK.raw")
dosage <- dosage[,-c(2:6)]
write.table(dosage, "/home/zhuyunqing/drugtarget_lipid/CAPEC_replication/CAPEC.cis.rep.APOC3.GCK.txt", row.names=FALSE, col.names=TRUE, quote=FALSE)
q()
n

cd home/zhuyunqing/drugtarget_lipid/CAPEC_replication/
/home/zhuyunqing/software/plink20/plink2 \
--bfile /home/share1/check/drugtargetmr_check/data/CAPEC/imputedSnpsDataset \
--pca \
--out /home/zhuyunqing/drugtarget_lipid/CAPEC_replication/CAPEC.imputed



