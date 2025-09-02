######
#replicate tsmr

####
#34.1 34.2是在CAPOC中跑panss减分率gwas的过程，我把生成的summary data放在check文件夹里了，复核中可不用再跑gwas


#generate phenotype panss.redu.rate
####
#replication using two-sample mr


library(broom)
library(data.table)
library(dplyr)
library(ggplot2)
library(forestploter)
library(TwoSampleMR)
library(mr.raps)
library(LDlinkR)
library(MRPRESSO)
library(corrplot)
library(reshape2)
library(grid)
library(ggrepel)


library("rio")

####
#生成GWAS分析的表型数据

#生成panss减分率 作为phenotype
capoc <- fread("/Users/zhuyunqing/Desktop/study/drugMR_lipid/data/onesamplemr/CAPOC.cis.capoc_cases.grs.pheno.new.txt", h=T)
attach(capoc)
capoc$FID = capoc$dn
capoc$IID = capoc$dn


model <- lm(PANSS_reduce_rate ~ age + age2 + as.factor(sex) + as.factor(centers) + PC1 + PC2 + PC3 + PC4 + PC5 + course + drugp07 + as.factor(medication), data = capoc)
capoc$PANSS_reduce_rate_res <- residuals(model)
model <- lm(PANSS_G_reduce_rate ~ age + age2 + as.factor(sex) + as.factor(centers) + PC1 + PC2 + PC3 + PC4 + PC5 + course + drugp07 + as.factor(medication), data = capoc)
capoc$PANSS_G_reduce_rate_res <- residuals(model)
model <- lm(PANSS_N_reduce_rate ~ age + age2 + as.factor(sex) + as.factor(centers) + PC1 + PC2 + PC3 + PC4 + PC5 + course + drugp07 + as.factor(medication), data = capoc)
capoc$PANSS_N_reduce_rate_res <- residuals(model)
model <- lm(PANSS_P_reduce_rate ~ age + age2 + as.factor(sex) + as.factor(centers) + PC1 + PC2 + PC3 + PC4 + PC5 + course + drugp07 + as.factor(medication), data = capoc)
capoc$PANSS_P_reduce_rate_res <- residuals(model)

model <- lm(ldl1_mg ~ age + age2 + sex + as.factor(centers) + PC1 + PC2 + PC3 + PC4 + PC5 + course + drugp07, data = capoc)
capoc$LDL_res <- residuals(model)
model <- lm(hdl1_mg ~ age + age2 + sex + as.factor(centers) + PC1 + PC2 + PC3 + PC4 + PC5 + course + drugp07, data = capoc)
capoc$HDL_res <- residuals(model)
model <- lm(tc1_mg ~ age + age2 + sex + as.factor(centers) + PC1 + PC2 + PC3 + PC4 + PC5 + course + drugp07, data = capoc)
capoc$TC_res <- residuals(model)
model <- lm(tg1_mg ~ age + age2 + sex + as.factor(centers) + PC1 + PC2 + PC3 + PC4 + PC5 + course + drugp07, data = capoc)
capoc$TG_res <- residuals(model)
model <- lm(glu1_mg ~ age + age2 + sex + as.factor(centers) + PC1 + PC2 + PC3 + PC4 + PC5 + course + drugp07, data = capoc)
capoc$glu_res <- residuals(model)


selected_vars <- c("FID", "IID", "PANSS_reduce_rate_res", "PANSS_G_reduce_rate_res", "PANSS_N_reduce_rate_res", "PANSS_P_reduce_rate_res",
                   "LDL_res", "HDL_res", "TC_res", "TG_res", "glu_res","PC1", "PC2", "PC3", "PC4", "PC5")

new_data <- dplyr::select(as_tibble(capoc), all_of(selected_vars))

fwrite(new_data, "/Users/zhuyunqing/Desktop/study/drugMR_lipid/data/onesamplemr/CAPOC.cis.capoc_cases.grs.panss_redu.lipid.glu.to_gwas.txt",
       sep = "\t")


