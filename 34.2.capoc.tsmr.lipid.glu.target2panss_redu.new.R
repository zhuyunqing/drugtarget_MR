#####
#在linux中运行

conda activate r_431

R

##双样本mr验证分析
#apoc3 tg1, gck glu1 ~ panss redu.rate

library(broom)
library(data.table)
library(dplyr)
library(ggplot2)
library(TwoSampleMR)
library(mr.raps)
library(MRPRESSO)
library(reshape2)
library(grid)
library(ggrepel)



####################################
#antidiabetic target 2 panss_redu
bim <- fread("/home/share1/check/863_qc/imputed_clean_qc.bim", header = FALSE)
bim <- bim[,2]
colnames(bim) <- "SNP"

glu.gwas <- fread("/home/share1/check/drugtargetmr_check/data/deduplicated_glucose_taiwan_summ_hg19.txt", h=T )
glu.gwas <- merge(bim, glu.gwas, by.x = "SNP", by.y = "name", all = F)
glu.gwas <- as.data.frame(glu.gwas)
#降糖作用 beta*-1
glu.gwas$BETA <- -1 * glu.gwas$BETA
glu.gwas.sig <- subset(glu.gwas, P < 1e-05)

glu.gwas.sig <- format_data(
  glu.gwas.sig,
  type = "exposure",
  snps = NULL,
  header = TRUE,
  #phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "EAF",
  effect_allele_col = "EA",
  other_allele_col = "OA",
  pval_col = "P",
  samplesize_col = "N",
  chr_col = "chr",
  pos_col = "bp19",
  log_pval = FALSE
)

fwrite(glu.gwas.sig, "/home/zhuyunqing/drugtarget_lipid/tsmr/glu.taiwan_hg19.sig.1e-05.exp.txt",
       sep = "\t")

glu.GCK.sig <- subset(glu.gwas.sig, chr.exposure == 7 & pos.exposure >= 44182812 - 100000 & pos.exposure <= 44229038 + 100000 )
fwrite(glu.GCK.sig, "/home/zhuyunqing/drugtarget_lipid/tsmr/GCK.glu.taiwan_hg19.sig.exp.txt",
       sep = "\t")






#########
#tg 药靶APOC3
tg.gwas <- fread("/home/share1/check/drugtargetmr_check/data/deduplicated_logTG_GLGC_EAS_2021_hg19.txt", h=T )
tg.gwas <- merge(bim, tg.gwas, by.x = "SNP", by.y = "name", all = F)
tg.gwas <- as.data.frame(tg.gwas)
tg.gwas$p <- as.numeric(tg.gwas$p)
#降脂作用 beta*-1
tg.gwas$b <- -1 * tg.gwas$b
tg.gwas.sig <- subset(tg.gwas, p < 1e-05)

tg.gwas.sig <- format_data(
  tg.gwas.sig,
  type = "exposure",
  snps = NULL,
  header = TRUE,
  #phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "b",
  se_col = "se",
  eaf_col = "EAF",
  effect_allele_col = "EA",
  other_allele_col = "OA",
  pval_col = "p",
  samplesize_col = "N",
  chr_col = "chr",
  pos_col = "bp19",
  log_pval = FALSE
)

fwrite(tg.gwas.sig, "/home/zhuyunqing/drugtarget_lipid/tsmr/tg.glgc_hg19.sig.1e-05.exp.txt",
       sep = "\t")


tg.APOC3.sig <- subset(tg.gwas.sig, chr.exposure == 11 & 
                         pos.exposure >= 116700623 - 100000 & 
                         pos.exposure <= 116703788 + 100000 )
fwrite(tg.APOC3.sig, "/home/zhuyunqing/drugtarget_lipid/tsmr/APOC3.tg.glgc_hg19.sig.exp.txt",
       sep = "\t")





#########
#tc 药靶APOC3
tc.gwas <- fread("/home/share1/check/drugtargetmr_check/data/deduplicated_TC_GLGC_EAS_2021_hg19.txt", h=T )
tc.gwas <- merge(bim, tc.gwas, by.x = "SNP", by.y = "name", all = F)
tc.gwas <- as.data.frame(tc.gwas)
tc.gwas$p <- as.numeric(tc.gwas$p)
#降脂作用 beta*-1
tc.gwas$b <- -1 * tc.gwas$b
tc.gwas.sig <- subset(tc.gwas, p < 1e-05)

tc.gwas.sig <- format_data(
  tc.gwas.sig,
  type = "exposure",
  snps = NULL,
  header = TRUE,
  #phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "b",
  se_col = "se",
  eaf_col = "EAF",
  effect_allele_col = "EA",
  other_allele_col = "OA",
  pval_col = "p",
  samplesize_col = "N",
  chr_col = "chr",
  pos_col = "bp19",
  log_pval = FALSE
)

fwrite(tc.gwas.sig, "/home/zhuyunqing/drugtarget_lipid/tsmr/tc.glgc_hg19.sig.1e-05.exp.txt",
       sep = "\t")


tc.APOC3.sig <- subset(tc.gwas.sig, chr.exposure == 11 & 
                         pos.exposure >= 116700623 - 100000 & 
                         pos.exposure <= 116703788 + 100000 )
fwrite(tc.APOC3.sig, "/home/zhuyunqing/drugtarget_lipid/tsmr/APOC3.tc.glgc_hg19.sig.exp.txt",
       sep = "\t")











#######
#drug target tsmr

glu.GCK.exp <- subset(glu.gwas.sig, SNP == "rs2284773" | 
                        SNP == "rs2908277" | 
                        SNP == "rs2908289" | 
                        SNP == "rs60463592" | 
                        SNP == "rs62459092" | 
                        SNP == "rs741038")
glu.GCK.exp$exposure <- "GCK_GLU"

tg.APOC3.exp <- subset(tg.APOC3.sig, SNP == "rs4938309" | 
                         SNP == "rs595049" | 
                         SNP == "rs11216190" | 
                         SNP == "rs17519093" | 
                         SNP == "rs888245" | 
                         SNP == "rs518181" | 
                         SNP == "rs116852373" | 
                         SNP == "rs79210031" | 
                         SNP == "rs139733873" | 
                         SNP == "rs80022746" | 
                         SNP == "rs2849169" | 
                         SNP == "rs583219" | 
                         SNP == "rs651821")

tg.APOC3.exp$exposure <- "APOC3_TG"

tc.APOC3.exp <- subset(tc.APOC3.sig, SNP == "rs2239013" | 
                         SNP == "rs595049" | 
                         SNP == "rs651821" | 
                         SNP == "rs918143")

tc.APOC3.exp$exposure <- "APOC3_TC"


######
#Outcome summ
#####
#outcome gwas - PANSS redu gwas
capoc.PANSS_reduce_rate.summ <- fread("/home/share1/check/drugtargetmr_check/data/CAPOC.PANSS_reduce_rate_res.gwas.PANSS_reduce_rate_res.glm.linear", h=T)
capoc.PANSS_reduce_rate.summ$other.allele <- ifelse(capoc.PANSS_reduce_rate.summ$A1 == capoc.PANSS_reduce_rate.summ$ALT, 
                                                    capoc.PANSS_reduce_rate.summ$REF, capoc.PANSS_reduce_rate.summ$ALT)
capoc.PANSS_reduce_rate.summ <- capoc.PANSS_reduce_rate.summ %>%
  distinct(ID, .keep_all = TRUE)
capoc.PANSS_reduce_rate.summ <- as.data.frame(capoc.PANSS_reduce_rate.summ)


capoc.PANSS_P_reduce_rate.summ <- fread("/home/share1/check/drugtargetmr_check/data/CAPOC.PANSS_P_reduce_rate_res.gwas.PANSS_P_reduce_rate_res.glm.linear", h=T)
capoc.PANSS_P_reduce_rate.summ$other.allele <- ifelse(capoc.PANSS_P_reduce_rate.summ$A1 == capoc.PANSS_P_reduce_rate.summ$ALT, 
                                                      capoc.PANSS_P_reduce_rate.summ$REF, capoc.PANSS_P_reduce_rate.summ$ALT)
capoc.PANSS_P_reduce_rate.summ <- capoc.PANSS_P_reduce_rate.summ %>%
  distinct(ID, .keep_all = TRUE)
capoc.PANSS_P_reduce_rate.summ <- as.data.frame(capoc.PANSS_P_reduce_rate.summ)


capoc.PANSS_N_reduce_rate.summ <- fread("/home/share1/check/drugtargetmr_check/data/CAPOC.PANSS_N_reduce_rate_res.gwas.PANSS_N_reduce_rate_res.glm.linear", h=T)
capoc.PANSS_N_reduce_rate.summ$other.allele <- ifelse(capoc.PANSS_N_reduce_rate.summ$A1 == capoc.PANSS_N_reduce_rate.summ$ALT, 
                                                      capoc.PANSS_N_reduce_rate.summ$REF, capoc.PANSS_N_reduce_rate.summ$ALT)
capoc.PANSS_N_reduce_rate.summ <- capoc.PANSS_N_reduce_rate.summ %>%
  distinct(ID, .keep_all = TRUE)
capoc.PANSS_N_reduce_rate.summ <- as.data.frame(capoc.PANSS_N_reduce_rate.summ)


capoc.PANSS_G_reduce_rate.summ <- fread("/home/share1/check/drugtargetmr_check/data/CAPOC.PANSS_G_reduce_rate_res.gwas.PANSS_G_reduce_rate_res.glm.linear", h=T)
capoc.PANSS_G_reduce_rate.summ$other.allele <- ifelse(capoc.PANSS_G_reduce_rate.summ$A1 == capoc.PANSS_G_reduce_rate.summ$ALT, 
                                                      capoc.PANSS_G_reduce_rate.summ$REF, capoc.PANSS_G_reduce_rate.summ$ALT)
capoc.PANSS_G_reduce_rate.summ <- capoc.PANSS_G_reduce_rate.summ %>%
  distinct(ID, .keep_all = TRUE)
capoc.PANSS_G_reduce_rate.summ <- as.data.frame(capoc.PANSS_G_reduce_rate.summ)



capoc.PANSS_reduce_rate.outcome <- format_data(
  capoc.PANSS_reduce_rate.summ,
  type = "outcome",
  snps = NULL,
  header = TRUE,
  #phenotype_col = "Phenotype",
  snp_col = "ID",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "A1_FREQ",
  effect_allele_col = "A1",
  other_allele_col = "other.allele",
  pval_col = "P",
  samplesize_col = "OBS_CT",
  chr_col = "#CHROM",
  pos_col = "POS",
  log_pval = FALSE
)
capoc.PANSS_reduce_rate.outcome$outcome <- "PANSS_reduce_rate"


capoc.PANSS_P_reduce_rate.outcome <- format_data(
  capoc.PANSS_P_reduce_rate.summ,
  type = "outcome",
  snps = NULL,
  header = TRUE,
  #phenotype_col = "Phenotype",
  snp_col = "ID",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "A1_FREQ",
  effect_allele_col = "A1",
  other_allele_col = "other.allele",
  pval_col = "P",
  samplesize_col = "OBS_CT",
  chr_col = "#CHROM",
  pos_col = "POS",
  log_pval = FALSE
)
capoc.PANSS_P_reduce_rate.outcome$outcome <- "PANSS_P_reduce_rate"


capoc.PANSS_N_reduce_rate.outcome <- format_data(
  capoc.PANSS_N_reduce_rate.summ,
  type = "outcome",
  snps = NULL,
  header = TRUE,
  #phenotype_col = "Phenotype",
  snp_col = "ID",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "A1_FREQ",
  effect_allele_col = "A1",
  other_allele_col = "other.allele",
  pval_col = "P",
  samplesize_col = "OBS_CT",
  chr_col = "#CHROM",
  pos_col = "POS",
  log_pval = FALSE
)
capoc.PANSS_N_reduce_rate.outcome$outcome <- "PANSS_N_reduce_rate"


capoc.PANSS_G_reduce_rate.outcome <- format_data(
  capoc.PANSS_G_reduce_rate.summ,
  type = "outcome",
  snps = NULL,
  header = TRUE,
  #phenotype_col = "Phenotype",
  snp_col = "ID",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "A1_FREQ",
  effect_allele_col = "A1",
  other_allele_col = "other.allele",
  pval_col = "P",
  samplesize_col = "OBS_CT",
  chr_col = "#CHROM",
  pos_col = "POS",
  log_pval = FALSE
)
capoc.PANSS_G_reduce_rate.outcome$outcome <- "PANSS_G_reduce_rate"











#####
#tsmr association

exposures <- c("tg.APOC3.exp", "tc.APOC3.exp", "glu.GCK.exp")
outcomes <- c("capoc.PANSS_reduce_rate.outcome", "capoc.PANSS_P_reduce_rate.outcome", 
              "capoc.PANSS_N_reduce_rate.outcome", "capoc.PANSS_G_reduce_rate.outcome")

results1 <- data.frame()
test1 <- data.frame()

for (exposure_name in exposures) {
  for (outcome_name in outcomes) {
    exposure_dat <- get(exposure_name)
    outcome_dat <- get(outcome_name)
    
    dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat, action = 1)
    
    res <- mr(dat)
    
    
    res$b <- as.numeric(res$b)
    res$se <- as.numeric(res$se)
    res$nsnp <- as.numeric(res$nsnp)
    # res$OR <- exp(res$b)
    # res$OR_L <- exp(res$b - 1.96 * res$se)
    # res$OR_U <- exp(res$b + 1.96 * res$se)
    
    # 添加 exposure 和 outcome 标签
    res$exposure <- exposure_name
    res$outcome <- outcome_name
    
    # 去除无用列
    res <- res[, c("exposure", "outcome", "method", "nsnp", "b", "se",  "pval")]
    
    # 合并结果
    results1 <- rbind(results1, res)
    
    
    
    #pleiotropy test
    pleio <- mr_pleiotropy_test(dat)
    pleio

    #heterogeneity test
    het <- mr_heterogeneity(dat)
    het
    
    #MR Steiger
    out <- directionality_test(dat)
    out
    
    #F-Statistics
    R2=out$snp_r2.exposure
    n=dat$samplesize.exposure[1] #注意修改样本???
    k=res$nsnp[1]
    F = R2/(1-R2) * (n-k-1)/k
    F
    
    
    test <- cbind(out$exposure, out$outcome, F, 
                  pleio$egger_intercept, pleio$se, pleio$pval, 
                  het$Q[2], het$Q_pval[2], 
                  out$snp_r2.exposure,out$snp_r2.outcome, 
                  out$correct_causal_direction, out$steiger_pval)
    
    colnames(test) <- c("Exposure","Outcome","F","Egger.inter",
                        "Egger.SE","P.inter","Q","P.heterogeneity",
                        "snp_r2.exposure","snp_r2.outcome","Direction","P.Steiger")
    
    test1 <- rbind(test1, test)
    
    
    
  }
}

results1
test1

write.csv(results1, "/home/zhuyunqing/drugtarget_lipid/tsmr/glu.panss_redu.APOC3.GCK.result.csv")
write.csv(test1, "/home/zhuyunqing/drugtarget_lipid/tsmr/glu.panss_redu.tsmr.APOC3.GCK.csv")

