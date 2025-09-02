library(data.table)
library(broom)
library(dplyr)

dat <- fread("/Users/zhuyunqing/Desktop/study/drugMR_lipid/data/CAPEC/CAPEC.grs.pheno.new.txt")

#first epi

outcomes <- c("PANSS_reduce_rate", "PANSS_P_reduce_rate", "PANSS_N_reduce_rate", "PANSS_G_reduce_rate")
results <- list()

for (outcome in outcomes) {
  for (epi_group in unique(dat$FirstEpi)) {
    dat_sub <- dat %>% filter(FirstEpi == epi_group)
    
    fml <- as.formula(paste0(outcome, " ~ grs_GCK_w_taiwan_glu + age + age2 + factor(Sex) + factor(center) + PC1 + PC2 + PC3 + PC4 + PC5 + drug"))
    
    fit <- lm(fml, data = dat_sub)
    
    tidy_fit <- tidy(fit, conf.int = TRUE)
    res <- tidy_fit %>%
      filter(term == "grs_GCK_w_taiwan_glu") %>%
      mutate(outcome = outcome,
             FirstEpi = epi_group) %>%
      select(outcome, FirstEpi, beta = estimate, conf.low, conf.high, p.value)
    
    results[[paste0(outcome, "_epi", epi_group)]] <- res
  }
}

final_result_firstepi_stra <- bind_rows(results)
colnames(final_result_firstepi_stra)[2] <- "group"
final_result_firstepi_stra$stra = "firstepi"


##logistic
outcomes <- c("PANSS_reduce_rate_2g", "PANSS_P_reduce_rate_2g", 
              "PANSS_N_reduce_rate_2g", "PANSS_G_reduce_rate_2g",
              "PANSS_reduce_rate_2g_med",
              "PANSS_P_reduce_rate_2g_med",
              "PANSS_N_reduce_rate_2g_med",
              "PANSS_G_reduce_rate_2g_med") 
results <- list()

for (outcome in outcomes) {
  for (epi_group in unique(dat$FirstEpi)) {
    dat_sub <- dat %>% filter(FirstEpi == epi_group)
    
    fml <- as.formula(paste0(outcome, " ~ grs_GCK_w_taiwan_glu + age + age2 + factor(Sex) + factor(center) + PC1 + PC2 + PC3 + PC4 + PC5 + drug"))
    
    fit <- glm(fml, data = dat_sub, family = binomial)
    
    tidy_fit <- tidy(fit, conf.int = TRUE, exponentiate = TRUE)  
    res <- tidy_fit %>%
      filter(term == "grs_GCK_w_taiwan_glu") %>%
      mutate(outcome = outcome,
             FirstEpi = epi_group) %>%
      select(outcome, FirstEpi, OR = estimate, CI_low = conf.low, CI_high = conf.high, p.value)
    
    results[[paste0(outcome, "_epi", epi_group)]] <- res
  }
}

final_result_logit_firstepi_stra <- bind_rows(results)
colnames(final_result_logit_firstepi_stra)[2] <- "group"
final_result_logit_firstepi_stra$stra = "firstepi"


######
#drug with mets effect
dat$drug_mets <- ifelse((dat$drug == "喹硫平" | dat$drug == "奥氮平" | dat$drug == "氯氮平"), 1, 0 )

outcomes <- c("PANSS_reduce_rate", "PANSS_P_reduce_rate", "PANSS_N_reduce_rate", "PANSS_G_reduce_rate")
results <- list()

for (outcome in outcomes) {
  for (epi_group in unique(dat$drug_mets)) {
    dat_sub <- dat %>% filter(drug_mets == epi_group)
    
    fml <- as.formula(paste0(outcome, " ~ grs_GCK_w_taiwan_glu + age + age2 + factor(Sex) + factor(center) + PC1 + PC2 + PC3 + PC4 + PC5 + FirstEpi"))
    
    fit <- lm(fml, data = dat_sub)
    
    tidy_fit <- tidy(fit, conf.int = TRUE)
    res <- tidy_fit %>%
      filter(term == "grs_GCK_w_taiwan_glu") %>%
      mutate(outcome = outcome,
             drug_mets = epi_group) %>%
      select(outcome, drug_mets, beta = estimate, conf.low, conf.high, p.value)
    
    results[[paste0(outcome, "_epi", epi_group)]] <- res
  }
}

final_result_drugmets_stra <- bind_rows(results)
colnames(final_result_drugmets_stra)[2] <- "group"
final_result_drugmets_stra$stra = "drugmets"




##logistic
outcomes <- c("PANSS_reduce_rate_2g", "PANSS_P_reduce_rate_2g", 
              "PANSS_N_reduce_rate_2g", "PANSS_G_reduce_rate_2g",
              "PANSS_reduce_rate_2g_med",
              "PANSS_P_reduce_rate_2g_med",
              "PANSS_N_reduce_rate_2g_med",
              "PANSS_G_reduce_rate_2g_med") 
results <- list()

for (outcome in outcomes) {
  for (epi_group in unique(dat$drug_mets)) {
    dat_sub <- dat %>% filter(drug_mets == epi_group)
    
    fml <- as.formula(paste0(outcome, " ~ grs_GCK_w_taiwan_glu + age + age2 + factor(Sex) + factor(center) + PC1 + PC2 + PC3 + PC4 + PC5 + drug_mets"))
    
    fit <- glm(fml, data = dat_sub, family = binomial)
    
    tidy_fit <- tidy(fit, conf.int = TRUE, exponentiate = TRUE)  
    res <- tidy_fit %>%
      filter(term == "grs_GCK_w_taiwan_glu") %>%
      mutate(outcome = outcome,
             drug_mets = epi_group) %>%
      select(outcome, drug_mets, OR = estimate, CI_low = conf.low, CI_high = conf.high, p.value)
    
    results[[paste0(outcome, "_epi", epi_group)]] <- res
  }
}

final_result_logit_drugmets_stra <- bind_rows(results)
colnames(final_result_logit_drugmets_stra)[2] <- "group"
final_result_logit_drugmets_stra$stra = "drugmets"



final_resuls_linear <- rbind(final_result_firstepi_stra, final_result_drugmets_stra)
final_results_log <- rbind(final_result_logit_firstepi_stra, final_result_logit_drugmets_stra)


write.csv(final_resuls_linear, "/Users/zhuyunqing/Desktop/study/drugMR_lipid/result/CAPEC.APOC3.GCK.GRS2glu.PANSSredu.stra.linear.csv")
#STable 33

write.csv(final_results_log, "/Users/zhuyunqing/Desktop/study/drugMR_lipid/result/CAPEC.APOC3.GCK.GRS2glu.PANSSredu.stra.log.csv")
#STable 34

