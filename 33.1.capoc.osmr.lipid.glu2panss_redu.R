########
#GLU 2 PANSS.redu

library(data.table)
library(broom)
library(dplyr)

dat <- fread("/Users/zhuyunqing/Desktop/study/drugMR_lipid/data/onesamplemr/CAPOC.cis.capoc_cases.grs.pheno.new.txt", h=T)


###########################
#observational analysis

####linear
outcomes <- c("PANSS_reduce_rate", "PANSS_P_reduce_rate", "PANSS_N_reduce_rate", "PANSS_G_reduce_rate")


res_list <- list()

for (outcome in outcomes) {
  formula_str <- paste0(outcome, " ~ glu1_mg + age + age2 + as.factor(sex) + as.factor(centers) + course + drugp07 + as.factor(medication)")
    fit <- try(lm(as.formula(formula_str), data = dat), silent = TRUE)
  
  if (!inherits(fit, "try-error")) {
    coef_summary <- summary(fit)$coefficients
    ci <- confint(fit)
    
    if ("glu1_mg" %in% rownames(coef_summary)) {
      beta <- coef_summary["glu1_mg", "Estimate"]
      pval <- coef_summary["glu1_mg", "Pr(>|t|)"]
      lower <- ci["glu1_mg", 1]
      upper <- ci["glu1_mg", 2]
      
      res_list[[outcome]] <- data.frame(
        Trait = "glu1_mg",
        Outcome = outcome,
        Beta = round(beta, 3),
        CI = paste0(sprintf("%.3f", lower), " ~ ", sprintf("%.3f", upper)),
        P = signif(pval, 3)
      )
    }
  }
}

res_df_glu1_lm <- bind_rows(res_list)
print(res_df_glu1_lm)






######logistic
outcomes <- c("PANSS_reduce_rate_2g", "PANSS_P_reduce_rate_2g", "PANSS_N_reduce_rate_2g", "PANSS_G_reduce_rate_2g",
              "PANSS_reduce_rate_2g_med", "PANSS_P_reduce_rate_2g_med", "PANSS_N_reduce_rate_2g_med", "PANSS_G_reduce_rate_2g_med")

res_list <- list()

for (outcome in outcomes) {
  
  formula_str <- paste0(outcome, " ~ glu1_mg + age + age2 + as.factor(sex) + as.factor(centers) + course + drugp07 + as.factor(medication)")
  
  fit <- try(glm(as.formula(formula_str), data = dat, family = binomial), silent = TRUE)
  
  if (!inherits(fit, "try-error")) {
    coef_summary <- summary(fit)$coefficients
    ci <- suppressMessages(confint(fit))  
    
    if ("glu1_mg" %in% rownames(coef_summary)) {
      beta <- coef_summary["glu1_mg", "Estimate"]
      or <- exp(beta)
      pval <- coef_summary["glu1_mg", "Pr(>|z|)"]
      lower <- exp(ci["glu1_mg", 1])
      upper <- exp(ci["glu1_mg", 2])
      
      res_list[[outcome]] <- data.frame(
        Trait = "glu1_mg",
        Outcome = outcome,
        OR = round(or, 3),
        CI = paste0(sprintf("%.3f", lower), " ~ ", sprintf("%.3f", upper)),
        P = signif(pval, 3)
      )
    }
  }
}

res_df_glu1_logistic <- bind_rows(res_list)
print(res_df_glu1_logistic)






#####
#observe tg to panss
outcomes <- c("PANSS_reduce_rate", "PANSS_P_reduce_rate", "PANSS_N_reduce_rate", "PANSS_G_reduce_rate")


res_list <- list()

for (outcome in outcomes) {
  formula_str <- paste0(outcome, " ~ tg1_mg + age + age2 + as.factor(sex) + as.factor(centers) + course + drugp07 + as.factor(medication)")
  fit <- try(lm(as.formula(formula_str), data = dat), silent = TRUE)
  
  if (!inherits(fit, "try-error")) {
    coef_summary <- summary(fit)$coefficients
    ci <- confint(fit)
    
    if ("tg1_mg" %in% rownames(coef_summary)) {
      beta <- coef_summary["tg1_mg", "Estimate"]
      pval <- coef_summary["tg1_mg", "Pr(>|t|)"]
      lower <- ci["tg1_mg", 1]
      upper <- ci["tg1_mg", 2]
      
      res_list[[outcome]] <- data.frame(
        Trait = "tg1_mg",
        Outcome = outcome,
        Beta = round(beta, 3),
        CI = paste0(sprintf("%.3f", lower), " ~ ", sprintf("%.3f", upper)),
        P = signif(pval, 3)
      )
    }
  }
}

res_df_tg1_lm <- bind_rows(res_list)
print(res_df_tg1_lm)






######logistic
outcomes <- c("PANSS_reduce_rate_2g", "PANSS_P_reduce_rate_2g", "PANSS_N_reduce_rate_2g", "PANSS_G_reduce_rate_2g",
              "PANSS_reduce_rate_2g_med", "PANSS_P_reduce_rate_2g_med", "PANSS_N_reduce_rate_2g_med", "PANSS_G_reduce_rate_2g_med")

res_list <- list()

for (outcome in outcomes) {
  
  formula_str <- paste0(outcome, " ~ tg1_mg + age + age2 + as.factor(sex) + as.factor(centers) + course + drugp07 + as.factor(medication)")
  
  fit <- try(glm(as.formula(formula_str), data = dat, family = binomial), silent = TRUE)
  
  if (!inherits(fit, "try-error")) {
    coef_summary <- summary(fit)$coefficients
    ci <- suppressMessages(confint(fit))  
    
    if ("tg1_mg" %in% rownames(coef_summary)) {
      beta <- coef_summary["tg1_mg", "Estimate"]
      or <- exp(beta)
      pval <- coef_summary["tg1_mg", "Pr(>|z|)"]
      lower <- exp(ci["tg1_mg", 1])
      upper <- exp(ci["tg1_mg", 2])
      
      res_list[[outcome]] <- data.frame(
        Trait = "tg1_mg",
        Outcome = outcome,
        OR = round(or, 3),
        CI = paste0(sprintf("%.3f", lower), " ~ ", sprintf("%.3f", upper)),
        P = signif(pval, 3)
      )
    }
  }
}

res_df_tg1_logistic <- bind_rows(res_list)
print(res_df_tg1_logistic)




#####
#observe tc to panss
outcomes <- c("PANSS_reduce_rate", "PANSS_P_reduce_rate", "PANSS_N_reduce_rate", "PANSS_G_reduce_rate")


res_list <- list()

for (outcome in outcomes) {
  formula_str <- paste0(outcome, " ~ tc1_mg + age + age2 + as.factor(sex) + as.factor(centers) + course + drugp07 + as.factor(medication)")
  fit <- try(lm(as.formula(formula_str), data = dat), silent = TRUE)
  
  if (!inherits(fit, "try-error")) {
    coef_summary <- summary(fit)$coefficients
    ci <- confint(fit)
    
    if ("tc1_mg" %in% rownames(coef_summary)) {
      beta <- coef_summary["tc1_mg", "Estimate"]
      pval <- coef_summary["tc1_mg", "Pr(>|t|)"]
      lower <- ci["tc1_mg", 1]
      upper <- ci["tc1_mg", 2]
      
      res_list[[outcome]] <- data.frame(
        Trait = "tc1_mg",
        Outcome = outcome,
        Beta = round(beta, 3),
        CI = paste0(sprintf("%.3f", lower), " ~ ", sprintf("%.3f", upper)),
        P = signif(pval, 3)
      )
    }
  }
}

res_df_tc1_lm <- bind_rows(res_list)
print(res_df_tc1_lm)






######logistic
outcomes <- c("PANSS_reduce_rate_2g", "PANSS_P_reduce_rate_2g", "PANSS_N_reduce_rate_2g", "PANSS_G_reduce_rate_2g",
              "PANSS_reduce_rate_2g_med", "PANSS_P_reduce_rate_2g_med", "PANSS_N_reduce_rate_2g_med", "PANSS_G_reduce_rate_2g_med")

res_list <- list()

for (outcome in outcomes) {
  
  formula_str <- paste0(outcome, " ~ tc1_mg + age + age2 + as.factor(sex) + as.factor(centers) + course + drugp07 + as.factor(medication)")
  
  fit <- try(glm(as.formula(formula_str), data = dat, family = binomial), silent = TRUE)
  
  if (!inherits(fit, "try-error")) {
    coef_summary <- summary(fit)$coefficients
    ci <- suppressMessages(confint(fit))  
    
    if ("tc1_mg" %in% rownames(coef_summary)) {
      beta <- coef_summary["tc1_mg", "Estimate"]
      or <- exp(beta)
      pval <- coef_summary["tc1_mg", "Pr(>|z|)"]
      lower <- exp(ci["tc1_mg", 1])
      upper <- exp(ci["tc1_mg", 2])
      
      res_list[[outcome]] <- data.frame(
        Trait = "tc1_mg",
        Outcome = outcome,
        OR = round(or, 3),
        CI = paste0(sprintf("%.3f", lower), " ~ ", sprintf("%.3f", upper)),
        P = signif(pval, 3)
      )
    }
  }
}

res_df_tc1_logistic <- bind_rows(res_list)
print(res_df_tc1_logistic)



res_lm <- rbind(res_df_glu1_lm, res_df_tg1_lm, res_df_tc1_lm)
res_log <- rbind(res_df_glu1_logistic, res_df_tg1_logistic, res_df_tc1_logistic)


write.csv(res_lm, 
          "/Users/zhuyunqing/Desktop/study/drugMR_lipid/data/onesamplemr/obser_glutctg2panss_redu.linear.csv", 
          row.names = FALSE)

########
#stable 15


write.csv(res_log, 
          "/Users/zhuyunqing/Desktop/study/drugMR_lipid/data/onesamplemr/obser_glutctg2panss_redu.logistic.csv", 
          row.names = FALSE)

########
#stable 16





########
#one sample mr analysis

#iv selection
#14.25

#gen.grs 
#14.5

#asso

library(data.table)
library(broom)
library(dplyr)

dat <- fread("/Users/zhuyunqing/Desktop/study/drugMR_lipid/data/onesamplemr/CAPOC.cis.capoc_cases.grs.pheno.new.txt", h=T)
fit <- lm(glu1_mg ~ grs_w_taiwan_glu_higher + age + age2 + as.factor(sex) + as.factor(centers) + course + drugp07 + PC1 + PC2 + PC3 + PC4 + PC5,
          data = dat)

dat$grs_pred_glu <- predict(fit)


fit <- lm(tg1_mg ~ grs_w_glgc_logTG_higher + age + age2 + as.factor(sex) + as.factor(centers) + course + drugp07 + PC1 + PC2 + PC3 + PC4 + PC5,
          data = dat)

dat$grs_pred_tg <- predict(fit)


fit <- lm(tc1_mg ~ grs_w_glgc_TC_higher + age + age2 + as.factor(sex) + as.factor(centers) + course + drugp07 + PC1 + PC2 + PC3 + PC4 + PC5,
          data = dat)

dat$grs_pred_tc <- predict(fit)



####linear
outcomes <- c("PANSS_reduce_rate", "PANSS_P_reduce_rate", "PANSS_N_reduce_rate", "PANSS_G_reduce_rate")


res_list <- list()

for (outcome in outcomes) {
  formula_str <- paste0(outcome, " ~ grs_pred_glu + age + age2 + as.factor(sex) + as.factor(centers) + course + drugp07 + as.factor(medication) + PC1 + PC2 + PC3 + PC4 + PC5")
  fit <- try(lm(as.formula(formula_str), data = dat), silent = TRUE)
  
  if (!inherits(fit, "try-error")) {
    coef_summary <- summary(fit)$coefficients
    ci <- confint(fit)
    
    if ("grs_pred_glu" %in% rownames(coef_summary)) {
      beta <- coef_summary["grs_pred_glu", "Estimate"]
      pval <- coef_summary["grs_pred_glu", "Pr(>|t|)"]
      lower <- ci["grs_pred_glu", 1]
      upper <- ci["grs_pred_glu", 2]
      
      res_list[[outcome]] <- data.frame(
        Trait = "grs_pred_glu",
        Outcome = outcome,
        Beta = round(beta, 4),
        CI = paste0(sprintf("%.3f", lower), " ~ ", sprintf("%.3f", upper)),
        P = signif(pval, 4)
      )
    }
  }
}

res_df_glu1_lm <- bind_rows(res_list)
print(res_df_glu1_lm)






######logistic
outcomes <- c("PANSS_reduce_rate_2g", "PANSS_P_reduce_rate_2g", "PANSS_N_reduce_rate_2g", "PANSS_G_reduce_rate_2g",
              "PANSS_reduce_rate_2g_med", "PANSS_P_reduce_rate_2g_med", "PANSS_N_reduce_rate_2g_med", "PANSS_G_reduce_rate_2g_med")

res_list <- list()

for (outcome in outcomes) {
  
  formula_str <- paste0(outcome, " ~ grs_pred_glu + age + age2 + as.factor(sex) + as.factor(centers) + course + drugp07 + as.factor(medication) + PC1 + PC2 + PC3 + PC4 + PC5")
  
  fit <- try(glm(as.formula(formula_str), data = dat, family = binomial), silent = TRUE)
  
  if (!inherits(fit, "try-error")) {
    coef_summary <- summary(fit)$coefficients
    ci <- suppressMessages(confint(fit))  
    
    if ("grs_pred_glu" %in% rownames(coef_summary)) {
      beta <- coef_summary["grs_pred_glu", "Estimate"]
      or <- exp(beta)
      pval <- coef_summary["grs_pred_glu", "Pr(>|z|)"]
      lower <- exp(ci["grs_pred_glu", 1])
      upper <- exp(ci["grs_pred_glu", 2])
      
      res_list[[outcome]] <- data.frame(
        Trait = "grs_pred_glu",
        Outcome = outcome,
        OR = round(or, 4),
        CI = paste0(sprintf("%.3f", lower), " ~ ", sprintf("%.3f", upper)),
        P = signif(pval, 4)
      )
    }
  }
}

res_df_glu1_logistic <- bind_rows(res_list)
print(res_df_glu1_logistic)









#######
#grs pred tg 2 panss

####linear
outcomes <- c("PANSS_reduce_rate", "PANSS_P_reduce_rate", "PANSS_N_reduce_rate", "PANSS_G_reduce_rate")


res_list <- list()

for (outcome in outcomes) {
  formula_str <- paste0(outcome, " ~ grs_pred_tg + age + age2 + as.factor(sex) + as.factor(centers) + course + drugp07 + as.factor(medication) + PC1 + PC2 + PC3 + PC4 + PC5")
  fit <- try(lm(as.formula(formula_str), data = dat), silent = TRUE)
  
  if (!inherits(fit, "try-error")) {
    coef_summary <- summary(fit)$coefficients
    ci <- confint(fit)
    
    if ("grs_pred_tg" %in% rownames(coef_summary)) {
      beta <- coef_summary["grs_pred_tg", "Estimate"]
      pval <- coef_summary["grs_pred_tg", "Pr(>|t|)"]
      lower <- ci["grs_pred_tg", 1]
      upper <- ci["grs_pred_tg", 2]
      
      res_list[[outcome]] <- data.frame(
        Trait = "grs_pred_tg",
        Outcome = outcome,
        Beta = round(beta, 4),
        CI = paste0(sprintf("%.3f", lower), " ~ ", sprintf("%.3f", upper)),
        P = signif(pval, 4)
      )
    }
  }
}

res_df_tg1_lm <- bind_rows(res_list)
print(res_df_tg1_lm)



######logistic
outcomes <- c("PANSS_reduce_rate_2g", "PANSS_P_reduce_rate_2g", "PANSS_N_reduce_rate_2g", "PANSS_G_reduce_rate_2g",
              "PANSS_reduce_rate_2g_med", "PANSS_P_reduce_rate_2g_med", "PANSS_N_reduce_rate_2g_med", "PANSS_G_reduce_rate_2g_med")

res_list <- list()

for (outcome in outcomes) {
  
  formula_str <- paste0(outcome, " ~ grs_pred_tg + age + age2 + as.factor(sex) + as.factor(centers) + course + drugp07 + as.factor(medication) + PC1 + PC2 + PC3 + PC4 + PC5")
  
  fit <- try(glm(as.formula(formula_str), data = dat, family = binomial), silent = TRUE)
  
  if (!inherits(fit, "try-error")) {
    coef_summary <- summary(fit)$coefficients
    ci <- suppressMessages(confint(fit))  
    
    if ("grs_pred_tg" %in% rownames(coef_summary)) {
      beta <- coef_summary["grs_pred_tg", "Estimate"]
      or <- exp(beta)
      pval <- coef_summary["grs_pred_tg", "Pr(>|z|)"]
      lower <- exp(ci["grs_pred_tg", 1])
      upper <- exp(ci["grs_pred_tg", 2])
      
      res_list[[outcome]] <- data.frame(
        Trait = "grs_pred_tg",
        Outcome = outcome,
        OR = round(or, 4),
        CI = paste0(sprintf("%.3f", lower), " ~ ", sprintf("%.3f", upper)),
        P = signif(pval, 4)
      )
    }
  }
}

res_df_tg1_logistic <- bind_rows(res_list)
print(res_df_tg1_logistic)









#######
#grs pred tc 2 panss

####linear
outcomes <- c("PANSS_reduce_rate", "PANSS_P_reduce_rate", "PANSS_N_reduce_rate", "PANSS_G_reduce_rate")


res_list <- list()

for (outcome in outcomes) {
  formula_str <- paste0(outcome, " ~ grs_pred_tc + age + age2 + as.factor(sex) + as.factor(centers) + course + drugp07 + as.factor(medication) + PC1 + PC2 + PC3 + PC4 + PC5")
  fit <- try(lm(as.formula(formula_str), data = dat), silent = TRUE)
  
  if (!inherits(fit, "try-error")) {
    coef_summary <- summary(fit)$coefficients
    ci <- confint(fit)
    
    if ("grs_pred_tc" %in% rownames(coef_summary)) {
      beta <- coef_summary["grs_pred_tc", "Estimate"]
      pval <- coef_summary["grs_pred_tc", "Pr(>|t|)"]
      lower <- ci["grs_pred_tc", 1]
      upper <- ci["grs_pred_tc", 2]
      
      res_list[[outcome]] <- data.frame(
        Trait = "grs_pred_tc",
        Outcome = outcome,
        Beta = round(beta, 4),
        CI = paste0(sprintf("%.3f", lower), " ~ ", sprintf("%.3f", upper)),
        P = signif(pval, 4)
      )
    }
  }
}

res_df_tc1_lm <- bind_rows(res_list)
print(res_df_tc1_lm)



######logistic
outcomes <- c("PANSS_reduce_rate_2g", "PANSS_P_reduce_rate_2g", "PANSS_N_reduce_rate_2g", "PANSS_G_reduce_rate_2g",
              "PANSS_reduce_rate_2g_med", "PANSS_P_reduce_rate_2g_med", "PANSS_N_reduce_rate_2g_med", "PANSS_G_reduce_rate_2g_med")

res_list <- list()

for (outcome in outcomes) {
  
  formula_str <- paste0(outcome, " ~ grs_pred_tc + age + age2 + as.factor(sex) + as.factor(centers) + course + drugp07 + as.factor(medication) + PC1 + PC2 + PC3 + PC4 + PC5")
  
  fit <- try(glm(as.formula(formula_str), data = dat, family = binomial), silent = TRUE)
  
  if (!inherits(fit, "try-error")) {
    coef_summary <- summary(fit)$coefficients
    ci <- suppressMessages(confint(fit))  
    
    if ("grs_pred_tc" %in% rownames(coef_summary)) {
      beta <- coef_summary["grs_pred_tc", "Estimate"]
      or <- exp(beta)
      pval <- coef_summary["grs_pred_tc", "Pr(>|z|)"]
      lower <- exp(ci["grs_pred_tc", 1])
      upper <- exp(ci["grs_pred_tc", 2])
      
      res_list[[outcome]] <- data.frame(
        Trait = "grs_pred_tc",
        Outcome = outcome,
        OR = round(or, 4),
        CI = paste0(sprintf("%.3f", lower), " ~ ", sprintf("%.3f", upper)),
        P = signif(pval, 4)
      )
    }
  }
}

res_df_tc1_logistic <- bind_rows(res_list)
print(res_df_tc1_logistic)


res_df_lm <- rbind(res_df_glu1_lm, res_df_tg1_lm, res_df_tc1_lm)
res_df_logistic <- rbind(res_df_glu1_logistic, res_df_tg1_logistic, res_df_tc1_logistic)




write.csv(res_df_lm, 
          "/Users/zhuyunqing/Desktop/study/drugMR_lipid/data/onesamplemr/capoc.osmr.res_grs.pred.lg2panss_redu.linear.csv", 
          row.names = FALSE)

########
#stable 17, figure 3


write.csv(res_df_logistic, 
          "/Users/zhuyunqing/Desktop/study/drugMR_lipid/data/onesamplemr/capoc.osmr.res_grs.pred.lg2panss_redu.logistic.csv", 
          row.names = FALSE)

########
#stable 18



