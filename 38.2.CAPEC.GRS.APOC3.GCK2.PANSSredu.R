library(data.table)
library(broom)
library(dplyr)

dat <- fread("/Users/zhuyunqing/Desktop/study/drugMR_lipid/data/CAPEC/CAPEC.grs.pheno.new.txt")

# #######
# #glu
genes <- c("GCK")  # 你可以换成你实际的 GRS 基因名

res_list <- list()

for (gene in genes) {
  grs_var <- paste0("grs_", gene, "_w_taiwan_glu")
  
  if (!(grs_var %in% names(dat))) {
    message("Skipping: ", grs_var, " not in dat")
    next
  }
  
  formula_str <- paste0("v1_glu_mg ~ ", grs_var,
                        " + dat$age + dat$age2 + as.factor(dat$Sex) + as.factor(dat$center) + ",
                        "dat$PC1 + dat$PC2 + dat$PC3 + dat$PC4 + dat$PC5 + dat$FirstEpi")
  
  # 构建数据框（包含 dat 和 m 变量）
  dat_model <- cbind(dat[, ..grs_var], dat[, .(age, age2, Sex, center, PC1, PC2, PC3, PC4, PC5, FirstEpi)], v1_glu_mg = dat$v1_glu_mg)
  colnames(dat_model)[1] <- grs_var  # 重命名为原始名，便于模型使用
  
  # 模型拟合
  fit <- try(lm(as.formula(formula_str), data = dat_model), silent = TRUE)
  if (inherits(fit, "try-error")) next
  
  # 提取第一项（即 GRS）系数结果
  s <- summary(fit)
  ci <- confint(fit)
  
  # GRS 是模型中的第一个变量
  beta <- coef(fit)[2]
  lower <- ci[2, 1]
  upper <- ci[2, 2]
  pval <- coef(s)[2, 4]
  
  res_list[[gene]] <- data.frame(
    Gene = gene,
    Trait = "glu",
    Beta = beta,
    CI = paste0(sprintf("%.3f", lower), " ~ ", sprintf("%.3f", upper)),
    P = signif(pval, 3)
  )
}

res_df_glu <- bind_rows(res_list)
print(res_df_glu)






#######
#降脂药靶分析
genes   <- c("APOC3")
traits  <- c("v1_tg_mg", "v1_hdl_mg", "v1_ldl_mg", "v1_tc_mg")

res_list <- list()

for (gene in genes) {
  grs_var <- paste0("grs_", gene, "_w_glgc_logTG")
  if (!(grs_var %in% names(dat))) next
  
  for (trait in traits) {
    if (!(trait %in% names(dat))) next
    
    dat_model <- data.table(
      GRS        = dat[[grs_var]],
      age        = dat$age,
      age2       = dat$age2,
      Sex        = as.factor(dat$Sex),
      center     = as.factor(dat$center),
      PC1        = dat$PC1, PC2 = dat$PC2, PC3 = dat$PC3,
      PC4        = dat$PC4, PC5 = dat$PC5,
      FirstEpi   = dat$FirstEpi,
      Outcome    = dat[[trait]]
    )
    
    fit <- try(
      lm(Outcome ~ GRS + age + age2 + Sex + center +
           PC1 + PC2 + PC3 + PC4 + PC5 + FirstEpi,
         data = dat_model),
      silent = TRUE
    )
    if (inherits(fit, "try-error")) next
    
    s   <- summary(fit)
    ci  <- confint(fit)
    
    beta   <- coef(fit)["GRS"]
    lower  <- ci["GRS", 1]
    upper  <- ci["GRS", 2]
    pval   <- s$coefficients["GRS", 4]
    
    res_list[[paste(gene, trait, sep = "_")]] <-
      data.frame(
        Gene  = gene,
        Trait = trait,
        Beta  = beta,
        CI    = sprintf("%.3f ~ %.3f", lower, upper),
        P     = signif(pval, 3)
      )
  }
}

res_df_TG <- bind_rows(res_list)
print(res_df_TG)



#######
#TC
res_list <- list()

for (gene in genes) {
  grs_var <- paste0("grs_", gene, "_w_glgc_TC")
  if (!(grs_var %in% names(dat))) next
  
  for (trait in traits) {
    if (!(trait %in% names(dat))) next
    
    dat_model <- data.table(
      GRS        = dat[[grs_var]],
      age        = dat$age,
      age2       = dat$age2,
      Sex        = as.factor(dat$Sex),
      center     = as.factor(dat$center),
      PC1        = dat$PC1, PC2 = dat$PC2, PC3 = dat$PC3,
      PC4        = dat$PC4, PC5 = dat$PC5,
      FirstEpi   = dat$FirstEpi,
      Outcome    = dat[[trait]]
    )
    
    fit <- try(
      lm(Outcome ~ GRS + age + age2 + Sex + center +
           PC1 + PC2 + PC3 + PC4 + PC5 + FirstEpi,
         data = dat_model),
      silent = TRUE
    )
    if (inherits(fit, "try-error")) next
    
    s   <- summary(fit)
    ci  <- confint(fit)
    
    beta   <- coef(fit)["GRS"]
    lower  <- ci["GRS", 1]
    upper  <- ci["GRS", 2]
    pval   <- s$coefficients["GRS", 4]
    
    res_list[[paste(gene, trait, sep = "_")]] <-
      data.frame(
        Gene  = gene,
        Trait = trait,
        Beta  = beta,
        CI    = sprintf("%.3f ~ %.3f", lower, upper),
        P     = signif(pval, 3)
      )
  }
}

res_df_TC <- bind_rows(res_list)
print(res_df_TC)



final_resuls <- rbind(res_df_glu, res_df_TG, res_df_TC)
write.csv(final_resuls, "/Users/zhuyunqing/Desktop/study/drugMR_lipid/result/CAPEC.APOC3.GCK.GRS2glu.lipid.csv")

########
#stable 30




########
#grs 2 panss

# 假设 dat 是你的数据框，包含以下列：
# PANSS_reduce_rate, PANSS_P_reduce_rate, PANSS_N_reduce_rate, PANSS_G_reduce_rate, …
# grs_GCK_w_taiwan_glu, age, age2, Sex, center, PC1–PC5, FirstEpi, drug

outcomes <- c(
  "PANSS_reduce_rate",
  "PANSS_P_reduce_rate",
  "PANSS_N_reduce_rate",
  "PANSS_G_reduce_rate"
)



results_list <- lapply(outcomes, function(outcome) {
  fmla <- as.formula(paste0(
    outcome, " ~ grs_GCK_w_taiwan_glu + age + age2 + ",
    "factor(Sex) + factor(center) + ",
    "PC1 + PC2 + PC3 + PC4 + PC5 + ",
    "FirstEpi + factor(drug)"
  ))
  
  fit <- lm(fmla, data = dat)
  
  tidy(fit, conf.int = TRUE) %>%
    filter(term == "grs_GCK_w_taiwan_glu") %>%
    mutate(outcome = outcome)
})


final_results_glu <- bind_rows(results_list)

print(final_results_glu)





results_list <- lapply(outcomes, function(outcome) {
  fmla <- as.formula(paste0(
    outcome, " ~ grs_APOC3_w_glgc_logTG + age + age2 + ",
    "factor(Sex) + factor(center) + ",
    "PC1 + PC2 + PC3 + PC4 + PC5 + ",
    "FirstEpi + factor(drug)"
  ))
  
  fit <- lm(fmla, data = dat)
  
  tidy(fit, conf.int = TRUE) %>%
    filter(term == "grs_APOC3_w_glgc_logTG") %>%
    mutate(outcome = outcome)
})


final_results_tg <- bind_rows(results_list)

print(final_results_tg)



results_list <- lapply(outcomes, function(outcome) {
  fmla <- as.formula(paste0(
    outcome, " ~ grs_APOC3_w_glgc_TC + age + age2 + ",
    "factor(Sex) + factor(center) + ",
    "PC1 + PC2 + PC3 + PC4 + PC5 + ",
    "FirstEpi + factor(drug)"
  ))
  
  fit <- lm(fmla, data = dat)
  
  tidy(fit, conf.int = TRUE) %>%
    filter(term == "grs_APOC3_w_glgc_TC") %>%
    mutate(outcome = outcome)
})


final_results_tc <- bind_rows(results_list)

print(final_results_tc)


final_resuls <- rbind(final_results_glu, final_results_tg, final_results_tc)
write.csv(final_resuls, "/Users/zhuyunqing/Desktop/study/drugMR_lipid/result/CAPEC.APOC3.GCK.GRS2glu.PANSSredu.csv")


########
#stable 31






#######
#logistic
outcomes <- c(
  "PANSS_reduce_rate_2g",
  "PANSS_P_reduce_rate_2g",
  "PANSS_N_reduce_rate_2g",
  "PANSS_G_reduce_rate_2g",
  "PANSS_reduce_rate_2g_med",
  "PANSS_P_reduce_rate_2g_med",
  "PANSS_N_reduce_rate_2g_med",
  "PANSS_G_reduce_rate_2g_med"
)

# GCK GRS
results_list <- lapply(outcomes, function(outcome) {
  
  fmla <- as.formula(paste0(
    outcome, " ~ grs_GCK_w_taiwan_glu + age + age2 + ",
    "factor(Sex) + factor(center) + ",
    "PC1 + PC2 + PC3 + PC4 + PC5 + ",
    "FirstEpi + factor(drug)"
  ))
  
  fit <- tryCatch(
    glm(fmla, data = dat, family = binomial()),
    error = function(e) NULL
  )
  
  if (is.null(fit)) return(NULL)
  
  tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term == "grs_GCK_w_taiwan_glu") %>%
    mutate(outcome = outcome) %>%
    rename(
      OR = estimate,
      cilower = conf.low,
      ciupper = conf.high,
      p = p.value
    ) %>%
    select(outcome, OR, cilower, ciupper, p)
})

final_results_glu <- bind_rows(results_list)
print(final_results_glu)



# APOC3
results_list <- lapply(outcomes, function(outcome) {
  
  fmla <- as.formula(paste0(
    outcome, " ~ grs_APOC3_w_glgc_logTG + age + age2 + ",
    "factor(Sex) + factor(center) + ",
    "PC1 + PC2 + PC3 + PC4 + PC5 + ",
    "FirstEpi + factor(drug)"
  ))
  
  fit <- tryCatch(
    glm(fmla, data = dat, family = binomial()),
    error = function(e) NULL
  )
  
  if (is.null(fit)) return(NULL)
  
  tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term == "grs_APOC3_w_glgc_logTG") %>%
    mutate(outcome = outcome) %>%
    rename(
      OR = estimate,
      cilower = conf.low,
      ciupper = conf.high,
      p = p.value
    ) %>%
    select(outcome, OR, cilower, ciupper, p)
})

final_results_tg <- bind_rows(results_list)
print(final_results_tg)




# APOC3
results_list <- lapply(outcomes, function(outcome) {
  
  fmla <- as.formula(paste0(
    outcome, " ~ grs_APOC3_w_glgc_TC + age + age2 + ",
    "factor(Sex) + factor(center) + ",
    "PC1 + PC2 + PC3 + PC4 + PC5 + ",
    "FirstEpi + factor(drug)"
  ))
  
  fit <- tryCatch(
    glm(fmla, data = dat, family = binomial()),
    error = function(e) NULL
  )
  
  if (is.null(fit)) return(NULL)
  
  tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term == "grs_APOC3_w_glgc_TC") %>%
    mutate(outcome = outcome) %>%
    rename(
      OR = estimate,
      cilower = conf.low,
      ciupper = conf.high,
      p = p.value
    ) %>%
    select(outcome, OR, cilower, ciupper, p)
})

final_results_tc <- bind_rows(results_list)
print(final_results_tc)


final_resuls <- rbind(final_results_glu, final_results_tg, final_results_tc)
write.csv(final_resuls, "/Users/zhuyunqing/Desktop/study/drugMR_lipid/result/CAPEC.APOC3.GCK.GRS2glu.PANSSredu.logi.csv")


########
#stable 32


