library(data.table)
library(broom)
library(dplyr)

dat <- fread("/Users/zhuyunqing/Desktop/study/drugMR_lipid/data/onesamplemr/CAPOC.cis.capoc_cases.grs.pheno.new.txt", h=T)

# #######
# #glu
genes <- c("ABCB11", "GCK", "GLP1R")  # 你可以换成你实际的 GRS 基因名

res_list <- list()

for (gene in genes) {
  grs_var <- paste0("grs_", gene, "_w_taiwan_glu")

  if (!(grs_var %in% names(dat))) {
    message("Skipping: ", grs_var, " not in dat")
    next
  }

  formula_str <- paste0("glu1_mg ~ ", grs_var,
                        " + dat$age + dat$age2 + as.factor(dat$sex) + as.factor(dat$centers) + ",
                        "dat$PC1 + dat$PC2 + dat$PC3 + dat$PC4 + dat$PC5 + dat$course + dat$drugp07")

  # 构建数据框（包含 dat 和 m 变量）
  dat_model <- cbind(dat[, ..grs_var], dat[, .(age, age2, sex, centers, PC1, PC2, PC3, PC4, PC5, course, drugp07)], glu1_mg = dat$glu1_mg)
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

#######
#LDL

genes <- c("ABCA1", "PCSK9", "APOB", "LDLR", "LPA", "HMGCR")  # 你可以换成你实际的 GRS 基因名

res_list <- list()

for (gene in genes) {
  grs_var <- paste0("grs_", gene, "_w_glgc_ldl")

  if (!(grs_var %in% names(dat))) {
    message("Skipping: ", grs_var, " not in dat")
    next
  }

  formula_str <- paste0("ldl1_mg ~ ", grs_var,
                        " + dat$age + dat$age2 + as.factor(dat$sex) + as.factor(dat$centers) + ",
                        "dat$PC1 + dat$PC2 + dat$PC3 + dat$PC4 + dat$PC5 + dat$course + dat$drugp07")

  # 构建数据框（包含 dat 和 m 变量）
  dat_model <- cbind(dat[, ..grs_var], dat[, .(age, age2, sex, centers, PC1, PC2, PC3, PC4, PC5, course, drugp07)], ldl1_mg = dat$ldl1_mg)
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
    Trait = "LDL",
    Beta = beta,
    CI = paste0(sprintf("%.3f", lower), " ~ ", sprintf("%.3f", upper)),
    P = signif(pval, 3)
  )
}

res_df_LDL <- bind_rows(res_list)
print(res_df_LDL)


#######
#HDL
genes <- c("PPARG", "HCAR3", "HCAR2", "ABCA1", "MTTP", "APOC3", "CETP", 
           "APOB", "LDLR", "LPL", "LPA")  # 你可以换成你实际的 GRS 基因名

res_list <- list()

for (gene in genes) {
  grs_var <- paste0("grs_", gene, "_w_glgc_hdl")

  if (!(grs_var %in% names(dat))) {
    message("Skipping: ", grs_var, " not in dat")
    next
  }

  formula_str <- paste0("hdl1_mg ~ ", grs_var,
                        " + dat$age + dat$age2 + as.factor(dat$sex) + as.factor(dat$centers) + ",
                        "dat$PC1 + dat$PC2 + dat$PC3 + dat$PC4 + dat$PC5 + dat$course + dat$drugp07")

  # 构建数据框（包含 dat 和 m 变量）
  dat_model <- cbind(dat[, ..grs_var], dat[, .(age, age2, sex, centers, PC1, PC2, PC3, PC4, PC5, course, drugp07)], hdl1_mg = dat$hdl1_mg)
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
    Trait = "HDL",
    Beta = beta,
    CI = paste0(sprintf("%.3f", lower), " ~ ", sprintf("%.3f", upper)),
    P = signif(pval, 3)
  )
}

res_df_HDL <- bind_rows(res_list)
print(res_df_HDL)




#######
#logTG
genes <- c("ANGPTL3", "APOC3", "CETP", "APOB", "LPL")  # 你可以换成你实际的 GRS 基因名

res_list <- list()

for (gene in genes) {
  grs_var <- paste0("grs_", gene, "_w_glgc_logTG")

  if (!(grs_var %in% names(dat))) {
    message("Skipping: ", grs_var, " not in dat")
    next
  }

  formula_str <- paste0("tg1_mg ~ ", grs_var,
                        " + dat$age + dat$age2 + as.factor(dat$sex) + as.factor(dat$centers) + ",
                        "dat$PC1 + dat$PC2 + dat$PC3 + dat$PC4 + dat$PC5 + dat$course + dat$drugp07")

  # 构建数据框（包含 dat 和 m 变量）
  dat_model <- cbind(dat[, ..grs_var], dat[, .(age, age2, sex, centers, PC1, PC2, PC3, PC4, PC5, course, drugp07)], tg1_mg = dat$tg1_mg)
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
    Trait = "TG",
    Beta = beta,
    CI = paste0(sprintf("%.3f", lower), " ~ ", sprintf("%.3f", upper)),
    P = signif(pval, 3)
  )
}

res_df_TG <- bind_rows(res_list)
print(res_df_TG)






#######
#TC
genes <- c("ABCA1", "ANGPTL3", "APOB", "APOC3", "CETP", "HMGCR", "LDLR", 
           "LPA", "MTTP", "PCSK9")  

res_list <- list()

for (gene in genes) {
  grs_var <- paste0("grs_", gene, "_w_glgc_TC")

  if (!(grs_var %in% names(dat))) {
    message("Skipping: ", grs_var, " not in dat")
    next
  }

  formula_str <- paste0("tc1_mg ~ ", grs_var,
                        " + dat$age + dat$age2 + as.factor(dat$sex) + as.factor(dat$centers) + ",
                        "dat$PC1 + dat$PC2 + dat$PC3 + dat$PC4 + dat$PC5 + dat$course + dat$drugp07")

  # 构建数据框（包含 dat 和 m 变量）
  dat_model <- cbind(dat[, ..grs_var], dat[, .(age, age2, sex, centers, PC1, PC2, PC3, PC4, PC5, course, drugp07)], tc1_mg = dat$tc1_mg)
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
    Trait = "TC",
    Beta = beta,
    CI = paste0(sprintf("%.3f", lower), " ~ ", sprintf("%.3f", upper)),
    P = signif(pval, 3)
  )
}

res_df_TC <- bind_rows(res_list)
print(res_df_TC)

# 
res_df <- rbind(res_df_LDL, res_df_HDL, res_df_TG, res_df_TC, res_df_glu)

write.csv(res_df, "/Users/zhuyunqing/Desktop/study/drugMR_lipid/data/onesamplemr/res_target_weight_lgs_2lipid.glu.sbp.new.csv", row.names = FALSE)

##########
# STable 9, Figure 2
# 



