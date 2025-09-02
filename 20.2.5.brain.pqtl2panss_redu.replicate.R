#PQTL-脑

library(data.table)
library(dplyr)
library(ggplot2)
library(forestploter)

library(corrplot)
library(reshape2)
library(grid)
library(ggrepel)

library("rio")
#########
#在863人群中算分
pqtl <- fread("/Users/zhuyunqing/Desktop/study/drugMR_lipid/data/onesamplemr/CAPOC.cis.brain.pqtl1.txt",h=T)

fill_missing_with_mode <- function(df) {
  get_mode <- function(x) {
    mode_val <- names(sort(table(x), decreasing = TRUE))[1]
    return(mode_val)
  }
  for (col in colnames(df)) {
    mode_value <- get_mode(df[[col]][!is.na(df[[col]])])
    df[[col]][is.na(df[[col]])] <- as.numeric(mode_value)
  }
  return(df)
}

pqtl <- fill_missing_with_mode(pqtl)
str(pqtl)






pqtl$SERPINE1.wustl.grs <- -1 * scale((pqtl$rs1006507_T * (-0.098) + 
                                         pqtl$rs1050955_G * (-1*0.19) + 
                                         pqtl$rs11546151_G * (-1*0.36) + 
                                         pqtl$rs12334006_C * (-1*0.11) + 
                                         pqtl$rs7792190_C * (-1))/5)


pqtl$SERPINE1.simp.wustl.grs <- -1 * scale((
                                         pqtl$rs1050955_G * (-1*0.19))/1)


pqtl$GCK.wustl.grs <- scale(pqtl$rs3757840_T * (-1*0.073))
pqtl$APOC3.wustl.grs <- -1 * scale(pqtl$rs11827682_C * 0.18)

write.table(pqtl, 
            "/Users/zhuyunqing/Desktop/study/drugMR_lipid/data/pQTL/capoc/CAPOC.cis.brain.pqtl.grs.txt",
            row.names=FALSE, col.names=TRUE, quote=FALSE)





########## 
#与panss减分率关联分析


library(data.table)
library(broom)
library(dplyr)


capoc <- fread("/Users/zhuyunqing/Desktop/study/drugMR_lipid/data/onesamplemr/CAPOC.cis.capoc_cases.grs.pheno.new.txt", h=T)
csf.pqtl <- fread("/Users/zhuyunqing/Desktop/study/drugMR_lipid/data/pQTL/capoc/CAPOC.cis.brain.pqtl.grs.txt", h=T)
dat <- merge(capoc, csf.pqtl, by.x = "dn", by.y = "FID", all=F)


genes <- c("GCK", "APOC3")  # 后续可扩展为多个基因
trait <- "wustl.grs"
outcomes <- c("PANSS_reduce_rate", "PANSS_P_reduce_rate", "PANSS_N_reduce_rate", "PANSS_G_reduce_rate"
)

res_list <- list()

for (gene in genes) {
  pred_var <- paste0(gene, ".", trait)
  
  if (!(pred_var %in% names(dat))) {
    message("Skipping: ", pred_var, " not found")
    next
  }
  
  for (outcome in outcomes) {
    if (!(outcome %in% names(dat))) {
      message("Skipping: ", outcome, " not found in m")
      next
    }
    
    model_data <- cbind(
      y = dat[[outcome]],
      dat[, ..pred_var],
      dat[, .(age, age2, sex, centers, PC1, PC2, PC3, PC4, PC5, course, drugp07, medication)]
    )
    colnames(model_data)[2] <- "pred"
    
    fit <- try(lm(y ~ pred + age + age2 + as.factor(sex) + as.factor(centers) + PC1 + PC2 + PC3 + PC4 + PC5 + course + drugp07 + as.factor(medication), data = model_data), silent = TRUE)
    if (inherits(fit, "try-error")) next
    
    ci <- confint(fit)
    beta <- coef(fit)["pred"]
    lower <- ci["pred", 1]
    upper <- ci["pred", 2]
    pval <- summary(fit)$coefficients["pred", "Pr(>|t|)"]
    
    res_list[[paste(gene, outcome, sep = "_")]] <- data.frame(
      Gene = gene,
      Trait = trait,
      Outcome = outcome,
      Beta = beta,
      CI = paste0(sprintf("%.3f", lower), " ~ ", sprintf("%.3f", upper)),
      P = signif(pval, 3)
    )
  }
}

res_df_linear <- bind_rows(res_list)
print(res_df_linear)




outcomes <- c("PANSS_reduce_rate_2g", "PANSS_P_reduce_rate_2g", "PANSS_N_reduce_rate_2g", "PANSS_G_reduce_rate_2g",
              "PANSS_reduce_rate_2g_med", "PANSS_P_reduce_rate_2g_med", "PANSS_N_reduce_rate_2g_med", "PANSS_G_reduce_rate_2g_med")

res_list <- list()

for (gene in genes) {
  pred_var <- paste0(gene, ".", trait)
  
  if (!(pred_var %in% names(dat))) {
    message("Skipping: ", pred_var, " not found")
    next
  }
  
  for (outcome in outcomes) {
    if (!(outcome %in% names(dat))) {
      message("Skipping: ", outcome, " not found in m")
      next
    }
    
    model_data <- cbind(
      y = dat[[outcome]],
      dat[, ..pred_var],
      dat[, .(age, age2, sex, centers, PC1, PC2, PC3, PC4, PC5, course, drugp07, medication)]
    )
    colnames(model_data)[2] <- "pred"
    
    fit <- try(glm(y ~ pred + age + age2 + as.factor(sex) + as.factor(centers) + PC1 + PC2 + PC3 + PC4 + PC5 + course + drugp07 + as.factor(medication),
                   data = model_data, family = binomial), silent = TRUE)
    if (inherits(fit, "try-error")) next
    
    beta <- coef(fit)["pred"]
    or <- exp(beta)
    ci <- suppressMessages(confint(fit))
    lower <- exp(ci["pred", 1])
    upper <- exp(ci["pred", 2])
    pval <- summary(fit)$coefficients["pred", "Pr(>|z|)"]
    
    res_list[[paste(gene, outcome, sep = "_")]] <- data.frame(
      Gene = gene,
      Trait = trait,
      Outcome = outcome,
      OR = round(or, 3),
      CI = paste0(sprintf("%.3f", lower), " ~ ", sprintf("%.3f", upper)),
      P = signif(pval, 3)
    )
  }
}

res_df_logistic <- bind_rows(res_list)
print(res_df_logistic)


write.csv(res_df_linear, "/Users/zhuyunqing/Desktop/study/drugMR_lipid/result/brain.pQTL2panss/res_brain.pqtl.grs2panss_redu.linear.csv")
write.csv(res_df_logistic, "/Users/zhuyunqing/Desktop/study/drugMR_lipid/result/brain.pQTL2panss/res_brain.pqtl.grs2panss_redu.logistic.csv")


########
#stable 27-28,Figure 3

