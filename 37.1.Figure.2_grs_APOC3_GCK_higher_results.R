library(data.table)
library(ggplot2)
library(forestploter)
library(broom)
library(dplyr)
library(corrplot)
library(reshape2)
library(grid)
library(ggrepel)
library(scales)
library("rio")
library(stringr)
library(patchwork)
library(tidyr)

setwd("/Users/zhuyunqing/Desktop/study/drugMR_lipid/data/onesamplemr")
dat <- fread("/Users/zhuyunqing/Desktop/study/drugMR_lipid/data/onesamplemr/CAPOC.cis.capoc_cases.grs.pheno.new.txt", h=T)
res <- fread("res_target_weight_lgs_2lipid.glu.sbp.new.csv")
result_APOC3 <- subset(res, Gene == "APOC3")

dat_glu <- dat[, .(
  grs_APOC3_w_glgc_logTG,
  age, age2, sex, centers,
  PC1, PC2, PC3, PC4, PC5,
  course, drugp07,
  glu1_mg
)]

fit_glu <- lm(
  glu1_mg ~ grs_APOC3_w_glgc_logTG +
    age + age2 + factor(sex) + factor(centers) +
    PC1 + PC2 + PC3 + PC4 + PC5 +
    course + drugp07,
  data = dat_glu
)

beta_glu  <- coef(fit_glu)[2]
ci_glu    <- confint(fit_glu)[2, ]
p_glu     <- summary(fit_glu)$coeff[2, 4]

new_row <- data.frame(
  Gene  = "APOC3",
  Trait = "GLU",
  Beta  = beta_glu,
  CI    = sprintf("%.3f ~ %.3f", ci_glu[1], ci_glu[2]),
  P     = signif(p_glu, 3)
)

result_APOC3 <- dplyr::bind_rows(result_APOC3, new_row)
print(result_APOC3)

result_APOC3 <- result_APOC3 %>%
  separate(CI, into = c("CI_Lower", "CI_Upper"), sep = "\\s*~\\s*", convert = TRUE)



fit_ldl <- lm(ldl1_mg ~ grs_APOC3_w_glgc_logTG + age + I(age^2) +
                factor(sex) + factor(centers) +
                PC1 + PC2 + PC3 + PC4 + PC5 +
                course + drugp07,
              data = dat)

# 提取该项的结果
ldl_result <- tidy(fit_ldl, conf.int = TRUE) %>%
  filter(term == "grs_APOC3_w_glgc_logTG") %>%
  transmute(
    Gene     = "APOC3",
    Trait    = "LDL",
    Beta     = estimate,
    CI_Lower = conf.low,
    CI_Upper = conf.high,
    P        = p.value
  )

# 合并进已有结果表
result_APOC3 <- bind_rows(result_APOC3, ldl_result)
result_APOC3$Trait <- c("HDL", "TG", "TC", "Glucose", "LDL")




#####
#GCK
result_GCK <- subset(res, Gene == "GCK")

library(dplyr)

models <- c(
  TC       = "tc1_mg",
  TG       = "tg1_mg",
  LDL      = "ldl1_mg",
  HDL      = "hdl1_mg"
)

grs_var <- "grs_GCK_w_taiwan_glu"
add_res <- list()

for (mdl in names(models)) {
  out_var <- models[[mdl]]
  df <- cbind(
    dat[, ..grs_var],
    dat[, .(age, age2, sex, centers)],
    dat[,  .(PC1, PC2, PC3, PC4, PC5, course, drugp07)],
    outcome = dat[[out_var]]
  )
  colnames(df)[1] <- grs_var
  colnames(df)[ncol(df)] <- out_var
  
  fit <- lm(
    as.formula(paste0(out_var, " ~ ", grs_var,
                      " + age + age2 + factor(sex) + factor(centers) + ",
                      "PC1 + PC2 + PC3 + PC4 + PC5 + course + drugp07")),
    data = df
  )
  
  ci <- confint(fit)[2, ]
  add_res[[mdl]] <- data.frame(
    Gene = "GCK",
    Trait    = mdl,
    Beta     = coef(fit)[2],
    CI_Lower = ci[1],
    CI_Upper = ci[2],
    P  = summary(fit)$coeff[2, 4]
  )
}

result_GCK <- bind_rows(result_GCK, bind_rows(add_res))

result_GCK <- result_GCK %>% 
  mutate(
    CI_Lower = ifelse(is.na(CI_Lower) & !is.na(CI),
                      as.numeric(str_extract(CI, "[-0-9.]+")), CI_Lower),
    CI_Upper = ifelse(is.na(CI_Upper) & !is.na(CI),
                      as.numeric(str_extract(CI, "(?<=~ )[-0-9.]+")), CI_Upper)
  ) %>% 
  select(-CI)

print(result_GCK)

result_GCK$Trait <- c("Glucose", "TC", "TG", "LDL", "HDL")




# dat <- fread("/Users/zhuyunqing/Desktop/study/drugMR_lipid/data/onesamplemr/CAPOC.cis.capoc_cases.grs.pheno.new.txt", h=T)
# 
# # 需要分析的连续结局变量（可按需增删、改名）
# response_vars <- c("tc1_mg", "tg1_mg", "ldl1_mg", "hdl1_mg", "glu1_mg")  # 示例
# 
# results_list <- list()
# 
# for (resp in response_vars) {
#   formula <- as.formula(
#     paste(resp, "~ grs_GCK_w_taiwan_glu + age + I(age^2) +",
#           "as.factor(sex) + as.factor(centers) +",
#           "PC1 + PC2 + PC3 + PC4 + PC5 + course + drugp07")
#   )
#   model <- lm(formula, data = dat)
#   
#   res <- tidy(model, conf.int = TRUE) %>%
#     filter(term == "grs_GCK_w_taiwan_glu") %>%
#     mutate(Model = resp,
#            Beta = estimate,
#            CI_Lower = conf.low,
#            CI_Upper = conf.high,
#            P_value = p.value) %>%
#     select(Model, Beta, CI_Lower, CI_Upper, P_value)
#   
#   results_list[[resp]] <- res
# }
# 
# results_GCK <- bind_rows(results_list)
# results_GCK$Model <- c("TC", "TG", "LDLC", "HDLC", "Glucose")
# print(results_GCK)
# 
# 
# #####
# #APOC3
# results_list <- list()
# 
# for (resp in response_vars) {
#   formula <- as.formula(
#     paste(resp, "~ grs_APOC3_w_glgc_logTG + age + I(age^2) +",
#           "as.factor(sex) + as.factor(centers) +",
#           "PC1 + PC2 + PC3 + PC4 + PC5 + course + drugp07")
#   )
#   model <- lm(formula, data = dat)
#   
#   res <- tidy(model, conf.int = TRUE) %>%
#     filter(term == "grs_APOC3_w_glgc_logTG") %>%
#     mutate(Model = resp,
#            Beta = estimate,
#            CI_Lower = conf.low,
#            CI_Upper = conf.high,
#            P_value = p.value) %>%
#     select(Model, Beta, CI_Lower, CI_Upper, P_value)
#   
#   results_list[[resp]] <- res
# }
# 
# results_APOC3 <- bind_rows(results_list)
# results_APOC3$Model <- c("TC", "TG", "LDLC", "HDLC", "Glucose")
# print(results_APOC3)





######
#画图

# 图 p1：APOC3
p1 <- ggplot(result_APOC3, aes(x = Trait, y = Beta)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(size = 3, color = "#2b5461") +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper),
                width = 0.2, size = 0.7, color = "#2b5461") +
  scale_y_continuous(limits = c(-9, 1), breaks = seq(-8, 0, 2)) +
  labs(x = NULL, y = "Difference in concentration, mg/dL") +
  theme_minimal(base_size = 14) +
  theme(
    panel.background   = element_blank(),
    panel.grid         = element_blank(),
    panel.border       = element_rect(color = "black", fill = NA),
    axis.title.y       = element_text(face = "bold", size = 11.5),
    axis.text.x        = element_text(face = "bold", size = 11.5),
    axis.ticks.y       = element_line(color = "black"),
    axis.ticks.length  = unit(4, "pt"),
    legend.position    = "none"
  )

p1

# 图 p2：GCK
p2 <- ggplot(result_GCK, aes(x = Trait, y = Beta)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(size = 3, color = "#8b8474") +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper),
                width = 0.2, size = 0.7, color = "#8b8474") +
  scale_y_continuous(limits = c(-3.5, 2.5), breaks = seq(-3, 2, 1)) +
  #scale_x_discrete(limits = c("TC", "TG", "LDL", "HDL", "Glucose")) +
  labs(x = NULL, y = "Difference in concentration, mg/dL") +
  theme_minimal(base_size = 14) +
  theme(
    panel.background   = element_blank(),
    panel.grid         = element_blank(),
    panel.border       = element_rect(color = "black", fill = NA),
    axis.title.y       = element_text(face = "bold", size = 11.5),
    axis.text.x        = element_text(face = "bold", size = 11.5),
    axis.ticks.y       = element_line(color = "black"),
    axis.ticks.length  = unit(4, "pt"),
    legend.position    = "none"
  )

p2

# 4. 拼接并导出 PDF
combined <- p1 / p2       # 上下排列
combined
ggsave("/Users/zhuyunqing/Desktop/study/drugMR_lipid/result/Figure.2_grs_APOC3_GCK_higher_results.pdf",
       combined,
       width  = 7.54,        # 总宽度（英寸）
       height = 6.59,       # 总高度（英寸）
       dpi    = 300)


