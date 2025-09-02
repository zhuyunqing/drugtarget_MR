######
#apoc3 gck的药靶交互效应

library(data.table)
library(broom)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)

dat <- fread("/Users/zhuyunqing/Desktop/study/drugMR_lipid/data/onesamplemr/CAPOC.cis.capoc_cases.grs.pheno.new.txt", h=T)

dat$grs_APOC3_GCK_4g <- ifelse(
  dat$grs_APOC3_w_glgc_logTG >= median(dat$grs_APOC3_w_glgc_logTG) & 
    dat$grs_GCK_w_taiwan_glu >= median(dat$grs_GCK_w_taiwan_glu), 4, #服两种药
  ifelse(
    dat$grs_APOC3_w_glgc_logTG >= median(dat$grs_APOC3_w_glgc_logTG) & 
      dat$grs_GCK_w_taiwan_glu < median(dat$grs_GCK_w_taiwan_glu), 2, #服APOC3
    ifelse(
      dat$grs_APOC3_w_glgc_logTG < median(dat$grs_APOC3_w_glgc_logTG) & 
        dat$grs_GCK_w_taiwan_glu >= median(dat$grs_GCK_w_taiwan_glu), 3, #服GCK
      1
    )
  )
)

table(dat$grs_APOC3_GCK_4g)

# table(dat$grs_APOC3_GCK_4g)
# 
# 1   2-apoc3   3-gck   4 
# 499 554       555     503 






vars <- c("glu1_mg", "tg1_mg", 
          "PANSS_reduce_rate", "PANSS_P_reduce_rate", 
          "PANSS_N_reduce_rate", "PANSS_G_reduce_rate")

summary_tbl <- dat %>%                               # dat 为原数据框
  group_by(grs_APOC3_GCK_4g) %>%                     # 分四组
  summarise(
    across(all_of(vars),
           list(mean = ~ mean(., na.rm = TRUE),
                q1   = ~ quantile(., 0.25, na.rm = TRUE),
                q3   = ~ quantile(., 0.75, na.rm = TRUE)),
           .names = "{.col}.{.fn}"),
    .groups = "drop"
  ) %>%
  # ── 将 “变量.统计量” 拆成两列 ─────────────────────────────────────────────
  pivot_longer(
    -grs_APOC3_GCK_4g,
    names_to      = c("variable", "stat"),
    names_pattern = "^(.*)\\.(mean|q1|q3)$"           # 只匹配 .mean / .q1 / .q3
  ) %>%
  pivot_wider(names_from = stat, values_from = value) %>%  # 回到宽格式：mean / q1 / q3 三列
  mutate(fmt = sprintf("%.1f (%.1f, %.1f)", mean, q1, q3)) %>%
  select(grs_APOC3_GCK_4g, variable, fmt) %>%
  pivot_wider(
    names_from  = grs_APOC3_GCK_4g,
    values_from = fmt,
    names_sort  = TRUE
  ) %>%
  arrange(match(variable, vars))                      # 保持原始变量顺序

print(summary_tbl, row.names = FALSE)

write.csv(summary_tbl, "/Users/zhuyunqing/Desktop/study/drugMR_lipid/result/des.apoc3.gck.4g.capoc.csv")

########
#stable 35. figure 5






########
#GRS4G TO PANSS


dat <- dat %>% mutate(
  grs_APOC3_GCK_4g = factor(grs_APOC3_GCK_4g, levels = 1:4)
)

outcomes <- c("glu1_mg", "tg1_mg", 
  "PANSS_reduce_rate",
  "PANSS_P_reduce_rate",
  "PANSS_N_reduce_rate",
  "PANSS_G_reduce_rate"
)

covars <- paste(
  "age", "age2",
  "as.factor(sex)", "as.factor(centers)",
  "PC1", "PC2", "PC3", "PC4", "PC5",
  "course", "drugp07",
  "as.factor(medication)",
  sep = " + "
)

res_all <- map_dfr(outcomes, function(outcome) {
  
  fml <- as.formula(
    paste(outcome, "~ as.factor(grs_APOC3_GCK_4g) +", covars)
  )
  
  fit <- lm(fml, data = dat)
  
  tidy(fit, conf.int = TRUE) %>%
    filter(grepl("^as\\.factor\\(grs_APOC3_GCK_4g\\)", term)) %>%
    transmute(
      outcome  = outcome,       
      level    = term,         
      beta     = estimate,
      lower95  = conf.low,
      upper95  = conf.high,
      p_value  = p.value
    )
})

print(res_all)
write.csv(res_all, 
          "/Users/zhuyunqing/Desktop/study/drugMR_lipid/result/GRS4g_PANSS_betas_95CI_p.csv", row.names = FALSE)



########
#stable 36,38, figure 5





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

res_all <- map_dfr(outcomes, function(outcome) {
  
  fml <- as.formula(
    paste(outcome, "~ as.factor(grs_APOC3_GCK_4g) +", covars)
  )
  
  fit <- glm(fml, data = dat, family = binomial)
  
  tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(grepl("^as\\.factor\\(grs_APOC3_GCK_4g\\)", term)) %>%
    transmute(
      outcome  = outcome,
      level    = term,
      OR       = estimate,
      lower95  = conf.low,
      upper95  = conf.high,
      p_value  = p.value
    )
})

# 输出结果
print(res_all)
write.csv(
  res_all,
  "/Users/zhuyunqing/Desktop/study/drugMR_lipid/result/GRS4g_PANSS_OR_95CI_p.csv",
  row.names = FALSE
)

########
#stable 39





#########
#interactive analysis
outcomes <- c("glu1_mg", "tg1_mg",
  "PANSS_reduce_rate",
  "PANSS_P_reduce_rate",
  "PANSS_N_reduce_rate",
  "PANSS_G_reduce_rate"
)

covars <- paste(
  "age", "age2",
  "as.factor(sex)", "as.factor(centers)",
  "PC1", "PC2", "PC3", "PC4", "PC5",
  "course", "drugp07",
  "as.factor(medication)",
  sep = " + "
)

res_all <- map_dfr(outcomes, function(outc) {
  
  fml <- as.formula(
    paste(outc,
          "~ grs_APOC3_w_glgc_logTG * grs_GCK_w_taiwan_glu +",
          covars)
  )
  
  fit <- lm(fml, data = dat)
  
  # 把三条系数 β / 95% CI / P 摘出来
  tidy(fit, conf.int = TRUE) %>%
    filter(term %in% c(
      "grs_APOC3_w_glgc_logTG",
      "grs_GCK_w_taiwan_glu",
      "grs_APOC3_w_glgc_logTG:grs_GCK_w_taiwan_glu"
    )) %>%
    transmute(
      outcome  = outc,                      # ← 新增列
      term     = term,
      beta     = estimate,
      lower95  = conf.low,
      upper95  = conf.high,
      p_value  = p.value
    )
})

print(res_all)

write.csv(res_all, 
          "/Users/zhuyunqing/Desktop/study/drugMR_lipid/result/apoc.gck.grs_PANSS_interaction.csv", 
          row.names = FALSE)

########
#stable 37,40, figure 5






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

res_all <- map_dfr(outcomes, function(outc) {
  
  fml <- as.formula(
    paste(
      outc,
      "~ grs_APOC3_w_glgc_logTG * grs_GCK_w_taiwan_glu +", 
      covars
    )
  )
  
  # 逻辑回归；估计量直接指数化为 OR
  fit <- glm(fml, data = dat, family = binomial)
  
  tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%          # exponentiate=TRUE → OR
    filter(term %in% c(
      "grs_APOC3_w_glgc_logTG",
      "grs_GCK_w_taiwan_glu",
      "grs_APOC3_w_glgc_logTG:grs_GCK_w_taiwan_glu"
    )) %>%
    transmute(
      outcome = outc,
      term    = term,
      OR      = estimate,
      lower95 = conf.low,
      upper95 = conf.high,
      p_value = p.value
    )
})

print(res_all)

write.csv(
  res_all,
  "/Users/zhuyunqing/Desktop/study/drugMR_lipid/result/GRS_interaction_OR_95CI_p.csv",
  row.names = FALSE
)

########
#stable 41







######
#apoc3_TC gck的药靶交互效应

library(data.table)
library(broom)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)

dat <- fread("/Users/zhuyunqing/Desktop/study/drugMR_lipid/data/onesamplemr/CAPOC.cis.capoc_cases.grs.pheno.new.txt", h=T)


#######
#APOC3_TC * GCK_TG
dat$grs_APOC3_GCK_4g <- ifelse(
  dat$grs_APOC3_w_glgc_TC >= median(dat$grs_APOC3_w_glgc_TC) & 
    dat$grs_GCK_w_taiwan_glu >= median(dat$grs_GCK_w_taiwan_glu), 4, #服两种药
  ifelse(
    dat$grs_APOC3_w_glgc_TC >= median(dat$grs_APOC3_w_glgc_TC) & 
      dat$grs_GCK_w_taiwan_glu < median(dat$grs_GCK_w_taiwan_glu), 2, #服APOC3
    ifelse(
      dat$grs_APOC3_w_glgc_TC < median(dat$grs_APOC3_w_glgc_TC) & 
        dat$grs_GCK_w_taiwan_glu >= median(dat$grs_GCK_w_taiwan_glu), 3, #服GCK
      1
    )
  )
)

table(dat$grs_APOC3_GCK_4g)

# > table(dat$grs_APOC3_GCK_4g)
# 
# 1   2   3   4 
# 428 625 484 574 



vars <- c("glu1_mg", "tc1_mg", 
          "PANSS_reduce_rate", "PANSS_P_reduce_rate", 
          "PANSS_N_reduce_rate", "PANSS_G_reduce_rate")

summary_tbl <- dat %>%                               
  group_by(grs_APOC3_GCK_4g) %>%                     
  summarise(
    across(all_of(vars),
           list(mean = ~ mean(., na.rm = TRUE),
                q1   = ~ quantile(., 0.25, na.rm = TRUE),
                q3   = ~ quantile(., 0.75, na.rm = TRUE)),
           .names = "{.col}.{.fn}"),
    .groups = "drop"
  ) %>%
  
  pivot_longer(
    -grs_APOC3_GCK_4g,
    names_to      = c("variable", "stat"),
    names_pattern = "^(.*)\\.(mean|q1|q3)$"           
  ) %>%
  pivot_wider(names_from = stat, values_from = value) %>%  
  mutate(fmt = sprintf("%.1f (%.1f, %.1f)", mean, q1, q3)) %>%
  select(grs_APOC3_GCK_4g, variable, fmt) %>%
  pivot_wider(
    names_from  = grs_APOC3_GCK_4g,
    values_from = fmt,
    names_sort  = TRUE
  ) %>%
  arrange(match(variable, vars))                     

print(summary_tbl, row.names = FALSE)

write.csv(summary_tbl, "/Users/zhuyunqing/Desktop/study/drugMR_lipid/result/des.apoc3_tc.gck.4g.capoc.csv")

########
#stable 35






########
#GRS4G TO PANSS

library(dplyr)
library(broom)
library(purrr)



dat <- dat %>% mutate(
  grs_APOC3_GCK_4g = factor(grs_APOC3_GCK_4g, levels = 1:4)
)

outcomes <- c("glu1_mg", "tc1_mg", 
              "PANSS_reduce_rate",
              "PANSS_P_reduce_rate",
              "PANSS_N_reduce_rate",
              "PANSS_G_reduce_rate"
)

covars <- paste(
  "age", "age2",
  "as.factor(sex)", "as.factor(centers)",
  "PC1", "PC2", "PC3", "PC4", "PC5",
  "course", "drugp07",
  "as.factor(medication)",
  sep = " + "
)

res_all <- map_dfr(outcomes, function(outcome) {
  
  fml <- as.formula(
    paste(outcome, "~ as.factor(grs_APOC3_GCK_4g) +", covars)
  )
  
  fit <- lm(fml, data = dat)
  
  tidy(fit, conf.int = TRUE) %>%
    filter(grepl("^as\\.factor\\(grs_APOC3_GCK_4g\\)", term)) %>%
    transmute(
      outcome  = outcome,       
      level    = term,         
      beta     = estimate,
      lower95  = conf.low,
      upper95  = conf.high,
      p_value  = p.value
    )
})

print(res_all)
write.csv(res_all, 
          "/Users/zhuyunqing/Desktop/study/drugMR_lipid/result/GRS4g_TC_PANSS_betas_95CI_p.csv", row.names = FALSE)

########
#stable 36,38





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

res_all <- map_dfr(outcomes, function(outcome) {
  
  fml <- as.formula(
    paste(outcome, "~ as.factor(grs_APOC3_GCK_4g) +", covars)
  )
  
  fit <- glm(fml, data = dat, family = binomial)
  
  tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(grepl("^as\\.factor\\(grs_APOC3_GCK_4g\\)", term)) %>%
    transmute(
      outcome  = outcome,
      level    = term,
      OR       = estimate,
      lower95  = conf.low,
      upper95  = conf.high,
      p_value  = p.value
    )
})

# 输出结果
print(res_all)
write.csv(
  res_all,
  "/Users/zhuyunqing/Desktop/study/drugMR_lipid/result/GRS4g_TC_PANSS_OR_95CI_p.csv",
  row.names = FALSE
)

########
#stable 39







#########
#interactive analysis
outcomes <- c("glu1_mg", "tc1_mg",
              "PANSS_reduce_rate",
              "PANSS_P_reduce_rate",
              "PANSS_N_reduce_rate",
              "PANSS_G_reduce_rate"
)

covars <- paste(
  "age", "age2",
  "as.factor(sex)", "as.factor(centers)",
  "PC1", "PC2", "PC3", "PC4", "PC5",
  "course", "drugp07",
  "as.factor(medication)",
  sep = " + "
)

res_all <- map_dfr(outcomes, function(outc) {
  
  fml <- as.formula(
    paste(outc,
          "~ grs_APOC3_w_glgc_TC * grs_GCK_w_taiwan_glu +",
          covars)
  )
  
  fit <- lm(fml, data = dat)
  
  # 把三条系数 β / 95% CI / P 摘出来
  tidy(fit, conf.int = TRUE) %>%
    filter(term %in% c(
      "grs_APOC3_w_glgc_TC",
      "grs_GCK_w_taiwan_glu",
      "grs_APOC3_w_glgc_TC:grs_GCK_w_taiwan_glu"
    )) %>%
    transmute(
      outcome  = outc,                      # ← 新增列
      term     = term,
      beta     = estimate,
      lower95  = conf.low,
      upper95  = conf.high,
      p_value  = p.value
    )
})

print(res_all)

write.csv(res_all, 
          "/Users/zhuyunqing/Desktop/study/drugMR_lipid/result/apoc3_TC.gck.grs_PANSS_interaction.csv", 
          row.names = FALSE)

########
#stable 37,40







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

res_all <- map_dfr(outcomes, function(outc) {
  
  fml <- as.formula(
    paste(
      outc,
      "~ grs_APOC3_w_glgc_TC * grs_GCK_w_taiwan_glu +", 
      covars
    )
  )
  
  # 逻辑回归；估计量直接指数化为 OR
  fit <- glm(fml, data = dat, family = binomial)
  
  tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%          # exponentiate=TRUE → OR
    filter(term %in% c(
      "grs_APOC3_w_glgc_TC",
      "grs_GCK_w_taiwan_glu",
      "grs_APOC3_w_glgc_TC:grs_GCK_w_taiwan_glu"
    )) %>%
    transmute(
      outcome = outc,
      term    = term,
      OR      = estimate,
      lower95 = conf.low,
      upper95 = conf.high,
      p_value = p.value
    )
})

print(res_all)

write.csv(
  res_all,
  "/Users/zhuyunqing/Desktop/study/drugMR_lipid/result/GRS_TC_interaction_OR_95CI_p.csv",
  row.names = FALSE
)

########
#stable 41
