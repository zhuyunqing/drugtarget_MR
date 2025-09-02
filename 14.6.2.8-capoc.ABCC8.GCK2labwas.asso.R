#进行labwas分析
library(data.table)
library(broom)
library(dplyr)

dat <- fread("/Users/zhuyunqing/Desktop/study/drugMR_lipid/data/onesamplemr/CAPOC.cis.capoc_cases.grs.pheno.new.txt", h=T)

# 先定义所有涉及的变量名（包括 sbp、dbp、bmi 等）
base_vars <- c("sbp0601", "sbp3101", "dbp0601", "dbp3101",
               "pls0601", "pls3101", "bmi0601", "bmi3101",
               "fw0601", "fw3101", "qtc3201", "qtc0801")

markers <- c("alt", "ast", "bun", "cre", "tc", 
             "tg", "hdl", "ldl", "glu", "prl", "hba")

marker_vars <- c(paste0(markers, "1"), paste0(markers, "3"))

all_vars <- c(base_vars, marker_vars)

for (var in all_vars) {
  dat[[var]] <- as.numeric(dat[[var]])
}

dat$sbp_change <- (dat$sbp3101 - dat$sbp0601) / dat$sbp0601
dat$dbp_change <- (dat$dbp3101 - dat$dbp0601) / dat$dbp0601
dat$pls_change <- (dat$pls3101 - dat$pls0601) / dat$pls0601
dat$bmi_change <- (dat$bmi3101 - dat$bmi0601) / dat$bmi0601
dat$fw_change  <- (dat$fw3101  - dat$fw0601)  / dat$fw0601
dat$qtc_change  <- (dat$qtc3201  - dat$qtc0801)  / dat$qtc0801

for (marker in markers) {
  var1 <- paste0(marker, "1")
  var3 <- paste0(marker, "3")
  var_change <- paste0(marker, "_change")
  
  dat[[var_change]] <- (dat[[var3]] - dat[[var1]]) / dat[[var1]]
}


#########
#关联分析
library(data.table)
library(broom)  # 用于提取模型结果 tidy()

outcomes <- c("sbp0601", "dbp0601", "pls0601", "bmi0601", "alt1", "ast1", 
              "bun1", "cre1", "tc1", "tg1", "hdl1", "ldl1", "prl1", "hba1", 
              "qtc0801")

genes <- c("grs_APOC3_w_glgc_logTG", "grs_GCK_w_taiwan_glu")


# 初始化结果列表
res_list <- list()

# 主循环：遍历 outcome × gene
for (y in outcomes) {
  for (gene in genes) {
    grs_var <- paste0(gene)
    
    # 确保变量存在
    if (!all(c(y, grs_var) %in% names(dat))) next
    
    # 所有模型变量
    vars_needed <- c(y, grs_var, "age", "age2", "sex", "centers",  "PC1", "PC2", "PC3", "PC4", "PC5", "course", "drugp07")
    
    # 保留所有变量都是 finite 的行（防 NA / NaN / Inf）
    is_valid <- dat[, lapply(.SD, is.finite), .SDcols = vars_needed]
    valid_rows <- apply(is_valid, 1, all)
    dat_sub <- dat[valid_rows, ..vars_needed]
    
    # 若样本太少则跳过
    if (nrow(dat_sub) < 10) next
    
    # 构建回归公式
    formula_str <- paste0(
      y, " ~ ", grs_var,
      " + age + age2 + as.factor(sex) + as.factor(centers) + PC1 + PC2 + PC3 + PC4 + PC5 + course + drugp07"
    )
    
    # 安全建模
    tryCatch({
      model <- lm(as.formula(formula_str), data = dat_sub)
      coef_info <- tidy(model)
      
      row_grs <- coef_info[coef_info$term == grs_var, ]
      if (nrow(row_grs) == 1) {
        beta <- row_grs$estimate
        se <- row_grs$std.error
        p <- row_grs$p.value
        ci_low <- beta - 1.96 * se
        ci_high <- beta + 1.96 * se
        
        res_list[[length(res_list) + 1]] <- data.table(
          outcome = y,
          gene = gene,
          beta = round(beta, 4),
          ci = sprintf("%.3f ~ %.3f", ci_low, ci_high),
          pval = signif(p, 3),
          n = nrow(dat_sub)
        )
      }
    }, error = function(e) {
      message("Model failed for ", y, " ~ ", grs_var, ": ", e$message)
    })
  }
}

# 合并输出结果
res_df <- rbindlist(res_list)

# 打印结果或保存
print(res_df)
write.csv(res_df, "/Users/zhuyunqing/Desktop/study/drugMR_lipid/result/apoc3_gck_labwas_regression_results.csv")

########
#STable 10




