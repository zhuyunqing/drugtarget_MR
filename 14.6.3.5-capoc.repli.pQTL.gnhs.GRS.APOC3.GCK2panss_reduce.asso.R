#重新清洗数据
#863数据预处理

library(data.table)
library(dplyr)
library(ggplot2)
library(forestploter)

library(corrplot)
library(reshape2)
library(grid)
library(ggrepel)

library("rio")



#############
#grs 2 panss
###############
m <- fread("/Users/zhuyunqing/Desktop/study/drugMR_lipid/data/onesamplemr/CAPOC.cis.capoc_cases.grs.pheno.new.txt", h=T)


##############
#grs-panss associations
###############
# 定义变量名列表
panss_vars <- c("PANSS_reduce_rate", "PANSS_P_reduce_rate", "PANSS_N_reduce_rate", "PANSS_G_reduce_rate")  # PANSS变量名称
genes <- c("apoc3_w_gnhs")          


# 遍历每个基因和PANSS变量的组合
for (gene in genes) {
  for (panss in panss_vars) {
    
    tryCatch({
      # 构建公式和模型名称
      formula1 <- as.formula(paste0("m$", panss, " ~ m$grs_", gene))  # 未调整模型
      model_name1 <- paste0("reg_", gene, "_", panss)                       # 例如：reg_apoa_PANSS1
      
      # 创建未调整模型
      assign(model_name1, lm(formula1, data = m))
      print(paste("Model:", model_name1))
      print(summary(get(model_name1)))
      
      # 构建调整后的公式和模型名称
      formula2 <- as.formula(paste0("m$", panss, " ~ m$grs_", gene, " + m$age + m$age2 + as.factor(m$sex) + m$course + m$drugp07 + m$medication + m$PC1 + m$PC2 + m$PC3 + m$PC4 + m$PC5"))
      model_name2 <- paste0("reg_", gene, "_", panss, "_adj")               # 例如：reg_apoa_PANSS1_adj
      
      # 创建调整后的模型
      assign(model_name2, lm(formula2, data = m))
      print(paste("Model:", model_name2))
      print(summary(get(model_name2)))
      
      
      cat(paste0("Successfully processed: ", gene, "\n"))
    }, error = function(e) {
      # 错误处理
      cat(paste0("Error in processing: ", gene, " - ", e$message, "\n"))
    })
    
  }
}





# 创建一个空列表来保存结果
results_asso <- list()

# 遍历每个基因和 PANSS 变量的组合
for (gene in genes) {
  for (panss in panss_vars) {
    
    tryCatch({
      
      # 构建模型名称
      model_name1 <- paste0("reg_", gene, "_", panss)            # 例如：reg_apoa_PANSS1
      model_name1_adj <- paste0("reg_", gene, "_", panss, "_adj") # 例如：reg_apoa_PANSS1_adj
      
      # 提取未调整模型的结果
      res_name1 <- paste0("res_", model_name1)  # 结果变量名称
      model1 <- get(model_name1)                # 获取模型对象
      res1 <- c(model1$coefficients[2], confint(model1)[2,], 
                summary(model1)$coefficients[paste0("m$grs_", gene), "Pr(>|t|)"])
      assign(res_name1, res1)
      print(paste("Results for", res_name1))
      print(res1)
      
      # 提取调整后模型的结果
      res_name1_adj <- paste0("res_", model_name1_adj)  # 结果变量名称
      model1_adj <- get(model_name1_adj)                # 获取调整后模型对象
      res1_adj <- c(model1_adj$coefficients[2], confint(model1_adj)[2,], 
                    summary(model1_adj)$coefficients[paste0("m$grs_", gene), "Pr(>|t|)"])
      assign(res_name1_adj, res1_adj)
      print(paste("Results for", res_name1_adj))
      print(res1_adj)
      
      # 将结果保存到列表中
      results_asso[[res_name1]] <- res1
      results_asso[[res_name1_adj]] <- res1_adj
      
      
      cat(paste0("Successfully processed: ", gene, "\n"))
    }, error = function(e) {
      # 错误处理
      cat(paste0("Error in processing: ", gene, " - ", e$message, "\n"))
    })
    
    
  }
}







# 定义变量
adjustments <- c("", "_adj")  # 空代表未调整，_adj代表调整后

# 初始化结果列表
res_list <- list()

# 循环组合对象名并提取
for (gene in genes) {
  for (ptype in panss_vars) {
    for (adj in adjustments) {
      obj_name <- paste0("res_reg_", gene, "_", ptype, adj)
      if (exists(obj_name)) {
        res <- get(obj_name)
        res$Gene <- gene
        res$PANSS <- ptype
        res$Adj <- ifelse(adj == "_adj", "Adjusted", "Unadjusted")
        res_list[[length(res_list) + 1]] <- res
      } else {
        message("Object not found: ", obj_name)
      }
    }
  }
}

# 合并所有结果为一个 data.frame
res_reg_gene_panss <- do.call(rbind, res_list)


write.csv(res_reg_gene_panss, "/Users/zhuyunqing/Desktop/study/drugMR_lipid/data/onesamplemr/replicate.APOC3.pqtlGRS2panss.redu.csv")

########
#stable 25











# 定义变量名列表（二分类）
panss_vars <- c("PANSS_reduce_rate_2g", "PANSS_P_reduce_rate_2g", "PANSS_N_reduce_rate_2g", "PANSS_G_reduce_rate_2g",
                "PANSS_reduce_rate_2g_med", "PANSS_P_reduce_rate_2g_med", 
                "PANSS_N_reduce_rate_2g_med", "PANSS_G_reduce_rate_2g_med"
                )
genes <- c("apoc3_w_gnhs")  

# 保存结果列表
res_list <- list()

# 遍历每个基因和 PANSS 变量的组合
for (gene in genes) {
  for (panss in panss_vars) {
    for (adj in c(FALSE, TRUE)) {
      tryCatch({
        # 构建变量名
        grs_gene <- paste0("grs_", gene)
        formula <- if (!adj) {
          as.formula(paste0(panss, " ~ ", grs_gene))
        } else {
          as.formula(paste0(panss, " ~ ", grs_gene, " + age + age2 + as.factor(sex) + course + drugp07 + medication + PC1 + PC2 + PC3 + PC4 + PC5"))
        }
        
        # 运行 logistic 回归
        model <- glm(formula, data = m, family = binomial())
        
        # 提取 GRS 的估计系数
        coef_name <- grs_gene
        est <- coef(summary(model))[coef_name, "Estimate"]
        se <- coef(summary(model))[coef_name, "Std. Error"]
        pval <- coef(summary(model))[coef_name, "Pr(>|z|)"]
        
        # OR 和置信区间
        OR <- exp(est)
        OR_L <- exp(est - 1.96 * se)
        OR_U <- exp(est + 1.96 * se)
        
        # 存储结构化结果
        res_df <- data.frame(
          Gene = gene,
          PANSS = panss,
          Adj = ifelse(adj, "Adjusted", "Unadjusted"),
          OR = OR,
          OR_lower = OR_L,
          OR_upper = OR_U,
          p_value = pval,
          stringsAsFactors = FALSE
        )
        
        res_list[[length(res_list) + 1]] <- res_df
        cat(paste0("Processed: ", gene, " ~ ", panss, if (adj) " (adj)" else "", "\n"))
      }, error = function(e) {
        cat(paste0("Error in: ", gene, " ~ ", panss, " - ", e$message, "\n"))
      })
    }
  }
}

# 合并所有结果
res_reg_gene_panss <- do.call(rbind, res_list)


write.csv(res_reg_gene_panss,"/Users/zhuyunqing/Desktop/study/drugMR_lipid/data/onesamplemr/replicate.APOC3.pqtlGRS2panss.redu_2g.csv" )

########
#stable 26, figure 3




