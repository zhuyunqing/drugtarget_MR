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


m <- fread("/Users/zhuyunqing/Desktop/study/drugMR_lipid/data/CAPEC/CAPEC.grs.pheno.new.txt")

################
# 计算中位数
median_GCK <- median(m$grs_GCK_w_taiwan_glu)
median_APOC3 <- median(m$grs_APOC3_w_glgc_logTG)
m$GCK_GRS_higher <- ifelse(m$grs_GCK_w_taiwan_glu >= median_GCK, 1, 0)
m$APOC3_GRS_higher <- ifelse(m$grs_APOC3_w_glgc_logTG >= median_APOC3, 1, 0)


# 计算所有研究对象中的男性人数及其占比
results <- m %>%
  summarise(
    male_count = sum(Sex == 1),
    total_count = n(),
    male_ratio = male_count / total_count * 100,
    male_count_ratio = paste0(male_count, " (", sprintf("%0.1f", male_ratio), ")")
  )

print(results)

male_count_ratio <- as.data.table(t(as.data.table(results$male_count_ratio)))
male_count_ratio$variable <- "male_count_ratio"




# 计算亚组中的女性人数及其占比
results <- m %>%
  summarise(
    female_count = sum(Sex == 2),
    total_count = n(),
    female_ratio = female_count / total_count * 100,
    female_count_ratio = paste0(female_count, " (", sprintf("%0.1f", female_ratio),")")
  )

# 查看结果
print(results)
female_count_ratio <- as.data.table(t(as.data.table(results$female_count_ratio)))
female_count_ratio$variable <- "female_count_ratio"



# 计算中位数和四分位数（P25, P75）
results <- m %>%
  summarise(
    median_age = median(age, na.rm = TRUE),
    q25 = quantile(age, 0.25, na.rm = TRUE),
    q75 = quantile(age, 0.75, na.rm = TRUE),
    age_median_iqr = paste0(
      sprintf("%0.1f", median_age), " (",
      sprintf("%0.1f", q25), ", ",
      sprintf("%0.1f", q75), ")"
    )
  )

# 查看结果
print(results)

# 转换为一行格式
age_median_iqr <- as.data.table(t(as.data.table(results$age_median_iqr)))
age_median_iqr$variable <- "age_median_iqr"







# 定义变量列表
variables <- c("center", "FirstEpi") 

results_cate_list <- list()

for (var in variables) {
  result <- m %>%
    summarise(
      count = sum(get(var) == 1), 
      total_count = n(),
      ratio = count / total_count * 100,
      count_ratio = paste0(count, " (", sprintf("%0.1f", ratio), ")")
    )
  
  result_t <- as.data.table(t(as.data.table(result$count_ratio)))
  result_t[, variable := var]  # 添加变量名列
  
  results_cate_list[[var]] <- result_t
}

final_results_cate <- rbindlist(results_cate_list, use.names = TRUE, fill = TRUE)
print(final_results_cate)








# 连续变量列表
continuous_variables <- c("v1_tc_mg", "v1_tg_mg", "v1_hdl_mg", "v1_ldl_mg", "v1_glu_mg")

# 初始化结果列表
results_continuous_list <- list()

for (var in continuous_variables) {
  # 计算中位数和四分位数
  result <- m %>%
    summarise(
      median_value = median(get(var), na.rm = TRUE),
      q25 = quantile(get(var), 0.25, na.rm = TRUE),
      q75 = quantile(get(var), 0.75, na.rm = TRUE),
      median_iqr = paste0(
        sprintf("%0.1f", median_value), " (",
        sprintf("%0.1f", q25), ", ",
        sprintf("%0.1f", q75), ")"
      )
    )
  
  # 转置结果并添加变量名
  result_t <- as.data.table(t(as.data.table(result$median_iqr)))
  result_t[, variable := var]
  
  # 保存到列表
  results_continuous_list[[var]] <- result_t
}

# 合并为最终表格
final_results_continuous2 <- rbindlist(results_continuous_list, use.names = TRUE, fill = TRUE)

# 查看结果
print(final_results_continuous2)



#combine为table1
table1_all <- rbind(male_count_ratio, female_count_ratio, age_median_iqr, final_results_cate, final_results_continuous2)
colnames(table1_all) <- c("all", "var")








#######################
# 计算各类别人数
count <- table(m$APOC3_GRS_higher)
print(count)

# 计算各类别占比
proportion <- prop.table(count)
print(proportion)

# 合并结果为数据框
result <- data.frame(Category = names(count), Count = as.vector(count), Proportion = as.vector(proportion * 100))
print(result)
N_prop <- as.data.table(t(as.data.table(paste0(result$Count, 
                                               " (", sprintf("%0.1f", result$Proportion), ")"))))

N_prop





######
#性别
# 计算每个打鼾亚组中的男性人数及其占比
results <- m %>%
  group_by(APOC3_GRS_higher) %>%
  summarise(
    male_count = sum(Sex == 1),
    total_count = n(),
    male_ratio = male_count / total_count * 100,
    male_count_ratio = paste0(male_count, " (", sprintf("%0.1f", male_ratio),")")
  )

# 查看结果
print(results)
male_count_ratio <- as.data.table(t(as.data.table(results$male_count_ratio)))
male_count_ratio$variable <- "male_count_ratio"




# 计算每个打鼾亚组中的女性人数及其占比
results <- m %>%
  group_by(APOC3_GRS_higher) %>%
  summarise(
    female_count = sum(Sex == 2),
    total_count = n(),
    female_ratio = female_count / total_count * 100,
    female_count_ratio = paste0(female_count, " (", sprintf("%0.1f", female_ratio),")")
  )

# 查看结果
print(results)
female_count_ratio <- as.data.table(t(as.data.table(results$female_count_ratio)))
female_count_ratio$variable <- "female_count_ratio"






# 计算中位数和四分位数（P25, P75）
results <- m %>%
  group_by(APOC3_GRS_higher) %>%
  summarise(
    median_age = median(age, na.rm = TRUE),
    q25 = quantile(age, 0.25, na.rm = TRUE),
    q75 = quantile(age, 0.75, na.rm = TRUE),
    age_median_iqr = paste0(
      sprintf("%0.1f", median_age), " (",
      sprintf("%0.1f", q25), ", ",
      sprintf("%0.1f", q75), ")"
    )
  )

# 查看结果
print(results)

# 转换为一行格式
age_median_iqr <- as.data.table(t(as.data.table(results$age_median_iqr)))
age_median_iqr$variable <- "age_median_iqr"





# 定义变量列表
variables <- c("center", "FirstEpi") 

results_cate_list <- list()

for (var in variables) {
  result <- m %>%
    group_by(APOC3_GRS_higher) %>%
    summarise(
      count = sum(get(var) == 1), 
      total_count = n(),
      ratio = count / total_count * 100,
      count_ratio = paste0(count, " (", sprintf("%0.1f", ratio), ")")
    )
  
  result_t <- as.data.table(t(as.data.table(result$count_ratio)))
  result_t[, variable := var]  # 添加变量名列
  
  results_cate_list[[var]] <- result_t
}

final_results_cate <- rbindlist(results_cate_list, use.names = TRUE, fill = TRUE)
print(final_results_cate)







# 连续变量列表
continuous_variables <- c("v1_tc_mg", "v1_tg_mg", "v1_hdl_mg", "v1_ldl_mg", "v1_glu_mg")

# 初始化结果列表
results_continuous_list <- list()

for (var in continuous_variables) {
  # 计算中位数和四分位数
  result <- m %>%
    group_by(APOC3_GRS_higher) %>%
    summarise(
      median_value = median(get(var), na.rm = TRUE),
      q25 = quantile(get(var), 0.25, na.rm = TRUE),
      q75 = quantile(get(var), 0.75, na.rm = TRUE),
      median_iqr = paste0(
        sprintf("%0.1f", median_value), " (",
        sprintf("%0.1f", q25), ", ",
        sprintf("%0.1f", q75), ")"
      )
    )
  
  # 转置结果并添加变量名
  result_t <- as.data.table(t(as.data.table(result$median_iqr)))
  result_t[, variable := var]
  
  # 保存到列表
  results_continuous_list[[var]] <- result_t
}

# 合并为最终表格
final_results_continuous2 <- rbindlist(results_continuous_list, use.names = TRUE, fill = TRUE)

# 查看结果
print(final_results_continuous2)




#combine为table1
table1_APOC3 <- rbind(male_count_ratio, female_count_ratio, age_median_iqr, final_results_cate, final_results_continuous2)
colnames(table1_APOC3) <- c("APOC3_GRS_lower", "APOC3_GRS_higher", "var")












#######################
# 计算各类别人数
count <- table(m$GCK_GRS_higher)
print(count)

# 计算各类别占比
proportion <- prop.table(count)
print(proportion)

# 合并结果为数据框
result <- data.frame(Category = names(count), Count = as.vector(count), Proportion = as.vector(proportion * 100))
print(result)
N_prop <- as.data.table(t(as.data.table(paste0(result$Count, 
                                               " (", sprintf("%0.1f", result$Proportion), ")"))))

N_prop





######
#性别
# 计算每个打鼾亚组中的男性人数及其占比
results <- m %>%
  group_by(GCK_GRS_higher) %>%
  summarise(
    male_count = sum(Sex == 1),
    total_count = n(),
    male_ratio = male_count / total_count * 100,
    male_count_ratio = paste0(male_count, " (", sprintf("%0.1f", male_ratio),")")
  )

# 查看结果
print(results)
male_count_ratio <- as.data.table(t(as.data.table(results$male_count_ratio)))
male_count_ratio$variable <- "male_count_ratio"




# 计算每个打鼾亚组中的女性人数及其占比
results <- m %>%
  group_by(GCK_GRS_higher) %>%
  summarise(
    female_count = sum(Sex == 2),
    total_count = n(),
    female_ratio = female_count / total_count * 100,
    female_count_ratio = paste0(female_count, " (", sprintf("%0.1f", female_ratio),")")
  )

# 查看结果
print(results)
female_count_ratio <- as.data.table(t(as.data.table(results$female_count_ratio)))
female_count_ratio$variable <- "female_count_ratio"






# 计算中位数和四分位数（P25, P75）
results <- m %>%
  group_by(GCK_GRS_higher) %>%
  summarise(
    median_age = median(age, na.rm = TRUE),
    q25 = quantile(age, 0.25, na.rm = TRUE),
    q75 = quantile(age, 0.75, na.rm = TRUE),
    age_median_iqr = paste0(
      sprintf("%0.1f", median_age), " (",
      sprintf("%0.1f", q25), ", ",
      sprintf("%0.1f", q75), ")"
    )
  )

# 查看结果
print(results)

# 转换为一行格式
age_median_iqr <- as.data.table(t(as.data.table(results$age_median_iqr)))
age_median_iqr$variable <- "age_median_iqr"





# 定义变量列表
variables <- c("center", "FirstEpi") 

results_cate_list <- list()

for (var in variables) {
  result <- m %>%
    group_by(GCK_GRS_higher) %>%
    summarise(
      count = sum(get(var) == 1), 
      total_count = n(),
      ratio = count / total_count * 100,
      count_ratio = paste0(count, " (", sprintf("%0.1f", ratio), ")")
    )
  
  result_t <- as.data.table(t(as.data.table(result$count_ratio)))
  result_t[, variable := var]  # 添加变量名列
  
  results_cate_list[[var]] <- result_t
}

final_results_cate <- rbindlist(results_cate_list, use.names = TRUE, fill = TRUE)
print(final_results_cate)







# 连续变量列表
continuous_variables <- c("v1_tc_mg", "v1_tg_mg", "v1_hdl_mg", "v1_ldl_mg", "v1_glu_mg")

# 初始化结果列表
results_continuous_list <- list()

for (var in continuous_variables) {
  # 计算中位数和四分位数
  result <- m %>%
    group_by(GCK_GRS_higher) %>%
    summarise(
      median_value = median(get(var), na.rm = TRUE),
      q25 = quantile(get(var), 0.25, na.rm = TRUE),
      q75 = quantile(get(var), 0.75, na.rm = TRUE),
      median_iqr = paste0(
        sprintf("%0.1f", median_value), " (",
        sprintf("%0.1f", q25), ", ",
        sprintf("%0.1f", q75), ")"
      )
    )
  
  # 转置结果并添加变量名
  result_t <- as.data.table(t(as.data.table(result$median_iqr)))
  result_t[, variable := var]
  
  # 保存到列表
  results_continuous_list[[var]] <- result_t
}

# 合并为最终表格
final_results_continuous2 <- rbindlist(results_continuous_list, use.names = TRUE, fill = TRUE)

# 查看结果
print(final_results_continuous2)




#combine为table1
table1_GCK <- rbind(male_count_ratio, female_count_ratio, age_median_iqr, final_results_cate, final_results_continuous2)
colnames(table1_GCK) <- c("GCK_GRS_lower", "GCK_GRS_higher", "var")





tables <- left_join(table1_all, table1_APOC3, by = "var")
tables <- left_join(tables, table1_GCK, by = "var")



############

model <- lm(m$grs_APOC3_w_glgc_logTG ~ m$Sex + m$age + m$center + m$PC1 + m$PC2 + m$PC3 + m$PC4 + m$PC5)
summary(model)
#仅地区显著

# 定义要替换的变量
variables <- c("FirstEpi", 
               "v1_tc_mg", "v1_tg_mg", "v1_hdl_mg", "v1_ldl_mg", "v1_glu_mg"
)
p_values_df <- data.frame(variable = character(), p_value = numeric(), stringsAsFactors = FALSE)

for (var in variables) {
  formula <- as.formula(paste("m$", var, " ~ m$grs_APOC3_w_glgc_logTG + m$Sex + m$age + m$age2 + m$center + m$PC1 + m$PC2 + m$PC3 + m$PC4 + m$PC5", sep = ""))
  model <- lm(formula, data = m)
  p_value <- summary(model)$coefficients[2, "Pr(>|t|)"]
  p_values_df <- rbind(p_values_df, data.frame(variable = var, p_value = p_value))
}


p_values_df$p_value <- ifelse(p_values_df$p_value < 0.001, "<0.001", sprintf("%0.3e", p_values_df$p_value))
p_values_APOC3 <- p_values_df
colnames(p_values_APOC3) <- c("var", "APOC3_p")




############

model <- lm(m$grs_GCK_w_taiwan_glu ~ m$Sex + m$age + m$center + m$PC1 + m$PC2 + m$PC3 + m$PC4 + m$PC5)
summary(model)
#仅地区显著

# 定义要替换的变量
variables <- c("FirstEpi", 
               "v1_tc_mg", "v1_tg_mg", "v1_hdl_mg", "v1_ldl_mg", "v1_glu_mg"
)
p_values_df <- data.frame(variable = character(), p_value = numeric(), stringsAsFactors = FALSE)

for (var in variables) {
  formula <- as.formula(paste("m$", var, " ~ m$grs_GCK_w_taiwan_glu + m$Sex + m$age + m$age2 + m$center + m$PC1 + m$PC2 + m$PC3 + m$PC4 + m$PC5", sep = ""))
  model <- lm(formula, data = m)
  p_value <- summary(model)$coefficients[2, "Pr(>|t|)"]
  p_values_df <- rbind(p_values_df, data.frame(variable = var, p_value = p_value))
}


p_values_df$p_value <- ifelse(p_values_df$p_value < 0.001, "<0.001", sprintf("%0.3e", p_values_df$p_value))
p_values_GCK <- p_values_df
colnames(p_values_GCK) <- c("var", "GCK_p")

p_values <- left_join(p_values_APOC3, p_values_GCK, by = "var")

fwrite(tables, "/Users/zhuyunqing/Desktop/study/drugMR_lipid/result/CAPEC.grs_4g.des.csv")
fwrite(p_values, "/Users/zhuyunqing/Desktop/study/drugMR_lipid/result/CAPEC.grs_4g.des.p.csv")


########
#stable 29


