library(data.table)
library(broom)
library(dplyr)

dat <- fread("/Users/zhuyunqing/Desktop/study/drugMR_lipid/data/onesamplemr/CAPOC.cis.capoc_cases.grs.pheno.new.txt", h=T)


#降糖药靶grs～panss减分率
#减分率连续变量

genes <- c("GCK", "ABCB11")
panss_traits <- c("PANSS_reduce_rate", "PANSS_P_reduce_rate", 
                  "PANSS_N_reduce_rate", "PANSS_G_reduce_rate")
trait <- "glu"

results <- list()

for (gene in genes) {
  grs_var <- paste0("grs_", gene, "_w_taiwan_glu")
  if (!(grs_var %in% names(dat))) next
  
  for (panss in panss_traits) {
    if (!(panss %in% names(dat))) next
    
    formula <- as.formula(
      paste(panss, "~", grs_var,
            "+ age + age2 + as.factor(sex) + as.factor(centers) +",
            "PC1 + PC2 + PC3 + PC4 + PC5 + course + drugp07 + as.factor(medication)")
    )
    
    model <- lm(formula, data = dat)
    
    tidy_model <- tidy(model, conf.int = TRUE) %>%
      filter(term == grs_var) %>%
      mutate(
        gene = gene,
        panss_trait = panss,
        trait = trait
      ) %>%
      select(gene, panss_trait, trait, estimate, conf.low, conf.high, p.value)
    
    results[[paste(gene, panss, sep = "_")]] <- tidy_model
  }
}

final_result_glu_linear <- bind_rows(results)
colnames(final_result_glu_linear) <- c("Gene", "PANSS_trait", "trait", 
                                "Beta", "CI_lower", "CI_upper", "P")

print(final_result_glu_linear)






#panss减分率2分类

genes <- c("GCK", "ABCB11")
panss_traits <- c("PANSS_reduce_rate_2g", "PANSS_P_reduce_rate_2g",
                  "PANSS_N_reduce_rate_2g", "PANSS_G_reduce_rate_2g",
                  "PANSS_reduce_rate_2g_med", "PANSS_P_reduce_rate_2g_med", 
                  "PANSS_N_reduce_rate_2g_med", "PANSS_G_reduce_rate_2g_med")
trait <- "glu"
results <- list()

for (gene in genes) {
  grs_var <- paste0("grs_", gene, "_w_taiwan_glu")
  if (!(grs_var %in% names(dat))) next
  
  for (panss in panss_traits) {
    if (!(panss %in% names(dat))) next
    
    formula <- as.formula(
      paste(panss, "~", grs_var,
            "+ age + age2 + as.factor(sex) + as.factor(centers) +",
            "PC1 + PC2 + PC3 + PC4 + PC5 + course + drugp07 + as.factor(medication)")
    )
    
    model <- glm(formula, data = dat, family = binomial)
    
    tidy_model <- tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
      filter(term == grs_var) %>%
      mutate(
        gene = gene,
        panss_trait = panss,
        trait = trait
      ) %>%
      select(gene, panss_trait, trait, estimate, conf.low, conf.high, p.value)
    
    results[[paste(gene, panss, sep = "_")]] <- tidy_model
  }
}

final_result_glu_logistic <- bind_rows(results)
colnames(final_result_glu_logistic) <- c("Gene", "PANSS_trait", "trait", 
                                         "OR", "CI_lower", "CI_upper", "P")
print(final_result_glu_logistic)







#####
#降脂药靶grs～panss减分率
#logTG


genes <- c("APOC3", "CETP", "APOB")
panss_traits <- c("PANSS_reduce_rate", "PANSS_P_reduce_rate", 
                  "PANSS_N_reduce_rate", "PANSS_G_reduce_rate")
trait <- "tg"

results <- list()

for (gene in genes) {
  grs_var <- paste0("grs_", gene, "_w_glgc_logTG")
  if (!(grs_var %in% names(dat))) next
  
  for (panss in panss_traits) {
    if (!(panss %in% names(dat))) next
    
    formula <- as.formula(
      paste(panss, "~", grs_var,
            "+ age + age2 + as.factor(sex) + as.factor(centers) +",
            "PC1 + PC2 + PC3 + PC4 + PC5 + course + drugp07 + as.factor(medication)")
    )
    
    model <- lm(formula, data = dat)
    
    tidy_model <- tidy(model, conf.int = TRUE) %>%
      filter(term == grs_var) %>%
      mutate(
        gene = gene,
        panss_trait = panss,
        trait = trait
      ) %>%
      select(gene, panss_trait, trait, estimate, conf.low, conf.high, p.value)
    
    results[[paste(gene, panss, sep = "_")]] <- tidy_model
  }
}

final_result_tg_linear <- bind_rows(results)
colnames(final_result_tg_linear) <- c("Gene", "PANSS_trait", "trait", 
                                      "Beta", "CI_lower", "CI_upper", "P")

print(final_result_tg_linear)






#panss减分率2分类

genes <- c("APOC3", "CETP", "APOB")
panss_traits <- c("PANSS_reduce_rate_2g", "PANSS_P_reduce_rate_2g",
                  "PANSS_N_reduce_rate_2g", "PANSS_G_reduce_rate_2g",
                  "PANSS_reduce_rate_2g_med", "PANSS_P_reduce_rate_2g_med", 
                  "PANSS_N_reduce_rate_2g_med", "PANSS_G_reduce_rate_2g_med")
trait <- "tg"
results <- list()

for (gene in genes) {
  grs_var <- paste0("grs_", gene, "_w_glgc_logTG")
  if (!(grs_var %in% names(dat))) next
  
  for (panss in panss_traits) {
    if (!(panss %in% names(dat))) next
    
    formula <- as.formula(
      paste(panss, "~", grs_var,
            "+ age + age2 + as.factor(sex) + as.factor(centers) +",
            "PC1 + PC2 + PC3 + PC4 + PC5 + course + drugp07 + as.factor(medication)")
    )
    
    model <- glm(formula, data = dat, family = binomial)
    
    tidy_model <- tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
      filter(term == grs_var) %>%
      mutate(
        gene = gene,
        panss_trait = panss,
        trait = trait
      ) %>%
      select(gene, panss_trait, trait, estimate, conf.low, conf.high, p.value)
    
    results[[paste(gene, panss, sep = "_")]] <- tidy_model
  }
}

final_result_tg_logistic <- bind_rows(results)
colnames(final_result_tg_logistic) <- c("Gene", "PANSS_trait", "trait", 
                                        "OR", "CI_lower", "CI_upper", "P")
print(final_result_tg_logistic)










#LDL


genes <- c("APOB", "HMGCR")
panss_traits <- c("PANSS_reduce_rate", "PANSS_P_reduce_rate", 
                  "PANSS_N_reduce_rate", "PANSS_G_reduce_rate")
trait <- "LDL"

results <- list()

for (gene in genes) {
  grs_var <- paste0("grs_", gene, "_w_glgc_ldl")
  if (!(grs_var %in% names(dat))) next
  
  for (panss in panss_traits) {
    if (!(panss %in% names(dat))) next
    
    formula <- as.formula(
      paste(panss, "~", grs_var,
            "+ age + age2 + as.factor(sex) + as.factor(centers) +",
            "PC1 + PC2 + PC3 + PC4 + PC5 + course + drugp07 + as.factor(medication)")
    )
    
    model <- lm(formula, data = dat)
    
    tidy_model <- tidy(model, conf.int = TRUE) %>%
      filter(term == grs_var) %>%
      mutate(
        gene = gene,
        panss_trait = panss,
        trait = trait
      ) %>%
      select(gene, panss_trait, trait, estimate, conf.low, conf.high, p.value)
    
    results[[paste(gene, panss, sep = "_")]] <- tidy_model
  }
}

final_result_LDL_linear <- bind_rows(results)
colnames(final_result_LDL_linear) <- c("Gene", "PANSS_trait", "trait", 
                                       "Beta", "CI_lower", "CI_upper", "P")

print(final_result_LDL_linear)






#panss减分率2分类

genes <- c("APOB", "HMGCR")
panss_traits <- c("PANSS_reduce_rate_2g", "PANSS_P_reduce_rate_2g",
                  "PANSS_N_reduce_rate_2g", "PANSS_G_reduce_rate_2g", 
                  "PANSS_reduce_rate_2g_med", "PANSS_P_reduce_rate_2g_med", 
                  "PANSS_N_reduce_rate_2g_med", "PANSS_G_reduce_rate_2g_med")
trait <- "LDL"
results <- list()

for (gene in genes) {
  grs_var <- paste0("grs_", gene, "_w_glgc_ldl")
  if (!(grs_var %in% names(dat))) next
  
  for (panss in panss_traits) {
    if (!(panss %in% names(dat))) next
    
    formula <- as.formula(
      paste(panss, "~", grs_var,
            "+ age + age2 + as.factor(sex) + as.factor(centers) +",
            "PC1 + PC2 + PC3 + PC4 + PC5 + course + drugp07 + as.factor(medication)")
    )
    
    model <- glm(formula, data = dat, family = binomial)
    
    tidy_model <- tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
      filter(term == grs_var) %>%
      mutate(
        gene = gene,
        panss_trait = panss,
        trait = trait
      ) %>%
      select(gene, panss_trait, trait, estimate, conf.low, conf.high, p.value)
    
    results[[paste(gene, panss, sep = "_")]] <- tidy_model
  }
}

final_result_LDL_logistic <- bind_rows(results)
colnames(final_result_LDL_logistic) <- c("Gene", "PANSS_trait", "trait", 
                                         "OR", "CI_lower", "CI_upper", "P")
print(final_result_LDL_logistic)










#HDL


genes <- c("ABCA1", "CETP", "LDLR", "LPL")
panss_traits <- c("PANSS_reduce_rate", "PANSS_P_reduce_rate", 
                  "PANSS_N_reduce_rate", "PANSS_G_reduce_rate")
trait <- "HDL"

results <- list()

for (gene in genes) {
  grs_var <- paste0("grs_", gene, "_w_glgc_hdl")
  if (!(grs_var %in% names(dat))) next
  
  for (panss in panss_traits) {
    if (!(panss %in% names(dat))) next
    
    formula <- as.formula(
      paste(panss, "~", grs_var,
            "+ age + age2 + as.factor(sex) + as.factor(centers) +",
            "PC1 + PC2 + PC3 + PC4 + PC5 + course + drugp07 + as.factor(medication)")
    )
    
    model <- lm(formula, data = dat)
    
    tidy_model <- tidy(model, conf.int = TRUE) %>%
      filter(term == grs_var) %>%
      mutate(
        gene = gene,
        panss_trait = panss,
        trait = trait
      ) %>%
      select(gene, panss_trait, trait, estimate, conf.low, conf.high, p.value)
    
    results[[paste(gene, panss, sep = "_")]] <- tidy_model
  }
}

final_result_HDL_linear <- bind_rows(results)
colnames(final_result_HDL_linear) <- c("Gene", "PANSS_trait", "trait", 
                                       "Beta", "CI_lower", "CI_upper", "P")

print(final_result_HDL_linear)






#panss减分率2分类

genes <- c("ABCA1", "CETP", "LDLR", "LPL")
panss_traits <- c("PANSS_reduce_rate_2g", "PANSS_P_reduce_rate_2g",
                  "PANSS_N_reduce_rate_2g", "PANSS_G_reduce_rate_2g", 
                  "PANSS_reduce_rate_2g_med", "PANSS_P_reduce_rate_2g_med", 
                  "PANSS_N_reduce_rate_2g_med", "PANSS_G_reduce_rate_2g_med")
trait <- "HDL"
results <- list()

for (gene in genes) {
  grs_var <- paste0("grs_", gene, "_w_glgc_hdl")
  if (!(grs_var %in% names(dat))) next
  
  for (panss in panss_traits) {
    if (!(panss %in% names(dat))) next
    
    formula <- as.formula(
      paste(panss, "~", grs_var,
            "+ age + age2 + as.factor(sex) + as.factor(centers) +",
            "PC1 + PC2 + PC3 + PC4 + PC5 + course + drugp07 + as.factor(medication)")
    )
    
    model <- glm(formula, data = dat, family = binomial)
    
    tidy_model <- tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
      filter(term == grs_var) %>%
      mutate(
        gene = gene,
        panss_trait = panss,
        trait = trait
      ) %>%
      select(gene, panss_trait, trait, estimate, conf.low, conf.high, p.value)
    
    results[[paste(gene, panss, sep = "_")]] <- tidy_model
  }
}

final_result_HDL_logistic <- bind_rows(results)
colnames(final_result_HDL_logistic) <- c("Gene", "PANSS_trait", "trait", 
                                         "OR", "CI_lower", "CI_upper", "P")
print(final_result_HDL_logistic)







#TC


genes <- c("APOB", "APOC3", "LDLR")
panss_traits <- c("PANSS_reduce_rate", "PANSS_P_reduce_rate", 
                  "PANSS_N_reduce_rate", "PANSS_G_reduce_rate")
trait <- "TC"

results <- list()

for (gene in genes) {
  grs_var <- paste0("grs_", gene, "_w_glgc_TC")
  if (!(grs_var %in% names(dat))) next
  
  for (panss in panss_traits) {
    if (!(panss %in% names(dat))) next
    
    formula <- as.formula(
      paste(panss, "~", grs_var,
            "+ age + age2 + as.factor(sex) + as.factor(centers) +",
            "PC1 + PC2 + PC3 + PC4 + PC5 + course + drugp07 + as.factor(medication)")
    )
    
    model <- lm(formula, data = dat)
    
    tidy_model <- tidy(model, conf.int = TRUE) %>%
      filter(term == grs_var) %>%
      mutate(
        gene = gene,
        panss_trait = panss,
        trait = trait
      ) %>%
      select(gene, panss_trait, trait, estimate, conf.low, conf.high, p.value)
    
    results[[paste(gene, panss, sep = "_")]] <- tidy_model
  }
}

final_result_TC_linear <- bind_rows(results)
colnames(final_result_TC_linear) <- c("Gene", "PANSS_trait", "trait", 
                                      "Beta", "CI_lower", "CI_upper", "P")

print(final_result_TC_linear)






#panss减分率2分类

genes <- c("APOB", "APOC3", "LDLR")
panss_traits <- c("PANSS_reduce_rate_2g", "PANSS_P_reduce_rate_2g",
                  "PANSS_N_reduce_rate_2g", "PANSS_G_reduce_rate_2g", 
                  "PANSS_reduce_rate_2g_med", "PANSS_P_reduce_rate_2g_med", 
                  "PANSS_N_reduce_rate_2g_med", "PANSS_G_reduce_rate_2g_med")
trait <- "TC"
results <- list()

for (gene in genes) {
  grs_var <- paste0("grs_", gene, "_w_glgc_TC")
  if (!(grs_var %in% names(dat))) next
  
  for (panss in panss_traits) {
    if (!(panss %in% names(dat))) next
    
    formula <- as.formula(
      paste(panss, "~", grs_var,
            "+ age + age2 + as.factor(sex) + as.factor(centers) +",
            "PC1 + PC2 + PC3 + PC4 + PC5 + course + drugp07 + as.factor(medication)")
    )
    
    model <- glm(formula, data = dat, family = binomial)
    
    tidy_model <- tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
      filter(term == grs_var) %>%
      mutate(
        gene = gene,
        panss_trait = panss,
        trait = trait
      ) %>%
      select(gene, panss_trait, trait, estimate, conf.low, conf.high, p.value)
    
    results[[paste(gene, panss, sep = "_")]] <- tidy_model
  }
}

final_result_TC_logistic <- bind_rows(results)
colnames(final_result_TC_logistic) <- c("Gene", "PANSS_trait", "trait", 
                                        "OR", "CI_lower", "CI_upper", "P")
print(final_result_TC_logistic)




##########
#combined results
res_df_linear <- rbind(final_result_glu_linear, final_result_tg_linear, final_result_LDL_linear, final_result_HDL_linear, final_result_TC_linear)
res_df_logistic <- rbind(final_result_glu_logistic, final_result_tg_logistic, final_result_LDL_logistic, final_result_HDL_logistic, final_result_TC_logistic)

write.csv(res_df_linear, "/Users/zhuyunqing/Desktop/study/drugMR_lipid/result/res_grs_lg_2panss_redu.linear.new.3SD.csv", row.names = FALSE)
write.csv(res_df_logistic, "/Users/zhuyunqing/Desktop/study/drugMR_lipid/result/res_grs_lg_2panss_redu.logistic.new.3SD.csv", row.names = FALSE)


########
#figure 3, stable 13-14



