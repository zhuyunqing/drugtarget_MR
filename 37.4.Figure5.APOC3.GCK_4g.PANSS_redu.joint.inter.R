########
#GLU 2 PANSS.redu

library(data.table)
library(broom)
library(dplyr)
library(ggplot2)


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

# > table(dat$grs_APOC3_GCK_4g)
# 
# 1   2   3   4 
# 499 554 555 503


dat <- dat %>%
  mutate(grs_APOC3_GCK_4g = factor(grs_APOC3_GCK_4g,
                                   levels = 1:4,
                                   labels  = paste0("Group ", 1:4)))



# 要分析的变量
vars <- c("glu1_mg", "tg1_mg")

# 分组描述性分析
desc_stats <- dat %>%
  group_by(grs_APOC3_GCK_4g) %>%
  summarise(
    across(all_of(vars), list(
      mean  = ~ mean(., na.rm = TRUE),
      median = ~ median(., na.rm = TRUE),
      Q1    = ~ quantile(., 0.25, na.rm = TRUE),
      Q3    = ~ quantile(., 0.75, na.rm = TRUE)
    ), .names = "{.col}_{.fn}")
  )

print(desc_stats)

desc_stats <- as.data.frame(desc_stats)
write.csv(desc_stats, "/Users/zhuyunqing/Desktop/study/drugMR_lipid/result/capoc.des.apoc3.gck.4g.TG.glucose.csv", row.names = FALSE)



######
#glucose

######
#regression
fml <- glu1_mg ~ as.factor(grs_APOC3_GCK_4g) + age + age2 + 
  as.factor(sex) + as.factor(centers) + PC1 + PC2 + PC3 + PC4 + PC5 +
  course + drugp07

fit <- lm(fml, data = dat)

res_glu <- tidy(fit) %>%
  filter(grepl("^as\\.factor\\(grs_APOC3_GCK_4g\\)", term)) %>%
  select(
    level   = term,
    beta    = estimate,
    se      = std.error,
    p_value = p.value
  )
res_glu$trait <- "glucose"


#######
#plot
my_cols <- c("grey", "#2b5461", "#8b8474", "#6a3d9a")

p <- ggplot(dat, aes(x = grs_APOC3_GCK_4g,
                y = glu1_mg,
                fill = grs_APOC3_GCK_4g)) +
  
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.8) +
  
  # 均值：空心菱形 (shape = 23)
  stat_summary(fun     = mean,
               geom    = "point",
               shape   = 23,
               size    = 3,
               stroke  = 0.8,
               fill    = "white",
               colour  = "black") +
  
  scale_fill_manual(values = my_cols, guide = "none") +
  labs(x = NULL,
       y = "Glucose (mg/dL)",
       title = NULL) +
  
  theme_bw(base_size = 13) +
  theme(
    plot.title      = element_text(hjust = 0.5),
    panel.grid.major.x = element_blank()
  )

pdf("/Users/zhuyunqing/Desktop/study/drugMR_lipid/result/Figure5.2.grs4g2glu.new.pdf", 
    width = 4.6, height = 4.4)

p
dev.off()




######
#tg

######
#regression
fml <- tg1_mg ~ as.factor(grs_APOC3_GCK_4g) + age + age2 + 
  as.factor(sex) + as.factor(centers) + PC1 + PC2 + PC3 + PC4 + PC5 +
  course + drugp07

fit <- lm(fml, data = dat)

res_tg <- tidy(fit) %>%
  filter(grepl("^as\\.factor\\(grs_APOC3_GCK_4g\\)", term)) %>%
  select(
    level   = term,
    beta    = estimate,
    se      = std.error,
    p_value = p.value
  )

res_tg$trait = "TG"


#######
#plot
my_cols <- c("grey", "#2b5461", "#8b8474", "#6a3d9a")

p <- ggplot(dat, aes(x = grs_APOC3_GCK_4g,
                y = tg1_mg,
                fill = grs_APOC3_GCK_4g)) +
  
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.8) +
  
  # 均值：空心菱形 (shape = 23)
  stat_summary(fun     = mean,
               geom    = "point",
               shape   = 23,
               size    = 3,
               stroke  = 0.8,
               fill    = "white",
               colour  = "black") +
  
  scale_fill_manual(values = my_cols, guide = "none") +
  labs(x = NULL,
       y = "TG (mg/dL)",
       title = NULL) +
  
  theme_bw(base_size = 13) +
  theme(
    plot.title      = element_text(hjust = 0.5),
    panel.grid.major.x = element_blank()
  ) + 
  coord_cartesian(ylim = c(0, 250))


pdf("/Users/zhuyunqing/Desktop/study/drugMR_lipid/result/Figure5.1.grs4g2tg.new.pdf", 
    width = 4.6, height = 4.4)

p
dev.off()



res <- rbind(res_glu, res_tg)
write.csv(res, "/Users/zhuyunqing/Desktop/study/drugMR_lipid/result/capoc.grs4g.apoc3.gck.glu.tg.csv", row.names = FALSE)









#########
#APOC3_TC * GCK

dat <- fread("/Users/zhuyunqing/Desktop/study/drugMR_lipid/data/onesamplemr/CAPOC.cis.capoc_cases.grs.pheno.new.txt", h=T)

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


dat <- dat %>%
  mutate(grs_APOC3_GCK_4g = factor(grs_APOC3_GCK_4g,
                                   levels = 1:4,
                                   labels  = paste0("Group ", 1:4)))



# 要分析的变量
vars <- c("glu1_mg", "tc1_mg")

# 分组描述性分析
desc_stats <- dat %>%
  group_by(grs_APOC3_GCK_4g) %>%
  summarise(
    across(all_of(vars), list(
      mean  = ~ mean(., na.rm = TRUE),
      median = ~ median(., na.rm = TRUE),
      Q1    = ~ quantile(., 0.25, na.rm = TRUE),
      Q3    = ~ quantile(., 0.75, na.rm = TRUE)
    ), .names = "{.col}_{.fn}")
  )

print(desc_stats)

desc_stats <- as.data.frame(desc_stats)
write.csv(desc_stats, "/Users/zhuyunqing/Desktop/study/drugMR_lipid/result/capoc.des.apoc3.gck.4g.TC.glucose.csv", row.names = FALSE)



######
#glucose

######
#regression
fml <- glu1_mg ~ as.factor(grs_APOC3_GCK_4g) + age + age2 + 
  as.factor(sex) + as.factor(centers) + PC1 + PC2 + PC3 + PC4 + PC5 +
  course + drugp07

fit <- lm(fml, data = dat)

res_glu <- tidy(fit) %>%
  filter(grepl("^as\\.factor\\(grs_APOC3_GCK_4g\\)", term)) %>%
  select(
    level   = term,
    beta    = estimate,
    se      = std.error,
    p_value = p.value
  )
res_glu$trait <- "glucose"





######
#tc

######
#regression
fml <- tc1_mg ~ as.factor(grs_APOC3_GCK_4g) + age + age2 + 
  as.factor(sex) + as.factor(centers) + PC1 + PC2 + PC3 + PC4 + PC5 +
  course + drugp07

fit <- lm(fml, data = dat)

res_tc <- tidy(fit) %>%
  filter(grepl("^as\\.factor\\(grs_APOC3_GCK_4g\\)", term)) %>%
  select(
    level   = term,
    beta    = estimate,
    se      = std.error,
    p_value = p.value
  )

res_tc$trait = "TC"





res <- rbind(res_glu, res_tc)
write.csv(res, "/Users/zhuyunqing/Desktop/study/drugMR_lipid/result/capoc.grs4g.apoc3_TC.gck_glu.csv", row.names = FALSE)












#######
#plot.grs4g.2.panss_redu_rate

#画图 grs_glu 2 glu

library(ggplot2)
library(dplyr)
library(data.table)
library(broom)
library(tidyr)
library(forestploter)

setwd("/Users/zhuyunqing/Desktop/study/drugMR_lipid/result")
df <- read.csv("grs.apoc3.gck2panss_redu.res.csv")

df <- df %>%
  mutate(
    Outcome = factor(Outcome, levels = c("PANSS", "PANSS_P", 
                                         "PANSS_N", "PANSS_G")),
    Group = factor(Group, levels = c("Both_lower", "APOC3_GRS_higher", "GCK_GRS_higher", "Both_higher"))
  ) %>%
  arrange(Outcome, Group)

Group_colors <- c("Both_lower" = "black", "APOC3_GRS_higher" = "#2b5461", 
                  "GCK_GRS_higher" = "#8b8474", "Both_higher" = "#6a3d9a")


#加上画图行
df$`plot` <- paste(rep(" ", 20), collapse = " ")
df$`s` <- paste(rep(" ", 3), collapse = " ")
#加上beta (95%ci)列
df$`Beta (95%CI)` <- paste(sprintf("%.2f (%.2f, %.2f)", 
                                   as.numeric(df$beta), 
                                   as.numeric(df$lower95), 
                                   as.numeric(df$upper95)))
df$Interaction <- sprintf("%.3f", df$Interaction)


tm <- forest_theme(base_size = 10,
                   ci_lwd = 1.2,
                   ci_pch = 16,
                   ci_Theight = 0.4,
                   vertline_lwd = 1,
                   vertline_col = "grey"
)




p <- forest(df[,c(1, 3, 2, 11, 9, 8)],
            est = df$beta,
            lower = df$lower95, 
            upper = df$upper95,
            xlim = c(-7.1, 4.1),
            ticks_at = c(-6 ,-3, 0, 3),
            ci_column = c(5),
            ref_line = 0,
            size=0.7,
            nudge_y = 0,
            base_size = 6,
            theme = tm
)


p


#改背景
p1 <- edit_plot(p, row = c(1:4, 9:12), which = "background",
                gp = gpar(fill="white") )
p1


p1 <- edit_plot(p1, row = c(5:8, 13:16), which = "background",
                gp = gpar(fill="#f0f0f0") )
p1



#改ci
p1 <- edit_plot(p1, row = c(1,5,9,13), col = 5, which = "ci",
                gp = gpar(col="black") )
p1

#改ci
p1 <- edit_plot(p1, row = c(2,6,10,14), col = 5, which = "ci",
                gp = gpar(col="#2b5461") )
p1

p1 <- edit_plot(p1, row = c(3,7,11,15), col = 5, which = "ci",
                gp = gpar(col="#8b8474") )
p1

p1 <- edit_plot(p1, row = c(4,8,12,16), col = 5, which = "ci",
                gp = gpar(col="#6a3d9a") )
p1

p1 <- add_border(p1, part = "header", where = "bottom")
p1


pdf("/Users/zhuyunqing/Desktop/study/drugMR_lipid/result/Figure5.3.APOC3.GCK_4g2panss_redu.new.pdf", 
    width = 8.45, height = 5.7)

p1
dev.off()



