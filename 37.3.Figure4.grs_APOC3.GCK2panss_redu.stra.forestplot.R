#画图 grs_glu 2 glu

library(ggplot2)
library(dplyr)
library(data.table)
library(broom)
library(tidyr)
library(forestploter)

setwd("/Users/zhuyunqing/Desktop/study/drugMR_lipid/data/onesamplemr")
dat <- read.csv("res_grs.APOC3.GCK_2panss_redu.linear.stra_medication_metabo_2g.csv")
dat <- dat[-c(9:16),]
dat$Outcome <- c("PANSS", "PANSS", "PANSS_P", "PANSS_P", 
                 "PANSS_N", "PANSS_N", "PANSS_G", "PANSS_G",
                 "PANSS", "PANSS", "PANSS_P", "PANSS_P", 
                 "PANSS_N", "PANSS_N", "PANSS_G", "PANSS_G")


dat <- dat %>%
  mutate(
    Outcome = factor(Outcome, levels = c("PANSS", "PANSS_P", 
                                         "PANSS_N", "PANSS_G")),
    gene = factor(gene, levels = c("APOC3", "GCK")),
    med_group = factor(med_group, levels = c(1,2))
  ) %>%
  arrange(med_group, gene, Outcome)

dat <- dat[,-1]

colnames(dat) <- c("Gene", "Medication", "Beta", "CI_lower", "CI_upper", "P", "Outcome")


blank_row <- dat[1, ]  # 复制一个结构一致的空行
blank_row[,] <- NA     # 全部设为 NA
dat <- bind_rows(
  dat[1:8, ],
  blank_row,
  dat[9:nrow(dat), ]
)




blank_row_APOC3 <- data.frame(
  Gene = "APOC3", Medication = NA,
  Beta = NA, CI_lower = NA, CI_upper = NA, P = NA, Outcome = NA
)

blank_row_GCK <- data.frame(
  Gene = "GCK", Medication = NA,
  Beta = NA, CI_lower = NA, CI_upper = NA, P = NA, Outcome = NA
)

dat_new <- bind_rows(
  blank_row_APOC3,
  dat[1:4, ],
  blank_row_GCK,
  dat[5:8, ],
  dat[9, ],  # 原空行
  blank_row_APOC3,
  dat[10:13, ],
  blank_row_GCK,
  dat[14:17, ]
)

# 重置行号
rownames(dat_new) <- NULL
dat = dat_new





##
#加上画图行
dat$`plot` <- paste(rep(" ", 18), collapse = " ")
dat$`s` <- paste(rep(" ", 0.001), collapse = " ")
#加上beta (95%ci)列
dat$`Beta (95%CI)` <- paste(sprintf("%.2f (%.2f, %.2f)", 
                                    as.numeric(dat$Beta), 
                                    as.numeric(dat$CI_lower), 
                                    as.numeric(dat$CI_upper)))
dat$P <- sprintf("%.3f", dat$P)



tm <- forest_theme(base_size = 10,
                   ci_lwd = 1.2,
                   ci_pch = 16,
                   ci_Theight = 0.4,
                   vertline_lwd = 1,
                   vertline_col = "grey"
)



p <- forest(dat[,c(1, 7, 10, 8, 6)],
            est = dat$Beta,
            lower = dat$CI_lower, 
            upper = dat$CI_upper,
            xlim = c(-3.5, 3),
            ticks_at = c( -2,  0, 2),
            ci_column = c(4),
            ref_line = 0,
            size=0.6,
            nudge_y = 0,
            base_size = 10,
            theme = tm
)

p



#background
rows_list <- list(1:nrow(dat))
fills <- c("white", "white", "white", "white", "white")

p1 <- p
for (i in seq_along(rows_list)) {
  p1 <- edit_plot(p1, row = rows_list[[i]], which = "background",
                  gp = gpar(fill = fills[i]))
}


#ci
rows_list <- list(c(2:5, 13:16), c(7:10, 18:21))
fills <- c("#2b5461", "#8b8474")

p1 <- p1  # 如果前面 edit_plot 过，继续用
for (i in seq_along(rows_list)) {
  p1 <- edit_plot(p1, row = rows_list[[i]], col = 4, which = "ci",
                  gp = gpar(col = fills[i]))
}

p1

p1 <- add_border(p1, part = "header", where = "bottom")
p1


pdf("/Users/zhuyunqing/Desktop/study/drugMR_lipid/result/res_grs_APOC3_GCK_2panss_redu.linear.stra_medication_2g.forest.pdf", 
    width = 6, height = 6)

p1
dev.off()








dat <- read.csv("res_grs.APOC3.GCK_2panss_redu.linear.stra_firstepi_nodrug_2g.csv")
dat <- dat[-c(9:16),]
dat$Outcome <- c("PANSS", "PANSS", "PANSS_P", "PANSS_P", 
                 "PANSS_N", "PANSS_N", "PANSS_G", "PANSS_G",
                 "PANSS", "PANSS", "PANSS_P", "PANSS_P", 
                 "PANSS_N", "PANSS_N", "PANSS_G", "PANSS_G")


dat <- dat %>%
  mutate(
    Outcome = factor(Outcome, levels = c("PANSS", "PANSS_P", 
                                         "PANSS_N", "PANSS_G")),
    gene = factor(gene, levels = c("APOC3", "GCK")),
    med_group = factor(med_group, levels = c(0,1))
  ) %>%
  arrange(med_group, gene, Outcome)

dat <- dat[,-1]

colnames(dat) <- c("Gene", "Medication", "Beta", "CI_lower", "CI_upper", "P", "Outcome")


blank_row <- dat[1, ]  # 复制一个结构一致的空行
blank_row[,] <- NA     # 全部设为 NA
dat <- bind_rows(
  dat[1:8, ],
  blank_row,
  dat[9:nrow(dat), ]
)




blank_row_APOC3 <- data.frame(
  Gene = "APOC3", Medication = NA,
  Beta = NA, CI_lower = NA, CI_upper = NA, P = NA, Outcome = NA
)

blank_row_GCK <- data.frame(
  Gene = "GCK", Medication = NA,
  Beta = NA, CI_lower = NA, CI_upper = NA, P = NA, Outcome = NA
)

dat_new <- bind_rows(
  blank_row_APOC3,
  dat[1:4, ],
  blank_row_GCK,
  dat[5:8, ],
  dat[9, ],  # 原空行
  blank_row_APOC3,
  dat[10:13, ],
  blank_row_GCK,
  dat[14:17, ]
)

# 重置行号
rownames(dat_new) <- NULL
dat = dat_new





##
#加上画图行
dat$`plot` <- paste(rep(" ", 18), collapse = " ")
dat$`s` <- paste(rep(" ", 0.001), collapse = " ")
#加上beta (95%ci)列
dat$`Beta (95%CI)` <- paste(sprintf("%.2f (%.2f, %.2f)", 
                                    as.numeric(dat$Beta), 
                                    as.numeric(dat$CI_lower), 
                                    as.numeric(dat$CI_upper)))
dat$P <- sprintf("%.3f", dat$P)



tm <- forest_theme(base_size = 10,
                   ci_lwd = 1.2,
                   ci_pch = 16,
                   ci_Theight = 0.4,
                   vertline_lwd = 1,
                   vertline_col = "grey"
)



p <- forest(dat[,c(1, 7, 10, 8, 6)],
            est = dat$Beta,
            lower = dat$CI_lower, 
            upper = dat$CI_upper,
            xlim = c(-5, 3),
            ticks_at = c( -4, -2, 0, 2),
            ci_column = c(4),
            ref_line = 0,
            size=0.6,
            nudge_y = 0,
            base_size = 10,
            theme = tm
)

p



#background
rows_list <- list(1:nrow(dat))
fills <- c("white", "white", "white", "white", "white")

p1 <- p
for (i in seq_along(rows_list)) {
  p1 <- edit_plot(p1, row = rows_list[[i]], which = "background",
                  gp = gpar(fill = fills[i]))
}


#ci
rows_list <- list(c(2:5, 13:16), c(7:10, 18:21))
fills <- c("#2b5461", "#8b8474")

p1 <- p1  # 如果前面 edit_plot 过，继续用
for (i in seq_along(rows_list)) {
  p1 <- edit_plot(p1, row = rows_list[[i]], col = 4, which = "ci",
                  gp = gpar(col = fills[i]))
}

p1

p1 <- add_border(p1, part = "header", where = "bottom")
p1


pdf("/Users/zhuyunqing/Desktop/study/drugMR_lipid/result/res_grs_APOC3_GCK_2panss_redu.linear.stra_firstepi_nodrug_2g.forest.pdf", 
    width = 6, height = 6)

p1
dev.off()


