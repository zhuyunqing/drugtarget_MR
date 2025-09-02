#画图 grs_glu 2 glu

library(ggplot2)
library(dplyr)
library(data.table)
library(broom)
library(tidyr)
library(forestploter)

dat <- fread("/Users/zhuyunqing/Desktop/study/drugMR_lipid/result/res_grs_lg_2panss_redu.linear.new.3SD.csv", h=T)
dat <- subset(dat, Gene=="GCK" | (Gene=="APOC3" & trait == "tg"))
dat$PANSS_trait <- c("PANSS", "PANSS_P", "PANSS_N", "PANSS_G",
                     "PANSS", "PANSS_P", "PANSS_N", "PANSS_G")
dat <- dat[,-3]

dat1 <- fread("/Users/zhuyunqing/Desktop/study/drugMR_lipid/data/onesamplemr/capoc.osmr.res_grs.pred.lg2panss_redu.linear.csv")
dat1 <- dat1[-c(9:12),]
dat1$Gene = c("Glucose", "Glucose", "Glucose", "Glucose", 
              "TG", "TG", "TG", "TG")
dat1$PANSS_trait <- c("PANSS", "PANSS_P", "PANSS_N", "PANSS_G",
                      "PANSS", "PANSS_P", "PANSS_N", "PANSS_G")
dat1 <- dat1 %>%
  separate(CI, into = c("CI_lower", "CI_upper"), sep = "\\s*~\\s*", convert = TRUE)
dat1 <- dat1[,c(7,8,3,4,5,6)]


dat <- rbind(dat, dat1)



dat <- dat %>%
  mutate(
    PANSS_trait = factor(PANSS_trait, levels = c("PANSS", "PANSS_P", 
                                                         "PANSS_N", "PANSS_G")),
    Gene = factor(Gene, levels = c("TG","APOC3", "Glucose","GCK"))
    
  ) %>%
  arrange(PANSS_trait, Gene)


trait_colors <- c("TG" = "#2b5461",  "APOC3" = "#2b5461", "Glucose" = "#8b8474", "GCK" = "#8b8474")

#加上画图行
dat$`plot` <- paste(rep(" ", 20), collapse = " ")
dat$`s` <- paste(rep(" ", 3), collapse = " ")
#加上beta (95%ci)列
dat$`Beta (95%CI)` <- paste(sprintf("%.2f (%.2f, %.2f)", 
                                    as.numeric(dat$Beta), 
                                    as.numeric(dat$CI_lower), 
                                    as.numeric(dat$CI_upper)))
dat$P <- sprintf("%.3f", dat$P)






#####
#lipid glu分开
dat_APOC3 <- subset(dat, Gene == "TG" | Gene == "APOC3")
dat_GCK <- subset(dat, Gene == "Glucose" | Gene == "GCK")



tm <- forest_theme(base_size = 10,
                   ci_lwd = 1.2,
                   ci_pch = 16,
                   ci_Theight = 0.4,
                   vertline_lwd = 1,
                   vertline_col = "grey"
)


p <- forest(dat_APOC3[,c(2, 8, 1,8, 9, 8, 7, 8, 6)],
            est = dat_APOC3$Beta,
            lower = dat_APOC3$CI_lower, 
            upper = dat_APOC3$CI_upper,
            xlim = c(-1, 2.2),
            ticks_at = c(-1,  0, 1, 2),
            ci_column = c(7),
            ref_line = 0,
            size=0.6,
            nudge_y = 0,
            base_size = 6,
            theme = tm
)


p


#改背景
p1 <- edit_plot(p, row = c(1:2), which = "background",
                gp = gpar(fill="white") )
p1


p1 <- edit_plot(p1, row = c(3:4), which = "background",
                gp = gpar(fill="#f0f0f0") )
p1

p1 <- edit_plot(p1, row = c(5:6), which = "background",
                gp = gpar(fill="white") )
p1

p1 <- edit_plot(p1, row = c(7:8), which = "background",
                gp = gpar(fill="#f0f0f0") )
p1




#改ci
p1 <- edit_plot(p1, row = c(1:8), col = 7, which = "ci",
                gp = gpar(col="#2b5461") )
p1

p1 <- add_border(p1, part = "header", where = "bottom")
p1



pdf("/Users/zhuyunqing/Desktop/study/drugMR_lipid/result/Figure3.1.grs_grs.apoc32panss_redu.new.pdf", 
     width = 8.45, height = 5.7)

p1
dev.off()




p <- forest(dat_GCK[,c(2, 8, 1,8, 9, 8, 7, 8, 6)],
            est = dat_GCK$Beta,
            lower = dat_GCK$CI_lower, 
            upper = dat_GCK$CI_upper,
            xlim = c(-3, 2),
            ticks_at = c(-3,  -2, -1, 0, 1, 2),
            ci_column = c(7),
            ref_line = 0,
            size=0.6,
            nudge_y = 0,
            base_size = 6,
            theme = tm
)


p


#改背景
p1 <- edit_plot(p, row = c(1:2), which = "background",
                gp = gpar(fill="white") )
p1


p1 <- edit_plot(p1, row = c(3:4), which = "background",
                gp = gpar(fill="#f0f0f0") )
p1

p1 <- edit_plot(p1, row = c(5:6), which = "background",
                gp = gpar(fill="white") )
p1

p1 <- edit_plot(p1, row = c(7:8), which = "background",
                gp = gpar(fill="#f0f0f0") )
p1




#改ci
p1 <- edit_plot(p1, row = c(1:8), col = 7, which = "ci",
                gp = gpar(col="#8b8474") )
p1

p1 <- add_border(p1, part = "header", where = "bottom")
p1



pdf("/Users/zhuyunqing/Desktop/study/drugMR_lipid/result/Figure3.2.grs_grs.GCK32panss_redu.new.pdf", 
    width = 8.45, height = 5.7)

p1
dev.off()





