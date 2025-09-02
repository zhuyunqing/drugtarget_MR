######
#capec验证
library(data.table)

setwd("/Users/zhuyunqing/Desktop/study/drugMR_lipid/data/CAPEC")


####baseline
capec.bl <- readxl::read_xlsx("cov_863_568.xlsx")
capec.bl$age <- as.numeric(capec.bl$age)
capec.bl$age2 <- capec.bl$age * capec.bl$age


######pcs
capec.pc <- fread("CAPEC.imputed.eigenvec")
capec.pc <- capec.pc[,-1]

######
##grs of apoc3 gck

dat <- fread("/Users/zhuyunqing/Desktop/study/drugMR_lipid/data/onesamplemr/CAPEC.cis.rep.APOC3.GCK.txt")

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

dat <- fill_missing_with_mode(dat)
str(dat)



commands <- c("dat$rs4938309_T = 2 - dat$rs4938309_C", 
              "dat$rs595049_C = 2 - dat$rs595049_A", 
              "dat$rs11216190_C = 2 - dat$rs11216190_T", 
              "dat$rs17519093_A = 2 - dat$rs17519093_G", 
              "dat$rs888245_G = 2 - dat$rs888245_C", 
              "dat$rs518181_A = 2 - dat$rs518181_C", 
              "dat$rs116852373_A = 2 - dat$rs116852373_G", 
              "dat$rs79210031_T = 2 - dat$rs79210031_A", 
              "dat$rs139733873_G = 2 - dat$rs139733873_A", 
              "dat$rs80022746_A = 2 - dat$rs80022746_G", 
              "dat$rs2849169_A = 2 - dat$rs2849169_C", 
              "dat$rs583219_C = 2 - dat$rs583219_T", 
              "dat$rs651821_C = 2 - dat$rs651821_T", 
              "dat$rs4938309_C = 2 - dat$rs4938309_T", 
              "dat$rs595049_A = 2 - dat$rs595049_C", 
              "dat$rs11216190_T = 2 - dat$rs11216190_C", 
              "dat$rs17519093_G = 2 - dat$rs17519093_A", 
              "dat$rs888245_C = 2 - dat$rs888245_G", 
              "dat$rs518181_C = 2 - dat$rs518181_A", 
              "dat$rs116852373_G = 2 - dat$rs116852373_A", 
              "dat$rs79210031_A = 2 - dat$rs79210031_T", 
              "dat$rs139733873_A = 2 - dat$rs139733873_G", 
              "dat$rs80022746_G = 2 - dat$rs80022746_A", 
              "dat$rs2849169_C = 2 - dat$rs2849169_A", 
              "dat$rs583219_T = 2 - dat$rs583219_C", 
              "dat$rs651821_T = 2 - dat$rs651821_C",
              
              "dat$rs2284773_A = 2 - dat$rs2284773_G", 
              "dat$rs2908277_A = 2 - dat$rs2908277_G", 
              "dat$rs2908289_A = 2 - dat$rs2908289_G", 
              "dat$rs60463592_T = 2 - dat$rs60463592_C", 
              "dat$rs62459092_A = 2 - dat$rs62459092_G", 
              "dat$rs741038_A = 2 - dat$rs741038_G", 
              "dat$rs2284773_G = 2 - dat$rs2284773_A", 
              "dat$rs2908277_G = 2 - dat$rs2908277_A", 
              "dat$rs2908289_G = 2 - dat$rs2908289_A", 
              "dat$rs60463592_C = 2 - dat$rs60463592_T", 
              "dat$rs62459092_G = 2 - dat$rs62459092_A", 
              "dat$rs741038_G = 2 - dat$rs741038_A", 
              
              "dat$rs2239013_T = 2 - dat$rs2239013_C", 
              "dat$rs595049_C = 2 - dat$rs595049_A", 
              "dat$rs651821_C = 2 - dat$rs651821_T", 
              "dat$rs918143_T = 2 - dat$rs918143_C", 
              "dat$rs2239013_C = 2 - dat$rs2239013_T", 
              "dat$rs595049_A = 2 - dat$rs595049_C", 
              "dat$rs651821_T = 2 - dat$rs651821_C", 
              "dat$rs918143_C = 2 - dat$rs918143_T")

# 自动执行所有语句，报错不中断
for (cmd in commands) {
  try(eval(parse(text = cmd)), silent = TRUE)
}




dat$grs_APOC3_w_glgc_logTG <- -1 * scale((dat$rs11216190_C * (-0.102214) + 
                                            dat$rs116852373_A * (-0.131741) + 
                                            dat$rs17519093_A * (-0.10379) + 
                                            dat$rs4938309_T * (-0.139702) + 
                                            dat$rs518181_A * (-0.0577882) + 
                                            dat$rs595049_C * (0.100009) + 
                                            dat$rs651821_C * (0.290942) + 
                                            dat$rs79210031_T * (-0.182449) + 
                                            dat$rs888245_G * (0.110682) )/9)

dat$grs_APOC3_w_glgc_TC <- -1 * scale(dat$rs2239013_T * (-0.0815331) + 
                                           dat$rs595049_C * (0.0418515) + 
                                           dat$rs651821_C * (0.0571283) + 
                                           dat$rs918143_T * (-0.0237512)/4)


dat$grs_GCK_w_taiwan_glu <- -1 * scale((dat$rs2284773_G * (-1) * (-0.0817197453915558) + 
                                          dat$rs2908289_G * (-1) * (0.129313068319187) + 
                                          dat$rs60463592_C * (-1) * (-0.0636837426427698) + 
                                          dat$rs62459092_G * (-1) * (0.0345144635541533) + 
                                          dat$rs741038_G * (-1) * (0.0448525318601902))/5)


write.table(dat, "/Users/zhuyunqing/Desktop/study/drugMR_lipid/data/CAPEC/CAPEC.grs.APOC3.GCK.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")



######
#phenotype
capec.panss <- readxl::read_xlsx("CAPEC_PANSS.xlsx")
#排除两名完全无数据者
capec.panss <- subset(capec.panss, Subjid!=30051 & Subjid!=40257 & Subjid!=40254)
#对于panss个别字条目缺失者，用1填补
variables <- c("v1_P1", "v1_P2", "v1_P3", "v1_P4", "v1_P5", "v1_P6", "v1_P7", 
               "v1_N1", "v1_N2", "v1_N3", "v1_N4", "v1_N5", "v1_N6", "v1_N7", 
               "v1_G1", "v1_G2", "v1_G3", "v1_G4", "v1_G5", "v1_G6", "v1_G7", 
               "v1_G8", "v1_G9", "v1_G10", "v1_G11", "v1_G12", "v1_G13", "v1_G14", 
               "v1_G15", "v1_G16")

for (var in variables) {
  capec.panss[[var]] <- ifelse(is.na(capec.panss[[var]]), 1, capec.panss[[var]])
}

#重新加和算总panss分
capec.panss$PANSS1 <- capec.panss$v1_P1 + capec.panss$v1_P2 + capec.panss$v1_P3 + 
  capec.panss$v1_P4 + capec.panss$v1_P5 + capec.panss$v1_P6 + capec.panss$v1_P7 + 
  capec.panss$v1_N1 + capec.panss$v1_N2 + capec.panss$v1_N3 + capec.panss$v1_N4 + 
  capec.panss$v1_N5 + capec.panss$v1_N6 + capec.panss$v1_N7 + 
  capec.panss$v1_G1 + capec.panss$v1_G2 + capec.panss$v1_G3 + capec.panss$v1_G4 + 
  capec.panss$v1_G5 + capec.panss$v1_G6 + capec.panss$v1_G7 + capec.panss$v1_G8 + 
  capec.panss$v1_G9 + capec.panss$v1_G10 + capec.panss$v1_G11 + capec.panss$v1_G12 + 
  capec.panss$v1_G13 + capec.panss$v1_G14 + capec.panss$v1_G15 + capec.panss$v1_G16

capec.panss$PANSS_P1 <- capec.panss$v1_P1 + capec.panss$v1_P2 + capec.panss$v1_P3 + 
  capec.panss$v1_P4 + capec.panss$v1_P5 + capec.panss$v1_P6 + capec.panss$v1_P7

capec.panss$PANSS_N1 <- capec.panss$v1_N1 + capec.panss$v1_N2 + capec.panss$v1_N3 + capec.panss$v1_N4 + 
  capec.panss$v1_N5 + capec.panss$v1_N6 + capec.panss$v1_N7

capec.panss$PANSS_G1 <- capec.panss$v1_G1 + capec.panss$v1_G2 + capec.panss$v1_G3 + capec.panss$v1_G4 + 
  capec.panss$v1_G5 + capec.panss$v1_G6 + capec.panss$v1_G7 + capec.panss$v1_G8 + 
  capec.panss$v1_G9 + capec.panss$v1_G10 + capec.panss$v1_G11 + capec.panss$v1_G12 + 
  capec.panss$v1_G13 + capec.panss$v1_G14 + capec.panss$v1_G15 + capec.panss$v1_G16



##########
#panss v2
#######V2-PANSS
variables <- c("v2_P1", "v2_P2", "v2_P3", "v2_P4", "v2_P5", "v2_P6", "v2_P7",
               "v2_N1", "v2_N2", "v2_N3", "v2_N4", "v2_N5", "v2_N6", "v2_N7", 
               "v2_G1", "v2_G2", "v2_G3", "v2_G4", "v2_G5", "v2_G6", "v2_G7", 
               "v2_G8", "v2_G9", "v2_G10", "v2_G11", "v2_G12", "v2_G13", "v2_G14", 
               "v2_G15", "v2_G16", "PANSS_total_2")

for (var in variables) {
  capec.panss[[var]] <- ifelse(is.na(capec.panss$PANSS_total_2), -9, capec.panss[[var]])
}


#对于panss个别字条目缺失者，用1填补
variables <- c("v2_P1", "v2_P2", "v2_P3", "v2_P4", "v2_P5", "v2_P6", "v2_P7",
               "v2_N1", "v2_N2", "v2_N3", "v2_N4", "v2_N5", "v2_N6", "v2_N7", 
               "v2_G1", "v2_G2", "v2_G3", "v2_G4", "v2_G5", "v2_G6", "v2_G7", 
               "v2_G8", "v2_G9", "v2_G10", "v2_G11", "v2_G12", "v2_G13", "v2_G14", 
               "v2_G15", "v2_G16")

for (var in variables) {
  capec.panss[[var]] <- as.numeric(capec.panss[[var]])
  capec.panss[[var]] <- ifelse(is.na(capec.panss[[var]]), 1, capec.panss[[var]])
}


for (var in variables) {
  print(summary(capec.panss[[var]]))
}


capec.panss$PANSS2 <- ifelse(capec.panss$PANSS_total_2 == -9, -9, 
                             capec.panss$v2_P1 + capec.panss$v2_P2 + capec.panss$v2_P3 + capec.panss$v2_P4 + 
                               capec.panss$v2_P5 + capec.panss$v2_P6 + capec.panss$v2_P7 + 
                               capec.panss$v2_N1 + capec.panss$v2_N2 + capec.panss$v2_N3 + capec.panss$v2_N4 + 
                               capec.panss$v2_N5 + capec.panss$v2_N6 + capec.panss$v2_N7 + 
                               capec.panss$v2_G1 + capec.panss$v2_G2 + capec.panss$v2_G3 + capec.panss$v2_G4 + 
                               capec.panss$v2_G5 + capec.panss$v2_G6 + capec.panss$v2_G7 + capec.panss$v2_G8 + 
                               capec.panss$v2_G9 + capec.panss$v2_G10 + capec.panss$v2_G11 + capec.panss$v2_G12 + 
                               capec.panss$v2_G13 + capec.panss$v2_G14 + capec.panss$v2_G15 + capec.panss$v2_G16)

capec.panss$PANSS_P2 <- ifelse(capec.panss$PANSS_total_2 == -9, -9, 
                               capec.panss$v2_P1 + capec.panss$v2_P2 + capec.panss$v2_P3 + capec.panss$v2_P4 + 
                                 capec.panss$v2_P5 + capec.panss$v2_P6 + capec.panss$v2_P7)

capec.panss$PANSS_N2 <- ifelse(capec.panss$PANSS_total_2 == -9, -9, 
                               capec.panss$v2_N1 + capec.panss$v2_N2 + capec.panss$v2_N3 + capec.panss$v2_N4 + 
                                 capec.panss$v2_N5 + capec.panss$v2_N6 + capec.panss$v2_N7)

capec.panss$PANSS_G2 <- ifelse(capec.panss$PANSS_total_2 == -9, -9, 
                               capec.panss$v2_G1 + capec.panss$v2_G2 + capec.panss$v2_G3 + capec.panss$v2_G4 + 
                                 capec.panss$v2_G5 + capec.panss$v2_G6 + capec.panss$v2_G7 + capec.panss$v2_G8 + 
                                 capec.panss$v2_G9 + capec.panss$v2_G10 + capec.panss$v2_G11 + capec.panss$v2_G12 + 
                                 capec.panss$v2_G13 + capec.panss$v2_G14 + capec.panss$v2_G15 + capec.panss$v2_G16)



#########
#panss v3
#######V3-PANSS
variables <- c("v3_P1", "v3_P2", "v3_P3", "v3_P4", "v3_P5", "v3_P6", "v3_P7", 
               "v3_N1", "v3_N2", "v3_N3", "v3_N4", "v3_N5", "v3_N6", "v3_N7", 
               "v3_G1", "v3_G2", "v3_G3", "v3_G4", "v3_G5", "v3_G6", "v3_G7", "v3_G8", 
               "v3_G9", "v3_G10", "v3_G11", "v3_G12", "v3_G13", "v3_G14", "v3_G15", "v3_G16", "PANSS_total_3")

for (var in variables) {
  capec.panss[[var]] <- ifelse(is.na(capec.panss$PANSS_total_3), -9, capec.panss[[var]])
}


#对于panss个别字条目缺失者，用1填补
variables <- c("v3_P1", "v3_P2", "v3_P3", "v3_P4", "v3_P5", "v3_P6", "v3_P7", 
               "v3_N1", "v3_N2", "v3_N3", "v3_N4", "v3_N5", "v3_N6", "v3_N7", 
               "v3_G1", "v3_G2", "v3_G3", "v3_G4", "v3_G5", "v3_G6", "v3_G7", 
               "v3_G8", "v3_G9", "v3_G10", "v3_G11", "v3_G12", "v3_G13", "v3_G14", "v3_G15", "v3_G16")

for (var in variables) {
  capec.panss[[var]] <- as.numeric(capec.panss[[var]])
  capec.panss[[var]] <- ifelse(is.na(capec.panss[[var]]), 1, capec.panss[[var]])
}


for (var in variables) {
  print(summary(capec.panss[[var]]))
}


capec.panss$PANSS3 <- ifelse(capec.panss$PANSS_total_3 == -9, -9, 
                             capec.panss$v3_P1 + capec.panss$v3_P2 + capec.panss$v3_P3 + capec.panss$v3_P4 + 
                               capec.panss$v3_P5 + capec.panss$v3_P6 + capec.panss$v3_P7 + 
                               capec.panss$v3_N1 + capec.panss$v3_N2 + capec.panss$v3_N3 + capec.panss$v3_N4 + 
                               capec.panss$v3_N5 + capec.panss$v3_N6 + capec.panss$v3_N7 + 
                               capec.panss$v3_G1 + capec.panss$v3_G2 + capec.panss$v3_G3 + capec.panss$v3_G4 + 
                               capec.panss$v3_G5 + capec.panss$v3_G6 + capec.panss$v3_G7 + capec.panss$v3_G8 + 
                               capec.panss$v3_G9 + capec.panss$v3_G10 + capec.panss$v3_G11 + capec.panss$v3_G12 + 
                               capec.panss$v3_G13 + capec.panss$v3_G14 + capec.panss$v3_G15 + capec.panss$v3_G16)

capec.panss$PANSS_P3 <- ifelse(capec.panss$PANSS_total_3 == -9, -9, 
                               capec.panss$v3_P1 + capec.panss$v3_P2 + capec.panss$v3_P3 + capec.panss$v3_P4 + 
                                 capec.panss$v3_P5 + capec.panss$v3_P6 + capec.panss$v3_P7)

capec.panss$PANSS_N3 <- ifelse(capec.panss$PANSS_total_3 == -9, -9, 
                               capec.panss$v3_N1 + capec.panss$v3_N2 + capec.panss$v3_N3 + capec.panss$v3_N4 + 
                                 capec.panss$v3_N5 + capec.panss$v3_N6 + capec.panss$v3_N7)

capec.panss$PANSS_G3 <- ifelse(capec.panss$PANSS_total_3 == -9, -9, 
                               capec.panss$v3_G1 + capec.panss$v3_G2 + capec.panss$v3_G3 + capec.panss$v3_G4 + 
                                 capec.panss$v3_G5 + capec.panss$v3_G6 + capec.panss$v3_G7 + capec.panss$v3_G8 + 
                                 capec.panss$v3_G9 + capec.panss$v3_G10 + capec.panss$v3_G11 + capec.panss$v3_G12 + 
                                 capec.panss$v3_G13 + capec.panss$v3_G14 + capec.panss$v3_G15 + capec.panss$v3_G16)



#没有随访信息者，用baseline填补
capec.panss$PANSS_LOCF <- capec.panss$PANSS3
capec.panss$PANSS_P_LOCF <- capec.panss$PANSS_P3
capec.panss$PANSS_G_LOCF <- capec.panss$PANSS_G3
capec.panss$PANSS_N_LOCF <- capec.panss$PANSS_N3



capec.panss$PANSS_LOCF <- ifelse((capec.panss$PANSS2 == -9 & capec.panss$PANSS3 == -9), 
                                 capec.panss$PANSS1, 
                                 capec.panss$PANSS_LOCF)
capec.panss$PANSS_P_LOCF <- ifelse((capec.panss$PANSS_P2 == -9 & capec.panss$PANSS_P3 == -9), 
                                   capec.panss$PANSS_P1, 
                                   capec.panss$PANSS_P_LOCF)
capec.panss$PANSS_N_LOCF <- ifelse((capec.panss$PANSS_N2 == -9 & capec.panss$PANSS_N3 == -9), 
                                   capec.panss$PANSS_N1, 
                                   capec.panss$PANSS_N_LOCF)
capec.panss$PANSS_G_LOCF <- ifelse((capec.panss$PANSS_G2 == -9 & capec.panss$PANSS_G3 == -9), 
                                   capec.panss$PANSS_G1, 
                                   capec.panss$PANSS_G_LOCF)




#capec.panss_logi <- subset(capec.panss, !(capec.panss$PANSS2 == -9 & capec.panss$PANSS3 == -9 & capec.panss$PANSS3 == -9)) #剩2107人

#剩余者以末位结转求panss随访最终值
capec.panss$PANSS_LOCF <- ifelse(capec.panss$PANSS_LOCF != -9, capec.panss$PANSS_LOCF,
                                 ifelse(capec.panss$PANSS3 != -9, capec.panss$PANSS3,
                                        ifelse(capec.panss$PANSS2 != -9, capec.panss$PANSS2, capec.panss$PANSS_LOCF)))

capec.panss$PANSS_P_LOCF <- ifelse(capec.panss$PANSS_P_LOCF != -9, capec.panss$PANSS_P_LOCF,
                                   ifelse(capec.panss$PANSS_P3 != -9, capec.panss$PANSS_P3,
                                          ifelse(capec.panss$PANSS_P2 != -9, capec.panss$PANSS_P2, capec.panss$PANSS_P_LOCF)))

capec.panss$PANSS_G_LOCF <- ifelse(capec.panss$PANSS_G_LOCF != -9, capec.panss$PANSS_G_LOCF,
                                   ifelse(capec.panss$PANSS_G3 != -9, capec.panss$PANSS_G3,
                                          ifelse(capec.panss$PANSS_G2 != -9, capec.panss$PANSS_G2, capec.panss$PANSS_G_LOCF)))

capec.panss$PANSS_N_LOCF <- ifelse(capec.panss$PANSS_N_LOCF != -9, capec.panss$PANSS_N_LOCF,
                                   ifelse(capec.panss$PANSS_N3 != -9, capec.panss$PANSS_N3,
                                          ifelse(capec.panss$PANSS_N2 != -9, capec.panss$PANSS_N2, capec.panss$PANSS_N_LOCF)))



#计算panss减分率
capec.panss$PANSS_reduce_rate <- (capec.panss$PANSS1 - capec.panss$PANSS_LOCF) / (capec.panss$PANSS1 - 30) *100
capec.panss$PANSS_P_reduce_rate <- (capec.panss$PANSS_P1 - capec.panss$PANSS_P_LOCF) / (capec.panss$PANSS_P1 ) *100
capec.panss$PANSS_G_reduce_rate <- (capec.panss$PANSS_G1 - capec.panss$PANSS_G_LOCF) / (capec.panss$PANSS_G1 ) *100
capec.panss$PANSS_N_reduce_rate <- (capec.panss$PANSS_N1 - capec.panss$PANSS_N_LOCF) / (capec.panss$PANSS_N1 ) *100

#panss绝对差值
capec.panss$PANSS_reduce <- capec.panss$PANSS1 - capec.panss$PANSS_LOCF
capec.panss$PANSS_P_reduce <- capec.panss$PANSS_P1 - capec.panss$PANSS_P_LOCF
capec.panss$PANSS_G_reduce <- capec.panss$PANSS_G1 - capec.panss$PANSS_G_LOCF
capec.panss$PANSS_N_reduce <- capec.panss$PANSS_N1 - capec.panss$PANSS_N_LOCF


#子量表panss的二分类
capec.panss$PANSS_reduce_rate_2g <- ifelse(capec.panss$PANSS_reduce_rate >= 50, 1, 0)
capec.panss$PANSS_reduce_rate_2g_med <- ifelse(capec.panss$PANSS_reduce_rate >= median(capec.panss$PANSS_reduce_rate), 1, 0)
capec.panss$PANSS_N_reduce_rate_2g <- ifelse(capec.panss$PANSS_N_reduce_rate >= 50, 1, 0)
capec.panss$PANSS_N_reduce_rate_2g_med <- ifelse(capec.panss$PANSS_N_reduce_rate >= median(capec.panss$PANSS_N_reduce_rate), 1, 0)
capec.panss$PANSS_P_reduce_rate_2g <- ifelse(capec.panss$PANSS_P_reduce_rate >= 50, 1, 0)
capec.panss$PANSS_P_reduce_rate_2g_med <- ifelse(capec.panss$PANSS_P_reduce_rate >= median(capec.panss$PANSS_P_reduce_rate), 1, 0)
capec.panss$PANSS_G_reduce_rate_2g <- ifelse(capec.panss$PANSS_G_reduce_rate >= 50, 1, 0)
capec.panss$PANSS_G_reduce_rate_2g_med <- ifelse(capec.panss$PANSS_G_reduce_rate >= median(capec.panss$PANSS_G_reduce_rate), 1, 0)


#####
#glu lipid
capec.lipid.glu <- read.csv("CAPEC.lipid.glu.csv")
capec.lipid.glu <- subset(capec.lipid.glu, !is.na(capec.lipid.glu$v1_tg) & !is.na(capec.lipid.glu$v1_glu))

#lipid_glu_trans
capec.lipid.glu$v1_tc_mg <- capec.lipid.glu$v1_tc * 38.67
capec.lipid.glu$v1_tg_mg <- capec.lipid.glu$v1_tg * 88.57
capec.lipid.glu$v1_hdl_mg <- capec.lipid.glu$v1_hdl * 38.67
capec.lipid.glu$v1_ldl_mg <- capec.lipid.glu$v1_ldl * 38.67
capec.lipid.glu$v1_glu_mg <- capec.lipid.glu$v1_glu * 18.02



#######
#merge
m <- merge(capec.bl, dat, by.x = "Subjid", by.y = "FID", all = F)
m <- merge(m, capec.pc, by.x = "Subjid", by.y = "IID", all = F)
m <- merge(m, capec.panss, by = "Subjid", all = F)
m <- merge(m, capec.lipid.glu, by = "Subjid", all = F)

fwrite(m, "/Users/zhuyunqing/Desktop/study/drugMR_lipid/data/CAPEC/CAPEC.grs.pheno.new.txt",
       sep = "\t")



