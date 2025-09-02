#重新清洗数据
#863数据预处理

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forestploter)

library(corrplot)
library(reshape2)
library(grid)
library(ggrepel)

library("rio")

setwd("/Users/zhuyunqing/Desktop/study/SCZ_resilience/data")
main_3030 <- fread("/Users/zhuyunqing/Desktop/study/863/cleaned.mainbank_863_n3030.csv",h=T)

#重复id
main_3030_dup <- main_3030[duplicated(main_3030$dn),]
#重复列名
main_3030_dup_columns <- duplicated(names(main_3030))
names(main_3030)[main_3030_dup_columns]

# diff_vd1 <- main_3030$vd1 != main_3030$vd1_1
# diff_vd2 <- main_3030$vd2 != main_3030$vd2_1
# diff_vd3 <- main_3030$vd3 != main_3030$vd3_1
# diff_vd4 <- main_3030$vd4 != main_3030$vd4_1
# 
# diff_PANSS1 <- main_3030$PANSS1 != main_3030$PANSS1_1
# diff_PANSS2 <- main_3030$PANSS2 != main_3030$PANSS2_1
# diff_PANSS3 <- main_3030$PANSS3 != main_3030$PANSS3_1
# diff_PANSS4 <- main_3030$PANSS4 != main_3030$PANSS4_1
# 
# 
# main_3030[diff_vd1,]
# main_3030[diff_vd2,]
# main_3030[diff_vd3,]
# main_3030[diff_vd4,]
# 
# main_3030[diff_PANSS1,]
# main_3030[diff_PANSS2,]
# main_3030[diff_PANSS3,]
# main_3030[diff_PANSS4,]






#合并pc数据
pc <- fread("/Users/zhuyunqing/Desktop/study/863/imputed_clean_qc.eigenvec",header=FALSE)
colnames(pc) <- c("FID", "IID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7",
                  "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15",
                  "PC16", "PC17", "PC18", "PC19", "PC20")

main <- merge(main_3030, pc, by.x="dn", by.y="FID", all=F)
main <- subset(main, IID!="")







#基线变量纳排
#性别 年龄 center
main[is.na(main$sex),]
main[is.na(main$age),]
main[is.na(main$centers),]
main[is.na(main$p041011),] #教育

main[is.na(main$PANSS1),] #10人缺失基线panss
main[is.na(main$p041021),] #首发年龄
main[is.na(main$p041024),] #首发

#保留2556人


nrow(main[is.na(main$course),]) #病程
nrow(main[is.na(main$bmi0601),]) #bmi缺失84人
nrow(main[is.na(main$weg0601),]) #体重缺1人
nrow(main[is.na(main$fw0601),]) #腹围缺3人



nrow(main[is.na(main$drugp07),]) #合并用药 缺6人
nrow(main[is.na(main$qtc0801),]) #qtc缺175人
nrow(main[is.na(main$qrs0801),]) #QRS缺82人
nrow(main[is.na(main$qtc3201),])
nrow(main[is.na(main$pls0601),]) #心率缺失1
nrow(main[is.na(main$sbp0601),]) #血压无缺失
nrow(main[is.na(main$dbp0601),])
nrow(main[is.na(main$tc1),]) #缺失24人
nrow(main[is.na(main$tg1),]) #缺失25人
nrow(main[is.na(main$hdl1),]) #缺失158人
nrow(main[is.na(main$ldl1),]) #缺失165人
nrow(main[is.na(main$glu1),]) #缺失11人
nrow(main[is.na(main$hba1)]) #缺失1922-不用了！

#血细胞
nrow(main[is.na(main$rbc1),]) #红细胞缺7人
nrow(main[is.na(main$hgb1),]) #血红蛋白缺8人
nrow(main[is.na(main$wbc1),]) #白细胞缺8人
nrow(main[is.na(main$plt1),]) #血小板缺8人

#代谢
nrow(main[is.na(main$alt1),]) #alt缺7人
nrow(main[is.na(main$ast1),]) #ast缺9人

nrow(main[is.na(main$bun1),]) #尿素氮缺12人
nrow(main[is.na(main$cre1),]) #肌酐缺14人


#######
#排除缺失值

main_exc_na <- main[!(is.na(main$course)),] 
main_exc_na <- main_exc_na[!(is.na(main_exc_na$drugp07)),] #2554

main_exc_na <- main_exc_na[!(is.na(main_exc_na$bmi0601)),] 
main_exc_na <- main_exc_na[!(is.na(main_exc_na$fw0601)),] 
main_exc_na <- main_exc_na[!(is.na(main_exc_na$weg0601)),] #2467


main_exc_na <- main_exc_na[!(is.na(main_exc_na$sbp0601)),] #2467
main_exc_na <- main_exc_na[!(is.na(main_exc_na$dbp0601)),] #2467

# main_exc_na <- main_exc_na[!(is.na(main_exc_na$pls0601)),]
# main_exc_na <- main_exc_na[!(is.na(main_exc_na$qrs0801)),]
# main_exc_na <- main_exc_na[!(is.na(main_exc_na$qtc0801)),] #2342

main_exc_na <- main_exc_na[!(is.na(main_exc_na$tc1)),] 
main_exc_na <- main_exc_na[!(is.na(main_exc_na$tg1)),] 
main_exc_na <- main_exc_na[!(is.na(main_exc_na$hdl1)),]  
main_exc_na <- main_exc_na[!(is.na(main_exc_na$ldl1)),]  
main_exc_na <- main_exc_na[!(is.na(main_exc_na$glu1)),] #2317


#剔除连续性变量的临床异常范围值
main_exc_na <- subset(main_exc_na, age>=18)
main_exc_na <- main_exc_na[abs(as.numeric(scale(bmi0601))) <= 3] #2288个 29


main_exc_na <- main_exc_na[abs(as.numeric(scale(glu1))) <= 3] #2307
main_exc_na <- main_exc_na[abs(as.numeric(scale(tc1))) <= 3] #2302
main_exc_na <- main_exc_na[abs(as.numeric(scale(tg1))) <= 3] #2291
main_exc_na <- main_exc_na[abs(as.numeric(scale(hdl1))) <= 3] #2285
main_exc_na <- main_exc_na[abs(as.numeric(scale(ldl1))) <= 3] #2158 130

main_exc_na <- main_exc_na[abs(as.numeric(scale(fw0601))) <= 3] #2282
main_exc_na <- main_exc_na[abs(as.numeric(scale(weg0601))) <= 3] #2124 34

main_exc_na <- main_exc_na[abs(as.numeric(scale(sbp0601))) <= 3]
main_exc_na <- main_exc_na[abs(as.numeric(scale(dbp0601))) <= 3] #2115 9


# main_exc_na <- subset(main_exc_na, bmi0601 >= 10 & bmi0601 <=50) #3
# main_exc_na <- subset(main_exc_na, glu1 >= 2.5 & glu1 <= 30) #4
# main_exc_na <- subset(main_exc_na, tc1 >= 1 & tc1 <= 15) #1
# main_exc_na <- subset(main_exc_na, tg1 >= 0.2 & tg1 <= 10) #8
# main_exc_na <- subset(main_exc_na, hdl1 >= 0.3 & hdl1 <= 5) #5
# main_exc_na <- subset(main_exc_na, ldl1 >= 0.5 & ldl1 <= 10) #8
# main_exc_na <- subset(main_exc_na, fw0601 >= 50 & fw0601 <= 150)
# main_exc_na <- subset(main_exc_na, weg0601 >= 30 & weg0601 <= 200)
# main_exc_na <- subset(main_exc_na, sbp0601 >= 60 & sbp0601 <= 250)
# main_exc_na <- subset(main_exc_na, dbp0601 >= 40 & dbp0601 <=150) #3

#余2288








#肝酶
# main_exc_na <- main_exc_na[!(is.na(main_exc_na$alt1)),] 
# main_exc_na <- main_exc_na[!(is.na(main_exc_na$ast1)),] #2197



#基线服用抗心血管病药物者-先不排除 算分
#120081-地奥心血康-扩张冠脉血管，改善心肌缺血
#120767-硝酸异山梨酯-用于治疗心绞痛，充血性心力衰竭和食管痉挛
#510227-阿替洛尔片-【治疗高血压 心律失常】


#随访期间服cvd药物者
#150102-普萘洛尔-高血压、劳力型心绞痛、心律失常(2周-结束)
#420430-普萘洛尔-高血压、劳力型心绞痛、心律失常（1周-结束）
#460499-贝他乐克-心肌缺血、快速性心律失常和胸痛（1周-结束）
#250365-洛伐他汀-高脂血症（1周-结束）

#和文章保持一致，排除基线患有心血管病者
id_anticvd_drug_prev <- c("120081", "120767", "510227")
main_exc_na <- subset(main_exc_na, dn!="120081" & dn!="120767" & dn!="510227")


#####
#血糖
# 210043-入院后三天确诊糖尿病
# 530159 糖尿病



# #血细胞
# main_exc_na <- main_exc_na[!(is.na(main_exc_na$rbc1)),]  
# main_exc_na <- main_exc_na[!(is.na(main_exc_na$hgb1)),]  
# main_exc_na <- main_exc_na[!(is.na(main_exc_na$wbc1)),]  
# main_exc_na <- main_exc_na[!(is.na(main_exc_na$plt1)),]
# 
# #尿素氮 肌酐
# main_exc_na <- main_exc_na[!(is.na(main_exc_na$bun1)),]
# main_exc_na <- main_exc_na[!(is.na(main_exc_na$cre1)),] #2189


#lipid_trans
main_exc_na$tc1_mg <- main_exc_na$tc1 * 38.67
main_exc_na$tg1_mg <- main_exc_na$tg1 * 88.57
main_exc_na$hdl1_mg <- main_exc_na$hdl1 * 38.67
main_exc_na$ldl1_mg <- main_exc_na$ldl1 * 38.67
main_exc_na$glu1_mg <- main_exc_na$glu1 * 18.02


#修改一下age2的定义
main_exc_na$age2 <- main_exc_na$age * main_exc_na$age


##############
#panss数据预处理


nrow(main_exc_na[is.na(main_exc_na$PANSS1),])
#main_exc_na <- subset(main_exc_na, PANSS1!=0)


#首先删掉基线panss完全无数据的一个人
nrow(main_exc_na[is.na(main_exc_na$PANSS1) & is.na(main_exc_na$PANSS1_N) & 
                   is.na(main_exc_na$PANSS1_P) & is.na(main_exc_na$PANSS1_G),])
main_exc_na <- subset(main_exc_na, !(is.na(main_exc_na$PANSS1) & is.na(main_exc_na$PANSS1_N) & 
                                       is.na(main_exc_na$PANSS1_P) & is.na(main_exc_na$PANSS1_G)))

#对于panss个别字条目缺失者，用1填补
variables <- c("p091011", "p091012", "p091013", "p091014", "p091015", "p091016", 
               "p091017", "p092011", "p092012", "p092013", "p092014", "p092015", 
               "p092016", "p092017", "p093011", "p093012", "p093013", "p093014", 
               "p093015", "p093016", "p093017", "p093018", "p093019", "p0930110", 
               "p0930111", "p0930112", "p0930113", "p0930114", "p0930115", "p0930116")

for (var in variables) {
  main_exc_na[[var]] <- ifelse(is.na(main_exc_na[[var]]), 1, main_exc_na[[var]])
}


#重新加和算总panss分
main_exc_na$PANSS1 <- main_exc_na$p091011+ main_exc_na$p091012+ main_exc_na$p091013+ 
  main_exc_na$p091014+ main_exc_na$p091015+ main_exc_na$p091016+ main_exc_na$p091017+ 
  main_exc_na$p092011+ main_exc_na$p092012+ main_exc_na$p092013+ main_exc_na$p092014+ 
  main_exc_na$p092015+ main_exc_na$p092016+ main_exc_na$p092017+ main_exc_na$p093011+ 
  main_exc_na$p093012+ main_exc_na$p093013+ main_exc_na$p093014+ main_exc_na$p093015+ 
  main_exc_na$p093016+ main_exc_na$p093017+ main_exc_na$p093018+ main_exc_na$p093019+ 
  main_exc_na$p0930110+ main_exc_na$p0930111+ main_exc_na$p0930112+ main_exc_na$p0930113+ 
  main_exc_na$p0930114+ main_exc_na$p0930115+ main_exc_na$p0930116


main_exc_na$PANSS1_P <- main_exc_na$p091011+ main_exc_na$p091012+ main_exc_na$p091013+ 
  main_exc_na$p091014+ main_exc_na$p091015+ main_exc_na$p091016+ main_exc_na$p091017

main_exc_na$PANSS1_N <- main_exc_na$p092011+ main_exc_na$p092012+ main_exc_na$p092013+ 
  main_exc_na$p092014+ main_exc_na$p092015+ main_exc_na$p092016+ main_exc_na$p092017

main_exc_na$PANSS1_G <- main_exc_na$p093011+ main_exc_na$p093012+ main_exc_na$p093013+ 
  main_exc_na$p093014+ main_exc_na$p093015+ main_exc_na$p093016+ main_exc_na$p093017+ 
  main_exc_na$p093018+ main_exc_na$p093019+ main_exc_na$p0930110+ main_exc_na$p0930111+ 
  main_exc_na$p0930112+ main_exc_na$p0930113+ main_exc_na$p0930114+ main_exc_na$p0930115+ main_exc_na$p0930116



########
#panss五因子划分
main_exc_na$PANSS1_F5_P <- main_exc_na$p091011 + main_exc_na$p091013 + main_exc_na$p093019 + 
  main_exc_na$p091016 + main_exc_na$p091015 + main_exc_na$p0930112
main_exc_na$PANSS1_F5_N <- main_exc_na$p092012 + main_exc_na$p092014 + main_exc_na$p092016 + 
  main_exc_na$p092011 + main_exc_na$p092013 + main_exc_na$p093017 + main_exc_na$p0930116
main_exc_na$PANSS1_F5_CD <- main_exc_na$p091012 + main_exc_na$p0930111 + main_exc_na$p092017 + 
  main_exc_na$p092015 + main_exc_na$p0930113 + main_exc_na$p093015 + main_exc_na$p0930110
main_exc_na$PANSS1_F5_DA <- main_exc_na$p093012 + main_exc_na$p093013 + main_exc_na$p093016
main_exc_na$PANSS1_F5_H <- main_exc_na$p091014 + main_exc_na$p091017 + main_exc_na$p093018 + 
  main_exc_na$p0930114 + main_exc_na$p093014




#######
#随访panss处理 减分率计算


#2week

nrow(main_exc_na[is.na(main_exc_na$PANSS2) & is.na(main_exc_na$PANSS2_N) & 
                   is.na(main_exc_na$PANSS2_P) & is.na(main_exc_na$PANSS2_G),])
missing_2nd <- (main_exc_na[is.na(main_exc_na$PANSS2) & is.na(main_exc_na$PANSS2_N) & 
                              is.na(main_exc_na$PANSS2_P) & is.na(main_exc_na$PANSS2_G),])
missing_2nd <- as.data.table(cbind(missing_2nd$dn, missing_2nd$p171011, 
                                   missing_2nd$p171012, missing_2nd$p171013, missing_2nd$p171014, 
                                   missing_2nd$p171015, missing_2nd$p171016, missing_2nd$p171017, 
                                   missing_2nd$p172011, missing_2nd$p172012, missing_2nd$p172013, 
                                   missing_2nd$p172014, missing_2nd$p172015, missing_2nd$p172016, 
                                   missing_2nd$p172017, missing_2nd$p173011, missing_2nd$p173012, 
                                   missing_2nd$p173013, missing_2nd$p173014, missing_2nd$p173015, 
                                   missing_2nd$p173016, missing_2nd$p173017, missing_2nd$p173018, 
                                   missing_2nd$p173019, missing_2nd$p1730110, missing_2nd$p1730111, 
                                   missing_2nd$p1730112, missing_2nd$p1730113, missing_2nd$p1730114, 
                                   missing_2nd$p1730115, missing_2nd$p1730116, missing_2nd$p174011))
sapply(missing_2nd, function(x) sum(is.na(x)))
#这173个人没有评估panss-标为-9


variables <- c("p171011", "p171012", "p171013", "p171014", "p171015", "p171016", 
               "p171017", "p172011", "p172012", "p172013", "p172014", "p172015", "p172016", 
               "p172017", "p173011", "p173012", "p173013", "p173014", "p173015", "p173016", 
               "p173017", "p173018", "p173019", "p1730110", "p1730111", "p1730112", "p1730113", 
               "p1730114", "p1730115", "p1730116", "p174011")

for (var in variables) {
  main_exc_na[[var]] <- ifelse((is.na(main_exc_na$PANSS2) & is.na(main_exc_na$PANSS2_N) & 
                                  is.na(main_exc_na$PANSS2_P) & is.na(main_exc_na$PANSS2_G)), 
                               -9, 
                               main_exc_na[[var]])
}





#对于panss个别字条目缺失者，用1填补
variables <- c("p171011", "p171012", "p171013", "p171014", "p171015", "p171016", 
               "p171017", "p172011", "p172012", "p172013", "p172014", "p172015", "p172016", 
               "p172017", "p173011", "p173012", "p173013", "p173014", "p173015", "p173016", 
               "p173017", "p173018", "p173019", "p1730110", "p1730111", "p1730112", "p1730113", 
               "p1730114", "p1730115", "p1730116")



for (var in variables) {
  main_exc_na[[var]] <- as.numeric(main_exc_na[[var]])
  main_exc_na[[var]] <- ifelse(is.na(main_exc_na[[var]]), 1, main_exc_na[[var]])
}

for (var in variables) {
  print(summary(main_exc_na[[var]]))
}

#计算2 week量表分
main_exc_na$PANSS2 <- ifelse(main_exc_na$p174011==-9, -9, 
                             main_exc_na$p171011 + main_exc_na$p171012 + main_exc_na$p171013 + 
                               main_exc_na$p171014 + main_exc_na$p171015 + main_exc_na$p171016 + 
                               main_exc_na$p171017 + main_exc_na$p172011 + main_exc_na$p172012 + 
                               main_exc_na$p172013 + main_exc_na$p172014 + main_exc_na$p172015 + 
                               main_exc_na$p172016 + main_exc_na$p172017 + main_exc_na$p173011 + 
                               main_exc_na$p173012 + main_exc_na$p173013 + main_exc_na$p173014 + 
                               main_exc_na$p173015 + main_exc_na$p173016 + main_exc_na$p173017 + 
                               main_exc_na$p173018 + main_exc_na$p173019 + main_exc_na$p1730110 + 
                               main_exc_na$p1730111 + main_exc_na$p1730112 + main_exc_na$p1730113 + 
                               main_exc_na$p1730114 + main_exc_na$p1730115 + main_exc_na$p1730116)

main_exc_na$PANSS2_P <- ifelse(main_exc_na$p174011==-9, -9, 
                               main_exc_na$p171011 + main_exc_na$p171012 + main_exc_na$p171013 + 
                                 main_exc_na$p171014 + main_exc_na$p171015 + main_exc_na$p171016 + main_exc_na$p171017)

main_exc_na$PANSS2_N <- ifelse(main_exc_na$p174011==-9, -9,
                               main_exc_na$p172011 + main_exc_na$p172012 + 
                                 main_exc_na$p172013 + main_exc_na$p172014 + main_exc_na$p172015 + 
                                 main_exc_na$p172016 + main_exc_na$p172017)

main_exc_na$PANSS2_G <- ifelse(main_exc_na$p174011==-9, -9,
                               main_exc_na$p173011 + 
                                 main_exc_na$p173012 + main_exc_na$p173013 + main_exc_na$p173014 + 
                                 main_exc_na$p173015 + main_exc_na$p173016 + main_exc_na$p173017 + 
                                 main_exc_na$p173018 + main_exc_na$p173019 + main_exc_na$p1730110 + 
                                 main_exc_na$p1730111 + main_exc_na$p1730112 + main_exc_na$p1730113 + 
                                 main_exc_na$p1730114 + main_exc_na$p1730115 + main_exc_na$p1730116)


########
#panss五因子划分
main_exc_na$PANSS2_F5_P <- main_exc_na$p171011 + main_exc_na$p171013 + main_exc_na$p173019 + 
  main_exc_na$p171016 + main_exc_na$p171015 + main_exc_na$p1730112
main_exc_na$PANSS2_F5_N <- main_exc_na$p172012 + main_exc_na$p172014 + main_exc_na$p172016 + 
  main_exc_na$p172011 + main_exc_na$p172013 + main_exc_na$p173017 + main_exc_na$p1730116
main_exc_na$PANSS2_F5_CD <- main_exc_na$p171012 + main_exc_na$p1730111 + main_exc_na$p172017 + 
  main_exc_na$p172015 + main_exc_na$p1730113 + main_exc_na$p173015 + main_exc_na$p1730110
main_exc_na$PANSS2_F5_DA <- main_exc_na$p173012 + main_exc_na$p173013 + main_exc_na$p173016
main_exc_na$PANSS2_F5_H <- main_exc_na$p171014 + main_exc_na$p171017 + main_exc_na$p173018 + 
  main_exc_na$p1730114 + main_exc_na$p173014



#######
#4 week
nrow(main_exc_na[is.na(main_exc_na$PANSS3) & is.na(main_exc_na$PANSS3_N) & 
                   is.na(main_exc_na$PANSS3_P) & is.na(main_exc_na$PANSS3_G),])
missing_3rd <- (main_exc_na[is.na(main_exc_na$PANSS3) & is.na(main_exc_na$PANSS3_N) & 
                              is.na(main_exc_na$PANSS3_P) & is.na(main_exc_na$PANSS3_G),])
missing_3rd <- as.data.table(cbind(missing_3rd$dn, missing_3rd$p251011, 
                                   missing_3rd$p251012, missing_3rd$p251013, missing_3rd$p251014, 
                                   missing_3rd$p251015, missing_3rd$p251016, missing_3rd$p251017, 
                                   missing_3rd$p252011, missing_3rd$p252012, missing_3rd$p252013, 
                                   missing_3rd$p252014, missing_3rd$p252015, missing_3rd$p252016, 
                                   missing_3rd$p252017, missing_3rd$p253011, missing_3rd$p253012, 
                                   missing_3rd$p253013, missing_3rd$p253014, missing_3rd$p253015, 
                                   missing_3rd$p253016, missing_3rd$p253017, missing_3rd$p253018, 
                                   missing_3rd$p253019, missing_3rd$p2530110, missing_3rd$p2530111, 
                                   missing_3rd$p2530112, missing_3rd$p2530113, missing_3rd$p2530114, 
                                   missing_3rd$p2530115, missing_3rd$p2530116, missing_3rd$p254011))
sapply(missing_3rd, function(x) sum(is.na(x)))



variables <- c("p251011", "p251012", "p251013", "p251014", "p251015", "p251016", 
               "p251017", "p252011", "p252012", "p252013", "p252014", "p252015", "p252016", 
               "p252017", "p253011", "p253012", "p253013", "p253014", "p253015", "p253016", 
               "p253017", "p253018", "p253019", "p2530110", "p2530111", "p2530112", "p2530113", 
               "p2530114", "p2530115", "p2530116", "p254011")

for (var in variables) {
  main_exc_na[[var]] <- ifelse((is.na(main_exc_na$PANSS3) & is.na(main_exc_na$PANSS3_N) & 
                                  is.na(main_exc_na$PANSS3_P) & is.na(main_exc_na$PANSS3_G)), -9, main_exc_na[[var]])
}


#对于panss个别字条目缺失者，用1填补
variables <- c("p251011", "p251012", "p251013", "p251014", "p251015", "p251016", 
               "p251017", "p252011", "p252012", "p252013", "p252014", "p252015", "p252016", 
               "p252017", "p253011", "p253012", "p253013", "p253014", "p253015", "p253016", 
               "p253017", "p253018", "p253019", "p2530110", "p2530111", "p2530112", "p2530113", 
               "p2530114", "p2530115", "p2530116")



for (var in variables) {
  main_exc_na[[var]] <- as.numeric(main_exc_na[[var]])
  main_exc_na[[var]] <- ifelse(is.na(main_exc_na[[var]]), 1, main_exc_na[[var]])
  
}


for (var in variables) {
  print(summary(main_exc_na[[var]]))
}


#计算2 week量表分
main_exc_na$PANSS3 <- ifelse(main_exc_na$p254011==-9, -9, 
                             main_exc_na$p251011 + main_exc_na$p251012 + main_exc_na$p251013 + 
                               main_exc_na$p251014 + main_exc_na$p251015 + main_exc_na$p251016 + 
                               main_exc_na$p251017 + main_exc_na$p252011 + main_exc_na$p252012 + 
                               main_exc_na$p252013 + main_exc_na$p252014 + main_exc_na$p252015 + 
                               main_exc_na$p252016 + main_exc_na$p252017 + main_exc_na$p253011 + 
                               main_exc_na$p253012 + main_exc_na$p253013 + main_exc_na$p253014 + 
                               main_exc_na$p253015 + main_exc_na$p253016 + main_exc_na$p253017 + 
                               main_exc_na$p253018 + main_exc_na$p253019 + main_exc_na$p2530110 + 
                               main_exc_na$p2530111 + main_exc_na$p2530112 + main_exc_na$p2530113 + 
                               main_exc_na$p2530114 + main_exc_na$p2530115 + main_exc_na$p2530116)

main_exc_na$PANSS3_P <- ifelse(main_exc_na$p254011==-9, -9, 
                               main_exc_na$p251011 + main_exc_na$p251012 + main_exc_na$p251013 + 
                                 main_exc_na$p251014 + main_exc_na$p251015 + main_exc_na$p251016 + main_exc_na$p251017)

main_exc_na$PANSS3_N <- ifelse(main_exc_na$p254011==-9, -9,
                               main_exc_na$p252011 + main_exc_na$p252012 + 
                                 main_exc_na$p252013 + main_exc_na$p252014 + main_exc_na$p252015 + 
                                 main_exc_na$p252016 + main_exc_na$p252017)

main_exc_na$PANSS3_G <- ifelse(main_exc_na$p254011==-9, -9,
                               main_exc_na$p253011 + 
                                 main_exc_na$p253012 + main_exc_na$p253013 + main_exc_na$p253014 + 
                                 main_exc_na$p253015 + main_exc_na$p253016 + main_exc_na$p253017 + 
                                 main_exc_na$p253018 + main_exc_na$p253019 + main_exc_na$p2530110 + 
                                 main_exc_na$p2530111 + main_exc_na$p2530112 + main_exc_na$p2530113 + 
                                 main_exc_na$p2530114 + main_exc_na$p2530115 + main_exc_na$p2530116)


########
#panss五因子划分
main_exc_na$PANSS3_F5_P <- main_exc_na$p251011 + main_exc_na$p251013 + main_exc_na$p253019 + 
  main_exc_na$p251016 + main_exc_na$p251015 + main_exc_na$p2530112
main_exc_na$PANSS3_F5_N <- main_exc_na$p252012 + main_exc_na$p252014 + main_exc_na$p252016 + 
  main_exc_na$p252011 + main_exc_na$p252013 + main_exc_na$p253017 + main_exc_na$p2530116
main_exc_na$PANSS3_F5_CD <- main_exc_na$p251012 + main_exc_na$p2530111 + main_exc_na$p252017 + 
  main_exc_na$p252015 + main_exc_na$p2530113 + main_exc_na$p253015 + main_exc_na$p2530110
main_exc_na$PANSS3_F5_DA <- main_exc_na$p253012 + main_exc_na$p253013 + main_exc_na$p253016
main_exc_na$PANSS3_F5_H <- main_exc_na$p251014 + main_exc_na$p251017 + main_exc_na$p253018 + 
  main_exc_na$p2530114 + main_exc_na$p253014







#######
#6 week
nrow(main_exc_na[is.na(main_exc_na$PANSS4) & is.na(main_exc_na$PANSS4_N) & 
                   is.na(main_exc_na$PANSS4_P) & is.na(main_exc_na$PANSS4_G),])
missing_4th <- (main_exc_na[is.na(main_exc_na$PANSS4) & is.na(main_exc_na$PANSS4_N) & 
                              is.na(main_exc_na$PANSS4_P) & is.na(main_exc_na$PANSS4_G),])
missing_4th <- as.data.table(cbind(missing_4th$dn, missing_4th$p331011, 
                                   missing_4th$p331012, missing_4th$p331013, missing_4th$p331014, 
                                   missing_4th$p331015, missing_4th$p331016, missing_4th$p331017, 
                                   missing_4th$p332011, missing_4th$p332012, missing_4th$p332013, 
                                   missing_4th$p332014, missing_4th$p332015, missing_4th$p332016, 
                                   missing_4th$p332017, missing_4th$p333011, missing_4th$p333012, 
                                   missing_4th$p333013, missing_4th$p333014, missing_4th$p333015, 
                                   missing_4th$p333016, missing_4th$p333017, missing_4th$p333018, 
                                   missing_4th$p333019, missing_4th$p3330110, missing_4th$p3330111, 
                                   missing_4th$p3330112, missing_4th$p3330113, missing_4th$p3330114, 
                                   missing_4th$p3330115, missing_4th$p3330116, missing_4th$p334011))
sapply(missing_4th, function(x) sum((x=="")))



variables <- c("p331011", "p331012", "p331013", "p331014", "p331015", "p331016", 
               "p331017", "p332011", "p332012", "p332013", "p332014", "p332015", "p332016", 
               "p332017", "p333011", "p333012", "p333013", "p333014", "p333015", "p333016", 
               "p333017", "p333018", "p333019", "p3330110", "p3330111", "p3330112", "p3330113", 
               "p3330114", "p3330115", "p3330116", "p334011")

for (var in variables) {
  main_exc_na[[var]] <- ifelse((is.na(main_exc_na$PANSS4) & is.na(main_exc_na$PANSS4_N) & 
                                  is.na(main_exc_na$PANSS4_P) & is.na(main_exc_na$PANSS4_G)), -9, main_exc_na[[var]])
}


#对于panss个别字条目缺失者，用1填补
variables <- c("p331011", "p331012", "p331013", "p331014", "p331015", "p331016", 
               "p331017", "p332011", "p332012", "p332013", "p332014", "p332015", "p332016", 
               "p332017", "p333011", "p333012", "p333013", "p333014", "p333015", "p333016", 
               "p333017", "p333018", "p333019", "p3330110", "p3330111", "p3330112", "p3330113", 
               "p3330114", "p3330115", "p3330116")

for (var in variables) {
  main_exc_na[[var]] <- as.numeric(main_exc_na[[var]])
  main_exc_na[[var]] <- ifelse(is.na(main_exc_na[[var]]), 1, main_exc_na[[var]])
}

for (var in variables) {
  print(summary(main_exc_na[[var]]))
}


#计算6 week量表分
main_exc_na$PANSS4 <- ifelse(main_exc_na$p334011==-9, -9, 
                             main_exc_na$p331011 + main_exc_na$p331012 + main_exc_na$p331013 + 
                               main_exc_na$p331014 + main_exc_na$p331015 + main_exc_na$p331016 + 
                               main_exc_na$p331017 + main_exc_na$p332011 + main_exc_na$p332012 + 
                               main_exc_na$p332013 + main_exc_na$p332014 + main_exc_na$p332015 + 
                               main_exc_na$p332016 + main_exc_na$p332017 + main_exc_na$p333011 + 
                               main_exc_na$p333012 + main_exc_na$p333013 + main_exc_na$p333014 + 
                               main_exc_na$p333015 + main_exc_na$p333016 + main_exc_na$p333017 + 
                               main_exc_na$p333018 + main_exc_na$p333019 + main_exc_na$p3330110 + 
                               main_exc_na$p3330111 + main_exc_na$p3330112 + main_exc_na$p3330113 + 
                               main_exc_na$p3330114 + main_exc_na$p3330115 + main_exc_na$p3330116)

main_exc_na$PANSS4_P <- ifelse(main_exc_na$p334011==-9, -9, 
                               main_exc_na$p331011 + main_exc_na$p331012 + main_exc_na$p331013 + 
                                 main_exc_na$p331014 + main_exc_na$p331015 + main_exc_na$p331016 + main_exc_na$p331017)

main_exc_na$PANSS4_N <- ifelse(main_exc_na$p334011==-9, -9,
                               main_exc_na$p332011 + main_exc_na$p332012 + 
                                 main_exc_na$p332013 + main_exc_na$p332014 + main_exc_na$p332015 + 
                                 main_exc_na$p332016 + main_exc_na$p332017)

main_exc_na$PANSS4_G <- ifelse(main_exc_na$p334011==-9, -9,
                               main_exc_na$p333011 + 
                                 main_exc_na$p333012 + main_exc_na$p333013 + main_exc_na$p333014 + 
                                 main_exc_na$p333015 + main_exc_na$p333016 + main_exc_na$p333017 + 
                                 main_exc_na$p333018 + main_exc_na$p333019 + main_exc_na$p3330110 + 
                                 main_exc_na$p3330111 + main_exc_na$p3330112 + main_exc_na$p3330113 + 
                                 main_exc_na$p3330114 + main_exc_na$p3330115 + main_exc_na$p3330116)



########
#panss五因子划分
main_exc_na$PANSS4_F5_P <- main_exc_na$p331011 + main_exc_na$p331013 + main_exc_na$p333019 + 
  main_exc_na$p331016 + main_exc_na$p331015 + main_exc_na$p3330112
main_exc_na$PANSS4_F5_N <- main_exc_na$p332012 + main_exc_na$p332014 + main_exc_na$p332016 + 
  main_exc_na$p332011 + main_exc_na$p332013 + main_exc_na$p333017 + main_exc_na$p3330116
main_exc_na$PANSS4_F5_CD <- main_exc_na$p331012 + main_exc_na$p3330111 + main_exc_na$p332017 + 
  main_exc_na$p332015 + main_exc_na$p3330113 + main_exc_na$p333015 + main_exc_na$p3330110
main_exc_na$PANSS4_F5_DA <- main_exc_na$p333012 + main_exc_na$p333013 + main_exc_na$p333016
main_exc_na$PANSS4_F5_H <- main_exc_na$p331014 + main_exc_na$p331017 + main_exc_na$p333018 + 
  main_exc_na$p3330114 + main_exc_na$p333014




#没有随访信息者，用baseline填补
main_exc_na$PANSS_LOCF <- main_exc_na$PANSS4
main_exc_na$PANSS_P_LOCF <- main_exc_na$PANSS4_P
main_exc_na$PANSS_G_LOCF <- main_exc_na$PANSS4_G
main_exc_na$PANSS_N_LOCF <- main_exc_na$PANSS4_N



main_exc_na$PANSS_LOCF <- ifelse((main_exc_na$PANSS2 == -9 & main_exc_na$PANSS3 == -9 & main_exc_na$PANSS4 == -9), 
                                 main_exc_na$PANSS1, 
                                 main_exc_na$PANSS_LOCF)
main_exc_na$PANSS_P_LOCF <- ifelse((main_exc_na$PANSS2_P == -9 & main_exc_na$PANSS3_P == -9 & main_exc_na$PANSS4_P == -9), 
                                   main_exc_na$PANSS1_P, 
                                   main_exc_na$PANSS_P_LOCF)
main_exc_na$PANSS_N_LOCF <- ifelse((main_exc_na$PANSS2_N == -9 & main_exc_na$PANSS3_N == -9 & main_exc_na$PANSS4_N == -9), 
                                   main_exc_na$PANSS1_N, 
                                   main_exc_na$PANSS_N_LOCF)
main_exc_na$PANSS_G_LOCF <- ifelse((main_exc_na$PANSS2_G == -9 & main_exc_na$PANSS3_G == -9 & main_exc_na$PANSS4_G == -9), 
                                   main_exc_na$PANSS1_G, 
                                   main_exc_na$PANSS_G_LOCF)




#main_exc_na_logi <- subset(main_exc_na, !(main_exc_na$PANSS2 == -9 & main_exc_na$PANSS3 == -9 & main_exc_na$PANSS4 == -9)) #剩2107人

#剩余者以末位结转求panss随访最终值
main_exc_na$PANSS_LOCF <- ifelse(main_exc_na$PANSS_LOCF != -9, main_exc_na$PANSS_LOCF,
                                 ifelse(main_exc_na$PANSS3 != -9, main_exc_na$PANSS3,
                                        ifelse(main_exc_na$PANSS2 != -9, main_exc_na$PANSS2, main_exc_na$PANSS_LOCF)))

main_exc_na$PANSS_P_LOCF <- ifelse(main_exc_na$PANSS_P_LOCF != -9, main_exc_na$PANSS_P_LOCF,
                                   ifelse(main_exc_na$PANSS3_P != -9, main_exc_na$PANSS3_P,
                                          ifelse(main_exc_na$PANSS2_P != -9, main_exc_na$PANSS2_P, main_exc_na$PANSS_P_LOCF)))

main_exc_na$PANSS_G_LOCF <- ifelse(main_exc_na$PANSS_G_LOCF != -9, main_exc_na$PANSS_G_LOCF,
                                   ifelse(main_exc_na$PANSS3_G != -9, main_exc_na$PANSS3_G,
                                          ifelse(main_exc_na$PANSS2_G != -9, main_exc_na$PANSS2_G, main_exc_na$PANSS_G_LOCF)))

main_exc_na$PANSS_N_LOCF <- ifelse(main_exc_na$PANSS_N_LOCF != -9, main_exc_na$PANSS_N_LOCF,
                                   ifelse(main_exc_na$PANSS3_N != -9, main_exc_na$PANSS3_N,
                                          ifelse(main_exc_na$PANSS2_N != -9, main_exc_na$PANSS2_N, main_exc_na$PANSS_N_LOCF)))



#计算panss减分率
main_exc_na$PANSS_reduce_rate <- (main_exc_na$PANSS1 - main_exc_na$PANSS_LOCF) / (main_exc_na$PANSS1 - 30) *100
main_exc_na$PANSS_P_reduce_rate <- (main_exc_na$PANSS1_P - main_exc_na$PANSS_P_LOCF) / (main_exc_na$PANSS1_P ) *100
main_exc_na$PANSS_G_reduce_rate <- (main_exc_na$PANSS1_G - main_exc_na$PANSS_G_LOCF) / (main_exc_na$PANSS1_G ) *100
main_exc_na$PANSS_N_reduce_rate <- (main_exc_na$PANSS1_N - main_exc_na$PANSS_N_LOCF) / (main_exc_na$PANSS1_N ) *100

#panss绝对差值
main_exc_na$PANSS_reduce <- main_exc_na$PANSS1 - main_exc_na$PANSS_LOCF
main_exc_na$PANSS_P_reduce <- main_exc_na$PANSS1_P - main_exc_na$PANSS_P_LOCF
main_exc_na$PANSS_G_reduce <- main_exc_na$PANSS1_G - main_exc_na$PANSS_G_LOCF
main_exc_na$PANSS_N_reduce <- main_exc_na$PANSS1_N - main_exc_na$PANSS_N_LOCF


#子量表panss的二分类
main_exc_na$PANSS_reduce_rate_2g <- ifelse(main_exc_na$PANSS_reduce_rate >= 50, 1, 0)
main_exc_na$PANSS_reduce_rate_2g_med <- ifelse(main_exc_na$PANSS_reduce_rate >= median(main_exc_na$PANSS_reduce_rate), 1, 0)
main_exc_na$PANSS_N_reduce_rate_2g <- ifelse(main_exc_na$PANSS_N_reduce_rate >= 50, 1, 0)
main_exc_na$PANSS_N_reduce_rate_2g_med <- ifelse(main_exc_na$PANSS_N_reduce_rate >= median(main_exc_na$PANSS_N_reduce_rate), 1, 0)
main_exc_na$PANSS_P_reduce_rate_2g <- ifelse(main_exc_na$PANSS_P_reduce_rate >= 50, 1, 0)
main_exc_na$PANSS_P_reduce_rate_2g_med <- ifelse(main_exc_na$PANSS_P_reduce_rate >= median(main_exc_na$PANSS_P_reduce_rate), 1, 0)
main_exc_na$PANSS_G_reduce_rate_2g <- ifelse(main_exc_na$PANSS_G_reduce_rate >= 50, 1, 0)
main_exc_na$PANSS_G_reduce_rate_2g_med <- ifelse(main_exc_na$PANSS_G_reduce_rate >= median(main_exc_na$PANSS_G_reduce_rate), 1, 0)

# main_exc_na$PANSS_reduce_rate_2g_20 <- ifelse(main_exc_na$PANSS_reduce_rate >= 0.2, 1, 0)
# main_exc_na$PANSS_N_reduce_rate_2g_20 <- ifelse(main_exc_na$PANSS_N_reduce_rate >= 0.2, 1, 0)
# main_exc_na$PANSS_P_reduce_rate_2g_20 <- ifelse(main_exc_na$PANSS_P_reduce_rate >= 0.2, 1, 0)
# main_exc_na$PANSS_G_reduce_rate_2g_20 <- ifelse(main_exc_na$PANSS_G_reduce_rate >= 0.2, 1, 0)
# 
# main_exc_na$PANSS_reduce_rate_2g_30 <- ifelse(main_exc_na$PANSS_reduce_rate >= 0.3, 1, 0)
# main_exc_na$PANSS_N_reduce_rate_2g_30 <- ifelse(main_exc_na$PANSS_N_reduce_rate >= 0.3, 1, 0)
# main_exc_na$PANSS_P_reduce_rate_2g_30 <- ifelse(main_exc_na$PANSS_P_reduce_rate >= 0.3, 1, 0)
# main_exc_na$PANSS_G_reduce_rate_2g_30 <- ifelse(main_exc_na$PANSS_G_reduce_rate >= 0.3, 1, 0)
# 
# main_exc_na$PANSS_reduce_rate_2g_40 <- ifelse(main_exc_na$PANSS_reduce_rate >= 0.4, 1, 0)
# main_exc_na$PANSS_N_reduce_rate_2g_40 <- ifelse(main_exc_na$PANSS_N_reduce_rate >= 0.4, 1, 0)
# main_exc_na$PANSS_P_reduce_rate_2g_40 <- ifelse(main_exc_na$PANSS_P_reduce_rate >= 0.4, 1, 0)
# main_exc_na$PANSS_G_reduce_rate_2g_40 <- ifelse(main_exc_na$PANSS_G_reduce_rate >= 0.4, 1, 0)
# 
# main_exc_na$PANSS_reduce_rate_2g_60 <- ifelse(main_exc_na$PANSS_reduce_rate >= 0.6, 1, 0)
# main_exc_na$PANSS_N_reduce_rate_2g_60 <- ifelse(main_exc_na$PANSS_N_reduce_rate >= 0.6, 1, 0)
# main_exc_na$PANSS_P_reduce_rate_2g_60 <- ifelse(main_exc_na$PANSS_P_reduce_rate >= 0.6, 1, 0)
# main_exc_na$PANSS_G_reduce_rate_2g_60 <- ifelse(main_exc_na$PANSS_G_reduce_rate >= 0.6, 1, 0)
# 
# main_exc_na$PANSS_reduce_rate_2g_70 <- ifelse(main_exc_na$PANSS_reduce_rate >= 0.7, 1, 0)
# main_exc_na$PANSS_N_reduce_rate_2g_70 <- ifelse(main_exc_na$PANSS_N_reduce_rate >= 0.7, 1, 0)
# main_exc_na$PANSS_P_reduce_rate_2g_70 <- ifelse(main_exc_na$PANSS_P_reduce_rate >= 0.7, 1, 0)
# main_exc_na$PANSS_G_reduce_rate_2g_70 <- ifelse(main_exc_na$PANSS_G_reduce_rate >= 0.7, 1, 0)
# 
# main_exc_na$PANSS_reduce_rate_2g_80 <- ifelse(main_exc_na$PANSS_reduce_rate >= 0.8, 1, 0)
# main_exc_na$PANSS_N_reduce_rate_2g_80 <- ifelse(main_exc_na$PANSS_N_reduce_rate >= 0.8, 1, 0)
# main_exc_na$PANSS_P_reduce_rate_2g_80 <- ifelse(main_exc_na$PANSS_P_reduce_rate >= 0.8, 1, 0)
# main_exc_na$PANSS_G_reduce_rate_2g_80 <- ifelse(main_exc_na$PANSS_G_reduce_rate >= 0.8, 1, 0)










main_exc_na$PANSS1_2g <- ifelse(main_exc_na$PANSS1>=70, 1, 0)
main_exc_na$PANSS_LOCF_2g <- ifelse(main_exc_na$PANSS_LOCF>=70, 1, 0)


main_exc_na <- main_exc_na %>%
  mutate(PANSS_trans = case_when(
    (main_exc_na$PANSS1_2g == 1 & main_exc_na$PANSS_LOCF_2g == 0) ~ "1",
    (main_exc_na$PANSS1_2g == 1 & main_exc_na$PANSS_LOCF_2g == 1 ) ~ "0",
    (main_exc_na$PANSS1_2g == 0 & main_exc_na$PANSS_LOCF_2g == 0 ) ~ "2",
    T~"3"
  ))

# table(main_exc_na$PANSS_trans)
# 0    1    2    3 
# 670 1346  176    8  

#绝大部分是从严重组变到不严重-case


main_exc_na <- main_exc_na %>%
  mutate(PANSS_trans_2g = case_when(
    (main_exc_na$PANSS_trans==1) ~ "1",
    T~"0"
  ))



########
#panss reduce rate-2w 4w
main_exc_na$PANSS_reduce_rate_2w = (main_exc_na$PANSS1 - main_exc_na$PANSS2) / (main_exc_na$PANSS1 - 30) * 100
main_exc_na$PANSS_reduce_rate_4w = (main_exc_na$PANSS1 - main_exc_na$PANSS3) / (main_exc_na$PANSS1 - 30) * 100
main_exc_na$PANSS_reduce_rate_2w_2g = ifelse(main_exc_na$PANSS_reduce_rate_2w >= 50, 1, 0)
main_exc_na$PANSS_reduce_rate_4w_2g = ifelse(main_exc_na$PANSS_reduce_rate_4w >= 50, 1, 0)
main_exc_na$PANSS_reduce_rate_2w_2g_med <- ifelse(main_exc_na$PANSS_reduce_rate_2w >= median(main_exc_na$PANSS_reduce_rate_2w, na.rm = TRUE),1, 0)
main_exc_na$PANSS_reduce_rate_4w_2g_med = ifelse(main_exc_na$PANSS_reduce_rate_4w >= median(main_exc_na$PANSS_reduce_rate_4w), 1, 0)

main_exc_na$PANSS_P_reduce_rate_2w = (main_exc_na$PANSS1_P - main_exc_na$PANSS2_P) / (main_exc_na$PANSS1_P) * 100
main_exc_na$PANSS_P_reduce_rate_4w = (main_exc_na$PANSS1_P - main_exc_na$PANSS3_P) / (main_exc_na$PANSS1_P) * 100
main_exc_na$PANSS_P_reduce_rate_2w_2g = ifelse(main_exc_na$PANSS_P_reduce_rate_2w >= 50, 1, 0)
main_exc_na$PANSS_P_reduce_rate_4w_2g = ifelse(main_exc_na$PANSS_P_reduce_rate_4w >= 50, 1, 0)
main_exc_na$PANSS_P_reduce_rate_2w_2g_med <- ifelse(main_exc_na$PANSS_P_reduce_rate_2w >= median(main_exc_na$PANSS_P_reduce_rate_2w, na.rm = TRUE),1, 0)
main_exc_na$PANSS_P_reduce_rate_4w_2g_med = ifelse(main_exc_na$PANSS_P_reduce_rate_4w >= median(main_exc_na$PANSS_P_reduce_rate_4w), 1, 0)

main_exc_na$PANSS_N_reduce_rate_2w = (main_exc_na$PANSS1_N - main_exc_na$PANSS2_N) / (main_exc_na$PANSS1_N) * 100
main_exc_na$PANSS_N_reduce_rate_4w = (main_exc_na$PANSS1_N - main_exc_na$PANSS3_N) / (main_exc_na$PANSS1_N) * 100
main_exc_na$PANSS_N_reduce_rate_2w_2g = ifelse(main_exc_na$PANSS_N_reduce_rate_2w >= 50, 1, 0)
main_exc_na$PANSS_N_reduce_rate_4w_2g = ifelse(main_exc_na$PANSS_N_reduce_rate_4w >= 50, 1, 0)
main_exc_na$PANSS_N_reduce_rate_2w_2g_med <- ifelse(main_exc_na$PANSS_N_reduce_rate_2w >= median(main_exc_na$PANSS_N_reduce_rate_2w, na.rm = TRUE),1, 0)
main_exc_na$PANSS_N_reduce_rate_4w_2g_med = ifelse(main_exc_na$PANSS_N_reduce_rate_4w >= median(main_exc_na$PANSS_N_reduce_rate_4w), 1, 0)

main_exc_na$PANSS_G_reduce_rate_2w = (main_exc_na$PANSS1_G - main_exc_na$PANSS2_G) / (main_exc_na$PANSS1_G) * 100
main_exc_na$PANSS_G_reduce_rate_4w = (main_exc_na$PANSS1_G - main_exc_na$PANSS3_G) / (main_exc_na$PANSS1_G) * 100
main_exc_na$PANSS_G_reduce_rate_2w_2g = ifelse(main_exc_na$PANSS_G_reduce_rate_2w >= 50, 1, 0)
main_exc_na$PANSS_G_reduce_rate_4w_2g = ifelse(main_exc_na$PANSS_G_reduce_rate_4w >= 50, 1, 0)
main_exc_na$PANSS_G_reduce_rate_2w_2g_med <- ifelse(main_exc_na$PANSS_G_reduce_rate_2w >= median(main_exc_na$PANSS_G_reduce_rate_2w, na.rm = TRUE),1, 0)
main_exc_na$PANSS_G_reduce_rate_4w_2g_med = ifelse(main_exc_na$PANSS_G_reduce_rate_4w >= median(main_exc_na$PANSS_G_reduce_rate_4w), 1, 0)







#########
#五因子减分率

#没有随访信息者，用baseline填补
main_exc_na$PANSS_F5_P_LOCF <- main_exc_na$PANSS4_F5_P
main_exc_na$PANSS_F5_CD_LOCF <- main_exc_na$PANSS4_F5_CD
main_exc_na$PANSS_F5_N_LOCF <- main_exc_na$PANSS4_F5_N
main_exc_na$PANSS_F5_DA_LOCF <- main_exc_na$PANSS4_F5_DA
main_exc_na$PANSS_F5_H_LOCF <- main_exc_na$PANSS4_F5_H


main_exc_na$PANSS_F5_P_LOCF <- ifelse((main_exc_na$PANSS2_F5_P == -9 & main_exc_na$PANSS3_F5_P == -9 & main_exc_na$PANSS4_F5_P == -9), 
                                      main_exc_na$PANSS1_F5_P, 
                                      main_exc_na$PANSS_F5_P_LOCF)
main_exc_na$PANSS_F5_N_LOCF <- ifelse((main_exc_na$PANSS2_F5_N == -9 & main_exc_na$PANSS3_F5_N == -9 & main_exc_na$PANSS4_F5_N == -9), 
                                      main_exc_na$PANSS1_F5_N, 
                                      main_exc_na$PANSS_F5_N_LOCF)
main_exc_na$PANSS_F5_CD_LOCF <- ifelse((main_exc_na$PANSS2_F5_CD == -9 & main_exc_na$PANSS3_F5_CD == -9 & main_exc_na$PANSS4_F5_CD == -9), 
                                       main_exc_na$PANSS1_F5_CD, 
                                       main_exc_na$PANSS_F5_CD_LOCF)

main_exc_na$PANSS_F5_DA_LOCF <- ifelse((main_exc_na$PANSS2_F5_DA == -9 & main_exc_na$PANSS3_F5_DA == -9 & main_exc_na$PANSS4_F5_DA == -9), 
                                       main_exc_na$PANSS1_F5_DA, 
                                       main_exc_na$PANSS_F5_DA_LOCF)

main_exc_na$PANSS_F5_H_LOCF <- ifelse((main_exc_na$PANSS2_F5_H == -9 & main_exc_na$PANSS3_F5_H == -9 & main_exc_na$PANSS4_F5_H == -9), 
                                      main_exc_na$PANSS1_F5_H, 
                                      main_exc_na$PANSS_F5_H_LOCF)





#main_exc_na_logi <- subset(main_exc_na, !(main_exc_na$PANSS2 == -9 & main_exc_na$PANSS3 == -9 & main_exc_na$PANSS4 == -9)) #剩2107人

#剩余者以末位结转求panss随访最终值

main_exc_na$PANSS_F5_P_LOCF <- ifelse(main_exc_na$PANSS_F5_P_LOCF != -9, main_exc_na$PANSS_F5_P_LOCF,
                                      ifelse(main_exc_na$PANSS3_F5_P != -9, main_exc_na$PANSS3_F5_P,
                                             ifelse(main_exc_na$PANSS2_F5_P != -9, main_exc_na$PANSS2_F5_P, main_exc_na$PANSS_F5_P_LOCF)))

main_exc_na$PANSS_F5_CD_LOCF <- ifelse(main_exc_na$PANSS_F5_CD_LOCF != -9, main_exc_na$PANSS_F5_CD_LOCF,
                                       ifelse(main_exc_na$PANSS3_F5_CD != -9, main_exc_na$PANSS3_F5_CD,
                                              ifelse(main_exc_na$PANSS2_F5_CD != -9, main_exc_na$PANSS2_F5_CD, main_exc_na$PANSS_F5_CD_LOCF)))

main_exc_na$PANSS_F5_N_LOCF <- ifelse(main_exc_na$PANSS_F5_N_LOCF != -9, main_exc_na$PANSS_F5_N_LOCF,
                                      ifelse(main_exc_na$PANSS3_F5_N != -9, main_exc_na$PANSS3_F5_N,
                                             ifelse(main_exc_na$PANSS2_F5_N != -9, main_exc_na$PANSS2_F5_N, main_exc_na$PANSS_F5_N_LOCF)))

main_exc_na$PANSS_F5_DA_LOCF <- ifelse(main_exc_na$PANSS_F5_DA_LOCF != -9, main_exc_na$PANSS_F5_DA_LOCF,
                                       ifelse(main_exc_na$PANSS3_F5_DA != -9, main_exc_na$PANSS3_F5_DA,
                                              ifelse(main_exc_na$PANSS2_F5_DA != -9, main_exc_na$PANSS2_F5_DA, main_exc_na$PANSS_F5_DA_LOCF)))

main_exc_na$PANSS_F5_H_LOCF <- ifelse(main_exc_na$PANSS_F5_H_LOCF != -9, main_exc_na$PANSS_F5_H_LOCF,
                                      ifelse(main_exc_na$PANSS3_F5_H != -9, main_exc_na$PANSS3_F5_H,
                                             ifelse(main_exc_na$PANSS2_F5_H != -9, main_exc_na$PANSS2_F5_H, main_exc_na$PANSS_F5_H_LOCF)))



#计算panss减分率
main_exc_na$PANSS_F5_P_reduce_rate <- (main_exc_na$PANSS1_F5_P - main_exc_na$PANSS_F5_P_LOCF) / (main_exc_na$PANSS1_F5_P ) * 100
main_exc_na$PANSS_F5_CD_reduce_rate <- (main_exc_na$PANSS1_F5_CD - main_exc_na$PANSS_F5_CD_LOCF) / (main_exc_na$PANSS1_F5_CD ) * 100
main_exc_na$PANSS_F5_N_reduce_rate <- (main_exc_na$PANSS1_F5_N - main_exc_na$PANSS_F5_N_LOCF) / (main_exc_na$PANSS1_F5_N ) * 100
main_exc_na$PANSS_F5_DA_reduce_rate <- (main_exc_na$PANSS1_F5_DA - main_exc_na$PANSS_F5_DA_LOCF) / (main_exc_na$PANSS1_F5_DA ) * 100
main_exc_na$PANSS_F5_H_reduce_rate <- (main_exc_na$PANSS1_F5_H - main_exc_na$PANSS_F5_H_LOCF) / (main_exc_na$PANSS1_F5_H ) * 100


#panss绝对差值
main_exc_na$PANSS_F5_P_reduce <- main_exc_na$PANSS1_F5_P - main_exc_na$PANSS_F5_P_LOCF
main_exc_na$PANSS_F5_CD_reduce <- main_exc_na$PANSS1_F5_CD - main_exc_na$PANSS_F5_CD_LOCF
main_exc_na$PANSS_F5_N_reduce <- main_exc_na$PANSS1_F5_N - main_exc_na$PANSS_F5_N_LOCF
main_exc_na$PANSS_F5_DA_reduce <- main_exc_na$PANSS1_F5_DA - main_exc_na$PANSS_F5_DA_LOCF
main_exc_na$PANSS_F5_H_reduce <- main_exc_na$PANSS1_F5_H - main_exc_na$PANSS_F5_H_LOCF



#子量表panss的二分类
main_exc_na$PANSS_F5_N_reduce_rate_2g <- ifelse(main_exc_na$PANSS_F5_N_reduce_rate >= 50, 1, 0)
main_exc_na$PANSS_F5_N_reduce_rate_2g_med <- ifelse(main_exc_na$PANSS_F5_N_reduce_rate >= median(main_exc_na$PANSS_F5_N_reduce_rate), 1, 0)
main_exc_na$PANSS_F5_P_reduce_rate_2g <- ifelse(main_exc_na$PANSS_F5_P_reduce_rate >= 50, 1, 0)
main_exc_na$PANSS_F5_P_reduce_rate_2g_med <- ifelse(main_exc_na$PANSS_F5_P_reduce_rate >= median(main_exc_na$PANSS_F5_P_reduce_rate), 1, 0)
main_exc_na$PANSS_F5_CD_reduce_rate_2g <- ifelse(main_exc_na$PANSS_F5_CD_reduce_rate >= 50, 1, 0)
main_exc_na$PANSS_F5_CD_reduce_rate_2g_med <- ifelse(main_exc_na$PANSS_F5_CD_reduce_rate >= median(main_exc_na$PANSS_F5_CD_reduce_rate), 1, 0)
main_exc_na$PANSS_F5_DA_reduce_rate_2g <- ifelse(main_exc_na$PANSS_F5_DA_reduce_rate >= 50, 1, 0)
main_exc_na$PANSS_F5_DA_reduce_rate_2g_med <- ifelse(main_exc_na$PANSS_F5_DA_reduce_rate >= median(main_exc_na$PANSS_F5_DA_reduce_rate), 1, 0)
main_exc_na$PANSS_F5_H_reduce_rate_2g <- ifelse(main_exc_na$PANSS_F5_H_reduce_rate >= 50, 1, 0)
main_exc_na$PANSS_F5_H_reduce_rate_2g_med <- ifelse(main_exc_na$PANSS_F5_H_reduce_rate >= median(main_exc_na$PANSS_F5_H_reduce_rate), 1, 0)







########
#panss reduce rate-2w 4w

main_exc_na$PANSS_F5_P_reduce_rate_2w = (main_exc_na$PANSS1_F5_P - main_exc_na$PANSS2_F5_P) / (main_exc_na$PANSS1_F5_P) * 100
main_exc_na$PANSS_F5_P_reduce_rate_4w = (main_exc_na$PANSS1_F5_P - main_exc_na$PANSS3_F5_P) / (main_exc_na$PANSS1_F5_P) * 100
main_exc_na$PANSS_F5_P_reduce_rate_2w_2g = ifelse(main_exc_na$PANSS_F5_P_reduce_rate_2w >= 50, 1, 0)
main_exc_na$PANSS_F5_P_reduce_rate_4w_2g = ifelse(main_exc_na$PANSS_F5_P_reduce_rate_4w >= 50, 1, 0)
main_exc_na$PANSS_F5_P_reduce_rate_2w_2g_med <- ifelse(main_exc_na$PANSS_F5_P_reduce_rate_2w >= median(main_exc_na$PANSS_F5_P_reduce_rate_2w, na.rm = TRUE),1, 0)
main_exc_na$PANSS_F5_P_reduce_rate_4w_2g_med = ifelse(main_exc_na$PANSS_F5_P_reduce_rate_4w >= median(main_exc_na$PANSS_F5_P_reduce_rate_4w), 1, 0)

main_exc_na$PANSS_F5_N_reduce_rate_2w = (main_exc_na$PANSS1_F5_N - main_exc_na$PANSS2_F5_N) / (main_exc_na$PANSS1_F5_N) * 100
main_exc_na$PANSS_F5_N_reduce_rate_4w = (main_exc_na$PANSS1_F5_N - main_exc_na$PANSS3_F5_N) / (main_exc_na$PANSS1_F5_N) * 100
main_exc_na$PANSS_F5_N_reduce_rate_2w_2g = ifelse(main_exc_na$PANSS_F5_N_reduce_rate_2w >= 50, 1, 0)
main_exc_na$PANSS_F5_N_reduce_rate_4w_2g = ifelse(main_exc_na$PANSS_F5_N_reduce_rate_4w >= 50, 1, 0)
main_exc_na$PANSS_F5_N_reduce_rate_2w_2g_med <- ifelse(main_exc_na$PANSS_F5_N_reduce_rate_2w >= median(main_exc_na$PANSS_F5_N_reduce_rate_2w, na.rm = TRUE),1, 0)
main_exc_na$PANSS_F5_N_reduce_rate_4w_2g_med = ifelse(main_exc_na$PANSS_F5_N_reduce_rate_4w >= median(main_exc_na$PANSS_F5_N_reduce_rate_4w), 1, 0)

main_exc_na$PANSS_F5_G_reduce_rate_2w = (main_exc_na$PANSS1_F5_G - main_exc_na$PANSS2_F5_G) / (main_exc_na$PANSS1_F5_G) * 100
main_exc_na$PANSS_F5_G_reduce_rate_4w = (main_exc_na$PANSS1_F5_G - main_exc_na$PANSS3_F5_G) / (main_exc_na$PANSS1_F5_G) * 100
main_exc_na$PANSS_F5_G_reduce_rate_2w_2g = ifelse(main_exc_na$PANSS_F5_G_reduce_rate_2w >= 50, 1, 0)
main_exc_na$PANSS_F5_G_reduce_rate_4w_2g = ifelse(main_exc_na$PANSS_F5_G_reduce_rate_4w >= 50, 1, 0)
main_exc_na$PANSS_F5_G_reduce_rate_2w_2g_med <- ifelse(main_exc_na$PANSS_F5_G_reduce_rate_2w >= median(main_exc_na$PANSS_F5_G_reduce_rate_2w, na.rm = TRUE),1, 0)
main_exc_na$PANSS_F5_G_reduce_rate_4w_2g_med = ifelse(main_exc_na$PANSS_F5_G_reduce_rate_4w >= median(main_exc_na$PANSS_F5_G_reduce_rate_4w), 1, 0)










#变量编码
table(main_exc_na$sex)
summary(main_exc_na$age)
table(main_exc_na$centers)

names(main_exc_na)[names(main_exc_na) == "p041011"] <- "edu"
names(main_exc_na)[names(main_exc_na) == "p041021"] <- "age_first_episode"
names(main_exc_na)[names(main_exc_na) == "p041024"] <- "first_episode_2g"


main_exc_na <- main_exc_na %>% 
  mutate(edu_4g = case_when(
    (main_exc_na$edu==7 | main_exc_na$edu==8) ~"primary_or_lower",
    main_exc_na$edu==6 ~"middle",
    main_exc_na$edu==5 ~"high",
    (main_exc_na$edu==4 | main_exc_na$edu==3 | 
       main_exc_na$edu==2 | main_exc_na$edu==1)~"college_or_higher",
    T~"NA"
  ))


# 
# main_exc_na <- main_exc_na %>%
#   mutate(qtc0801_2g = case_when(
#     (main_exc_na$qtc0801 > 450 & main_exc_na$sex == 1) ~ "1",
#     (main_exc_na$qtc0801 > 470 & main_exc_na$sex == 2) ~ "1",
#     T~"0"
#   ))


main_exc_na$PANSS1_2g <- ifelse(main_exc_na$PANSS1>=70, 1, 0)

#对血脂指标做标化
main_exc_na$tc1_z <- scale(main_exc_na$tc1)
main_exc_na$tg1_z <- scale(main_exc_na$tg1)
main_exc_na$hdl1_z <- scale(main_exc_na$hdl1)
main_exc_na$ldl1_z <- scale(main_exc_na$ldl1)
main_exc_na$glu1_z <- scale(main_exc_na$glu1)

main_exc_na$marry_2g <- ifelse(main_exc_na$marry=="2",1,0)


##########
#定义心血管危险因素者
library(dplyr)
main_exc_na <- main_exc_na %>%
  mutate(lipid_dys = ifelse(
    ldl1 >= 4.1 | 
      tg1 >= 1.7 | 
      (hdl1 < 1.04 & sex == 1) | 
      (hdl1 < 1.29 & sex == 2) | 
      (fw0601 >= 90 & sex == 1) | 
      (fw0601 >= 85 & sex == 2) | 
      glu1 >= 6.1 | 
      sbp0601 >= 130 | 
      dbp0601 >= 85, 1, 0
  ))

#没有cvd危险因素 且未服用有降压作用的抗精神病药物者
#排除基线服用降压药 降糖药 抗cvd药物者
main_exc_na$nocvrf_noantihyp <- ifelse(main_exc_na$lipid_dys==0 & main_exc_na$medication != 6 & 
                                         main_exc_na$medication != 7 & 
                                         main_exc_na$dn!="120081" & main_exc_na$dn!="120767" & main_exc_na$dn!="510227" & 
                                         main_exc_na$dn!="210043" & main_exc_na$dn!="530159", 
                                       1, 0)

main_exc_na$abdo_adiposity <- ifelse((main_exc_na$fw0601 >= 90 & main_exc_na$sex == 1) | 
                                       (main_exc_na$fw0601 >= 85 & main_exc_na$sex == 2), 1, 0)


main_exc_na$hyp <- ifelse(main_exc_na$sbp0601 >= 130 | main_exc_na$dbp0601 >= 85, 1, 0)

main_exc_na$glu_dys_2g <- ifelse(main_exc_na$glu1 >= 6.1, 1, 0)

main_exc_na$lipid_dys_2g <- ifelse(main_exc_na$ldl1 >= 4.1 | 
                                     main_exc_na$tg1 >= 1.7 | 
                                     (main_exc_na$hdl1 < 1.04 & main_exc_na$sex == 1) | 
                                     (main_exc_na$hdl1 < 1.29 & main_exc_na$sex == 2), 1, 0)


# 210043-入院后三天确诊糖尿病
# 530159 糖尿病

main_exc_na$fw0601_mets <- ifelse(((main_exc_na$fw0601 >= 90 & main_exc_na$sex == 1) | 
                                     (main_exc_na$fw0601 >= 85 & main_exc_na$sex == 2)),1, 0)
main_exc_na$glu1_mets <- ifelse((main_exc_na$glu1 >= 6.1 | main_exc_na$dn=="210043" | 
                                   main_exc_na$dn=="530159"), 1, 0)
main_exc_na$sbp0601_mets <- ifelse((main_exc_na$sbp0601 >= 130 | main_exc_na$dbp0601 >= 85), 1, 0)
main_exc_na$tg1_mets <- ifelse(main_exc_na$tg1 >= 1.7, 1, 0)
main_exc_na$hdl1_mets <- ifelse(((main_exc_na$hdl1 < 1.04 & main_exc_na$sex == 1) | 
                                   (main_exc_na$hdl1 < 1.29 & main_exc_na$sex == 2)),1,0)
main_exc_na$mets_score_bl <- main_exc_na$fw0601_mets + main_exc_na$glu1_mets + 
  main_exc_na$sbp0601_mets + main_exc_na$tg1_mets + main_exc_na$hdl1_mets

main_exc_na$mets_bl <- ifelse(main_exc_na$mets_score_bl >= 3, 1, 0)



main_exc_na <- main_exc_na %>%
  mutate(cvrh = ifelse(
    ldl1 >= 4.1 | mets_bl==1 | (hba1>= 5.7 & !is.na(hba1)), 1, 0
  ))


fwrite(main_exc_na, "/Users/zhuyunqing/Desktop/study/PhD/data/863/capoc_bl_fl_pheno_forpanss.txt", row.names = FALSE, sep = "\t")




#######
#merge with grs


#merge grs数据
grs <- fread("/Users/zhuyunqing/Desktop/study/drugMR_lipid/data/onesamplemr/CAPOC.cis.capoc_cases.grs.txt",
             h=T)
grs$FID <- as.character(grs$FID)
#grs <- grs %>% select(-grs_DPP4_uw_ukb, -grs_DPP4_w_ukb)
main_exc_na$dn <- as.character(main_exc_na$dn)

m <- merge(main_exc_na, grs, by.x="dn", by.y = "FID", all=F)
#m <- subset(m, IID!="")


fwrite(m, "/Users/zhuyunqing/Desktop/study/drugMR_lipid/data/onesamplemr/CAPOC.cis.capoc_cases.grs.pheno.new.txt",
       sep = "\t")






