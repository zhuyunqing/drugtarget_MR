###########
#weighted from gnhs
##########

library(data.table)
library(dplyr)
library(tidyr)

setwd("/Users/zhuyunqing/Desktop/study/drugMR_lipid/data/onesamplemr")
dat_gnhs <- fread("CAPOC.cis.capoc_cases.snp.raw", h=T)





############
#lipid glu sbp grs的snp
dat_lipid_glu_sbp <- fread("CAPOC.cis.lipid_glu_sbp.snp.txt", h=T)

var_names <- names(dat_lipid_glu_sbp)
core_names <- gsub("_[ACGT]+$", "", var_names)  # 去掉尾部等位基因信息
first_occurrence <- !duplicated(core_names)
keep_cols <- var_names[first_occurrence]
dat_clean <- dat_lipid_glu_sbp[, ..keep_cols]





#MERGE
dat <- merge(dat_gnhs, dat_clean, by="FID", all=F)




##############
#根据与血脂的关联，将其调成降脂的方向


#################
#weighted-grs weight来自gnhs


dat$grs_apoc3_w_gnhs <- -1 * scale((dat$`rs651821_T(/C)` * (-0.106423))/1)



##############
#构建以血脂为weight的血脂药靶GRS
commands <- c("dat$rs1236460_G = 2 - dat$rs1236460_A",
              "dat$rs17195054_A = 2 - dat$rs17195054_G",
              "dat$rs12233358_G = 2 - dat$rs12233358_C",
              "dat$rs17195054_A = 2 - dat$rs17195054_G",
              "dat$rs58561056_T = 2 - dat$rs58561056_G",
              "dat$rs12633551_T = 2 - dat$rs12633551_C",
              "dat$rs75157957_T = 2 - dat$rs75157957_C",
              "dat$rs12485478_G = 2 - dat$rs12485478_A",
              "dat$rs7628899_A = 2 - dat$rs7628899_T",
              "dat$rs6457713_T = 2 - dat$rs6457713_C",
              "dat$rs10177492_T = 2 - dat$rs10177492_C",
              "dat$rs35131127_C = 2 - dat$rs35131127_T",
              "dat$rs512535_T = 2 - dat$rs512535_C",
              "dat$rs57825321_A = 2 - dat$rs57825321_T",
              "dat$rs12992267_T = 2 - dat$rs12992267_C",
              "dat$rs77258720_A = 2 - dat$rs77258720_G",
              "dat$rs57825321_A = 2 - dat$rs57825321_T",
              "dat$rs6706783_A = 2 - dat$rs6706783_T",
              "dat$rs10177492_T = 2 - dat$rs10177492_C",
              "dat$rs17041688_A = 2 - dat$rs17041688_G",
              "dat$rs1998848_A = 2 - dat$rs1998848_G",
              "dat$rs2337382_A = 2 - dat$rs2337382_G",
              "dat$rs312026_C = 2 - dat$rs312026_G",
              "dat$rs57825321_A = 2 - dat$rs57825321_T",
              "dat$rs934197_A = 2 - dat$rs934197_G",
              "dat$rs16872670_A = 2 - dat$rs16872670_G",
              "dat$rs4604177_C = 2 - dat$rs4604177_T",
              "dat$rs4703676_C = 2 - dat$rs4703676_G",
              "dat$rs3797580_A = 2 - dat$rs3797580_G",
              "dat$rs6480685_T = 2 - dat$rs6480685_C",
              "dat$rs4604177_C = 2 - dat$rs4604177_T",
              "dat$rs4703676_C = 2 - dat$rs4703676_G",
              "dat$rs5744661_T = 2 - dat$rs5744661_C",
              "dat$rs1798192_T = 2 - dat$rs1798192_G",
              "dat$rs643285_G = 2 - dat$rs643285_A",
              "dat$rs1798192_T = 2 - dat$rs1798192_G",
              "dat$rs643285_G = 2 - dat$rs643285_A",
              "dat$rs6503531_G = 2 - dat$rs6503531_C",
              "dat$rs709591_A = 2 - dat$rs709591_T",
              "dat$rs1883025_T = 2 - dat$rs1883025_C",
              "dat$rs1800977_A = 2 - dat$rs1800977_G",
              "dat$rs186824556_A = 2 - dat$rs186824556_G",
              "dat$rs2472474_C = 2 - dat$rs2472474_T",
              "dat$rs2472508_A = 2 - dat$rs2472508_G",
              "dat$rs2740480_G = 2 - dat$rs2740480_A",
              "dat$rs2740487_A = 2 - dat$rs2740487_G",
              "dat$rs2740488_C = 2 - dat$rs2740488_A",
              "dat$rs4149307_C = 2 - dat$rs4149307_T",
              "dat$rs4743758_T = 2 - dat$rs4743758_C",
              "dat$rs7024300_T = 2 - dat$rs7024300_C",
              "dat$rs73664373_C = 2 - dat$rs73664373_G",
              "dat$rs10761090_C = 2 - dat$rs10761090_G",
              "dat$rs2515613_A = 2 - dat$rs2515613_G",
              "dat$rs2740480_G = 2 - dat$rs2740480_A",
              "dat$rs2740487_A = 2 - dat$rs2740487_G",
              "dat$rs2740488_C = 2 - dat$rs2740488_A",
              "dat$rs4743762_C = 2 - dat$rs4743762_A",
              "dat$rs73664373_C = 2 - dat$rs73664373_G",
              "dat$rs2306986_C = 2 - dat$rs2306986_G",
              "dat$rs1115721_G = 2 - dat$rs1115721_T",
              "dat$rs2306986_C = 2 - dat$rs2306986_G",
              "dat$rs10888897_T = 2 - dat$rs10888897_C",
              "dat$rs17111652_T = 2 - dat$rs17111652_C",
              "dat$rs2479404_A = 2 - dat$rs2479404_C",
              "dat$rs34232196_T = 2 - dat$rs34232196_C",
              "dat$rs565436_G = 2 - dat$rs565436_A",
              "dat$rs6663252_C = 2 - dat$rs6663252_T",
              "dat$rs2479404_A = 2 - dat$rs2479404_C",
              "dat$rs28615248_C = 2 - dat$rs28615248_T",
              "dat$rs34232196_T = 2 - dat$rs34232196_C",
              "dat$rs534473_G = 2 - dat$rs534473_T",
              "dat$rs565436_G = 2 - dat$rs565436_A",
              "dat$rs6663252_C = 2 - dat$rs6663252_T",
              "dat$rs67171713_C = 2 - dat$rs67171713_T",
              "dat$rs77498878_A = 2 - dat$rs77498878_G",
              "dat$rs6102358_T = 2 - dat$rs6102358_C",
              "dat$rs113879410_C = 2 - dat$rs113879410_G",
              "dat$rs4812493_C = 2 - dat$rs4812493_T",
              "dat$rs2539981_T = 2 - dat$rs2539981_C",
              "dat$rs4495740_G = 2 - dat$rs4495740_T",
              "dat$rs80148997_T = 2 - dat$rs80148997_C",
              "dat$rs10889348_T = 2 - dat$rs10889348_A",
              "dat$rs72807578_A = 2 - dat$rs72807578_G",
              "dat$rs2737265_G = 2 - dat$rs2737265_A",
              "dat$rs111884008_C = 2 - dat$rs111884008_G",
              "dat$rs116852373_A = 2 - dat$rs116852373_G",
              "dat$rs1263172_A = 2 - dat$rs1263172_G",
              "dat$rs17519093_A = 2 - dat$rs17519093_G",
              "dat$rs2737212_C = 2 - dat$rs2737212_T",
              "dat$rs4938309_T = 2 - dat$rs4938309_C",
              "dat$rs5094_A = 2 - dat$rs5094_G",
              "dat$rs651821_C = 2 - dat$rs651821_T",
              "dat$rs9297543_C = 2 - dat$rs9297543_T",
              "dat$rs11216190_C = 2 - dat$rs11216190_T",
              "dat$rs116852373_A = 2 - dat$rs116852373_G",
              "dat$rs139733873_G = 2 - dat$rs139733873_A",
              "dat$rs17519093_A = 2 - dat$rs17519093_G",
              "dat$rs2849169_A = 2 - dat$rs2849169_C",
              "dat$rs4938309_T = 2 - dat$rs4938309_C",
              "dat$rs518181_A = 2 - dat$rs518181_C",
              "dat$rs583219_C = 2 - dat$rs583219_T",
              "dat$rs595049_C = 2 - dat$rs595049_A",
              "dat$rs651821_C = 2 - dat$rs651821_T",
              "dat$rs79210031_T = 2 - dat$rs79210031_A",
              "dat$rs80022746_A = 2 - dat$rs80022746_G",
              "dat$rs888245_G = 2 - dat$rs888245_C",
              "dat$rs2239013_T = 2 - dat$rs2239013_C",
              "dat$rs2737265_G = 2 - dat$rs2737265_A",
              "dat$rs2849169_A = 2 - dat$rs2849169_C",
              "dat$rs595049_C = 2 - dat$rs595049_A",
              "dat$rs651821_C = 2 - dat$rs651821_T",
              "dat$rs918143_T = 2 - dat$rs918143_C",
              "dat$rs1138429_T = 2 - dat$rs1138429_A",
              "dat$rs12934632_T = 2 - dat$rs12934632_C",
              "dat$rs148808285_T = 2 - dat$rs148808285_C",
              "dat$rs148859408_C = 2 - dat$rs148859408_G",
              "dat$rs17369578_A = 2 - dat$rs17369578_G",
              "dat$rs185307787_C = 2 - dat$rs185307787_A",
              "dat$rs188149780_C = 2 - dat$rs188149780_G",
              "dat$rs201825234_A = 2 - dat$rs201825234_G",
              "dat$rs247617_A = 2 - dat$rs247617_C",
              "dat$rs28495885_T = 2 - dat$rs28495885_C",
              "dat$rs289716_T = 2 - dat$rs289716_A",
              "dat$rs289752_G = 2 - dat$rs289752_A",
              "dat$rs36229485_T = 2 - dat$rs36229485_C",
              "dat$rs37481_C = 2 - dat$rs37481_T",
              "dat$rs4783965_A = 2 - dat$rs4783965_G",
              "dat$rs62035923_A = 2 - dat$rs62035923_G",
              "dat$rs7187275_A = 2 - dat$rs7187275_C",
              "dat$rs79582366_A = 2 - dat$rs79582366_C",
              "dat$rs8059431_T = 2 - dat$rs8059431_C",
              "dat$rs9925265_A = 2 - dat$rs9925265_G",
              "dat$rs9926440_C = 2 - dat$rs9926440_G",
              "dat$rs821840_G = 2 - dat$rs821840_A",
              "dat$rs17231506_T = 2 - dat$rs17231506_C",
              "dat$rs188149780_C = 2 - dat$rs188149780_G",
              "dat$rs201825234_A = 2 - dat$rs201825234_G",
              "dat$rs36229485_T = 2 - dat$rs36229485_C",
              "dat$rs9926440_C = 2 - dat$rs9926440_G",
              "dat$rs35131127_C = 2 - dat$rs35131127_T",
              "dat$rs512535_T = 2 - dat$rs512535_C",
              "dat$rs57825321_A = 2 - dat$rs57825321_T",
              "dat$rs12992267_T = 2 - dat$rs12992267_C",
              "dat$rs34374768_G = 2 - dat$rs34374768_A",
              "dat$rs77258720_A = 2 - dat$rs77258720_G",
              "dat$rs57825321_A = 2 - dat$rs57825321_T",
              "dat$rs6706783_A = 2 - dat$rs6706783_T",
              "dat$rs144090867_T = 2 - dat$rs144090867_C",
              "dat$rs2337382_A = 2 - dat$rs2337382_G",
              "dat$rs57825321_A = 2 - dat$rs57825321_T",
              "dat$rs934197_A = 2 - dat$rs934197_G",
              "dat$rs1019194_C = 2 - dat$rs1019194_T",
              "dat$rs11085760_C = 2 - dat$rs11085760_T",
              "dat$rs17001200_C = 2 - dat$rs17001200_A",
              "dat$rs17243032_C = 2 - dat$rs17243032_A",
              "dat$rs17699089_G = 2 - dat$rs17699089_A",
              "dat$rs2738464_G = 2 - dat$rs2738464_C",
              "dat$rs442363_A = 2 - dat$rs442363_G",
              "dat$rs5927_A = 2 - dat$rs5927_G",
              "dat$rs9789328_T = 2 - dat$rs9789328_G",
              "dat$rs17699089_G = 2 - dat$rs17699089_A",
              "dat$rs4804154_T = 2 - dat$rs4804154_C",
              "dat$rs1019194_C = 2 - dat$rs1019194_T",
              "dat$rs11085760_C = 2 - dat$rs11085760_T",
              "dat$rs116886614_G = 2 - dat$rs116886614_A",
              "dat$rs149152494_T = 2 - dat$rs149152494_C",
              "dat$rs17243032_C = 2 - dat$rs17243032_A",
              "dat$rs2738446_G = 2 - dat$rs2738446_C",
              "dat$rs2738464_G = 2 - dat$rs2738464_C",
              "dat$rs3810308_C = 2 - dat$rs3810308_T",
              "dat$rs442363_A = 2 - dat$rs442363_G",
              "dat$rs5927_A = 2 - dat$rs5927_G",
              "dat$rs7250652_A = 2 - dat$rs7250652_G",
              "dat$rs8107532_T = 2 - dat$rs8107532_C",
              "dat$rs11998681_G = 2 - dat$rs11998681_T",
              "dat$rs1561748_C = 2 - dat$rs1561748_G",
              "dat$rs343_A = 2 - dat$rs343_C",
              "dat$rs4922108_C = 2 - dat$rs4922108_T",
              "dat$rs4922113_A = 2 - dat$rs4922113_G",
              "dat$rs77069344_G = 2 - dat$rs77069344_T",
              "dat$rs78299715_G = 2 - dat$rs78299715_C",
              "dat$rs9644636_G = 2 - dat$rs9644636_T",
              "dat$rs117038725_A = 2 - dat$rs117038725_G",
              "dat$rs1441778_C = 2 - dat$rs1441778_T",
              "dat$rs1561748_C = 2 - dat$rs1561748_G",
              "dat$rs343_A = 2 - dat$rs343_C",
              "dat$rs4922113_A = 2 - dat$rs4922113_G",
              "dat$rs73004966_T = 2 - dat$rs73004966_C",
              "dat$rs74855321_T = 2 - dat$rs74855321_C",
              "dat$rs78299715_G = 2 - dat$rs78299715_C",
              "dat$rs7770628_C = 2 - dat$rs7770628_T",
              "dat$rs9355296_A = 2 - dat$rs9355296_G",
              "dat$rs2457550_T = 2 - dat$rs2457550_C",
              "dat$rs7770628_C = 2 - dat$rs7770628_T",
              "dat$rs57176252_A = 2 - dat$rs57176252_C",
              "dat$rs6453131_T = 2 - dat$rs6453131_G",
              "dat$rs73763094_T = 2 - dat$rs73763094_C",
              "dat$rs7717415_G = 2 - dat$rs7717415_A",
              "dat$rs11816594_A = 2 - dat$rs11816594_G",
              "dat$rs2303152_A = 2 - dat$rs2303152_G",
              "dat$rs57176252_A = 2 - dat$rs57176252_C",
              "dat$rs6453131_T = 2 - dat$rs6453131_G",
              "dat$rs7717415_G = 2 - dat$rs7717415_A",
              "dat$rs80223115_A = 2 - dat$rs80223115_G",
              "dat$rs6545966_G = 2 - dat$rs6545966_A",
              "dat$rs1334806_A = 2 - dat$rs1334806_G",
              "dat$rs7528706_G = 2 - dat$rs7528706_A",
              "dat$rs11125941_A = 2 - dat$rs11125941_C",
              "dat$rs13002939_C = 2 - dat$rs13002939_A",
              "dat$rs4341893_A = 2 - dat$rs4341893_G",
              "dat$rs4426495_C = 2 - dat$rs4426495_T",
              "dat$rs6487176_A = 2 - dat$rs6487176_T",
              "dat$rs10201379_C = 2 - dat$rs10201379_A",
              "dat$rs74358761_A = 2 - dat$rs74358761_G",
              "dat$rs144090867_T = 2 - dat$rs144090867_C",
              "dat$rs4341893_A = 2 - dat$rs4341893_G",
              "dat$rs75509845_A = 2 - dat$rs75509845_G",
              "dat$rs12628618_T = 2 - dat$rs12628618_G",
              "dat$rs111859808_A = 2 - dat$rs111859808_G",
              "dat$rs6104410_G = 2 - dat$rs6104410_A",
              "dat$rs16990971_G = 2 - dat$rs16990971_A",
              "dat$rs1236460_A = 2 - dat$rs1236460_G",
              "dat$rs17195054_G = 2 - dat$rs17195054_A",
              "dat$rs12233358_C = 2 - dat$rs12233358_G",
              "dat$rs17195054_G = 2 - dat$rs17195054_A",
              "dat$rs58561056_G = 2 - dat$rs58561056_T",
              "dat$rs12633551_C = 2 - dat$rs12633551_T",
              "dat$rs75157957_C = 2 - dat$rs75157957_T",
              "dat$rs12485478_A = 2 - dat$rs12485478_G",
              "dat$rs7628899_T = 2 - dat$rs7628899_A",
              "dat$rs6457713_C = 2 - dat$rs6457713_T",
              "dat$rs10177492_C = 2 - dat$rs10177492_T",
              "dat$rs35131127_T = 2 - dat$rs35131127_C",
              "dat$rs512535_C = 2 - dat$rs512535_T",
              "dat$rs57825321_T = 2 - dat$rs57825321_A",
              "dat$rs12992267_C = 2 - dat$rs12992267_T",
              "dat$rs77258720_G = 2 - dat$rs77258720_A",
              "dat$rs57825321_T = 2 - dat$rs57825321_A",
              "dat$rs6706783_T = 2 - dat$rs6706783_A",
              "dat$rs10177492_C = 2 - dat$rs10177492_T",
              "dat$rs17041688_G = 2 - dat$rs17041688_A",
              "dat$rs1998848_G = 2 - dat$rs1998848_A",
              "dat$rs2337382_G = 2 - dat$rs2337382_A",
              "dat$rs312026_G = 2 - dat$rs312026_C",
              "dat$rs57825321_T = 2 - dat$rs57825321_A",
              "dat$rs934197_G = 2 - dat$rs934197_A",
              "dat$rs16872670_G = 2 - dat$rs16872670_A",
              "dat$rs4604177_T = 2 - dat$rs4604177_C",
              "dat$rs4703676_G = 2 - dat$rs4703676_C",
              "dat$rs3797580_G = 2 - dat$rs3797580_A",
              "dat$rs6480685_C = 2 - dat$rs6480685_T",
              "dat$rs4604177_T = 2 - dat$rs4604177_C",
              "dat$rs4703676_G = 2 - dat$rs4703676_C",
              "dat$rs5744661_C = 2 - dat$rs5744661_T",
              "dat$rs1798192_G = 2 - dat$rs1798192_T",
              "dat$rs643285_A = 2 - dat$rs643285_G",
              "dat$rs1798192_G = 2 - dat$rs1798192_T",
              "dat$rs643285_A = 2 - dat$rs643285_G",
              "dat$rs6503531_C = 2 - dat$rs6503531_G",
              "dat$rs709591_T = 2 - dat$rs709591_A",
              "dat$rs1883025_C = 2 - dat$rs1883025_T",
              "dat$rs1800977_G = 2 - dat$rs1800977_A",
              "dat$rs186824556_G = 2 - dat$rs186824556_A",
              "dat$rs2472474_T = 2 - dat$rs2472474_C",
              "dat$rs2472508_G = 2 - dat$rs2472508_A",
              "dat$rs2740480_A = 2 - dat$rs2740480_G",
              "dat$rs2740487_G = 2 - dat$rs2740487_A",
              "dat$rs2740488_A = 2 - dat$rs2740488_C",
              "dat$rs4149307_T = 2 - dat$rs4149307_C",
              "dat$rs4743758_C = 2 - dat$rs4743758_T",
              "dat$rs7024300_C = 2 - dat$rs7024300_T",
              "dat$rs73664373_G = 2 - dat$rs73664373_C",
              "dat$rs10761090_G = 2 - dat$rs10761090_C",
              "dat$rs2515613_G = 2 - dat$rs2515613_A",
              "dat$rs2740480_A = 2 - dat$rs2740480_G",
              "dat$rs2740487_G = 2 - dat$rs2740487_A",
              "dat$rs2740488_A = 2 - dat$rs2740488_C",
              "dat$rs4743762_A = 2 - dat$rs4743762_C",
              "dat$rs73664373_G = 2 - dat$rs73664373_C",
              "dat$rs2306986_G = 2 - dat$rs2306986_C",
              "dat$rs1115721_T = 2 - dat$rs1115721_G",
              "dat$rs2306986_G = 2 - dat$rs2306986_C",
              "dat$rs10888897_C = 2 - dat$rs10888897_T",
              "dat$rs17111652_C = 2 - dat$rs17111652_T",
              "dat$rs2479404_C = 2 - dat$rs2479404_A",
              "dat$rs34232196_C = 2 - dat$rs34232196_T",
              "dat$rs565436_A = 2 - dat$rs565436_G",
              "dat$rs6663252_T = 2 - dat$rs6663252_C",
              "dat$rs2479404_C = 2 - dat$rs2479404_A",
              "dat$rs28615248_T = 2 - dat$rs28615248_C",
              "dat$rs34232196_C = 2 - dat$rs34232196_T",
              "dat$rs534473_T = 2 - dat$rs534473_G",
              "dat$rs565436_A = 2 - dat$rs565436_G",
              "dat$rs6663252_T = 2 - dat$rs6663252_C",
              "dat$rs67171713_T = 2 - dat$rs67171713_C",
              "dat$rs77498878_G = 2 - dat$rs77498878_A",
              "dat$rs6102358_C = 2 - dat$rs6102358_T",
              "dat$rs113879410_G = 2 - dat$rs113879410_C",
              "dat$rs4812493_T = 2 - dat$rs4812493_C",
              "dat$rs2539981_C = 2 - dat$rs2539981_T",
              "dat$rs4495740_T = 2 - dat$rs4495740_G",
              "dat$rs80148997_C = 2 - dat$rs80148997_T",
              "dat$rs10889348_A = 2 - dat$rs10889348_T",
              "dat$rs72807578_G = 2 - dat$rs72807578_A",
              "dat$rs2737265_A = 2 - dat$rs2737265_G",
              "dat$rs111884008_G = 2 - dat$rs111884008_C",
              "dat$rs116852373_G = 2 - dat$rs116852373_A",
              "dat$rs1263172_G = 2 - dat$rs1263172_A",
              "dat$rs17519093_G = 2 - dat$rs17519093_A",
              "dat$rs2737212_T = 2 - dat$rs2737212_C",
              "dat$rs4938309_C = 2 - dat$rs4938309_T",
              "dat$rs5094_G = 2 - dat$rs5094_A",
              "dat$rs651821_T = 2 - dat$rs651821_C",
              "dat$rs9297543_T = 2 - dat$rs9297543_C",
              "dat$rs11216190_T = 2 - dat$rs11216190_C",
              "dat$rs116852373_G = 2 - dat$rs116852373_A",
              "dat$rs139733873_A = 2 - dat$rs139733873_G",
              "dat$rs17519093_G = 2 - dat$rs17519093_A",
              "dat$rs2849169_C = 2 - dat$rs2849169_A",
              "dat$rs4938309_C = 2 - dat$rs4938309_T",
              "dat$rs518181_C = 2 - dat$rs518181_A",
              "dat$rs583219_T = 2 - dat$rs583219_C",
              "dat$rs595049_A = 2 - dat$rs595049_C",
              "dat$rs651821_T = 2 - dat$rs651821_C",
              "dat$rs79210031_A = 2 - dat$rs79210031_T",
              "dat$rs80022746_G = 2 - dat$rs80022746_A",
              "dat$rs888245_C = 2 - dat$rs888245_G",
              "dat$rs2239013_C = 2 - dat$rs2239013_T",
              "dat$rs2737265_A = 2 - dat$rs2737265_G",
              "dat$rs2849169_C = 2 - dat$rs2849169_A",
              "dat$rs595049_A = 2 - dat$rs595049_C",
              "dat$rs651821_T = 2 - dat$rs651821_C",
              "dat$rs918143_C = 2 - dat$rs918143_T",
              "dat$rs1138429_A = 2 - dat$rs1138429_T",
              "dat$rs12934632_C = 2 - dat$rs12934632_T",
              "dat$rs148808285_C = 2 - dat$rs148808285_T",
              "dat$rs148859408_G = 2 - dat$rs148859408_C",
              "dat$rs17369578_G = 2 - dat$rs17369578_A",
              "dat$rs185307787_A = 2 - dat$rs185307787_C",
              "dat$rs188149780_G = 2 - dat$rs188149780_C",
              "dat$rs201825234_G = 2 - dat$rs201825234_A",
              "dat$rs247617_C = 2 - dat$rs247617_A",
              "dat$rs28495885_C = 2 - dat$rs28495885_T",
              "dat$rs289716_A = 2 - dat$rs289716_T",
              "dat$rs289752_A = 2 - dat$rs289752_G",
              "dat$rs36229485_C = 2 - dat$rs36229485_T",
              "dat$rs37481_T = 2 - dat$rs37481_C",
              "dat$rs4783965_G = 2 - dat$rs4783965_A",
              "dat$rs62035923_G = 2 - dat$rs62035923_A",
              "dat$rs7187275_C = 2 - dat$rs7187275_A",
              "dat$rs79582366_C = 2 - dat$rs79582366_A",
              "dat$rs8059431_C = 2 - dat$rs8059431_T",
              "dat$rs9925265_G = 2 - dat$rs9925265_A",
              "dat$rs9926440_G = 2 - dat$rs9926440_C",
              "dat$rs821840_A = 2 - dat$rs821840_G",
              "dat$rs17231506_C = 2 - dat$rs17231506_T",
              "dat$rs188149780_G = 2 - dat$rs188149780_C",
              "dat$rs201825234_G = 2 - dat$rs201825234_A",
              "dat$rs36229485_C = 2 - dat$rs36229485_T",
              "dat$rs9926440_G = 2 - dat$rs9926440_C",
              "dat$rs35131127_T = 2 - dat$rs35131127_C",
              "dat$rs512535_C = 2 - dat$rs512535_T",
              "dat$rs57825321_T = 2 - dat$rs57825321_A",
              "dat$rs12992267_C = 2 - dat$rs12992267_T",
              "dat$rs34374768_A = 2 - dat$rs34374768_G",
              "dat$rs77258720_G = 2 - dat$rs77258720_A",
              "dat$rs57825321_T = 2 - dat$rs57825321_A",
              "dat$rs6706783_T = 2 - dat$rs6706783_A",
              "dat$rs144090867_C = 2 - dat$rs144090867_T",
              "dat$rs2337382_G = 2 - dat$rs2337382_A",
              "dat$rs57825321_T = 2 - dat$rs57825321_A",
              "dat$rs934197_G = 2 - dat$rs934197_A",
              "dat$rs1019194_T = 2 - dat$rs1019194_C",
              "dat$rs11085760_T = 2 - dat$rs11085760_C",
              "dat$rs17001200_A = 2 - dat$rs17001200_C",
              "dat$rs17243032_A = 2 - dat$rs17243032_C",
              "dat$rs17699089_A = 2 - dat$rs17699089_G",
              "dat$rs2738464_C = 2 - dat$rs2738464_G",
              "dat$rs442363_G = 2 - dat$rs442363_A",
              "dat$rs5927_G = 2 - dat$rs5927_A",
              "dat$rs9789328_G = 2 - dat$rs9789328_T",
              "dat$rs17699089_A = 2 - dat$rs17699089_G",
              "dat$rs4804154_C = 2 - dat$rs4804154_T",
              "dat$rs1019194_T = 2 - dat$rs1019194_C",
              "dat$rs11085760_T = 2 - dat$rs11085760_C",
              "dat$rs116886614_A = 2 - dat$rs116886614_G",
              "dat$rs149152494_C = 2 - dat$rs149152494_T",
              "dat$rs17243032_A = 2 - dat$rs17243032_C",
              "dat$rs2738446_C = 2 - dat$rs2738446_G",
              "dat$rs2738464_C = 2 - dat$rs2738464_G",
              "dat$rs3810308_T = 2 - dat$rs3810308_C",
              "dat$rs442363_G = 2 - dat$rs442363_A",
              "dat$rs5927_G = 2 - dat$rs5927_A",
              "dat$rs7250652_G = 2 - dat$rs7250652_A",
              "dat$rs8107532_C = 2 - dat$rs8107532_T",
              "dat$rs11998681_T = 2 - dat$rs11998681_G",
              "dat$rs1561748_G = 2 - dat$rs1561748_C",
              "dat$rs343_C = 2 - dat$rs343_A",
              "dat$rs4922108_T = 2 - dat$rs4922108_C",
              "dat$rs4922113_G = 2 - dat$rs4922113_A",
              "dat$rs77069344_T = 2 - dat$rs77069344_G",
              "dat$rs78299715_C = 2 - dat$rs78299715_G",
              "dat$rs9644636_T = 2 - dat$rs9644636_G",
              "dat$rs117038725_G = 2 - dat$rs117038725_A",
              "dat$rs1441778_T = 2 - dat$rs1441778_C",
              "dat$rs1561748_G = 2 - dat$rs1561748_C",
              "dat$rs343_C = 2 - dat$rs343_A",
              "dat$rs4922113_G = 2 - dat$rs4922113_A",
              "dat$rs73004966_C = 2 - dat$rs73004966_T",
              "dat$rs74855321_C = 2 - dat$rs74855321_T",
              "dat$rs78299715_C = 2 - dat$rs78299715_G",
              "dat$rs7770628_T = 2 - dat$rs7770628_C",
              "dat$rs9355296_G = 2 - dat$rs9355296_A",
              "dat$rs2457550_C = 2 - dat$rs2457550_T",
              "dat$rs7770628_T = 2 - dat$rs7770628_C",
              "dat$rs57176252_C = 2 - dat$rs57176252_A",
              "dat$rs6453131_G = 2 - dat$rs6453131_T",
              "dat$rs73763094_C = 2 - dat$rs73763094_T",
              "dat$rs7717415_A = 2 - dat$rs7717415_G",
              "dat$rs11816594_G = 2 - dat$rs11816594_A",
              "dat$rs2303152_G = 2 - dat$rs2303152_A",
              "dat$rs57176252_C = 2 - dat$rs57176252_A",
              "dat$rs6453131_G = 2 - dat$rs6453131_T",
              "dat$rs7717415_A = 2 - dat$rs7717415_G",
              "dat$rs80223115_G = 2 - dat$rs80223115_A",
              "dat$rs6545966_A = 2 - dat$rs6545966_G",
              "dat$rs1334806_G = 2 - dat$rs1334806_A",
              "dat$rs7528706_A = 2 - dat$rs7528706_G",
              "dat$rs11125941_C = 2 - dat$rs11125941_A",
              "dat$rs13002939_A = 2 - dat$rs13002939_C",
              "dat$rs4341893_G = 2 - dat$rs4341893_A",
              "dat$rs4426495_T = 2 - dat$rs4426495_C",
              "dat$rs6487176_T = 2 - dat$rs6487176_A",
              "dat$rs10201379_A = 2 - dat$rs10201379_C",
              "dat$rs74358761_G = 2 - dat$rs74358761_A",
              "dat$rs144090867_C = 2 - dat$rs144090867_T",
              "dat$rs4341893_G = 2 - dat$rs4341893_A",
              "dat$rs75509845_G = 2 - dat$rs75509845_A",
              "dat$rs12628618_G = 2 - dat$rs12628618_T",
              "dat$rs111859808_G = 2 - dat$rs111859808_A",
              "dat$rs6104410_A = 2 - dat$rs6104410_G",
              "dat$rs16990971_A = 2 - dat$rs16990971_G",
              "dat$rs1402837_T = 2 - dat$rs1402837_C",
              "dat$rs151057583_T = 2 - dat$rs151057583_C",
              "dat$rs3770583_T = 2 - dat$rs3770583_G",
              "dat$rs486838_T = 2 - dat$rs486838_C",
              "dat$rs495074_T = 2 - dat$rs495074_C",
              "dat$rs560887_T = 2 - dat$rs560887_C",
              "dat$rs8076224_A = 2 - dat$rs8076224_T",
              "dat$rs1406981_T = 2 - dat$rs1406981_C",
              "dat$rs73119018_A = 2 - dat$rs73119018_G",
              "dat$rs146564328_A = 2 - dat$rs146564328_C",
              "dat$rs17205533_C = 2 - dat$rs17205533_G",
              "dat$rs7161785_C = 2 - dat$rs7161785_G",
              "dat$rs2284773_A = 2 - dat$rs2284773_G",
              "dat$rs2908277_A = 2 - dat$rs2908277_G",
              "dat$rs2908289_A = 2 - dat$rs2908289_G",
              "dat$rs60463592_T = 2 - dat$rs60463592_C",
              "dat$rs62459092_A = 2 - dat$rs62459092_G",
              "dat$rs741038_A = 2 - dat$rs741038_G",
              "dat$rs10023_A = 2 - dat$rs10023_G",
              "dat$rs7773695_T = 2 - dat$rs7773695_G",
              "dat$rs880347_A = 2 - dat$rs880347_G",
              "dat$rs9296285_A = 2 - dat$rs9296285_G",
              "dat$rs17390223_T = 2 - dat$rs17390223_C",
              "dat$rs8032477_T = 2 - dat$rs8032477_C",
              "dat$rs59165976_C = 2 - dat$rs59165976_G",
              "dat$rs8076224_A = 2 - dat$rs8076224_T",
              "dat$rs143129080_A = 2 - dat$rs143129080_G",
              "dat$rs35674932_A = 2 - dat$rs35674932_T",
              "dat$rs4507153_T = 2 - dat$rs4507153_C",
              "dat$rs4924458_C = 2 - dat$rs4924458_G",
              "dat$rs4507153_T = 2 - dat$rs4507153_C",
              "dat$rs117068364_T = 2 - dat$rs117068364_G",
              "dat$rs139500534_C = 2 - dat$rs139500534_G",
              "dat$rs1447105_T = 2 - dat$rs1447105_C",
              "dat$rs340515_T = 2 - dat$rs340515_G",
              "dat$rs4953153_A = 2 - dat$rs4953153_G",
              "dat$rs529438_T = 2 - dat$rs529438_C",
              "dat$rs561078_A = 2 - dat$rs561078_G",
              "dat$rs737447_T = 2 - dat$rs737447_G",
              "dat$rs2273800_A = 2 - dat$rs2273800_G",
              "dat$rs4905947_A = 2 - dat$rs4905947_G",
              "dat$rs2284773_A = 2 - dat$rs2284773_G",
              "dat$rs2908277_A = 2 - dat$rs2908277_G",
              "dat$rs2908289_A = 2 - dat$rs2908289_G",
              "dat$rs60463592_T = 2 - dat$rs60463592_C",
              "dat$rs62459092_A = 2 - dat$rs62459092_G",
              "dat$rs741038_A = 2 - dat$rs741038_G",
              "dat$rs3738002_A = 2 - dat$rs3738002_T",
              "dat$rs371747459_A = 2 - dat$rs371747459_G",
              "dat$rs1402837_C = 2 - dat$rs1402837_T",
              "dat$rs151057583_C = 2 - dat$rs151057583_T",
              "dat$rs3770583_G = 2 - dat$rs3770583_T",
              "dat$rs486838_C = 2 - dat$rs486838_T",
              "dat$rs495074_C = 2 - dat$rs495074_T",
              "dat$rs560887_C = 2 - dat$rs560887_T",
              "dat$rs8076224_T = 2 - dat$rs8076224_A",
              "dat$rs1406981_C = 2 - dat$rs1406981_T",
              "dat$rs73119018_G = 2 - dat$rs73119018_A",
              "dat$rs146564328_C = 2 - dat$rs146564328_A",
              "dat$rs17205533_G = 2 - dat$rs17205533_C",
              "dat$rs7161785_G = 2 - dat$rs7161785_C",
              "dat$rs2284773_G = 2 - dat$rs2284773_A",
              "dat$rs2908277_G = 2 - dat$rs2908277_A",
              "dat$rs2908289_G = 2 - dat$rs2908289_A",
              "dat$rs60463592_C = 2 - dat$rs60463592_T",
              "dat$rs62459092_G = 2 - dat$rs62459092_A",
              "dat$rs741038_G = 2 - dat$rs741038_A",
              "dat$rs10023_G = 2 - dat$rs10023_A",
              "dat$rs7773695_G = 2 - dat$rs7773695_T",
              "dat$rs880347_G = 2 - dat$rs880347_A",
              "dat$rs9296285_G = 2 - dat$rs9296285_A",
              "dat$rs17390223_C = 2 - dat$rs17390223_T",
              "dat$rs8032477_C = 2 - dat$rs8032477_T",
              "dat$rs59165976_G = 2 - dat$rs59165976_C",
              "dat$rs8076224_T = 2 - dat$rs8076224_A",
              "dat$rs143129080_G = 2 - dat$rs143129080_A",
              "dat$rs35674932_T = 2 - dat$rs35674932_A",
              "dat$rs4507153_C = 2 - dat$rs4507153_T",
              "dat$rs4924458_G = 2 - dat$rs4924458_C",
              "dat$rs4507153_C = 2 - dat$rs4507153_T",
              "dat$rs117068364_G = 2 - dat$rs117068364_T",
              "dat$rs139500534_G = 2 - dat$rs139500534_C",
              "dat$rs1447105_C = 2 - dat$rs1447105_T",
              "dat$rs340515_G = 2 - dat$rs340515_T",
              "dat$rs4953153_G = 2 - dat$rs4953153_A",
              "dat$rs529438_C = 2 - dat$rs529438_T",
              "dat$rs561078_G = 2 - dat$rs561078_A",
              "dat$rs737447_G = 2 - dat$rs737447_T",
              "dat$rs2273800_G = 2 - dat$rs2273800_A",
              "dat$rs4905947_G = 2 - dat$rs4905947_A",
              "dat$rs2284773_G = 2 - dat$rs2284773_A",
              "dat$rs2908277_G = 2 - dat$rs2908277_A",
              "dat$rs2908289_G = 2 - dat$rs2908289_A",
              "dat$rs60463592_C = 2 - dat$rs60463592_T",
              "dat$rs62459092_G = 2 - dat$rs62459092_A",
              "dat$rs741038_G = 2 - dat$rs741038_A",
              "dat$rs3738002_T = 2 - dat$rs3738002_A",
              "dat$rs371747459_G = 2 - dat$rs371747459_A",
              "dat$rs11196563_G = 2 - dat$rs11196563_A",
              "dat$rs57667171_T = 2 - dat$rs57667171_C")


# 自动执行所有语句，报错不中断
for (cmd in commands) {
  try(eval(parse(text = cmd)), silent = TRUE)
}


#################
#####logTG
dat$grs_ANGPTL3_w_glgc_logTG <- -1 * scale((dat$rs4495740_G * (-0.0817782))/1)

dat$grs_APOB_w_glgc_logTG <- -1 * scale((dat$rs57825321_A * (-0.0693446) + 
                                           dat$rs6706783_A * (0.0621935) )/2)


dat$grs_APOC3_w_glgc_logTG <- -1 * scale((dat$rs11216190_C * (-0.102214) + 
                                            dat$rs116852373_A * (-0.131741) + 
                                            dat$rs139733873_G * (0.243581) + 
                                            dat$rs17519093_A * (-0.10379) + 
                                            dat$rs4938309_T * (-0.139702) + 
                                            dat$rs518181_A * (-0.0577882) + 
                                            dat$rs595049_C * (0.100009) + 
                                            dat$rs651821_C * (0.290942) + 
                                            dat$rs79210031_T * (-0.182449) + 
                                            dat$rs888245_G * (0.110682) )/10)

dat$grs_APOC3_w_glgc_logTG_higher <- scale((dat$rs11216190_C * (-0.102214) + 
                                              dat$rs116852373_A * (-0.131741) + 
                                              dat$rs139733873_G * (0.243581) + 
                                              dat$rs17519093_A * (-0.10379) + 
                                              dat$rs4938309_T * (-0.139702) + 
                                              dat$rs518181_A * (-0.0577882) + 
                                              dat$rs595049_C * (0.100009) + 
                                              dat$rs651821_C * (0.290942) + 
                                              dat$rs79210031_T * (-0.182449) + 
                                              dat$rs888245_G * (0.110682) )/10)

dat$grs_CETP_w_glgc_logTG <- -1 * scale((dat$rs821840_G * (-0.0444137)))
# dat$grs_LDLR_w_glgc_logTG <- -1 * scale((dat$rs4804154_T * (-0.0277958) ))

dat$grs_LPL_w_glgc_logTG <- -1 * scale((dat$rs1561748_C * (-0.0618891) + 
                                          dat$rs4922113_A * (-0.0643988) + 
                                          dat$rs74855321_T * (-0.177221) + 
                                          dat$rs78299715_G * (-0.122596) )/4)



#ldl药靶GRS
dat$grs_ABCA1_w_glgc_ldl <- -1 * scale(dat$rs1883025_T * (-0.0397131))
dat$grs_APOB_w_glgc_ldl <- -1 * scale((dat$rs512535_T * (0.0635354) +
                                         dat$rs57825321_A * (-0.121726))/2)

dat$grs_HMGCR_w_glgc_ldl <- -1 * scale((
  dat$rs6453131_T * (-0.0815135) + 
    dat$rs7717415_G * (0.0407539))/2)

dat$grs_LDLR_w_glgc_ldl <- -1 * scale((dat$rs17243032_C * (0.0505803) + 
                                         dat$rs17699089_G * (-0.0620006) + 
                                         dat$rs2738464_G * (-0.10797) + 
                                         dat$rs442363_A * (0.0370919) + 
                                         dat$rs5927_A * (-0.0858549) + 
                                         dat$rs9789328_T * (0.0682901))/6)

dat$grs_LPA_w_glgc_ldl <- -1 * scale(dat$rs7770628_C * (0.0606087))

# dat$grs_MTTP_w_glgc_ldl <- -1 * scale(dat$rs2306986_C * (0.029681))

dat$grs_PCSK9_w_glgc_ldl <- -1 * scale((dat$rs10888897_T * (-0.0603903) + 
                                         dat$rs17111652_T * (0.0762897) + 
                                         dat$rs2479404_A * (-0.0352754) + 
                                         dat$rs34232196_T * (-0.0607008) + 
                                         dat$rs6663252_C * (-0.0750854))/5)









#hdl
dat$grs_ABCA1_w_glgc_hdl <- -1 * scale((dat$rs1800977_A * (0.0394685) + 
                                          dat$rs2472474_C * (0.0274939) + 
                                          dat$rs2740480_G * (-0.0492384) + 
                                          dat$rs2740487_A * (0.0396971) + 
                                          dat$rs2740488_C * (-0.0990417) + 
                                          dat$rs4149307_C * (-0.0745372) + 
                                          dat$rs73664373_C * (0.0721872))/7)


dat$grs_APOB_w_glgc_hdl <- -1 * scale((dat$rs12992267_T * (-0.0551877))/1)

dat$grs_APOC3_w_glgc_hdl <- -1 * scale((dat$rs111884008_C * (0.0741208) + 
                                         dat$rs1263172_A * (0.0296824) +
                                         dat$rs17519093_A * (0.103763) +  
                                         dat$rs4938309_T * (0.104774) + 
                                         dat$rs5094_A * (-0.0570135) + 
                                         dat$rs651821_C * (-0.168859))/6)

dat$grs_CETP_w_glgc_hdl <- -1 * scale((dat$rs1138429_T * (-0.0496921) +  
                                          dat$rs12934632_T * (-0.0967808) +  
                                          dat$rs148808285_T * (-0.134958) +  
                                          dat$rs148859408_C * (-0.17089) +  
                                          dat$rs17369578_A * (0.0532908) +  
                                          dat$rs185307787_C * (-0.150831) +  
                                          dat$rs188149780_C * (0.497995) +  
                                          dat$rs201825234_A * (-0.145084) +  
                                          dat$rs247617_A * (0.250976) +  
                                          dat$rs28495885_T * (0.0640403) +  
                                          dat$rs289716_T * (0.0406714) +  
                                          dat$rs289752_G * (-0.0356044) +  
                                          dat$rs36229485_T * (-0.129137) +  
                                          dat$rs62035923_A * (-0.0290072) +  
                                          dat$rs7187275_A * (0.0536496) +  
                                          dat$rs79582366_A * (0.0771163) +  
                                          dat$rs8059431_T * (-0.0751834) +  
                                          dat$rs9925265_A * (-0.056215) +  
                                          dat$rs9926440_C * (-0.150676))/21)


dat$grs_HCAR2_w_glgc_hdl <- -1 * scale(dat$rs1798192_T * (-0.032519))
dat$grs_HCAR3_w_glgc_hdl <- -1 * scale(dat$rs1798192_T * (-0.032519))
                                       
dat$grs_LDLR_w_glgc_hdl <- -1 * scale(dat$rs17699089_G * (-0.0628966))
dat$grs_LPA_w_glgc_hdl <- -1 * scale(dat$rs9355296_A * (-0.0386281))

dat$grs_LPL_w_glgc_hdl <- -1 * scale((dat$rs1561748_C * (0.0521463) + 
                                       dat$rs343_A * (0.0490636) + 
                                       dat$rs4922108_C * (-0.0297934) + 
                                       dat$rs4922113_A * (0.0542404) + 
                                       dat$rs77069344_G * (0.166403) + 
                                       dat$rs78299715_G * (0.126277) + 
                                       dat$rs9644636_G * (-0.0459704) )/8)


dat$grs_MTTP_w_glgc_hdl <- -1 * scale(dat$rs1115721_G * (0.0304242) )

dat$grs_PPARG_w_glgc_hdl <- -1 * scale((dat$rs12633551_T * (-0.040662) + 
                                           dat$rs75157957_T * (0.0342701) )/2)



##########
#TOTAL CHOL
dat$grs_ABCA1_w_glgc_TC <- -1 * scale((dat$rs2740480_G * (-0.0260511) + 
                                        dat$rs2740488_C * (-0.0765346) + 
                                        dat$rs4743762_C * (-0.0460668)) /3)

dat$grs_ANGPTL3_w_glgc_TC <- -1 * scale((dat$rs10889348_T * (-0.0497546))/1)

dat$grs_APOB_w_glgc_TC <- -1 * scale((dat$rs2337382_A * (0.0502508) + 
                                        dat$rs57825321_A * (-0.114637) + 
                                        dat$rs934197_A * (0.0639705))/3)

dat$grs_APOC3_w_glgc_TC <- -1 * scale((dat$rs2239013_T * (-0.0815331) + 
                                        dat$rs595049_C * (0.0418515) + 
                                        dat$rs651821_C * (0.0571283) + 
                                        dat$rs918143_T * (-0.0237512))/4)

dat$grs_CETP_w_glgc_TC <- -1 * scale((dat$rs17231506_T * (0.0716347) + 
                                        dat$rs188149780_C * (0.138118) + 
                                        dat$rs201825234_A * (-0.0440066) + 
                                        dat$rs9926440_C * (-0.0504504))/4)

dat$grs_HMGCR_w_glgc_TC <- -1 * scale((dat$rs2303152_A * (0.043022) + 
                                         dat$rs6453131_T * (-0.0677423) + 
                                         dat$rs7717415_G * (0.0359768))/3)

dat$grs_LDLR_w_glgc_TC <- -1 * scale((dat$rs17243032_C * (0.0366121) + 
                                        dat$rs2738446_G * (0.0690489) + 
                                        dat$rs2738464_G * (-0.0850829) + 
                                        dat$rs3810308_C * (-0.0750637) + 
                                        dat$rs442363_A * (0.0337199) + 
                                        dat$rs5927_A * (-0.0662167)  + 
                                        dat$rs8107532_T * (-0.075551))/7)

dat$grs_LPA_w_glgc_TC <- -1 * scale((dat$rs7770628_C * (0.0635399) )/1)

dat$grs_PCSK9_w_glgc_TC <- -1 * scale((dat$rs2479404_A * (-0.031361) + 
                                         dat$rs34232196_T * (-0.0597366) + 
                                         dat$rs534473_G * (-0.0533049) + 
                                         dat$rs565436_G * (-0.041905) + 
                                         dat$rs6663252_C * (-0.0743687) + 
                                         dat$rs67171713_C * (0.0630991))/6)




##########
#GLU-WEIGHT 降糖药靶GRS
dat$grs_ABCB11_w_taiwan_glu <- -1 * scale((dat$rs1402837_T * (0.105601901785905) + 
                                             dat$rs151057583_T * (0.0894239428855306) + 
                                             dat$rs3770583_T * (0.0753054488301792) + 
                                             dat$rs495074_T * (-0.0647278102067481) + 
                                             dat$rs560887_T * (-0.160159747452076))/6)

dat$grs_GCK_w_taiwan_glu <- -1 * scale((dat$rs2284773_G * (-1) * (-0.0817197453915558) + 
                                          dat$rs2908277_G * (-1) * (0.081647089637489) + 
                                          dat$rs2908289_G * (-1) * (0.129313068319187) + 
                                          dat$rs60463592_C * (-1) * (-0.0636837426427698) + 
                                          dat$rs62459092_G * (-1) * (0.0345144635541533) + 
                                          dat$rs741038_G * (-1) * (0.0448525318601902))/6)

dat$grs_GCK_w_taiwan_glu_higher <- scale((dat$rs2284773_G * (-1) * (-0.0817197453915558) + 
                                          dat$rs2908277_G * (-1) * (0.081647089637489) + 
                                          dat$rs2908289_G * (-1) * (0.129313068319187) + 
                                          dat$rs60463592_C * (-1) * (-0.0636837426427698) + 
                                          dat$rs62459092_G * (-1) * (0.0345144635541533) + 
                                          dat$rs741038_G * (-1) * (0.0448525318601902))/6)


dat$grs_GLP1R_w_taiwan_glu <- -1 * scale((dat$rs10023_A * (-0.0352008750721861) + 
                                            dat$rs880347_A * (-0.0296842304021061) + 
                                            dat$rs9296285_A * (0.0483843098910462))/3)








#############
#血糖GRS
glu.grs.snp <- fread("CAPOC.cis.glu_grs.snp.txt", h=T)

setDT(dat)
setDT(glu.grs.snp)

common_cols <- intersect(names(dat), names(glu.grs.snp))
common_cols <- setdiff(common_cols, "FID")  

# 从 glu.grs.snp 中删除重复列（除 FID）
glu.grs.snp_clean <- glu.grs.snp[, setdiff(names(glu.grs.snp), common_cols), with = FALSE]

dat1 <- merge(dat, glu.grs.snp_clean, by = "FID", all.x = TRUE)
dat = dat1




#########
#缺失值填上众数
##########

fill_missing_with_mode <- function(df) {
  get_mode <- function(x) {
    tab <- table(x)
    if (length(tab) == 0) return(NA)  
    mode_val <- names(which.max(tab))
    return(mode_val[1])  
  }
  
  for (col in colnames(df)) {
    non_na <- df[[col]][!is.na(df[[col]])]
    if (length(non_na) == 0) next  
    mode_value <- get_mode(non_na)
    
    if (is.numeric(df[[col]])) {
      mode_value <- as.numeric(mode_value)
    } else if (is.factor(df[[col]])) {
      mode_value <- as.factor(mode_value)
    }
    df[[col]][is.na(df[[col]])] <- mode_value
  }
  return(df)
}



dat <- fill_missing_with_mode(dat)
str(dat)






#填上snp
commands <- c("dat$rs10010131_A = 2 - dat$rs10010131_G",
              "dat$rs10158845_A = 2 - dat$rs10158845_G",
              "dat$rs10267664_C = 2 - dat$rs10267664_G",
              "dat$rs10441113_A = 2 - dat$rs10441113_G",
              "dat$rs10811658_A = 2 - dat$rs10811658_G",
              "dat$rs10830963_C = 2 - dat$rs10830963_G",
              "dat$rs10882083_C = 2 - dat$rs10882083_G",
              "dat$rs12429454_A = 2 - dat$rs12429454_G",
              "dat$rs12443160_T = 2 - dat$rs12443160_C",
              "dat$rs1260326_T = 2 - dat$rs1260326_C",
              "dat$rs13023591_T = 2 - dat$rs13023591_C",
              "dat$rs1402837_T = 2 - dat$rs1402837_C",
              "dat$rs1406981_T = 2 - dat$rs1406981_C",
              "dat$rs1409333_T = 2 - dat$rs1409333_C",
              "dat$rs1420566_A = 2 - dat$rs1420566_G",
              "dat$rs16922302_C = 2 - dat$rs16922302_G",
              "dat$rs17085675_A = 2 - dat$rs17085675_T",
              "dat$rs17168486_T = 2 - dat$rs17168486_C",
              "dat$rs174554_A = 2 - dat$rs174554_G",
              "dat$rs1974619_T = 2 - dat$rs1974619_C",
              "dat$rs204926_A = 2 - dat$rs204926_G",
              "dat$rs2237896_A = 2 - dat$rs2237896_G",
              "dat$rs2273800_A = 2 - dat$rs2273800_G",
              "dat$rs231840_T = 2 - dat$rs231840_C",
              "dat$rs243020_A = 2 - dat$rs243020_G",
              "dat$rs2908289_A = 2 - dat$rs2908289_G",
              "dat$rs340515_T = 2 - dat$rs340515_G",
              "dat$rs340878_C = 2 - dat$rs340878_G",
              "dat$rs35674932_A = 2 - dat$rs35674932_T",
              "dat$rs363404_T = 2 - dat$rs363404_C",
              "dat$rs3755934_T = 2 - dat$rs3755934_C",
              "dat$rs3787186_T = 2 - dat$rs3787186_C",
              "dat$rs3802177_A = 2 - dat$rs3802177_G",
              "dat$rs419145_A = 2 - dat$rs419145_G",
              "dat$rs4339696_T = 2 - dat$rs4339696_G",
              "dat$rs459193_A = 2 - dat$rs459193_G",
              "dat$rs519084_A = 2 - dat$rs519084_T",
              "dat$rs55804588_C = 2 - dat$rs55804588_G",
              "dat$rs56237081_T = 2 - dat$rs56237081_C",
              "dat$rs56252704_A = 2 - dat$rs56252704_G",
              "dat$rs59299413_C = 2 - dat$rs59299413_G",
              "dat$rs6048206_T = 2 - dat$rs6048206_C",
              "dat$rs61839365_A = 2 - dat$rs61839365_G",
              "dat$rs6495240_T = 2 - dat$rs6495240_C",
              "dat$rs67131976_T = 2 - dat$rs67131976_C",
              "dat$rs6780171_A = 2 - dat$rs6780171_T",
              "dat$rs6932473_A = 2 - dat$rs6932473_T",
              "dat$rs6972160_C = 2 - dat$rs6972160_G",
              "dat$rs7109575_A = 2 - dat$rs7109575_G",
              "dat$rs7161785_C = 2 - dat$rs7161785_G",
              "dat$rs7165887_A = 2 - dat$rs7165887_G",
              "dat$rs7209518_A = 2 - dat$rs7209518_C",
              "dat$rs883541_A = 2 - dat$rs883541_G",
              "dat$rs9263964_T = 2 - dat$rs9263964_C",
              "dat$rs9296285_A = 2 - dat$rs9296285_G",
              "dat$rs9788833_T = 2 - dat$rs9788833_C",
              "dat$rs9873341_T = 2 - dat$rs9873341_C",
              "dat$rs10010131_G = 2 - dat$rs10010131_A",
              "dat$rs10158845_G = 2 - dat$rs10158845_A",
              "dat$rs10267664_G = 2 - dat$rs10267664_C",
              "dat$rs10441113_G = 2 - dat$rs10441113_A",
              "dat$rs10811658_G = 2 - dat$rs10811658_A",
              "dat$rs10830963_G = 2 - dat$rs10830963_C",
              "dat$rs10882083_G = 2 - dat$rs10882083_C",
              "dat$rs12429454_G = 2 - dat$rs12429454_A",
              "dat$rs12443160_C = 2 - dat$rs12443160_T",
              "dat$rs1260326_C = 2 - dat$rs1260326_T",
              "dat$rs13023591_C = 2 - dat$rs13023591_T",
              "dat$rs1402837_C = 2 - dat$rs1402837_T",
              "dat$rs1406981_C = 2 - dat$rs1406981_T",
              "dat$rs1409333_C = 2 - dat$rs1409333_T",
              "dat$rs1420566_G = 2 - dat$rs1420566_A",
              "dat$rs16922302_G = 2 - dat$rs16922302_C",
              "dat$rs17085675_T = 2 - dat$rs17085675_A",
              "dat$rs17168486_C = 2 - dat$rs17168486_T",
              "dat$rs174554_G = 2 - dat$rs174554_A",
              "dat$rs1974619_C = 2 - dat$rs1974619_T",
              "dat$rs204926_G = 2 - dat$rs204926_A",
              "dat$rs2237896_G = 2 - dat$rs2237896_A",
              "dat$rs2273800_G = 2 - dat$rs2273800_A",
              "dat$rs231840_C = 2 - dat$rs231840_T",
              "dat$rs243020_G = 2 - dat$rs243020_A",
              "dat$rs2908289_G = 2 - dat$rs2908289_A",
              "dat$rs340515_G = 2 - dat$rs340515_T",
              "dat$rs340878_G = 2 - dat$rs340878_C",
              "dat$rs35674932_T = 2 - dat$rs35674932_A",
              "dat$rs363404_C = 2 - dat$rs363404_T",
              "dat$rs3755934_C = 2 - dat$rs3755934_T",
              "dat$rs3787186_C = 2 - dat$rs3787186_T",
              "dat$rs3802177_G = 2 - dat$rs3802177_A",
              "dat$rs419145_G = 2 - dat$rs419145_A",
              "dat$rs4339696_G = 2 - dat$rs4339696_T",
              "dat$rs459193_G = 2 - dat$rs459193_A",
              "dat$rs519084_T = 2 - dat$rs519084_A",
              "dat$rs55804588_G = 2 - dat$rs55804588_C",
              "dat$rs56237081_C = 2 - dat$rs56237081_T",
              "dat$rs56252704_G = 2 - dat$rs56252704_A",
              "dat$rs59299413_G = 2 - dat$rs59299413_C",
              "dat$rs6048206_C = 2 - dat$rs6048206_T",
              "dat$rs61839365_G = 2 - dat$rs61839365_A",
              "dat$rs6495240_C = 2 - dat$rs6495240_T",
              "dat$rs67131976_C = 2 - dat$rs67131976_T",
              "dat$rs6780171_T = 2 - dat$rs6780171_A",
              "dat$rs6932473_T = 2 - dat$rs6932473_A",
              "dat$rs6972160_G = 2 - dat$rs6972160_C",
              "dat$rs7109575_G = 2 - dat$rs7109575_A",
              "dat$rs7161785_G = 2 - dat$rs7161785_C",
              "dat$rs7165887_G = 2 - dat$rs7165887_A",
              "dat$rs7209518_C = 2 - dat$rs7209518_A",
              "dat$rs883541_G = 2 - dat$rs883541_A",
              "dat$rs9263964_C = 2 - dat$rs9263964_T",
              "dat$rs9296285_G = 2 - dat$rs9296285_A",
              "dat$rs9788833_C = 2 - dat$rs9788833_T",
              "dat$rs9873341_C = 2 - dat$rs9873341_T"
              )



# 
# snps <- c(
#   "rs10010131_A", "rs10158845_A", "rs10267664_C", "rs10441113_A",
#   "rs10811658_A", "rs10830963_C", "rs10882083_C", "rs12429454_A",
#   "rs12443160_T", "rs1260326_T", "rs13023591_T", "rs1402837_T",
#   "rs1406981_T", "rs1409333_T", "rs1420566_A", "rs16922302_C",
#   "rs17085675_A", "rs17168486_T", "rs174554_A", "rs1974619_T",
#   "rs204926_A",  "rs2237896_A", "rs2273800_A", "rs231840_T",
#   "rs243020_A",  "rs2908289_A", "rs340515_T",  "rs340878_C",
#   "rs35674932_A","rs363404_T",  "rs3755934_T", "rs3787186_T",
#   "rs3802177_A", "rs419145_A",  "rs4339696_T", "rs459193_A",
#   "rs519084_A",  "rs55804588_C","rs56237081_T","rs56252704_A",
#   "rs59299413_C","rs6048206_T", "rs61839365_A","rs6495240_T",
#   "rs67131976_T","rs6780171_A", "rs6932473_A","rs6972160_C",
#   "rs7109575_A", "rs7161785_C", "rs7165887_A","rs7209518_A",
#   "rs883541_A",  "rs9263964_T", "rs9296285_A","rs9788833_T",
#   "rs9873341_T"
# )
# 
# 
# miss_tbl <- dat %>%
#   summarise(across(all_of(snps),
#                    ~ sum(is.na(.)),
#                    .names = "{.col}_n_miss")) %>%
#   pivot_longer(everything(),
#                names_to = "SNP",
#                values_to = "n_missing") %>%
#   mutate(
#     SNP        = sub("_n_miss$", "", SNP),
#     pct_missing = round(100 * n_missing / nrow(dat), 2)
#   ) %>%
#   arrange(desc(n_missing))
# 
# print(miss_tbl)


# 自动执行所有语句，报错不中断
for (cmd in commands) {
  try(eval(parse(text = cmd)), silent = TRUE)
}
  
dat$grs_w_taiwan_glu <- -1 * scale((
  dat$rs10010131_A * (-0.0451670514129492) + 
    dat$rs10158845_A * (-0.0259456984212956) + 
    dat$rs10267664_C * (-0.0305392036558187) + 
    dat$rs10441113_A * (-0.0412104662653703) + 
    dat$rs10811658_A * (-0.042256593899443) + 
    dat$rs10830963_C * (-0.0842954850959991) + 
    dat$rs10882083_C * (0.0260495500577934) + 
    dat$rs12429454_A * (0.046686630181167) + 
    dat$rs12443160_T * (0.0243938711183734) + 
    dat$rs1260326_T * (-0.0594136979948764) + 
    dat$rs13023591_T * (0.0498786022734075) + 
    dat$rs1402837_T * (0.105601901785905) + 
    # dat$rs1406981_T * (0.0273556088532637) + 
    dat$rs1409333_T * (-0.0284358835109245) + 
    dat$rs1420566_A * (0.0238911815324376) + 
    dat$rs16922302_C * (0.0589596384204964) + 
    dat$rs17085675_A * (0.0428630219722681) + 
    dat$rs17168486_T * (0.0554472760385346) + 
    dat$rs174554_A * (0.0468679843667174) + 
    dat$rs1974619_T * (0.0689586889419305) + 
    dat$rs204926_A * (-0.038340226524706) + 
    dat$rs2237896_A * (-0.0450624737994237) + 
    # dat$rs2273800_A * (0.0272437235729481) + 
    dat$rs231840_T * (0.0321836680555029) + 
    dat$rs243020_A * (-0.0306622391550593) + 
    dat$rs2908289_A * (0.129313068319187) + 
    # dat$rs340515_T * (-0.0974012087552421) + 
    dat$rs340878_C * (0.0362209879034592) + 
    # dat$rs35674932_A * (-0.0369319081921918) + 
    dat$rs363404_T * (-0.0249720922345129) + 
    dat$rs3755934_T * (-0.0387763836916252) + 
    dat$rs3787186_T * (0.029649866198795) + 
    dat$rs3802177_A * (-0.0721402337753005) + 
    dat$rs419145_A * (0.0295237607135612) + 
    dat$rs4339696_T * (-0.0370588077769799) + 
    dat$rs459193_A * (-0.0286219384743933) + 
    dat$rs519084_A * (0.0264699724389698) + 
    dat$rs55804588_C * (-0.0413204829589097) + 
    dat$rs56237081_T * (0.0371573032310029) + 
    dat$rs56252704_A * (0.0471790644316816) + 
    dat$rs59299413_C * (0.0400771964396585) + 
    dat$rs6048206_T * (0.0963695476923827) + 
    dat$rs61839365_A * (-0.0363028517804961) + 
    dat$rs6495240_T * (-0.0292543869156661) + 
    dat$rs67131976_T * (0.0605120835069212) + 
    dat$rs6780171_A * (0.0293634698714103) + 
    dat$rs6932473_A * (-0.025278949222475) + 
    dat$rs6972160_C * (-0.0361365495259625) + 
    dat$rs7109575_A * (-0.0656517922061712) + 
    # dat$rs7161785_C * (-0.0544378939226819) + 
    dat$rs7165887_A * (-0.0267101818249656) + 
    dat$rs7209518_A * (0.0256995173290557) + 
    dat$rs883541_A * (-0.0258200700879679) + 
    dat$rs9263964_T * (0.0252471121107633) + 
    dat$rs9296285_A * (0.0483843098910462) + 
    dat$rs9788833_T * (0.0266577539431445) + 
    dat$rs9873341_T * (0.0315126762426803)) /52)



dat$grs_w_taiwan_glu_higher <- scale((
  dat$rs10010131_A * (-0.0451670514129492) + 
    dat$rs10158845_A * (-0.0259456984212956) + 
    dat$rs10267664_C * (-0.0305392036558187) + 
    dat$rs10441113_A * (-0.0412104662653703) + 
    dat$rs10811658_A * (-0.042256593899443) + 
    dat$rs10830963_C * (-0.0842954850959991) + 
    dat$rs10882083_C * (0.0260495500577934) + 
    dat$rs12429454_A * (0.046686630181167) + 
    dat$rs12443160_T * (0.0243938711183734) + 
    dat$rs1260326_T * (-0.0594136979948764) + 
    dat$rs13023591_T * (0.0498786022734075) + 
    dat$rs1402837_T * (0.105601901785905) + 
    # dat$rs1406981_T * (0.0273556088532637) + 
    dat$rs1409333_T * (-0.0284358835109245) + 
    dat$rs1420566_A * (0.0238911815324376) + 
    dat$rs16922302_C * (0.0589596384204964) + 
    dat$rs17085675_A * (0.0428630219722681) + 
    dat$rs17168486_T * (0.0554472760385346) + 
    dat$rs174554_A * (0.0468679843667174) + 
    dat$rs1974619_T * (0.0689586889419305) + 
    dat$rs204926_A * (-0.038340226524706) + 
    dat$rs2237896_A * (-0.0450624737994237) + 
    # dat$rs2273800_A * (0.0272437235729481) + 
    dat$rs231840_T * (0.0321836680555029) + 
    dat$rs243020_A * (-0.0306622391550593) + 
    dat$rs2908289_A * (0.129313068319187) + 
    # dat$rs340515_T * (-0.0974012087552421) + 
    dat$rs340878_C * (0.0362209879034592) + 
    # dat$rs35674932_A * (-0.0369319081921918) + 
    dat$rs363404_T * (-0.0249720922345129) + 
    dat$rs3755934_T * (-0.0387763836916252) + 
    dat$rs3787186_T * (0.029649866198795) + 
    dat$rs3802177_A * (-0.0721402337753005) + 
    dat$rs419145_A * (0.0295237607135612) + 
    dat$rs4339696_T * (-0.0370588077769799) + 
    dat$rs459193_A * (-0.0286219384743933) + 
    dat$rs519084_A * (0.0264699724389698) + 
    dat$rs55804588_C * (-0.0413204829589097) + 
    dat$rs56237081_T * (0.0371573032310029) + 
    dat$rs56252704_A * (0.0471790644316816) + 
    dat$rs59299413_C * (0.0400771964396585) + 
    dat$rs6048206_T * (0.0963695476923827) + 
    dat$rs61839365_A * (-0.0363028517804961) + 
    dat$rs6495240_T * (-0.0292543869156661) + 
    dat$rs67131976_T * (0.0605120835069212) + 
    dat$rs6780171_A * (0.0293634698714103) + 
    dat$rs6932473_A * (-0.025278949222475) + 
    dat$rs6972160_C * (-0.0361365495259625) + 
    dat$rs7109575_A * (-0.0656517922061712) + 
    # dat$rs7161785_C * (-0.0544378939226819) + 
    dat$rs7165887_A * (-0.0267101818249656) + 
    dat$rs7209518_A * (0.0256995173290557) + 
    dat$rs883541_A * (-0.0258200700879679) + 
    dat$rs9263964_T * (0.0252471121107633) + 
    dat$rs9296285_A * (0.0483843098910462) + 
    dat$rs9788833_T * (0.0266577539431445) + 
    dat$rs9873341_T * (0.0315126762426803)) /52)







#############
#tg.GRS
tg.grs.snp <- fread("CAPOC.cis.tc_tg_grs.snp.txt", h=T)

setDT(dat)
setDT(tg.grs.snp)

common_cols <- intersect(names(dat), names(tg.grs.snp))
common_cols <- setdiff(common_cols, "FID")  

# 从 tg.grs.snp 中删除重复列（除 FID）
tg.grs.snp_clean <- tg.grs.snp[, setdiff(names(tg.grs.snp), common_cols), with = FALSE]

dat1 <- merge(dat, tg.grs.snp_clean, by = "FID", all.x = TRUE)
dat = dat1



#########
#缺失值填上众数
##########

fill_missing_with_mode <- function(df) {
  get_mode <- function(x) {
    tab <- table(x)
    if (length(tab) == 0) return(NA)  
    mode_val <- names(which.max(tab))
    return(mode_val[1])  
  }
  
  for (col in colnames(df)) {
    non_na <- df[[col]][!is.na(df[[col]])]
    if (length(non_na) == 0) next  
    mode_value <- get_mode(non_na)
    
    if (is.numeric(df[[col]])) {
      mode_value <- as.numeric(mode_value)
    } else if (is.factor(df[[col]])) {
      mode_value <- as.factor(mode_value)
    }
    df[[col]][is.na(df[[col]])] <- mode_value
  }
  return(df)
}



dat <- fill_missing_with_mode(dat)
str(dat)


#填上snp
commands <- c("dat$rs1077835_G = 2 - dat$rs1077835_A",
              "dat$rs116111528_T = 2 - dat$rs116111528_C",
              "dat$rs1260326_C = 2 - dat$rs1260326_T",
              "dat$rs13246993_A = 2 - dat$rs13246993_G",
              "dat$rs139154032_A = 2 - dat$rs139154032_G",
              "dat$rs1475537_T = 2 - dat$rs1475537_C",
              "dat$rs16990971_G = 2 - dat$rs16990971_A",
              "dat$rs174559_A = 2 - dat$rs174559_G",
              "dat$rs2954021_A = 2 - dat$rs2954021_G",
              "dat$rs34828061_G = 2 - dat$rs34828061_A",
              "dat$rs35570672_C = 2 - dat$rs35570672_T",
              "dat$rs3752442_G = 2 - dat$rs3752442_A",
              "dat$rs3755980_T = 2 - dat$rs3755980_C",
              "dat$rs4418728_G = 2 - dat$rs4418728_T",
              "dat$rs4495740_G = 2 - dat$rs4495740_T",
              "dat$rs4704834_A = 2 - dat$rs4704834_G",
              "dat$rs483082_T = 2 - dat$rs483082_G",
              "dat$rs58542926_T = 2 - dat$rs58542926_C",
              "dat$rs651821_C = 2 - dat$rs651821_T",
              "dat$rs6706783_A = 2 - dat$rs6706783_T",
              "dat$rs7078456_C = 2 - dat$rs7078456_T",
              "dat$rs7165077_T = 2 - dat$rs7165077_C",
              "dat$rs74855321_T = 2 - dat$rs74855321_C",
              "dat$rs7928320_G = 2 - dat$rs7928320_C",
              "dat$rs7928320_G = 2 - dat$rs7928320_C",
              "dat$rs821840_G = 2 - dat$rs821840_A",
              "dat$rs9687832_A = 2 - dat$rs9687832_G",
              
              "dat$rs1077835_A = 2 - dat$rs1077835_G",
              "dat$rs116111528_C = 2 - dat$rs116111528_T",
              "dat$rs1260326_T = 2 - dat$rs1260326_C",
              "dat$rs13246993_G = 2 - dat$rs13246993_A",
              "dat$rs139154032_G = 2 - dat$rs139154032_A",
              "dat$rs1475537_C = 2 - dat$rs1475537_T",
              "dat$rs16990971_A = 2 - dat$rs16990971_G",
              "dat$rs174559_G = 2 - dat$rs174559_A",
              "dat$rs2954021_G = 2 - dat$rs2954021_A",
              "dat$rs34828061_A = 2 - dat$rs34828061_G",
              "dat$rs35570672_T = 2 - dat$rs35570672_C",
              "dat$rs3752442_A = 2 - dat$rs3752442_G",
              "dat$rs3755980_C = 2 - dat$rs3755980_T",
              "dat$rs4418728_T = 2 - dat$rs4418728_G",
              "dat$rs4495740_T = 2 - dat$rs4495740_G",
              "dat$rs4704834_G = 2 - dat$rs4704834_A",
              "dat$rs483082_G = 2 - dat$rs483082_T",
              "dat$rs58542926_C = 2 - dat$rs58542926_T",
              "dat$rs651821_T = 2 - dat$rs651821_C",
              "dat$rs6706783_T = 2 - dat$rs6706783_A",
              "dat$rs7078456_T = 2 - dat$rs7078456_C",
              "dat$rs7165077_C = 2 - dat$rs7165077_T",
              "dat$rs74855321_C = 2 - dat$rs74855321_T",
              "dat$rs7928320_C = 2 - dat$rs7928320_G",
              "dat$rs7928320_C = 2 - dat$rs7928320_G",
              "dat$rs821840_A = 2 - dat$rs821840_G",
              "dat$rs9687832_G = 2 - dat$rs9687832_A", 
              
              
              
              "dat$rs10438978_T = 2 - dat$rs10438978_C",
              "dat$rs1065853_T = 2 - dat$rs1065853_G",
              "dat$rs1077835_G = 2 - dat$rs1077835_A",
              "dat$rs10846744_G = 2 - dat$rs10846744_C",
              "dat$rs10860591_C = 2 - dat$rs10860591_T",
              "dat$rs10889348_T = 2 - dat$rs10889348_A",
              "dat$rs11102967_C = 2 - dat$rs11102967_T",
              "dat$rs112784971_T = 2 - dat$rs112784971_C",
              "dat$rs11320208_A = 2 - dat$rs11320208_C",
              "dat$rs11601507_A = 2 - dat$rs11601507_C",
              "dat$rs12162136_G = 2 - dat$rs12162136_A",
              "dat$rs12173764_T = 2 - dat$rs12173764_C",
              "dat$rs12233358_G = 2 - dat$rs12233358_C",
              "dat$rs12314392_G = 2 - dat$rs12314392_A",
              "dat$rs12484801_T = 2 - dat$rs12484801_C",
              "dat$rs12995972_A = 2 - dat$rs12995972_C",
              "dat$rs151021730_A = 2 - dat$rs151021730_G",
              "dat$rs17129638_C = 2 - dat$rs17129638_T",
              "dat$rs17231506_T = 2 - dat$rs17231506_C",
              "dat$rs174565_G = 2 - dat$rs174565_C",
              "dat$rs17660635_G = 2 - dat$rs17660635_A",
              "dat$rs184643955_A = 2 - dat$rs184643955_G",
              "dat$rs2206689_A = 2 - dat$rs2206689_G",
              "dat$rs2263608_T = 2 - dat$rs2263608_A",
              "dat$rs2306986_C = 2 - dat$rs2306986_G",
              "dat$rs2714875_A = 2 - dat$rs2714875_G",
              "dat$rs2737265_G = 2 - dat$rs2737265_A",
              "dat$rs2738464_G = 2 - dat$rs2738464_C",
              "dat$rs2740488_C = 2 - dat$rs2740488_A",
              "dat$rs2807834_T = 2 - dat$rs2807834_G",
              "dat$rs2954021_A = 2 - dat$rs2954021_G",
              "dat$rs3809868_G = 2 - dat$rs3809868_A",
              "dat$rs3810308_C = 2 - dat$rs3810308_T",
              "dat$rs4344974_T = 2 - dat$rs4344974_C",
              "dat$rs4454000_G = 2 - dat$rs4454000_A",
              "dat$rs4622073_A = 2 - dat$rs4622073_G",
              "dat$rs4690014_A = 2 - dat$rs4690014_G",
              "dat$rs532436_A = 2 - dat$rs532436_G",
              "dat$rs55686478_G = 2 - dat$rs55686478_A",
              "dat$rs56130071_C = 2 - dat$rs56130071_G",
              "dat$rs57176252_A = 2 - dat$rs57176252_C",
              "dat$rs5758361_C = 2 - dat$rs5758361_A",
              "dat$rs57825321_A = 2 - dat$rs57825321_T",
              "dat$rs61946969_A = 2 - dat$rs61946969_G",
              "dat$rs6437108_T = 2 - dat$rs6437108_C",
              "dat$rs6453131_T = 2 - dat$rs6453131_G",
              "dat$rs651821_C = 2 - dat$rs651821_T",
              "dat$rs6663252_C = 2 - dat$rs6663252_T",
              "dat$rs6874202_T = 2 - dat$rs6874202_C",
              "dat$rs6959252_A = 2 - dat$rs6959252_G",
              "dat$rs7165077_T = 2 - dat$rs7165077_C",
              "dat$rs72607108_G = 2 - dat$rs72607108_T",
              "dat$rs7312955_A = 2 - dat$rs7312955_C",
              "dat$rs7400614_A = 2 - dat$rs7400614_C",
              "dat$rs75352129_T = 2 - dat$rs75352129_C",
              "dat$rs77303550_T = 2 - dat$rs77303550_C",
              "dat$rs7770628_C = 2 - dat$rs7770628_T",
              "dat$rs78753419_T = 2 - dat$rs78753419_A",
              "dat$rs7954039_A = 2 - dat$rs7954039_C",
              "dat$rs814295_G = 2 - dat$rs814295_A",
              "dat$rs9376090_C = 2 - dat$rs9376090_T",
              "dat$rs9399668_T = 2 - dat$rs9399668_C",
              "dat$rs9953437_A = 2 - dat$rs9953437_G",
              "dat$rs10438978_C = 2 - dat$rs10438978_T",
              "dat$rs1065853_G = 2 - dat$rs1065853_T",
              "dat$rs1077835_A = 2 - dat$rs1077835_G",
              "dat$rs10846744_C = 2 - dat$rs10846744_G",
              "dat$rs10860591_T = 2 - dat$rs10860591_C",
              "dat$rs10889348_A = 2 - dat$rs10889348_T",
              "dat$rs11102967_T = 2 - dat$rs11102967_C",
              "dat$rs112784971_C = 2 - dat$rs112784971_T",
              "dat$rs11320208_C = 2 - dat$rs11320208_A",
              "dat$rs11601507_C = 2 - dat$rs11601507_A",
              "dat$rs12162136_A = 2 - dat$rs12162136_G",
              "dat$rs12173764_C = 2 - dat$rs12173764_T",
              "dat$rs12233358_C = 2 - dat$rs12233358_G",
              "dat$rs12314392_A = 2 - dat$rs12314392_G",
              "dat$rs12484801_C = 2 - dat$rs12484801_T",
              "dat$rs12995972_C = 2 - dat$rs12995972_A",
              "dat$rs151021730_G = 2 - dat$rs151021730_A",
              "dat$rs17129638_T = 2 - dat$rs17129638_C",
              "dat$rs17231506_C = 2 - dat$rs17231506_T",
              "dat$rs174565_C = 2 - dat$rs174565_G",
              "dat$rs17660635_A = 2 - dat$rs17660635_G",
              "dat$rs184643955_G = 2 - dat$rs184643955_A",
              "dat$rs2206689_G = 2 - dat$rs2206689_A",
              "dat$rs2263608_A = 2 - dat$rs2263608_T",
              "dat$rs2306986_G = 2 - dat$rs2306986_C",
              "dat$rs2714875_G = 2 - dat$rs2714875_A",
              "dat$rs2737265_A = 2 - dat$rs2737265_G",
              "dat$rs2738464_C = 2 - dat$rs2738464_G",
              "dat$rs2740488_A = 2 - dat$rs2740488_C",
              "dat$rs2807834_G = 2 - dat$rs2807834_T",
              "dat$rs2954021_G = 2 - dat$rs2954021_A",
              "dat$rs3809868_A = 2 - dat$rs3809868_G",
              "dat$rs3810308_T = 2 - dat$rs3810308_C",
              "dat$rs4344974_C = 2 - dat$rs4344974_T",
              "dat$rs4454000_A = 2 - dat$rs4454000_G",
              "dat$rs4622073_G = 2 - dat$rs4622073_A",
              "dat$rs4690014_G = 2 - dat$rs4690014_A",
              "dat$rs532436_G = 2 - dat$rs532436_A",
              "dat$rs55686478_A = 2 - dat$rs55686478_G",
              "dat$rs56130071_G = 2 - dat$rs56130071_C",
              "dat$rs57176252_C = 2 - dat$rs57176252_A",
              "dat$rs5758361_A = 2 - dat$rs5758361_C",
              "dat$rs57825321_T = 2 - dat$rs57825321_A",
              "dat$rs61946969_G = 2 - dat$rs61946969_A",
              "dat$rs6437108_C = 2 - dat$rs6437108_T",
              "dat$rs6453131_G = 2 - dat$rs6453131_T",
              "dat$rs651821_T = 2 - dat$rs651821_C",
              "dat$rs6663252_T = 2 - dat$rs6663252_C",
              "dat$rs6874202_C = 2 - dat$rs6874202_T",
              "dat$rs6959252_G = 2 - dat$rs6959252_A",
              "dat$rs7165077_C = 2 - dat$rs7165077_T",
              "dat$rs72607108_T = 2 - dat$rs72607108_G",
              "dat$rs7312955_C = 2 - dat$rs7312955_A",
              "dat$rs7400614_C = 2 - dat$rs7400614_A",
              "dat$rs75352129_C = 2 - dat$rs75352129_T",
              "dat$rs77303550_C = 2 - dat$rs77303550_T",
              "dat$rs7770628_T = 2 - dat$rs7770628_C",
              "dat$rs78753419_A = 2 - dat$rs78753419_T",
              "dat$rs7954039_C = 2 - dat$rs7954039_A",
              "dat$rs814295_A = 2 - dat$rs814295_G",
              "dat$rs9376090_T = 2 - dat$rs9376090_C",
              "dat$rs9399668_C = 2 - dat$rs9399668_T",
              "dat$rs9953437_G = 2 - dat$rs9953437_A"
)



# snps2 <- unique(c(
#   "rs1077835_G", "rs116111528_T", "rs1260326_C", "rs13246993_A",
#   "rs139154032_A", "rs1475537_T", "rs16990971_G", "rs174559_A",
#   "rs2954021_A", "rs34828061_G", "rs35570672_C", "rs3752442_G",
#   "rs3755980_T", "rs4418728_G", "rs4495740_G", "rs4704834_A",
#   "rs483082_T", "rs58542926_T", "rs651821_C", "rs6706783_A",
#   "rs7078456_C", "rs7165077_T", "rs74855321_T", "rs7928320_G",
#   "rs821840_G", "rs9687832_A"
# ))
# 
# ## 2. 计算每个 SNP 的缺失行数与缺失占比 ------------------------
# miss_tbl2 <- dat %>%
#   summarise(across(all_of(snps2),
#                    ~ sum(is.na(.)),
#                    .names = "{.col}_n_miss")) %>%
#   pivot_longer(everything(),
#                names_to  = "SNP",
#                values_to = "n_missing") %>%
#   mutate(
#     SNP          = sub("_n_miss$", "", SNP),
#     pct_missing  = round(100 * n_missing / nrow(dat), 2)
#   ) %>%
#   arrange(desc(n_missing))
# 
# print(miss_tbl2)



# 自动执行所有语句，报错不中断
for (cmd in commands) {
  try(eval(parse(text = cmd)), silent = TRUE)
}

dat$grs_w_glgc_logTG <- -1 * scale((
  dat$rs1077835_G * (0.0702666) + 
    dat$rs116111528_T * (0.0868355) + 
    dat$rs1260326_C * (-0.0982948) + 
    dat$rs13246993_A * (-0.134801) + 
    dat$rs139154032_A * (-0.137224) + 
    dat$rs1475537_T * (0.0369542) + 
    # dat$rs16990971_G * (0.0792017) + 
    dat$rs174559_A * (0.0413253) + 
    dat$rs2954021_A * (0.0707375) + 
    dat$rs34828061_G * (-0.0333604) + 
    dat$rs35570672_C * (-0.0407339) + 
    dat$rs3752442_G * (-0.0381573) + 
    dat$rs3755980_T * (0.0316024) + 
    dat$rs4418728_G * (0.0394342) + 
    dat$rs4495740_G * (-0.0817782) + 
    dat$rs4704834_A * (-0.0532378) + 
    dat$rs483082_T * (0.112978) + 
    dat$rs58542926_T * (-0.089976) + 
    dat$rs651821_C * (0.290942) + 
    dat$rs6706783_A * (0.0621935) + 
    dat$rs7078456_C * (0.0343025) + 
    dat$rs7165077_T * (-0.0496796) + 
    dat$rs74855321_T * (-0.177221) + 
    dat$rs7928320_G * (-0.116974) + 
    dat$rs7928320_G * (-0.116974) + 
    dat$rs821840_G * (-0.0444137) + 
    dat$rs9687832_A * (0.0451107)
) /25)



dat$grs_w_glgc_logTG_higher <- scale((
  dat$rs1077835_G * (0.0702666) + 
    dat$rs116111528_T * (0.0868355) + 
    dat$rs1260326_C * (-0.0982948) + 
    dat$rs13246993_A * (-0.134801) + 
    dat$rs139154032_A * (-0.137224) + 
    dat$rs1475537_T * (0.0369542) + 
  #  dat$rs16990971_G * (0.0792017) + 
    dat$rs174559_A * (0.0413253) + 
    dat$rs2954021_A * (0.0707375) + 
    dat$rs34828061_G * (-0.0333604) + 
    dat$rs35570672_C * (-0.0407339) + 
    dat$rs3752442_G * (-0.0381573) + 
    dat$rs3755980_T * (0.0316024) + 
    dat$rs4418728_G * (0.0394342) + 
    dat$rs4495740_G * (-0.0817782) + 
    dat$rs4704834_A * (-0.0532378) + 
    dat$rs483082_T * (0.112978) + 
    dat$rs58542926_T * (-0.089976) + 
    dat$rs651821_C * (0.290942) + 
    dat$rs6706783_A * (0.0621935) + 
    dat$rs7078456_C * (0.0343025) + 
    dat$rs7165077_T * (-0.0496796) + 
    dat$rs74855321_T * (-0.177221) + 
    dat$rs7928320_G * (-0.116974) + 
    dat$rs7928320_G * (-0.116974) + 
    dat$rs821840_G * (-0.0444137) + 
    dat$rs9687832_A * (0.0451107)
) /25)




dat$grs_w_glgc_TC <- -1 * scale((dat$rs10438978_T * (-0.028137) + 
                                      dat$rs1065853_T * (-0.370244) + 
                                      dat$rs1077835_G * (0.0633784) + 
                                      dat$rs10846744_G * (-0.0291246) + 
                                      dat$rs10860591_C * (0.0341989) + 
                                      dat$rs10889348_T * (-0.0497546) + 
                                      dat$rs11102967_C * (-0.137375) + 
                                      dat$rs112784971_T * (-0.0462477) + 
                                      dat$rs11320208_A * (-0.0412222) + 
                                      dat$rs11601507_A * (0.0590618) + 
                                      dat$rs12162136_G * (0.0340664) + 
                                      dat$rs12173764_T * (-0.0322456) + 
                                      dat$rs12314392_G * (0.0270593) + 
                                      dat$rs12484801_T * (-0.0275028) + 
                                      dat$rs12995972_A * (-0.034681) + 
                                      dat$rs151021730_A * (-0.0619466) + 
                                      dat$rs17129638_C * (-0.0291265) + 
                                      dat$rs17231506_T * (0.0716347) + 
                                      dat$rs174565_G * (-0.0397076) + 
                                      dat$rs17660635_G * (-0.10843) + 
                                      dat$rs184643955_A * (0.0986131) + 
                                      dat$rs2206689_A * (0.028864) + 
                                      dat$rs2263608_T * (0.038167) + 
                                      dat$rs2306986_C * (0.0402142) + 
                                      dat$rs2714875_A * (0.0250345) + 
                                      dat$rs2738464_G * (-0.0850829) + 
                                      dat$rs2740488_C * (-0.0765346) + 
                                      dat$rs2807834_T * (-0.0422582) + 
                                      dat$rs2954021_A * (0.0547932) + 
                                      dat$rs3809868_G * (0.0296915) + 
                                      dat$rs3810308_C * (-0.0750637) + 
                                      dat$rs4344974_T * (-0.0401339) + 
                                      dat$rs4454000_G * (-0.0272272) + 
                                      dat$rs4622073_A * (-0.0554232) + 
                                      dat$rs4690014_A * (-0.0331682) + 
                                      dat$rs532436_A * (0.0685843) + 
                                      dat$rs55686478_G * (-0.024068) + 
                                      dat$rs56130071_C * (0.0299187) + 
                                      dat$rs5758361_C * (0.0276585) + 
                                      dat$rs57825321_A * (-0.114637) + 
                                      dat$rs61946969_A * (-0.0315782) + 
                                      dat$rs6437108_T * (-0.0393568) + 
                                      dat$rs6453131_T * (-0.0677423) + 
                                      dat$rs651821_C * (0.0571283) + 
                                      dat$rs6663252_C * (-0.0743687) + 
                                      dat$rs6874202_T * (-0.062447) + 
                                      dat$rs6959252_A * (0.0259676) + 
                                      dat$rs7165077_T * (-0.0513008) + 
                                      dat$rs72607108_G * (-0.0494518) + 
                                      dat$rs7312955_A * (0.0254173) + 
                                      dat$rs7400614_A * (0.0304648) + 
                                      dat$rs75352129_T * (0.136357) + 
                                      dat$rs77303550_T * (-0.0881541) + 
                                      dat$rs7770628_C * (0.0635399) + 
                                      dat$rs78753419_T * (0.0345461) + 
                                      dat$rs7954039_A * (0.0318793) + 
                                      dat$rs814295_G * (-0.0510712) + 
                                      dat$rs9376090_C * (-0.0448819) + 
                                      dat$rs9399668_T * (-0.0257582) + 
                                      dat$rs9953437_A * (0.0594315) )/60)



dat$grs_w_glgc_TC_higher <- scale((dat$rs10438978_T * (-0.028137) + 
                                      dat$rs1065853_T * (-0.370244) + 
                                      dat$rs1077835_G * (0.0633784) + 
                                      dat$rs10846744_G * (-0.0291246) + 
                                      dat$rs10860591_C * (0.0341989) + 
                                      dat$rs10889348_T * (-0.0497546) + 
                                      dat$rs11102967_C * (-0.137375) + 
                                      dat$rs112784971_T * (-0.0462477) + 
                                      dat$rs11320208_A * (-0.0412222) + 
                                      dat$rs11601507_A * (0.0590618) + 
                                      dat$rs12162136_G * (0.0340664) + 
                                      dat$rs12173764_T * (-0.0322456) + 
                                      dat$rs12314392_G * (0.0270593) + 
                                      dat$rs12484801_T * (-0.0275028) + 
                                      dat$rs12995972_A * (-0.034681) + 
                                      dat$rs151021730_A * (-0.0619466) + 
                                      dat$rs17129638_C * (-0.0291265) + 
                                      dat$rs17231506_T * (0.0716347) + 
                                      dat$rs174565_G * (-0.0397076) + 
                                      dat$rs17660635_G * (-0.10843) + 
                                      dat$rs184643955_A * (0.0986131) + 
                                      dat$rs2206689_A * (0.028864) + 
                                      dat$rs2263608_T * (0.038167) + 
                                      dat$rs2306986_C * (0.0402142) + 
                                      dat$rs2714875_A * (0.0250345) + 
                                      dat$rs2738464_G * (-0.0850829) + 
                                      dat$rs2740488_C * (-0.0765346) + 
                                      dat$rs2807834_T * (-0.0422582) + 
                                      dat$rs2954021_A * (0.0547932) + 
                                      dat$rs3809868_G * (0.0296915) + 
                                      dat$rs3810308_C * (-0.0750637) + 
                                      dat$rs4344974_T * (-0.0401339) + 
                                      dat$rs4454000_G * (-0.0272272) + 
                                      dat$rs4622073_A * (-0.0554232) + 
                                      dat$rs4690014_A * (-0.0331682) + 
                                      dat$rs532436_A * (0.0685843) + 
                                      dat$rs55686478_G * (-0.024068) + 
                                      dat$rs56130071_C * (0.0299187) + 
                                      dat$rs5758361_C * (0.0276585) + 
                                      dat$rs57825321_A * (-0.114637) + 
                                      dat$rs61946969_A * (-0.0315782) + 
                                      dat$rs6437108_T * (-0.0393568) + 
                                      dat$rs6453131_T * (-0.0677423) + 
                                      dat$rs651821_C * (0.0571283) + 
                                      dat$rs6663252_C * (-0.0743687) + 
                                      dat$rs6874202_T * (-0.062447) + 
                                      dat$rs6959252_A * (0.0259676) + 
                                      dat$rs7165077_T * (-0.0513008) + 
                                      dat$rs72607108_G * (-0.0494518) + 
                                      dat$rs7312955_A * (0.0254173) + 
                                      dat$rs7400614_A * (0.0304648) + 
                                      dat$rs75352129_T * (0.136357) + 
                                      dat$rs77303550_T * (-0.0881541) + 
                                      dat$rs7770628_C * (0.0635399) + 
                                      dat$rs78753419_T * (0.0345461) + 
                                      dat$rs7954039_A * (0.0318793) + 
                                      dat$rs814295_G * (-0.0510712) + 
                                      dat$rs9376090_C * (-0.0448819) + 
                                      dat$rs9399668_T * (-0.0257582) + 
                                      dat$rs9953437_A * (0.0594315) )/60)




write.table(dat, "CAPOC.cis.capoc_cases.grs.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")



######
#检查一下 是否有构建的GRS是NA的 那可能是SNP未提取全造成的
