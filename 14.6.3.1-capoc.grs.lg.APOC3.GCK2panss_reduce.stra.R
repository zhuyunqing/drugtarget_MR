library(data.table)
library(broom)
library(dplyr)
library(purrr) 

#########
#stratified

######################
#根据药物类别



dat <- fread("/Users/zhuyunqing/Desktop/study/drugMR_lipid/data/onesamplemr/CAPOC.cis.capoc_cases.grs.pheno.new.txt", h=T)
dat$age2 = dat$age * dat$age
dat$medication_metabo <- ifelse((dat$medication==1 | dat$medication==2 | dat$medication==3), 1,2 )
dat$course_3g <- ifelse(dat$first_episode_2g == 1, 1, ifelse(dat$course > 60, 3, 2))


exposures <- c(
  "grs_APOC3_w_glgc_logTG",
  "grs_APOC3_w_glgc_TC",
  "grs_GCK_w_taiwan_glu"
)

outcomes <- c(
  "PANSS_reduce_rate",
  "PANSS_P_reduce_rate",
  "PANSS_N_reduce_rate",
  "PANSS_G_reduce_rate"
)

covars_common <- paste(
  "age", "age2",                  
  "as.factor(sex)", "as.factor(centers)",
  "PC1", "PC2", "PC3", "PC4", "PC5",
  "course", "drugp07",
  sep = " + "
)



res_all <- map_dfr(exposures, function(expo) {
  
  gene_name <- sub("^grs_([^_]+).*", "\\1", expo)
  map_dfr(outcomes, function(outcome) {
    
    map_dfr(1:2, function(med) {
      
      fml <- as.formula(
        paste(outcome, "~", expo, "+", covars_common)
      )
      
      fit <- tryCatch(
        lm(fml, data = filter(dat, medication_metabo == med)),
        error = function(e) NULL
      )
      if (is.null(fit)) return(NULL)
      
      tidy(fit, conf.int = TRUE) %>%
        filter(term == expo) %>%  
        transmute(
          
          panss_cat = outcome,
          gene      = gene_name,
          med_group = med,
          beta      = estimate,
          lower95   = conf.low,
          upper95   = conf.high,
          p         = p.value
          
        )
    })
  })
})

stra_medication_metabo2 <- res_all
write.csv(stra_medication_metabo2, 
          "/Users/zhuyunqing/Desktop/study/drugMR_lipid/data/onesamplemr/res_grs.APOC3.GCK_2panss_redu.linear.stra_medication_metabo_2g.csv", 
          row.names = FALSE)

########
#stable 21 figure 4



####
#logistic
outcomes <- c(
  "PANSS_reduce_rate_2g",
  "PANSS_P_reduce_rate_2g",
  "PANSS_N_reduce_rate_2g",
  "PANSS_G_reduce_rate_2g",
  "PANSS_reduce_rate_2g_med",
  "PANSS_P_reduce_rate_2g_med",
  "PANSS_N_reduce_rate_2g_med",
  "PANSS_G_reduce_rate_2g_med"
)


res_all <- map_dfr(exposures, function(expo) {
  gene_name <- sub("^grs_([^_]+).*", "\\1", expo)
  map_dfr(outcomes, function(outcome) {
    map_dfr(1:2, function(med) {
      
      fml <- as.formula(
        paste(outcome, "~", expo, "+", covars_common)
      )
      
      fit <- tryCatch(
        glm(fml, data = filter(dat, medication_metabo == med), family = binomial()),
        error = function(e) NULL
      )
      if (is.null(fit)) return(NULL)
      
      tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
        filter(term == expo) %>%
        transmute(
          panss_cat = outcome,
          gene = gene_name,
          med_group = med,
          OR = estimate,
          cilower = conf.low,
          ciupper = conf.high,
          p = p.value
        )
    })
  })
})


stra_medication_metabo2_log <- res_all
write.csv(stra_medication_metabo2_log, 
          "/Users/zhuyunqing/Desktop/study/drugMR_lipid/data/onesamplemr/res_grs.APOC3.GCK_2panss_redu.linear.stra_medication_metabo_2g_log.csv", 
          row.names = FALSE)


########
#stable 22




#######
#首发未用药

idlist <- c(
  140004, 510002, 410002, 420005, 140010, 420009, 420011, 420013, 250012, 250014, 
  210011, 410015, 210015, 330005, 520004, 140026, 520008, 530012, 120029, 140033, 
  240022, 330011, 240026, 320012, 120034, 320016, 210028, 210029, 220033, 140043, 
  320018, 320019, 230035, 530018, 130057, 140053, 150050, 140055, 140058, 210036, 
  140060, 140067, 140073, 240046, 320025, 210053, 310031, 530031, 210047, 310032, 
  520035, 130083, 230051, 440055, 440056, 520040, 130082, 320034, 340033, 530046, 
  310038, 320035, 530047, 230062, 230063, 410057, 520054, 520056, 140090, 410063, 
  130092, 230068, 220071, 140098, 310051, 410066, 450070, 510063, 220074, 220077, 
  230080, 530064, 220079, 230065, 330063, 450072, 520070, 330068, 230092, 330071, 
  350072, 510075, 120111, 530080, 420076, 120112, 130115, 310079, 210100, 330089, 
  530085, 310085, 140125, 330093, 410079, 130136, 220110, 150132, 330094, 410078, 
  520098, 240116, 410085, 530104, 140142, 140144, 230120, 230121, 210124, 230129, 
  220133, 230132, 520105, 120152, 130155, 140149, 230136, 310099, 310115, 210138, 
  140156, 330105, 330111, 340109, 340110, 410099, 140168, 230149, 310120, 440101, 
  310124, 140171, 230150, 320114, 120172, 120175, 320125, 230158, 330123, 330126, 
  330132, 310134, 320135, 120186, 430116, 520125, 120193, 220171, 140201, 140203, 
  520128, 230172, 320147, 120207, 340148, 120208, 320153, 330152, 240189, 220184, 
  220186, 230183, 150219, 230188, 220191, 330162, 130230, 310169, 230199, 320166, 
  130231, 120233, 220202, 230204, 140243, 140244, 440138, 410142, 220214, 410145, 
  510193, 130258, 120254, 330186, 230218, 320187, 340185, 510209, 120269, 220223, 
  320191, 230225, 120276, 220226, 510212, 510218, 230229, 510221, 430157, 120284, 
  140285, 430161, 510224, 310198, 230238, 130295, 450165, 120299, 310200, 220252, 
  250254, 310203, 430168, 220258, 130309, 220259, 220263, 220264, 120311, 220265, 
  320207, 340205, 450172, 220270, 220271, 220273, 320210, 220275, 340211, 420176, 
  140324, 510246, 310213, 130337, 120339, 130351, 220290, 130356, 240291, 320220, 
  440186, 240292, 220295, 120369, 310222, 140371, 310242, 320225, 130377, 110382, 
  140380, 230302, 340226, 230304, 230305, 150384, 230309, 140388, 230312, 510272, 
  520278, 230316, 340233, 150398, 130406, 510286, 340235, 150409, 140411, 510292, 
  150415, 520294, 140421, 510298, 130419, 520305, 520306, 210330, 210332, 340241, 
  150429, 340243, 220334, 320249, 130433, 130443, 230338, 440233, 160442, 160443, 
  340255, 340265, 160441, 230343, 380271, 130447, 160450, 330276, 520322, 310278, 
  310285, 340282, 510327, 130458, 160461, 340283, 160464, 370289, 380288, 520329, 
  250348, 240350, 320301, 340305, 380304, 130475, 340311, 390313, 130490, 160492, 
  130494, 470268, 160507, 320329, 520337, 140500, 320330, 340327, 340328, 230360, 
  250362, 120512, 130513, 130514, 390347, 130517, 380356, 470288, 130524, 370359, 
  240379, 250381, 250382, 140529, 140530, 250392, 140536, 320375, 130542, 150537, 
  160541, 320382, 430296, 250397, 340387, 340395, 130543, 230395, 140547, 160561, 
  160560, 210403, 340407, 340416, 160564, 410310, 140559, 340423, 370419, 410312, 
  430313, 140562, 150563, 140565, 120567, 260407, 390430, 420322, 160570, 410324, 
  210414, 340508, 120578, 230412, 370452, 540348, 140581, 160579, 370448, 320475, 
  370461, 250415, 540350, 390464, 210417, 210418, 510355, 240421, 370470, 120600, 
  250424, 160603, 370481, 390484, 160610, 430346, 430347, 430348, 450362, 130616, 
  160619, 160621, 340491, 370488, 230427, 320492, 130640, 160636, 160637, 210429, 
  210430, 450360, 140644, 460367, 130652, 140650, 230436, 450373, 540366, 340519, 
  470376, 230439, 340522, 390502, 140671, 230443, 510379, 130682, 540380, 140680, 
  340539, 340548, 370527, 440385, 140689, 320550, 460391, 470387, 130688, 140692, 
  440396, 470398, 540401, 450418, 340562, 540409, 320568, 370559, 540415, 130715, 
  510414, 140717, 160721, 540428, 460422, 540426, 130722, 160729, 110738, 120727, 
  390574, 390576, 390577, 540431, 540433, 540440, 150730, 210462, 370579, 460428, 
  140735, 140737, 210466, 220458, 230463, 420430, 510446, 210465, 210467, 230469, 
  450436, 130741, 130743, 390585, 540455, 540457, 540463, 160749, 540465, 540468, 
  340591, 340592, 340593, 340596, 370590, 540469, 540471, 210481, 220479, 220482, 
  310580, 380588, 450448, 470445, 230497, 340601, 260496, 460465, 120780, 140786, 
  220506, 230507, 140796, 140797, 340608, 510488, 510491, 370613, 440476, 340617, 
  130800, 220516, 520491, 160801, 210524, 410490, 540502, 540504, 540506, 540507, 
  540512, 130806, 160810, 160819, 450503, 460505, 470510, 130815, 140818, 450517, 
  140823, 160828, 260544, 380627, 160825, 340633
)

# 按 dn 筛选
dat$firstepi_nodrug <- ifelse(dat$dn %in% idlist, 1, 0)

outcomes <- c(
  "PANSS_reduce_rate",
  "PANSS_P_reduce_rate",
  "PANSS_N_reduce_rate",
  "PANSS_G_reduce_rate"
)

covars_common <- paste(
  "age", "age2",                  
  "as.factor(sex)", "as.factor(centers)",
  "PC1", "PC2", "PC3", "PC4", "PC5",
  "course", "drugp07", "as.factor(medication)", 
  sep = " + "
)


res_all <- map_dfr(exposures, function(expo) {
  
  gene_name <- sub("^grs_([^_]+).*", "\\1", expo)
  
  map_dfr(outcomes, function(outcome) {
    
    map_dfr(0:1, function(med) {
      
      fml <- as.formula(
        paste(outcome, "~", expo, "+", covars_common)
      )
      
      fit <- tryCatch(
        lm(fml, data = filter(dat, firstepi_nodrug == med)),
        error = function(e) NULL
      )
      if (is.null(fit)) return(NULL)
      
      tidy(fit, conf.int = TRUE) %>%
        filter(term == expo) %>%  
        transmute(
          
          panss_cat = outcome,
          gene      = gene_name,
          med_group = med,
          beta      = estimate,
          lower95   = conf.low,
          upper95   = conf.high,
          p         = p.value
          
        )
    })
  })
})

stra_firstepi_nodrug2 <- res_all
write.csv(stra_firstepi_nodrug2, 
          "/Users/zhuyunqing/Desktop/study/drugMR_lipid/data/onesamplemr/res_grs.APOC3.GCK_2panss_redu.linear.stra_firstepi_nodrug_2g.csv", 
          row.names = FALSE)

########
#stable 19 figure 4




#logistic
outcomes <- c(
  "PANSS_reduce_rate_2g",
  "PANSS_P_reduce_rate_2g",
  "PANSS_N_reduce_rate_2g",
  "PANSS_G_reduce_rate_2g",
  "PANSS_reduce_rate_2g_med",
  "PANSS_P_reduce_rate_2g_med",
  "PANSS_N_reduce_rate_2g_med",
  "PANSS_G_reduce_rate_2g_med"
)

# 执行logistic回归 + 分层 firstepi_nodrug
res_all <- map_dfr(exposures, function(expo) {
  gene_name <- sub("^grs_([^_]+).*", "\\1", expo)
  map_dfr(outcomes, function(outcome) {
    map_dfr(0:1, function(firstepi_val) {
      
      fml <- as.formula(
        paste(outcome, "~", expo, "+", covars_common)
      )
      
      fit <- tryCatch(
        glm(fml, data = filter(dat, firstepi_nodrug == firstepi_val), family = binomial()),
        error = function(e) NULL
      )
      if (is.null(fit)) return(NULL)
      
      tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
        filter(term == expo) %>%
        transmute(
          panss_cat = outcome,
          gene = gene_name,
          firstepi_nodrug = firstepi_val,
          OR = estimate,
          cilower = conf.low,
          ciupper = conf.high,
          p = p.value
        )
    })
  })
})

stra_firstepi_nodrug2_log <- res_all
write.csv(stra_firstepi_nodrug2_log, 
          "/Users/zhuyunqing/Desktop/study/drugMR_lipid/data/onesamplemr/res_grs.APOC3.GCK_2panss_redu.linear.stra_firstepi_nodrug_2g_logistic.csv", 
          row.names = FALSE)

########
#stable 20 

