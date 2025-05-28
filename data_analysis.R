######################### AC1 Compute #####################################
library(readxl)
library(irrCAC)
library(reshape2)

# Whole brain --------------------------------------------------------
preprocess_data <- function(df){
  rownames(df) <- df$name
  return(as.matrix(df[,-1]))
}

meg <- preprocess_data(read_excel("meg_loc.xlsx"))
seeg <- preprocess_data(read_excel("res_seeg_loc.xlsx"))

# AC1 Compute -----------------------------------------
calculate_gwet_full <- function(method1, method2){
  region_acu <- function(ratings){
    if(all(ratings[,1] == ratings[,2])){
      po <- mean(ratings[,1] == ratings[,2])
      pe <- 2 * (mean(ratings[,1]==0) + mean(ratings[,2]==0))/2 * 
        (1 - (mean(ratings[,1]==0) + mean(ratings[,2]==0))/2)
      return(c(po, pe, (po-pe)/(1-pe), NA, NA))
    }
    res <- gwet.ac1.raw(ratings)
    c(res$est$pa, res$est$pe, res$est$coeff.val, 
      res$est$coeff.se, 2*pnorm(-abs(res$est$coeff.val/res$est$coeff.se)))
  }
  
  all_ratings <- cbind(as.vector(method1), as.vector(method2))
  all_res <- region_acu(all_ratings)
  
  region_res <- t(sapply(1:ncol(method1), function(i){
    region_acu(cbind(method1[,i], method2[,i]))
  }))
  
  results <- data.frame(
    Region = c(colnames(method1), "Whole Brain"),
    Percent_Agreement = c(region_res[,1], all_res[1]),
    Chance_Agreement = c(region_res[,2], all_res[2]),
    AC1 = c(region_res[,3], all_res[3]),
    SE = c(region_res[,4], all_res[4]),
    p_value = c(region_res[,5], all_res[5]),
    stringsAsFactors = FALSE
  )
  
  return(results)
}

results <- calculate_gwet_full(meg, seeg)
results$FDR <- c(p.adjust(results$p_value[1:(nrow(results)-1)], "fdr"), NA)
results$CI_lower <- results$AC1 - 1.96*results$SE
results$CI_upper <- results$AC1 + 1.96*results$SE
print(results[c(nrow(results), 1:(nrow(results)-1)), ])  
write.csv(results, "full_brain_agreement.csv", row.names = FALSE)

# temporal ---------------------------------------------------------------
library(readxl)
library(irrCAC)
library(reshape2)
preprocess_data <- function(df){
  rownames(df) <- df$name
  return(as.matrix(df[,-1]))
}
t_meg <- preprocess_data(read_excel("meg_loc_t.xlsx"))
t_seeg <- preprocess_data(read_excel("res_seeg_loc_t.xlsx"))

# AC1 Compute -----------------------------------------
calculate_gwet_full <- function(method1, method2){
  region_acu <- function(ratings){
    if(all(ratings[,1] == ratings[,2])){
      po <- mean(ratings[,1] == ratings[,2])
      pe <- 2 * (mean(ratings[,1]==0) + mean(ratings[,2]==0))/2 * 
        (1 - (mean(ratings[,1]==0) + mean(ratings[,2]==0))/2)
      return(c(po, pe, (po-pe)/(1-pe), NA, NA))
    }
    res <- gwet.ac1.raw(ratings)
    c(res$est$pa, res$est$pe, res$est$coeff.val, 
      res$est$coeff.se, 2*pnorm(-abs(res$est$coeff.val/res$est$coeff.se)))
  }
  
  all_ratings <- cbind(as.vector(method1), as.vector(method2))
  all_res <- region_acu(all_ratings)

  region_res <- t(sapply(1:ncol(method1), function(i){
    region_acu(cbind(method1[,i], method2[,i]))
  }))
  
  results <- data.frame(
    Region = c(colnames(method1), "Temporal"),
    Percent_Agreement = c(region_res[,1], all_res[1]),
    Chance_Agreement = c(region_res[,2], all_res[2]),
    AC1 = c(region_res[,3], all_res[3]),
    SE = c(region_res[,4], all_res[4]),
    p_value = c(region_res[,5], all_res[5]),
    stringsAsFactors = FALSE
  )
  
  return(results)
}

results <- calculate_gwet_full(t_meg, t_seeg)
results$FDR <- c(p.adjust(results$p_value[1:(nrow(results)-1)], "fdr"), NA)
results$CI_lower <- results$AC1 - 1.96*results$SE
results$CI_upper <- results$AC1 + 1.96*results$SE
print(results[c(nrow(results), 1:(nrow(results)-1)), ])  
write.csv(results, "t_brain_agreement.csv", row.names = FALSE)

# Extemporal ---------------------------------------------------------------
library(readxl)
library(irrCAC)
library(reshape2)
preprocess_data <- function(df){
  rownames(df) <- df$name
  return(as.matrix(df[,-1]))
}
ext_meg <- preprocess_data(read_excel("meg_loc_ext.xlsx"))
ext_seeg <- preprocess_data(read_excel("res_seeg_loc_ext.xlsx"))

# AC1 Compute -----------------------------------------
calculate_gwet_full <- function(method1, method2){
  region_acu <- function(ratings){
    if(all(ratings[,1] == ratings[,2])){
      po <- mean(ratings[,1] == ratings[,2])
      pe <- 2 * (mean(ratings[,1]==0) + mean(ratings[,2]==0))/2 * 
        (1 - (mean(ratings[,1]==0) + mean(ratings[,2]==0))/2)
      return(c(po, pe, (po-pe)/(1-pe), NA, NA))
    }
    res <- gwet.ac1.raw(ratings)
    c(res$est$pa, res$est$pe, res$est$coeff.val, 
      res$est$coeff.se, 2*pnorm(-abs(res$est$coeff.val/res$est$coeff.se)))
  }

  all_ratings <- cbind(as.vector(method1), as.vector(method2))
  all_res <- region_acu(all_ratings)
  
  region_res <- t(sapply(1:ncol(method1), function(i){
    region_acu(cbind(method1[,i], method2[,i]))
  }))

  results <- data.frame(
    Region = c(colnames(method1), "exTemporal"),
    Percent_Agreement = c(region_res[,1], all_res[1]),
    Chance_Agreement = c(region_res[,2], all_res[2]),
    AC1 = c(region_res[,3], all_res[3]),
    SE = c(region_res[,4], all_res[4]),
    p_value = c(region_res[,5], all_res[5]),
    stringsAsFactors = FALSE
  )
  
  return(results)
}

results <- calculate_gwet_full(ext_meg, ext_seeg)
results$FDR <- c(p.adjust(results$p_value[1:(nrow(results)-1)], "fdr"), NA)
results$CI_lower <- results$AC1 - 1.96*results$SE
results$CI_upper <- results$AC1 + 1.96*results$SE
print(results[c(nrow(results), 1:(nrow(results)-1)), ]) 
write.csv(results, "ext_brain_agreement.csv", row.names = FALSE)

############################# Distance Compute ###################################
library(readxl)    
library(dplyr)     
library(ggplot2)   
library(ggpubr)   
library(rstatix)   
library(ggsci)
# distance and agreement--------------------------------------------------------
df_concordance <- read_excel("Patients.xlsx") %>%
  filter(!is.na(Distances)) %>%
  mutate(Concordance = factor(Agreement,   
                              levels = c(0, 1),
                              labels = c("Discordance", "Concordance")))

# statistical analysis------------------
shapiro_test <- df_concordance %>%
  group_by(Concordance) %>%
  summarise(p.value = shapiro.test(Distance)$p.value)
print(shapiro_test)

levene_test <- levene_test(Distance ~ Concordance, data = df_concordance)
print(levene_test)

if (all(shapiro_test$p.value > 0.05) & levene_test$p > 0.05) {
  test_result <- t_test(Distance ~ Concordance, data = df_concordance, detailed = TRUE)
  method_used <- "Independent samples t-test"
} else {
  test_result <- wilcox_test(Distance ~ Concordance, data = df_concordance, detailed = TRUE)
  method_used <- "Mann-Whitney U test"
}

# plot ----------------------------------
library(ggbeeswarm)

test_result <- test_result %>% 
  add_significance("p") %>% 
  add_xy_position(x = "Concordance")

ggplot(df_concordance, aes(x = Concordance, y = Distance, 
                           color = Concordance)) +
  geom_boxplot(width = 0.2, alpha = 0.7, outlier.shape = NA) +
  geom_quasirandom(size = 2.5, alpha = 0.8, width = 0.15) +  
  stat_summary(fun.data = mean_cl_normal, geom = "pointrange",
               color = "black", size = 0.8) +  
  scale_color_manual(values = c("#0072B2", "#E69F00")) +
  labs(title = "Centre-of-mass distance of lesions in the concordant and non-concordant groups") +
  stat_pvalue_manual(test_result, label = "p = {p}{p.signif}") +
  theme_bw(base_size = 22)

########################分区~距离#############################

# distance and region-------------------------------------------
df_region <- read_excel("Patients.xlsx") %>%
  filter(!is.na(Distances)) %>%    
  mutate(Region = factor(Localization,   
                         levels = c(0, 1),
                         labels = c("Extemporal", "Temporal")))


# statistical analysis------------------
shapiro_test <- df_region %>%
  group_by(Region) %>%
  summarise(p.value = shapiro.test(Distance)$p.value)
print(shapiro_test)

levene_test <- levene_test(Distance ~ Region, data = df_region)
print(levene_test)

if (all(shapiro_test$p.value > 0.05) & levene_test$p > 0.05) {
  test_result <- t_test(Distance ~ Region, data = df_region, detailed = TRUE)
  method_used <- "Independent samples t-test"
} else {
  test_result <- wilcox_test(Distance ~ Region, data = df_region, detailed = TRUE)
  method_used <- "Mann-Whitney U test"
}

#plot-----------------------------------
library(ggbeeswarm)

test_result <- test_result %>% 
  add_significance("p") %>% 
  add_xy_position(x = "Concordance")

ggplot(df_region, aes(x = Region, y = Distance, 
                      color = Region)) +
  geom_boxplot(width = 0.2, alpha = 0.7, outlier.shape = NA) +
  geom_quasirandom(size = 2.5, alpha = 0.8, width = 0.15) +  
  stat_summary(fun.data = mean_cl_normal, geom = "pointrange",
               color = "black", size = 0.8) +  #
  scale_color_manual(values = c("#0072B2", "#E69F00")) +
  labs(title = "Centre-of-mass distance of lesions in temporal and extemporal") +
  stat_pvalue_manual(test_result, label = "p = {p}{p.signif}") +
  theme_bw(base_size = 14)


########################Clinical Effectiveness Assessment#############################
library(pacman)
library(readxl)
library(dplyr)
library(epiR)
library(ggplot2)
library(tidyr)
library(broom)
library(binom)
library(kableExtra)
library(gt)
library(webshot2)
# ILAE Compute-------------------------------------------------------------
df_ILAE <- read_excel("Patients.xlsx") %>%
  filter(!is.na(ILAE)) %>%  
  mutate(across(c(Agreement, ILAE), as.factor))
contingency_table_ILAE <- table(df_ILAE$ILAE, df_ILAE$Agreement, 
                                dnn = c("Reference", "Test"))
# TP, TN, FP, FN
TP_ILAE <- contingency_table_ILAE[2, 2]
TN_ILAE <- contingency_table_ILAE[1, 1]
FP_ILAE <- contingency_table_ILAE[1, 2]
FN_ILAE <- contingency_table_ILAE[2, 1]

# 95%CI
sensitivity_ILAE <- TP_ILAE / (TP_ILAE + FN_ILAE)
specificity_ILAE <- TN_ILAE / (TN_ILAE + FP_ILAE)
ppv_ILAE <- TP_ILAE / (TP_ILAE + FP_ILAE)
npv_ILAE <- TN_ILAE / (TN_ILAE + FN_ILAE)
ACC_ILAE <- (TN_ILAE + TP_ILAE) / (TP_ILAE + FP_ILAE + TN_ILAE + FN_ILAE)
sensitivity_ci_ILAE <- binom.test(TP_ILAE, TP_ILAE + FN_ILAE)$conf.int
specificity_ci_ILAE <- binom.test(TN_ILAE, TN_ILAE + FP_ILAE)$conf.int
ppv_ci_ILAE <- binom.test(TP_ILAE, TP_ILAE + FP_ILAE)$conf.int
npv_ci_ILAE <- binom.test(TN_ILAE, TN_ILAE + FN_ILAE)$conf.int
ACC_ci_ILAE <- binom.test(TN_ILAE + TP_ILAE, TP_ILAE + FP_ILAE + TN_ILAE + FN_ILAE)$conf.int

# Engle Compute-----------------------------------------------------------------
df_Engle <- read_excel("Patients.xlsx") %>%
  filter(!is.na(Engle)) %>% 
  mutate(across(c(Agreement, Engle), as.factor))
contingency_table_Engle <- table(df_Engle$Engle, df_Engle$Agreement, 
                                 dnn = c("Reference", "Test"))
# TP, TN, FP, FN
TP_Engle <- contingency_table_Engle[2, 2]
TN_Engle <- contingency_table_Engle[1, 1]
FP_Engle <- contingency_table_Engle[1, 2]
FN_Engle <- contingency_table_Engle[2, 1]

# 95%CI
sensitivity_Engle <- TP_Engle / (TP_Engle + FN_Engle)
specificity_Engle <- TN_Engle / (TN_Engle + FP_Engle)
ppv_Engle <- TP_Engle / (TP_Engle + FP_Engle)
npv_Engle <- TN_Engle / (TN_Engle + FN_Engle)
ACC_Engle <- (TN_Engle + TP_Engle) / (TP_Engle + FP_Engle + TN_Engle + FN_Engle)
sensitivity_ci_Engle <- binom.test(TP_Engle, TP_Engle + FN_Engle)$conf.int
specificity_ci_Engle <- binom.test(TN_Engle, TN_Engle + FP_Engle)$conf.int
ppv_ci_Engle <- binom.test(TP_Engle, TP_Engle + FP_Engle)$conf.int
npv_ci_Engle <- binom.test(TN_Engle, TN_Engle + FN_Engle)$conf.int
ACC_ci_Engle <- binom.test(TN_Engle + TP_Engle, TP_Engle + FP_Engle + TN_Engle + FN_Engle)$conf.int
