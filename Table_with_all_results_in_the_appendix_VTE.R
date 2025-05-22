#############################
### TWO SAMPLE MR FOR PQTLS AND VTE
###########################
library(TwoSampleMR)
#read in the credible set data from sun etal 2023
prt_cs = rio::import("/proj/sens2017538/nobackup/Tijana/prots_cs_disc_ukbpp.xlsx")
vte <- vroom::vroom("/proj/sens2017538/nobackup/Tijana/VTE_DKonly_hg19.txt")

colnames(vte)
head(prt_cs)

#clean 
library(tidyr)
library(dplyr)
cpt= prt_cs %>% 
  select(1:10) %>% 
  separate('UKBPPP ProteinID', into = c("UKBPPP", "ProteinID", "Version", "vab"), sep = ":") %>% 
  separate(Variant_ID, into = c("Chromosome", "Position", "Ref_Allele", "Alt_Allele", "Type", "Version"), sep = ":")

#Retain biallelic pQTLs only
cpt_cln = cpt %>% 
  filter(nchar(Ref_Allele) == 1 & nchar(Alt_Allele) == 1)

#Rename variables
wfp = cpt_cln %>% 
  select(-ProteinID, -vab, -Type, -Version) %>% 
  rename(protein = UKBPPP,
         CHR = Chromosome,
         Pos = Position,
         NEA = Ref_Allele,
         EA = Alt_Allele,
         Beta = 'Beta(cond)',
         SE = 'SE(cond)',
         '-log10P' = '-log10(P)(cond)',
         EAF = 'ALT freq',
         cis_trans = 'Cis/trans') %>% 
  mutate(P = 10^(-(`-log10P`)))

head(wfp)

#count the number of trans or cis for each protein
#Gene expression differences between species are driven by both cis and trans effects. Whereas cis effects are caused by genetic variants located on the same DNA molecule as the target gene, trans effects are due to genetic variants that affect diffusible elements. 
prot_xics <-  wfp %>% 
  group_by(protein) %>% 
  count(cis_trans) %>%  
  pivot_wider(names_from = cis_trans, values_from = n) %>% 
  rename(prot_abbrv = protein)

head(prot_xics)

#Select proteins with at least one cis and at least 3 SNPs in total 
useful_prots = prot_xics %>% 
  filter(cis >= 1 & cis + trans >= 3) 

#get the proteins to use
proteins_to_use = wfp[wfp$protein %in% useful_prots$prot_abbrv,]


#Nest by protein to carry out modelling
nested_prots = proteins_to_use %>% 
  nest(.by = protein)

head(nested_prots)


# #prepreation for outcome data(venous thromboembolism)
# vte_new <- vte[,c("var_name","Effect.Meta","Baseline.Meta","Beta.meta", "sdE.meta", "chi.meta", "p.meta", "dir.meta")]
# 
# #Retain biallelic pQTLs only
# vte_cln = vte_new %>% 
#   filter(nchar(Effect.Meta) == 1 & nchar(Baseline.Meta) == 1)
# 
# 
# can <- vte_cln %>% 
#   select(1:8) %>% 
#   separate('var_name', into = c("Chromosome", "Position", "Ref_Allele", "Alt_Allele"), sep = "_")
# 
# 
# can <- can[,c(-5,-6)] 
# names(can) <- c("CHR", "Pos", "NEA" , "EA", "Beta", "SE", "CHi" , "P" , "dir")
# 
# can_new =  can %>%
#   unite(., c("CHR", "Pos"), col = "SNP" ,sep = ":")

#can_new <- merge(can, wfp[, c("CHR", "Pos", "NEA", "EA", "rsID")], by = c("CHR", "Pos", "NEA", "EA"))
#can_new1 <- distinct(can_new)
#Instead of 'cad_new1' and 'can_new', I used vte

df =  wfp %>%
  unite(., c("CHR", "Pos"), col = "SNP" ,sep = ":")

vte2 =  vte %>%
  unite(., c("CHR", "BP"), col = "SNP" ,sep = ":")

#get str data: outcome data also in protein file; warning message: unknown or uninitialised column: 'SNP'.
veins <- vte2[vte2$SNP %in% df$SNP,]



# Since we don't have standard deviation, we used this package to calculate it based on estimates, p-values and the total number of samples used (cases -18, 569; controls - 213,503; 18,569+ 213,503 = 232,072)
# We couldn't calculate it in Bianca since we had to install the package from GitHub, so we ran this part on our own computer
# if (!require("devtools")) {
#   install.packages("devtools")
# }
# devtools::install_github("MathiasHarrer/dmetar")
# library(dmetar)
# 

# # Calculate SE result
# se_result <- se.from.p(effect.size = vte2$BETA, p = vte2$PVAL, N = 232072, effect.size.type = 'difference', calculate.g = FALSE)
# 
# # Add rsID from vte2 to se_result
# se_result$rsID <- vte2$RS
# 
# # Add rsID from vte2 to se_result
# se_result$SNP <- vte2$SNP
# 
# # Save the SE result as a new file with rsID column
# write.csv(se_result, "C:/Users/tikas/Documents/Master's programme/Degree project/Coding/se_result_new_2.csv", row.names = FALSE)



# Read the SE result file
se_data <- read.csv("/proj/sens2017538/nobackup/Tijana/se_result_new_2.csv")

colnames(vte2)[colnames(vte2) == "RS"] <- "rsID"

vte2 <- merge(vte2, se_data[, c("SNP", "StandardError")], by = "SNP", all.x = TRUE)

head(vte2)

#model the nested data
#I changed EA in select() to A1
nested_models  = nested_prots %>% 
  mutate(mr_mods = lapply(data, function(df) {
    #format and harmonise
    
    fmtd = df %>% 
      mutate(phenotype = "protein") %>% 
      TwoSampleMR::format_data(., type = "exposure",
                               phenotype_col = "phenotype",
                               snp_col = "rsID",
                               effect_allele_col = "EA",
                               other_allele_col = "NEA",
                               beta_col = "Beta",
                               se_col = "SE",
                               pval_col = "P",
                               eaf_col = "EAF")
    #get vte data ==> REPLACE THIS WITH YOUR VTE DATA
    vte_prt <- vte2[vte2$rsID %in% df$rsID,]  %>%
      dplyr::select(rsID, A1, A0,BETA, StandardError, PVAL) %>%
      mutate(phenotype = "Venous thromboembolism") %>%
      TwoSampleMR::format_data(., type = "outcome",
                               phenotype_col = "phenotype",
                               snp_col = "rsID",
                               effect_allele_col = "A1", #changed from EA to A1
                               other_allele_col = "A0", #changed from NEA to A0
                               beta_col = "BETA", #changed from Beta to BETA
                               se_col = "StandardError",
                               pval_col = "PVAL")
    #harmonize
    hmz <- TwoSampleMR::harmonise_data(fmtd, vte_prt)
    
    #TSMR
    if(nrow(subset(hmz, mr_keep == TRUE)) > 3) {
      mr_res = TwoSampleMR::mr(hmz, method_list = c("mr_ivw","mr_egger_regression",
                                                    "mr_weighted_median"))%>% 
        data.frame()
      
      pleiotropy = TwoSampleMR::mr_pleiotropy_test(hmz) %>% 
        data.frame() %>% 
        mutate(method = "Egger Intercept", 
               b = egger_intercept,
               nsnp = mr_res$nsnp[1]) %>% 
        select(-egger_intercept)
      
      results = bind_rows(mr_res, pleiotropy)
      #return results
      return(results)
    } else {
      mr_res = TwoSampleMR::mr(hmz, method_list = c("mr_ivw","mr_egger_regression",
                                                    "mr_weighted_median"))%>% 
        data.frame()
      #return results
      return(mr_res)
    }
  }))

head(nested_models)
#save results
saveRDS(nested_models, "prot_vte_TSMR_atleast_1cis.rds")

#unnest results into a dataframe
newmods = nested_models %>% 
  select(1, mr_mods) %>% 
  unnest(., cols = mr_mods)


newmods$OR <- exp(newmods$b)
newmods$CI_lower <- exp(newmods$b - 1.96 * newmods$se)
newmods$CI_upper <- exp(newmods$b + 1.96 * newmods$se)

library(openxlsx)
write.xlsx(newmods, "/proj/sens2017538/nobackup/Tijana/mr_results_vte_proteins.xlsx", rowNames = FALSE)

inverse_vte <- newmods[which(newmods$method == "Inverse variance weighted"),]


#to check how many proteins we have in both oral combined and inverse vte datasets
library(readxl)
oralcombined = read_excel("/proj/sens2017538/nobackup/Tijana/HC_Routesofadmin/top_significant_oralcombined.xlsx")
colnames(oralcombined)[1] <- "protein"
#oralcombined$protein <- as.character(oralcombined$protein)
#inverse_vte$protein <- as.character(inverse_vte$protein)
mergedoralcombined_inverse_vte <- right_join(inverse_vte, oralcombined[,1])


#threshold for bonferoni
threshold_oralcombined_vte <- 0.05/nrow(mergedoralcombined_inverse_vte)
lower_significant_oralcombined_vte <- mergedoralcombined_inverse_vte[which(mergedoralcombined_inverse_vte$pval <= threshold_oralcombined_vte), ]



#implant/injection and inverse vte datasets
library(readxl)
implant_injection = read_excel("/proj/sens2017538/nobackup/Tijana/HC_Routesofadmin/top_significant_implant_injection_progestogen.xlsx")
colnames(implant_injection)[1] <- "protein"
#oralcombined$protein <- as.character(oralcombined$protein)
#inverse_vte$protein <- as.character(inverse_vte$protein)
merged_implant_injection_inverse_vte <- right_join(inverse_vte, implant_injection[,1])


#threshold for bonferoni
threshold_implant_injection_vte <- 0.05/nrow(merged_implant_injection_inverse_vte)
lower_significant_implant_injection_vte <- merged_implant_injection_inverse_vte[which(merged_implant_injection_inverse_vte$pval <= threshold_implant_injection_vte), ]




#local progestin (IUD) and inverse vte datasets
library(readxl)
local_progestin = read_excel("/proj/sens2017538/nobackup/Tijana/HC_Routesofadmin/top_significant_localprogestin.xlsx")
colnames(local_progestin)[1] <- "protein"
#oralcombined$protein <- as.character(oralcombined$protein)
#inverse_vte$protein <- as.character(inverse_vte$protein)
merged_local_progestin_inverse_vte <- right_join(inverse_vte, local_progestin[,1])


#threshold for bonferoni
threshold_local_progestin_vte <- 0.05/nrow(merged_local_progestin_inverse_vte)
lower_significant_local_progestin_vte <- merged_local_progestin_inverse_vte[which(merged_local_progestin_inverse_vte$pval <= threshold_local_progestin_vte), ]




#oral progestin and inverse vte datasets
library(readxl)
oral_progestin = read_excel("/proj/sens2017538/nobackup/Tijana/HC_Routesofadmin/top_significant_oralprogestin.xlsx")
colnames(oral_progestin)[1] <- "protein"
#oralcombined$protein <- as.character(oralcombined$protein)
#inverse_vte$protein <- as.character(inverse_vte$protein)
merged_oral_progestin_inverse_vte <- right_join(inverse_vte, oral_progestin[,1])


#threshold for bonferoni
threshold_oral_progestin_vte <- 0.05/nrow(merged_oral_progestin_inverse_vte)
lower_significant_oral_progestin_vte <- merged_oral_progestin_inverse_vte[which(merged_oral_progestin_inverse_vte$pval <= threshold_oral_progestin_vte), ]



##################################################################################################################################################################
pro_specific <- df[which(df$protein == "ANGPTL3"),]
vte3 <- distinct(vte2)
fmtd = df %>% 
  mutate(phenotype = "protein") %>% 
  TwoSampleMR::format_data(., type = "exposure",
                           phenotype_col = "phenotype",
                           snp_col = "rsID",
                           effect_allele_col = "EA",
                           other_allele_col = "NEA",
                           beta_col = "Beta",
                           se_col = "SE",
                           pval_col = "P",
                           eaf_col = "EAF")
#get vte data ==> REPLACE THIS WITH YOUR VTE DATA
vte_prt <- vte3[vte3$rsID %in% df$rsID,]  %>%
  dplyr::select(rsID, A1, A0,BETA, StandardError, PVAL) %>%
  mutate(phenotype = "Venous thromboembolism") %>%
  TwoSampleMR::format_data(., type = "outcome",
                           phenotype_col = "phenotype",
                           snp_col = "rsID",
                           effect_allele_col = "A1", #changed from EA to A1
                           other_allele_col = "A0", #changed from NEA to A0
                           beta_col = "BETA", #changed from Beta to BETA
                           se_col = "StandardError",
                           pval_col = "PVAL")
hmz <- TwoSampleMR::harmonise_data(fmtd, vte_prt)

#heterogeniety test
heterogeniet_result <- mr_heterogeneity(hmz)


hmz <- hmz[c("SNP","id.exposure","id.outcome","beta.exposure","beta.outcome","se.exposure", "se.outcome")]

newmods = rio::import("/proj/sens2017538/nobackup/Tijana/mr_results_vte_proteins.xlsx")
result_specific <- newmods[which(newmods$protein == "ANGPTL3"),]

intercept <- result_specific[which(result_specific$method == "Egger Intercept"), "b"]
# Extract results for specific methods
result_specific <- subset(result_specific, method %in% c("Inverse variance weighted", "MR Egger", "Weighted median"))

result_specific$lower <- result_specific$b - 1.96 * result_specific$se
result_specific$upper <- result_specific$b + 1.96 * result_specific$se
# Load necessary libraries
library(ggplot2)

# Create a new column for the modified y-values
hmz$y_modified <- sign(hmz$beta.exposure) * hmz$beta.outcome
hmz$x_modified <- abs(hmz$beta.exposure)

# Create the scatter plot with modified y-axis
new_plot6 <- ggplot(hmz, aes(x = x_modified, y = y_modified)) +
  geom_point() +
  geom_errorbar(aes(ymin = y_modified - 1.96 * se.outcome, ymax = y_modified + 1.96 * se.outcome), width = 0) +
  geom_errorbarh(aes(xmin = x_modified - 1.96 * se.exposure, xmax = x_modified + 1.96 * se.exposure), height = 0) +
  labs(
    x = "Effect of SNP on Exposure",
    y = "Effect of SNP on Outcome",
    title = "MR Scatter Plot for ANGPTL3 protein"
  ) +
  theme_minimal() +
  geom_abline(data = subset(result_specific, method == "Inverse variance weighted"), 
              aes(intercept = 0, slope = b, color = method), size = 1) +  # IVW through origin
  geom_abline(data = subset(result_specific, method == "Weighted median"), 
              aes(intercept = 0, slope = b, color = method), size = 1) +  # Weighted median through origin
  geom_abline(data = subset(result_specific, method == "MR Egger"), 
              aes(intercept = intercept, slope = b, color = method), size = 1) +  # Egger with intercept
  scale_color_manual(values = c("Inverse variance weighted" = "#FFB3E6", 
                                "MR Egger" = "tomato", 
                                "Weighted median" = "#00AFBB")) +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(size = 24, face = "bold"),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    legend.text = element_text(size = 20)  # Adjust the legend text font size
  ) +
  xlim(c(0, 0.1))  # Set x-axis limits

# Print the plot
print(new_plot6)

#FOLR1
pro_specific <- df[which(df$protein == "FOLR1"),]
fmtd = df %>% 
  mutate(phenotype = "protein") %>% 
  TwoSampleMR::format_data(., type = "exposure",
                           phenotype_col = "phenotype",
                           snp_col = "rsID",
                           effect_allele_col = "EA",
                           other_allele_col = "NEA",
                           beta_col = "Beta",
                           se_col = "SE",
                           pval_col = "P",
                           eaf_col = "EAF")
#get vte data ==> REPLACE THIS WITH YOUR VTE DATA
vte_prt <- vte2[vte2$rsID %in% df$rsID,]  %>%
  dplyr::select(rsID, A1, A0,BETA, StandardError, PVAL) %>%
  mutate(phenotype = "Venous thromboembolism") %>%
  TwoSampleMR::format_data(., type = "outcome",
                           phenotype_col = "phenotype",
                           snp_col = "rsID",
                           effect_allele_col = "A1", #changed from EA to A1
                           other_allele_col = "A0", #changed from NEA to A0
                           beta_col = "BETA", #changed from Beta to BETA
                           se_col = "StandardError",
                           pval_col = "PVAL")
hmz <- TwoSampleMR::harmonise_data(fmtd, vte_prt)

#heterogeniety test
heterogeniet_result <- mr_heterogeneity(hmz)


hmz <- hmz[c("SNP","id.exposure","id.outcome","beta.exposure","beta.outcome","se.exposure", "se.outcome")]

result_specific <- newmods[which(newmods$protein == "FOLR1"),]

intercept <- result_specific[which(result_specific$method == "Egger Intercept"), "b"]
# Extract results for specific methods
result_specific <- subset(result_specific, method %in% c("Inverse variance weighted", "MR Egger", "Weighted median"))

result_specific$lower <- result_specific$b - 1.96 * result_specific$se
result_specific$upper <- result_specific$b + 1.96 * result_specific$se
# Load necessary libraries
library(ggplot2)

# Create a new column for the modified y-values
hmz$y_modified <- sign(hmz$beta.exposure) * hmz$beta.outcome
hmz$x_modified <- abs(hmz$beta.exposure)

# Create the scatter plot with modified y-axis
new_plot2 <- ggplot(hmz, aes(x = x_modified, y = y_modified)) +
  geom_point() +
  geom_errorbar(aes(ymin = y_modified - 1.96 * se.outcome, ymax = y_modified + 1.96 * se.outcome), width = 0) +
  geom_errorbarh(aes(xmin = x_modified - 1.96 * se.exposure, xmax = x_modified + 1.96 * se.exposure), height = 0) +
  labs(
    x = "Effect of SNP on Exposure",
    y = "Effect of SNP on Outcome",
    title = "MR Scatter Plot for FOLR1 protein"
  ) +
  theme_minimal() +
  geom_abline(data = subset(result_specific, method == "Inverse variance weighted"), 
              aes(intercept = 0, slope = b, color = method), size = 1) +  # IVW through origin
  geom_abline(data = subset(result_specific, method == "Weighted median"), 
              aes(intercept = 0, slope = b, color = method), size = 1) +  # Weighted median through origin
  geom_abline(data = subset(result_specific, method == "MR Egger"), 
              aes(intercept = intercept, slope = b, color = method), size = 1) +  # Egger with intercept
  scale_color_manual(values = c("Inverse variance weighted" = "#FFB3E6", 
                                "MR Egger" = "tomato", 
                                "Weighted median" = "#00AFBB")) +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(size = 24, face = "bold"),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    legend.text = element_text(size = 20)  # Adjust the legend text font size
  ) +
  xlim(c(0, 0.1))  # Set x-axis limits

# Print the plot
print(new_plot2)

#L1CAM
pro_specific <- df[which(df$protein == "L1CAM"),]
fmtd = df %>% 
  mutate(phenotype = "protein") %>% 
  TwoSampleMR::format_data(., type = "exposure",
                           phenotype_col = "phenotype",
                           snp_col = "rsID",
                           effect_allele_col = "EA",
                           other_allele_col = "NEA",
                           beta_col = "Beta",
                           se_col = "SE",
                           pval_col = "P",
                           eaf_col = "EAF")
#get vte data ==> REPLACE THIS WITH YOUR VTE DATA
vte_prt <- vte2[vte2$rsID %in% df$rsID,]  %>%
  dplyr::select(rsID, A1, A0,BETA, StandardError, PVAL) %>%
  mutate(phenotype = "Venous thromboembolism") %>%
  TwoSampleMR::format_data(., type = "outcome",
                           phenotype_col = "phenotype",
                           snp_col = "rsID",
                           effect_allele_col = "A1", #changed from EA to A1
                           other_allele_col = "A0", #changed from NEA to A0
                           beta_col = "BETA", #changed from Beta to BETA
                           se_col = "StandardError",
                           pval_col = "PVAL")
hmz <- TwoSampleMR::harmonise_data(fmtd, vte_prt)

#heterogeniety test
heterogeniet_result <- mr_heterogeneity(hmz)


hmz <- hmz[c("SNP","id.exposure","id.outcome","beta.exposure","beta.outcome","se.exposure", "se.outcome")]

result_specific <- newmods[which(newmods$protein == "L1CAM"),]

intercept <- result_specific[which(result_specific$method == "Egger Intercept"), "b"]
# Extract results for specific methods
result_specific <- subset(result_specific, method %in% c("Inverse variance weighted", "MR Egger", "Weighted median"))

result_specific$lower <- result_specific$b - 1.96 * result_specific$se
result_specific$upper <- result_specific$b + 1.96 * result_specific$se
# Load necessary libraries
library(ggplot2)

# Create a new column for the modified y-values
hmz$y_modified <- sign(hmz$beta.exposure) * hmz$beta.outcome
hmz$x_modified <- abs(hmz$beta.exposure)

# Create the scatter plot with modified y-axis
new_plot2 <- ggplot(hmz, aes(x = x_modified, y = y_modified)) +
  geom_point() +
  geom_errorbar(aes(ymin = y_modified - 1.96 * se.outcome, ymax = y_modified + 1.96 * se.outcome), width = 0) +
  geom_errorbarh(aes(xmin = x_modified - 1.96 * se.exposure, xmax = x_modified + 1.96 * se.exposure), height = 0) +
  labs(
    x = "Effect of SNP on Exposure",
    y = "Effect of SNP on Outcome",
    title = "MR Scatter Plot for L1CAM protein"
  ) +
  theme_minimal() +
  geom_abline(data = subset(result_specific, method == "Inverse variance weighted"), 
              aes(intercept = 0, slope = b, color = method), size = 1) +  # IVW through origin
  geom_abline(data = subset(result_specific, method == "Weighted median"), 
              aes(intercept = 0, slope = b, color = method), size = 1) +  # Weighted median through origin
  geom_abline(data = subset(result_specific, method == "MR Egger"), 
              aes(intercept = intercept, slope = b, color = method), size = 1) +  # Egger with intercept
  scale_color_manual(values = c("Inverse variance weighted" = "#FFB3E6", 
                                "MR Egger" = "tomato", 
                                "Weighted median" = "#00AFBB")) +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(size = 24, face = "bold"),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    legend.text = element_text(size = 20)  # Adjust the legend text font size
  ) +
  xlim(c(0, 0.1))  # Set x-axis limits

# Print the plot
print(new_plot2)


#PDGFRA
pro_specific <- df[which(df$protein == "PDGFRA"),]
fmtd = df %>% 
  mutate(phenotype = "protein") %>% 
  TwoSampleMR::format_data(., type = "exposure",
                           phenotype_col = "phenotype",
                           snp_col = "rsID",
                           effect_allele_col = "EA",
                           other_allele_col = "NEA",
                           beta_col = "Beta",
                           se_col = "SE",
                           pval_col = "P",
                           eaf_col = "EAF")
#get vte data ==> REPLACE THIS WITH YOUR VTE DATA
vte_prt <- vte2[vte2$rsID %in% df$rsID,]  %>%
  dplyr::select(rsID, A1, A0,BETA, StandardError, PVAL) %>%
  mutate(phenotype = "Venous thromboembolism") %>%
  TwoSampleMR::format_data(., type = "outcome",
                           phenotype_col = "phenotype",
                           snp_col = "rsID",
                           effect_allele_col = "A1", #changed from EA to A1
                           other_allele_col = "A0", #changed from NEA to A0
                           beta_col = "BETA", #changed from Beta to BETA
                           se_col = "StandardError",
                           pval_col = "PVAL")
hmz <- TwoSampleMR::harmonise_data(fmtd, vte_prt)

#heterogeneity test
heterogeniet_result <- mr_heterogeneity(hmz)


hmz <- hmz[c("SNP","id.exposure","id.outcome","beta.exposure","beta.outcome","se.exposure", "se.outcome")]

result_specific <- newmods[which(newmods$protein == "PDGFRA"),]

intercept <- result_specific[which(result_specific$method == "Egger Intercept"), "b"]
# Extract results for specific methods
result_specific <- subset(result_specific, method %in% c("Inverse variance weighted", "MR Egger", "Weighted median"))

result_specific$lower <- result_specific$b - 1.96 * result_specific$se
result_specific$upper <- result_specific$b + 1.96 * result_specific$se
# Load necessary libraries
library(ggplot2)

# Create a new column for the modified y-values
hmz$y_modified <- sign(hmz$beta.exposure) * hmz$beta.outcome
hmz$x_modified <- abs(hmz$beta.exposure)

# Create the scatter plot with modified y-axis
new_plot2 <- ggplot(hmz, aes(x = x_modified, y = y_modified)) +
  geom_point() +
  geom_errorbar(aes(ymin = y_modified - 1.96 * se.outcome, ymax = y_modified + 1.96 * se.outcome), width = 0) +
  geom_errorbarh(aes(xmin = x_modified - 1.96 * se.exposure, xmax = x_modified + 1.96 * se.exposure), height = 0) +
  labs(
    x = "Effect of SNP on Exposure",
    y = "Effect of SNP on Outcome",
    title = "MR Scatter Plot for PDGFRA protein"
  ) +
  theme_minimal() +
  geom_abline(data = subset(result_specific, method == "Inverse variance weighted"), 
              aes(intercept = 0, slope = b, color = method), size = 1) +  # IVW through origin
  geom_abline(data = subset(result_specific, method == "Weighted median"), 
              aes(intercept = 0, slope = b, color = method), size = 1) +  # Weighted median through origin
  geom_abline(data = subset(result_specific, method == "MR Egger"), 
              aes(intercept = intercept, slope = b, color = method), size = 1) +  # Egger with intercept
  scale_color_manual(values = c("Inverse variance weighted" = "#FFB3E6", 
                                "MR Egger" = "tomato", 
                                "Weighted median" = "#00AFBB")) +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(size = 24, face = "bold"),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    legend.text = element_text(size = 20)  # Adjust the legend text font size
  ) +
  xlim(c(0, 0.1))  # Set x-axis limits

# Print the plot
print(new_plot2)

#######heterogeneity test for MR egger results
egger_rows <- which(newmods$method == "MR Egger")
egger_pvalues <- newmods[egger_rows,]
#unique_proteins <- data_frame(unique(egger_pvalues$protein))

#extract the MR egger results for the significant proteins
protein <- c("PDGFRA","FOLR1","L1CAM","ANGPTL3")
pro_specific <- egger_pvalues[egger_pvalues$protein %in% protein,]

heterogeniety_result1 <- list()  # Initialize as a list
for (i in 1:nrow(pro_specific)) {  # Loop over all unique proteins
  protein_name <- protein[[i]] # Extract protein name as a scalar
  pro_specific1 <- df[which(df$protein == protein_name), ]  # Filter by protein name
  # Format exposure data
  fmtd = pro_specific1 %>% 
    mutate(phenotype = "protein") %>% 
    TwoSampleMR::format_data(., type = "exposure",
                             phenotype_col = "phenotype",
                             snp_col = "rsID",
                             effect_allele_col = "EA",
                             other_allele_col = "NEA",
                             beta_col = "Beta",
                             se_col = "SE",
                             pval_col = "P",
                             eaf_col = "EAF")
  #get vte data ==> REPLACE THIS WITH YOUR VTE DATA
  vte_prt <- vte2[vte2$rsID %in% df$rsID,]  %>%
    dplyr::select(rsID, A1, A0,BETA, StandardError, PVAL) %>%
    mutate(phenotype = "Venous thromboembolism") %>%
    TwoSampleMR::format_data(., type = "outcome",
                             phenotype_col = "phenotype",
                             snp_col = "rsID",
                             effect_allele_col = "A1", #changed from EA to A1
                             other_allele_col = "A0", #changed from NEA to A0
                             beta_col = "BETA", #changed from Beta to BETA
                             se_col = "StandardError",
                             pval_col = "PVAL")
  # Harmonize data
  hmz <- TwoSampleMR::harmonise_data(fmtd, vte_prt)
  # Perform heterogeneity test and save result in the list
  heterogeniety_result1[[protein_name]] <- mr_heterogeneity(hmz)  # Use `protein_name` as the key
}
heterogeniry_finalresults1 <- do.call(rbind, heterogeniety_result1)

library(tibble)
heterogeniry_finalresults1 <- rownames_to_column(heterogeniry_finalresults1, var = "Protein")

library(writexl)
write_xlsx(heterogeniry_finalresults1, "/proj/sens2017538/nobackup/Tijana/Heterogeneity_vte.xlsx")
