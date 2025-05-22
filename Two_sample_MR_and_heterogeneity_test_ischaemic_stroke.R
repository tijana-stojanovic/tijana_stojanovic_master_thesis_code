library(TwoSampleMR)
#read in the credible set data from sun etal 2023
prt_cs = rio::import("/proj/sens2017538/nobackup/Tijana/prots_cs_disc_ukbpp.xlsx")
Strok <- vroom::vroom("/proj/sens2017538/nobackup/Tijana/Meta_IS_BothMandF_GIGA_EA_allfiles_metaAnalyseSEbased_wogc_2021June11_NeNs3rdDp1e6.gwascatalog.tsv_2")
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



df =  wfp %>%
  unite(., c("CHR", "Pos"), col = "SNP" ,sep = ":")

strok = Strok  %>%
  unite(., c("chromosome", "base_pair_location"), col = "SNP" ,sep = ":")

#get str data: outcome data also in protein file
Strok2 <-strok[strok$SNP %in% df$SNP,]
#add snps: for getting rs'ids to outcome file
dsnps = select(df, rsID, SNP)
#add rsids and format
d2 <- left_join(Strok2, dsnps, by = "SNP")
strok_new1 <- distinct(d2)

colnames(strok_new1) <- c("SNP", "EAF", "Beta", "SE", "P", "OR", "ci_lower", "ci_uppwe", "EA", "NEA", "rsID")
#model the nested data
strok_new2 <- strok_new1[,c("SNP", "NEA", "EA", "rsID", "Beta", "SE", "EAF", "P")]

nested_models = nested_prots %>%
  mutate(mr_mods = lapply(data, function(df) {
    
    # Ensure EAF is numeric in both datasets
    df$EAF <- as.numeric(df$EAF)
    strok_new2$EAF <- as.numeric(strok_new2$EAF)
    
    # Format exposure data
    fmtd = df %>%
      mutate(phenotype = "protein") %>%
      TwoSampleMR::format_data(type = "exposure",
                               phenotype_col = "phenotype",
                               snp_col = "rsID",
                               effect_allele_col = "EA",
                               other_allele_col = "NEA",
                               beta_col = "Beta",
                               se_col = "SE",
                               pval_col = "P",
                               eaf_col = "EAF")
    
    # Ensure formatted exposure data is not empty
    if (nrow(fmtd) == 0) return(NULL)
    
    # Filter outcome data
    stroke_prt = strok_new2 %>%
      filter(rsID %in% df$rsID) %>%
      dplyr::select(rsID, EA, NEA, Beta, SE, P, EAF) %>%
      mutate(phenotype = "Ischemic Stroke") %>%
      TwoSampleMR::format_data(type = "outcome",
                               phenotype_col = "phenotype",
                               snp_col = "rsID",
                               effect_allele_col = "EA",
                               other_allele_col = "NEA",
                               beta_col = "Beta",
                               se_col = "SE",
                               pval_col = "P",
                               eaf_col = "EAF")
    
    # Ensure outcome data is not empty
    if (nrow(stroke_prt) == 0) return(NULL)
    
    # Harmonize data
    hmz = TwoSampleMR::harmonise_data(fmtd, stroke_prt)
    
    # Ensure harmonization contains 'mr_keep' and SNPs were retained
    if (!"mr_keep" %in% colnames(hmz) || sum(hmz$mr_keep, na.rm = TRUE) <= 3) return(NULL)
    
    # Run MR analysis
    mr_res = TwoSampleMR::mr(hmz, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")) %>%
      data.frame()
    
    # Pleiotropy test
    pleiotropy = TwoSampleMR::mr_pleiotropy_test(hmz) %>%
      data.frame() %>%
      mutate(method = "Egger Intercept", 
             b = egger_intercept,
             nsnp = mr_res$nsnp[1]) %>%
      select(-egger_intercept)
    
    return(bind_rows(mr_res, pleiotropy))
  }))



head(nested_models)
#save results
#saveRDS(nested_models, "prot_breastcancer_TSMR_atleast_1cis.rds")

#unnest results into a dataframe
newmods = nested_models %>% 
  select(1, mr_mods) %>% 
  unnest(., cols = mr_mods)


newmods$OR <- exp(newmods$b)
newmods$CI_lower <- exp(newmods$b - 1.96 * newmods$se)
newmods$CI_upper <- exp(newmods$b + 1.96 * newmods$se)



inverse_strok <- newmods[which(newmods$method == "Inverse variance weighted"),]

#library(openxlsx)
#write.xlsx(newmods, "/proj/sens2017538/nobackup/Tijana/mr_results_stroke_proteins.xlsx", rowNames = FALSE)

#to check how many proteins we have in both oral combined and inverse stroke datasets
library(readxl)
oralcombined = read_excel("/proj/sens2017538/nobackup/Tijana/HC_Routesofadmin/top_significant_oralcombined.xlsx")
colnames(oralcombined)[1] <- "protein"
#oralcombined$protein <- as.character(oralcombined$protein)
#inverse_strok$protein <- as.character(inverse_strok$protein)
mergedoralcombined_inverse <- right_join(inverse_strok, oralcombined[,1])


#threshold for bonferoni
threshold_oralcombined <- 0.05/nrow(mergedoralcombined_inverse)

lower_significant_oralcombined <- mergedoralcombined_inverse[which(mergedoralcombined_inverse$pval <= threshold_oralcombined), ]



#implant/injection and inverse stroke datasets
library(readxl)
implant_injection = read_excel("/proj/sens2017538/nobackup/Tijana/HC_Routesofadmin/top_significant_implant_injection_progestogen.xlsx")
colnames(implant_injection)[1] <- "protein"
#oralcombined$protein <- as.character(oralcombined$protein)
#inverse_strok$protein <- as.character(inverse_strok$protein)
merged_implant_injection_inverse <- right_join(inverse_strok, implant_injection[,1])


#threshold for bonferoni
threshold_implant_injection <- 0.05/nrow(merged_implant_injection_inverse)

lower_significant_implant_injection <- merged_implant_injection_inverse[which(merged_implant_injection_inverse$pval <= threshold_implant_injection), ]



#local progestin (IUD) and inverse stroke datasets
library(readxl)
local_progestin = read_excel("/proj/sens2017538/nobackup/Tijana/HC_Routesofadmin/top_significant_localprogestin.xlsx")
colnames(local_progestin)[1] <- "protein"
#oralcombined$protein <- as.character(oralcombined$protein)
#inverse_strok$protein <- as.character(inverse_strok$protein)
merged_local_progestin_inverse <- right_join(inverse_strok, local_progestin[,1])


#threshold for bonferoni
threshold_local_progestin <- 0.05/nrow(merged_local_progestin_inverse)

lower_significant_local_progestin <- merged_local_progestin_inverse[which(merged_local_progestin_inverse$pval <= threshold_local_progestin), ]



#oral progestin and inverse stroke datasets
library(readxl)
oral_progestin = read_excel("/proj/sens2017538/nobackup/Tijana/HC_Routesofadmin/top_significant_oralprogestin.xlsx")
colnames(oral_progestin)[1] <- "protein"
#oralcombined$protein <- as.character(oralcombined$protein)
#inverse_strok$protein <- as.character(inverse_strok$protein)
merged_oral_progestin_inverse <- right_join(inverse_strok, oral_progestin[,1])


#threshold for bonferoni
threshold_oral_progestin <- 0.05/nrow(merged_oral_progestin_inverse)

lower_significant_oral_progestin <- merged_oral_progestin_inverse[which(merged_oral_progestin_inverse$pval <= threshold_oral_progestin), ]

##################################################################################################################################################################
pro_specific <- df[which(df$protein == "IL18BP"),]
#vte3 <- distinct(vte2)
fmtd = df %>%
  mutate(phenotype = "protein") %>%
  TwoSampleMR::format_data(type = "exposure",
                           phenotype_col = "phenotype",
                           snp_col = "rsID",
                           effect_allele_col = "EA",
                           other_allele_col = "NEA",
                           beta_col = "Beta",
                           se_col = "SE",
                           pval_col = "P",
                           eaf_col = "EAF")

# Ensure formatted exposure data is not empty
if (nrow(fmtd) == 0) return(NULL)

# Filter outcome data
stroke_prt = strok_new2 %>%
  filter(rsID %in% df$rsID) %>%
  dplyr::select(rsID, EA, NEA, Beta, SE, P, EAF) %>%
  mutate(phenotype = "Ischemic Stroke") %>%
  TwoSampleMR::format_data(type = "outcome",
                           phenotype_col = "phenotype",
                           snp_col = "rsID",
                           effect_allele_col = "EA",
                           other_allele_col = "NEA",
                           beta_col = "Beta",
                           se_col = "SE",
                           pval_col = "P",
                           eaf_col = "EAF")
hmz <- TwoSampleMR::harmonise_data(fmtd, stroke_prt)

#heterogeniety test
heterogeniet_result <- mr_heterogeneity(hmz)


hmz <- hmz[c("SNP","id.exposure","id.outcome","beta.exposure","beta.outcome","se.exposure", "se.outcome")]

newmods = rio::import("/proj/sens2017538/nobackup/Tijana/mr_results_vte_proteins.xlsx")
result_specific <- newmods[which(newmods$protein == "IL18BP"),]

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
new_plot1 <- ggplot(hmz, aes(x = x_modified, y = y_modified)) +
  geom_point() +
  geom_errorbar(aes(ymin = y_modified - 1.96 * se.outcome, ymax = y_modified + 1.96 * se.outcome), width = 0) +
  geom_errorbarh(aes(xmin = x_modified - 1.96 * se.exposure, xmax = x_modified + 1.96 * se.exposure), height = 0) +
  labs(
    x = "Effect of SNP on Exposure",
    y = "Effect of SNP on Outcome",
    title = "MR Scatter Plot for IL18BP protein"
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
print(new_plot1)

#SELE
pro_specific <- df[which(df$protein == "SELE"),]
fmtd = df %>%
  mutate(phenotype = "protein") %>%
  TwoSampleMR::format_data(type = "exposure",
                           phenotype_col = "phenotype",
                           snp_col = "rsID",
                           effect_allele_col = "EA",
                           other_allele_col = "NEA",
                           beta_col = "Beta",
                           se_col = "SE",
                           pval_col = "P",
                           eaf_col = "EAF")

# Ensure formatted exposure data is not empty
if (nrow(fmtd) == 0) return(NULL)

# Filter outcome data
stroke_prt = strok_new2 %>%
  filter(rsID %in% df$rsID) %>%
  dplyr::select(rsID, EA, NEA, Beta, SE, P, EAF) %>%
  mutate(phenotype = "Ischemic Stroke") %>%
  TwoSampleMR::format_data(type = "outcome",
                           phenotype_col = "phenotype",
                           snp_col = "rsID",
                           effect_allele_col = "EA",
                           other_allele_col = "NEA",
                           beta_col = "Beta",
                           se_col = "SE",
                           pval_col = "P",
                           eaf_col = "EAF")
hmz <- TwoSampleMR::harmonise_data(fmtd, stroke_prt)

#heterogeniety test
heterogeniet_result <- mr_heterogeneity(hmz)


hmz <- hmz[c("SNP","id.exposure","id.outcome","beta.exposure","beta.outcome","se.exposure", "se.outcome")]

result_specific <- newmods[which(newmods$protein == "SELE"),]

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
    title = "MR Scatter Plot for SELE protein"
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

#FASLG
pro_specific <- df[which(df$protein == "FASLG"),]
fmtd = df %>%
  mutate(phenotype = "protein") %>%
  TwoSampleMR::format_data(type = "exposure",
                           phenotype_col = "phenotype",
                           snp_col = "rsID",
                           effect_allele_col = "EA",
                           other_allele_col = "NEA",
                           beta_col = "Beta",
                           se_col = "SE",
                           pval_col = "P",
                           eaf_col = "EAF")

# Ensure formatted exposure data is not empty
if (nrow(fmtd) == 0) return(NULL)

# Filter outcome data
stroke_prt = strok_new2 %>%
  filter(rsID %in% df$rsID) %>%
  dplyr::select(rsID, EA, NEA, Beta, SE, P, EAF) %>%
  mutate(phenotype = "Ischemic Stroke") %>%
  TwoSampleMR::format_data(type = "outcome",
                           phenotype_col = "phenotype",
                           snp_col = "rsID",
                           effect_allele_col = "EA",
                           other_allele_col = "NEA",
                           beta_col = "Beta",
                           se_col = "SE",
                           pval_col = "P",
                           eaf_col = "EAF")
hmz <- TwoSampleMR::harmonise_data(fmtd, stroke_prt)

#heterogeniety test
heterogeniet_result <- mr_heterogeneity(hmz)


hmz <- hmz[c("SNP","id.exposure","id.outcome","beta.exposure","beta.outcome","se.exposure", "se.outcome")]

result_specific <- newmods[which(newmods$protein == "FASLG"),]

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
new_plot3 <- ggplot(hmz, aes(x = x_modified, y = y_modified)) +
  geom_point() +
  geom_errorbar(aes(ymin = y_modified - 1.96 * se.outcome, ymax = y_modified + 1.96 * se.outcome), width = 0) +
  geom_errorbarh(aes(xmin = x_modified - 1.96 * se.exposure, xmax = x_modified + 1.96 * se.exposure), height = 0) +
  labs(
    x = "Effect of SNP on Exposure",
    y = "Effect of SNP on Outcome",
    title = "MR Scatter Plot for FASLG protein"
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
print(new_plot3)

#######heterogeneity test for MR egger results
egger_rows <- which(newmods$method == "MR Egger")
egger_pvalues <- newmods[egger_rows,]
unique_proteins <- data_frame(unique(egger_pvalues$protein))

#extract the MR egger results for the significant proteins
protein <- c("SELE","FASLG","IL18BP")
pro_specific <- egger_pvalues[egger_pvalues$protein %in% protein,]

heterogeniety_result1 <- list()  # Initialize as a list
for (i in 1:nrow(pro_specific)) {  # Loop over all unique proteins
  protein_name <- protein[[i]] # Extract protein name as a scalar
  pro_specific1 <- df[which(df$protein == protein_name), ]  # Filter by protein name
  # Format exposure data
  fmtd = pro_specific1 %>%
    mutate(phenotype = "protein") %>%
    TwoSampleMR::format_data(type = "exposure",
                             phenotype_col = "phenotype",
                             snp_col = "rsID",
                             effect_allele_col = "EA",
                             other_allele_col = "NEA",
                             beta_col = "Beta",
                             se_col = "SE",
                             pval_col = "P",
                             eaf_col = "EAF")
  
  # Filter outcome data
  stroke_prt = strok_new2 %>%
    filter(rsID %in% df$rsID) %>%
    dplyr::select(rsID, EA, NEA, Beta, SE, P, EAF) %>%
    mutate(phenotype = "Ischemic Stroke") %>%
    TwoSampleMR::format_data(type = "outcome",
                             phenotype_col = "phenotype",
                             snp_col = "rsID",
                             effect_allele_col = "EA",
                             other_allele_col = "NEA",
                             beta_col = "Beta",
                             se_col = "SE",
                             pval_col = "P",
                             eaf_col = "EAF")
  # Harmonize data
  hmz <- TwoSampleMR::harmonise_data(fmtd, stroke_prt)
  # Perform heterogeneity test and save result in the list
  heterogeniety_result1[[protein_name]] <- mr_heterogeneity(hmz)  # Use `protein_name` as the key
}
heterogeniry_finalresults1 <- do.call(rbind, heterogeniety_result1)

library(tibble)
heterogeniry_finalresults1 <- rownames_to_column(heterogeniry_finalresults1, var = "Protein")


library(writexl)
write_xlsx(heterogeniry_finalresults1, "/proj/sens2017538/nobackup/Tijana/Heterogeneity_stroke.xlsx")
