library(TwoSampleMR)
#read in the credible set data from sun etal 2023
prt_cs = rio::import("/proj/sens2017538/nobackup/Tijana/prots_cs_disc_ukbpp.xlsx")
cad <- vroom::vroom("/proj/sens2017538/nobackup/Tijana/cad.add.160614.website.txt")
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

cad_snp = cad  %>%
  unite(., c("chr", "bp_hg19"), col = "SNP" ,sep = ":")

#get str data: outcome data also in protein file
cad2 <-cad_snp[cad_snp$SNP %in% df$SNP,]

# We had 8 SNPs that were duplicated, so we removed the ones starting with 'chr'
notunique_cad2 <- cad2 %>%
  group_by(SNP) %>%
  filter(n() > 1) %>%
  ungroup()

rm_ch_fromnotunique_cad2 <- notunique_cad2[grepl("^chr", notunique_cad2$markername),]

cad_clean <- anti_join(cad2, rm_ch_fromnotunique_cad2)

cad_new <- distinct(cad_clean)


colnames(cad_new) <- c("rsID", "SNP", "EA", "NEA", "EAF", "median_info", "model", "Beta", "SE", "P", "het_pvalue", "n_studies")

#model the nested data
#cad_new2 <- cad_new[,c("SNP", "NEA", "EA", "rsID", "Beta", "SE", "EAF", "P")]

#analysis
nested_models = nested_prots %>%
  mutate(mr_mods = lapply(data, function(df) {
    
    # Ensure EAF is numeric in both datasets
    df$EAF <- as.numeric(df$EAF)
    cad_new$EAF <- as.numeric(cad_new$EAF)
    
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
    cad_prt = cad_new %>%
      filter(rsID %in% df$rsID) %>%
      dplyr::select(rsID, EA, NEA, Beta, SE, P, EAF) %>%
      mutate(phenotype = "Coronary artery disease") %>%
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
    if (nrow(cad_prt) == 0) return(NULL)
    
    # Harmonize data
    hmz = TwoSampleMR::harmonise_data(fmtd, cad_prt)
    
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
saveRDS(nested_models, "prot_coronary_artery_disease_TSMR_atleast_1cis.rds")

#unnest results into a dataframe
newmods = nested_models %>% 
  select(1, mr_mods) %>% 
  unnest(., cols = mr_mods)


newmods$OR <- exp(newmods$b)
newmods$CI_lower <- exp(newmods$b - 1.96 * newmods$se)
newmods$CI_upper <- exp(newmods$b + 1.96 * newmods$se)



inverse_cad <- newmods[which(newmods$method == "Inverse variance weighted"),]

#library(openxlsx)
write.xlsx(newmods, "/proj/sens2017538/nobackup/Tijana/mr_results_cad_proteins.xlsx", rowNames = FALSE)



inverse_cad <- newmods[which(newmods$method == "Inverse variance weighted"),]


#to check how many proteins we have in both oral combined and inverse cad datasets
library(readxl)
oralcombined = read_excel("/proj/sens2017538/nobackup/Tijana/HC_Routesofadmin/top_significant_oralcombined.xlsx")
colnames(oralcombined)[1] <- "protein"
#oralcombined$protein <- as.character(oralcombined$protein)
#inverse_cad$protein <- as.character(inverse_cad$protein)
mergedoralcombined_inverse_cad <- right_join(inverse_cad, oralcombined[,1])


#threshold for bonferoni
threshold_oralcombined_cad <- 0.05/nrow(mergedoralcombined_inverse_cad)
lower_significant_oralcombined_cad <- mergedoralcombined_inverse_cad[which(mergedoralcombined_inverse_cad$pval <= threshold_oralcombined_cad), ]



#implant/injection and inverse cad datasets
library(readxl)
implant_injection = read_excel("/proj/sens2017538/nobackup/Tijana/HC_Routesofadmin/top_significant_implant_injection_progestogen.xlsx")
colnames(implant_injection)[1] <- "protein"
#oralcombined$protein <- as.character(oralcombined$protein)
#inverse_cad$protein <- as.character(inverse_cad$protein)
merged_implant_injection_inverse_cad <- right_join(inverse_cad, implant_injection[,1])


#threshold for bonferoni
threshold_implant_injection_cad <- 0.05/nrow(merged_implant_injection_inverse_cad)
lower_significant_implant_injection_cad <- merged_implant_injection_inverse_cad[which(merged_implant_injection_inverse_cad$pval <= threshold_implant_injection_cad), ]




#local progestin (IUD) and inverse cad datasets
library(readxl)
local_progestin = read_excel("/proj/sens2017538/nobackup/Tijana/HC_Routesofadmin/top_significant_localprogestin.xlsx")
colnames(local_progestin)[1] <- "protein"
#oralcombined$protein <- as.character(oralcombined$protein)
#inverse_cad$protein <- as.character(inverse_cad$protein)
merged_local_progestin_inverse_cad <- right_join(inverse_cad, local_progestin[,1])


#threshold for bonferoni
threshold_local_progestin_cad <- 0.05/nrow(merged_local_progestin_inverse_cad)
lower_significant_local_progestin_cad <- merged_local_progestin_inverse_cad[which(merged_local_progestin_inverse_cad$pval <= threshold_local_progestin_cad), ]




#oral progestin and inverse cad datasets
library(readxl)
oral_progestin = read_excel("/proj/sens2017538/nobackup/Tijana/HC_Routesofadmin/top_significant_oralprogestin.xlsx")
colnames(oral_progestin)[1] <- "protein"
#oralcombined$protein <- as.character(oralcombined$protein)
#inverse_cad$protein <- as.character(inverse_cad$protein)
merged_oral_progestin_inverse_cad <- right_join(inverse_cad, oral_progestin[,1])


#threshold for bonferoni
threshold_oral_progestin_cad <- 0.05/nrow(merged_oral_progestin_inverse_cad)
lower_significant_oral_progestin_cad <- merged_oral_progestin_inverse_cad[which(merged_oral_progestin_inverse_cad$pval <= threshold_oral_progestin_cad), ]



#######heterogeneity test for MR egger results
egger_rows <- which(newmods$method == "MR Egger")
egger_pvalues <- newmods[egger_rows,]
#unique_proteins <- data_frame(unique(egger_pvalues$protein))

#extract the MR egger results for the significant proteins
protein <- c("SUSD2","HJV","PODXL","CD34","PLA2G7","EPHA4")
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
  cad_prt = cad_new %>%
    filter(rsID %in% df$rsID) %>%
    dplyr::select(rsID, EA, NEA, Beta, SE, P, EAF) %>%
    mutate(phenotype = "Coronary artery disease") %>%
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
  hmz <- TwoSampleMR::harmonise_data(fmtd, cad_prt)
  # Perform heterogeneity test and save result in the list
  heterogeniety_result1[[protein_name]] <- mr_heterogeneity(hmz)  # Use `protein_name` as the key
}
heterogeniry_finalresults1 <- do.call(rbind, heterogeniety_result1)

library(tibble)
heterogeniry_finalresults1 <- rownames_to_column(heterogeniry_finalresults1, var = "Protein")

library(writexl)
write_xlsx(heterogeniry_finalresults1, "/proj/sens2017538/nobackup/Tijana/Heterogeneity_cad.xlsx")
