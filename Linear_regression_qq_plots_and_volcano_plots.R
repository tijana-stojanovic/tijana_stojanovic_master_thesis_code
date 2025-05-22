library(readxl)
fffinal_data  <- read_excel("/proj/sens2017538/nobackup/Tijana/prots_with_pt_batch.xlsx")

fffinal_data1  <- vroom::vroom("/proj/sens2017538/nobackup/Tijana/New_proteins_and_participants.tsv")

fffinal_data2 <- merge(fffinal_data1, fffinal_data , by = "FID")

###exlude Age<50 menopause==1 , Hysterectomy==1 , Bilateral_OOpherectomy==1
new <- subset(fffinal_data2, Age < 50)
new <- new[is.na(new$Menopause) | new$Menopause == 0,]
new <- new[is.na(new$Hysterectomy) | new$Hysterectomy == 0,]
new <- new[is.na(new$Bilateral_oophorectomy) | new$Bilateral_oophorectomy == 0,]
######################exlude medication except OCP and MHT
# Load the readxl package
library(readxl)

# Specify the path to your Excel file
excel_file_path <- "/proj/sens2017538/nobackup/Tijana/my_medication.xlsx"

# Read the Excel file into a data frame
my_medication <- read_excel(excel_file_path)

# Read the medication code from UK biobank
medication_coding  <- vroom::vroom("/proj/sens2017538/nobackup/Tijana/coding4_medication.tsv")
names(medication_coding)[1] <- "Coding_a"
names(my_medication)[2] <- "Coding_a"

OCP_MHT_medication <- merge(my_medication,medication_coding, by="Coding_a")

new_row <- c(Coding_a="99999",
             Category = "Free-text entry, unable to be coded",
             `Medication ATC code`="none",
             `Drug name` = "Free-text entry, unable to be coded",
             meaning = "Free-text entry, unable to be coded")

#warning message: in rbind(deparse.level,...) : number of columns of result, 6, is not a multiple of vector length 5 of arg 2 
OCP_MHT_medication <- rbind(OCP_MHT_medication, new_row) 

OCP_medication_only <- new
##############
med_col <- paste0("medication_", 1:48)
un <- data.frame(Unique_Values = unique(stack(OCP_medication_only[med_col])$values))
names(un) <- "Coding_a"
un2 <- merge(un, OCP_MHT_medication, by= "Coding_a")


########################
operation_AUD <- data.frame(Coding_a = c("P315","Q121","Q122","Q123","Q124","Q128","Q129","S625"),
                            Category = c("Removal of intrauterine contraceptive device from pouch of Douglas","Introduction of intrauterine contraceptive device","Replacement of intrauterine contraceptive device"
                                         ,"Removal of displaced intrauterine contraceptive device NEC","Removal of intrauterine contraceptive device NEC","Other specified intrauterine contraceptive device",
                                         "Unspecified intrauterine contraceptive device", "Removal of hormone implant from subcutaneous tissue"))




operation_frequencies <- numeric(nrow(operation_AUD))
h <- which(colnames(OCP_medication_only) == "operation_1")
c <- which(colnames(OCP_medication_only) == "operation_117")

# Loop through each medication ID to find the frequency for each operative procedure
#check the frequency of each operative procedure we are interested in
for (i in 1:nrow(operation_AUD)) {
  # Extract the medication ID
  op_id <- operation_AUD[i, "Coding_a"]
  
  # Count occurrences in each column from 11 to 60
  operation_frequency <- sum(OCP_medication_only[, h:c] == op_id, na.rm = TRUE)
  
  # Store the frequency
  operation_frequencies[i] <- operation_frequency
}

# Add the frequencies to stroke_meds
operation_AUD$frequency <- operation_frequencies


#loop in the operation columns to find the cells containing operation procedures I want and make then see if it happened in 5 years from the attending
r <- which(colnames(OCP_medication_only) == "operation_date_1")
rr <- which(colnames(OCP_medication_only) == "operation_date_117")


#operation_col <- MHT_medication_only_new[,h:c]
#Operation_dates <- MHT_medication_only_new[,r:rr]
##############################
#check if we have more than 1 operation of interest in each row

OCP_operation_unique <- unique(operation_AUD[,1])
common_values_list <- list()
for (i in 1:nrow(OCP_medication_only)){
  
  data_unique <- unique(OCP_medication_only[i,h:c])
  common_values <- data_unique[data_unique %in% OCP_operation_unique]
  common_values_list[[i]] <- common_values
}


# Create a new column in your dataset
OCP_medication_only$common_operation_column <- NA

# Iterate over each element of common_values_list
for (i in seq_along(common_values_list)) {
  # Check if common_values_list[[i]] is not empty
  if (length(common_values_list[[i]]) > 0) {
    # Convert the list of common values to a single string
    common_values_string <- toString(common_values_list[[i]])
    # Update the corresponding row in the dataset with the common values
    OCP_medication_only[i, "common_operation_column"] <- common_values_string
  }
}
#############################
operation_AUD <- operation_AUD[operation_AUD$frequency != 0 ,]

########
#here, I tried to find the participants who have done any operational procedure within 5 years before initial assessment and if so we put 1 in a column named the same as operational procedure type.
# Initialize columns for each operation
for (op_code in unique(operation_AUD$Coding_a)) {
  OCP_medication_only[paste0("operation_within_5years_", op_code)] <- NA
}

for (op_code in unique(operation_AUD$Coding_a)) {
  OCP_medication_only[paste0("operation_within_5years_date_", op_code)] <- NA
}
# Loop through each row in the operation_AUD data frame
for (i in 1:nrow(operation_AUD)) {
  # Find rows where the operation matches
  op_rows <- OCP_medication_only[, h:c] == operation_AUD[i, "Coding_a"]
  
  # Find row indices where op_rows is TRUE
  op_row_indices <- which(op_rows, arr.ind = TRUE)[, 1]
  
  op_row_col_spec <- data.frame(which(op_rows, arr.ind = TRUE))
  op_row_col_spec$operation_col <- op_row_col_spec$col + h -1
  op_row_col_spec$operation_date_col <- op_row_col_spec$col + r - 1
  
  # Extract dates of attending
  attending_dates <- data.frame(OCP_medication_only[op_row_col_spec$row, "Date_of_attending"])
  # Calculate the difference between operation date and attending date
  for (b in 1:nrow(op_row_col_spec)) {
    row_indi <- op_row_col_spec[b, 1]
    col_indi <- op_row_col_spec[b, 4]
    attending_dates[b, 2] <-  OCP_medication_only[row_indi, col_indi]
  }
  
  names(attending_dates)[1] <- "Date_of_attending"
  names(attending_dates)[2] <- "operation_date_1"
  attending_dates$time_since_op <- attending_dates$Date_of_attending - attending_dates$operation_date_1
  
  # Check if the time since operation is less than 5 years (assuming 365 days in a year)
  attending_dates$within_5_years <- attending_dates$time_since_op < (5 * 365)
  attending_dates$within_5_years <- ifelse(attending_dates$within_5_years == "TRUE", 1, 0)
  attending_dates$within_5_years <- ifelse(attending_dates$time_since_op < 0, 0, attending_dates$within_5_years)
  
  attending_dates <- cbind(op_row_col_spec[,c(1,3,4)], attending_dates)
  
  # Store information for each operation within the row in separate columns
  for (l in 1:nrow(attending_dates)) {
    row_indi <- attending_dates[l, "row"]
    if (attending_dates[l, 7] == 1) {
      OCP_medication_only[row_indi, paste0("operation_within_5years_", operation_AUD[i, "Coding_a"])] <- operation_AUD[i, "Coding_a"]
      OCP_medication_only[row_indi, paste0("operation_within_5years_date_", operation_AUD[i, "Coding_a"])] <- attending_dates[l,"operation_date_1"]
      
    }
  }
}

###########################factor medications based on combined or progestogens only
#medication part
un2$category <- ifelse(grepl("^G03AA|^G03AB|^G03AD", un2$`Medication ATC code`), "Combined", 
                       ifelse(grepl("^G03AC|^G03DA|^G03DC02|^G03DC06", un2$`Medication ATC code`), "progestogens","none"))


un2 <- subset(un2, !grepl("^G03F", un2$`Medication ATC code`))
un2 <- subset(un2, !grepl("^G03C", un2$`Medication ATC code`))
un2 <- subset(un2, !grepl("G03DA04", un2$`Medication ATC code`))


un2$category <- ifelse(grepl("implant", un2$meaning), "implant_progestin",un2$category)
un2$category <- ifelse(grepl("pessary", un2$meaning), "pessary_progestin",un2$category)
un2$category <- ifelse(grepl("injection", un2$meaning), "injection_progestin",un2$category)
un2$category <- ifelse(grepl("product", un2$meaning), "product_progestin",un2$category)
un2$category <- ifelse(grepl("etonogestrel", un2$meaning), "implant_progestin",un2$category)


# Specify the file path
#file_path <- "/proj/sens2017538/nobackup/Tijana/patientsmedication.xlsx"

# Save the data frame with row and column names
#write.xlsx(un2, file = file_path, colNames = TRUE)

matched_data_list <- vector("list", length = 48)

# Loop through each column
for (i in 1:48) {
  # Match the medication column in data1 with data2
  matches <- match(OCP_medication_only[[paste0("medication_", i)]], un2$Coding_a)
  
  # Find non-NA matches in data2$category
  non_na_matches <- !is.na(matches)
  matched_data2 <- un2[matches[non_na_matches], ]
  
  # Update values in data1$medication based on another column in data2
  OCP_medication_only[[paste0("medication_", i)]][non_na_matches] <- matched_data2$category
  
  # Store the matched data for this column
  matched_data_list[[i]] <- matched_data2
}


##################
OCP_medication_unique <- unique(un2$category[2:nrow(un2)])
common_values_list <- list()
for (i in 1:nrow(OCP_medication_only)) {
  
  data_unique <- unique(OCP_medication_only[i,med_col])
  common_values <- data_unique[data_unique %in% OCP_medication_unique]
  common_values_list[[i]] <- common_values
}


# Create a new column in your dataset
OCP_medication_only$common_medications_column <- NA

# Iterate over each element of common_values_list
for (i in seq_along(common_values_list)) {
  # Check if common_values_list[[i]] is not empty
  if (length(common_values_list[[i]]) > 0) {
    # Convert the list of common values to a single string
    common_values_string <- toString(common_values_list[[i]])
    # Update the corresponding row in the dataset with the common values
    OCP_medication_only[i, "common_medications_column"] <- common_values_string
  }
}
##############################
#we have to see if there are any estrogen only individuals that have our interested operations or not, if so change it to combined estrogen and progestogen
MHT_med <- unique(un2$category)
subset_OCP_medication_only <- subset(OCP_medication_only, !is.na(operation_within_5years_Q121) | !is.na(operation_within_5years_Q122))
new_ocp <- subset(subset_OCP_medication_only, select = c("FID","Age","Date_of_attending","pre_OCP","current_OCP","common_medications_column", "operation_within_5years_Q121",
                                                         "operation_within_5years_date_Q121","operation_within_5years_Q122","operation_within_5years_date_Q122"
                                                         ,"operation_within_5years_Q123","operation_within_5years_date_Q123" ,"operation_within_5years_Q124","operation_within_5years_date_Q124"))



new_ocp <- subset(new_ocp, !(!is.na(new_ocp$operation_within_5years_Q121) & !is.na(new_ocp$operation_within_5years_Q124)))
new_ocp <- subset(new_ocp, !(!is.na(new_ocp$operation_within_5years_Q122) & !is.na(new_ocp$operation_within_5years_Q124)))
new_ocp <- subset(new_ocp, !(!is.na(new_ocp$operation_within_5years_Q121) & !is.na(new_ocp$operation_within_5years_Q123)))
new_ocp <- subset(new_ocp, !(!is.na(new_ocp$operation_within_5years_Q122) & !is.na(new_ocp$operation_within_5years_Q123)))

#change those taking Estradiol and have operation Q121 and Q122 as combined estrogen and intrauterine device
new_ocp$common_medications_column[is.na(new_ocp$common_medications_column) & (!is.na(new_ocp$operation_within_5years_Q121) | !is.na(new_ocp$operation_within_5years_Q122))] <- "progestogens1"
new_ocp$common_medications_column[new_ocp$common_medications_column == "Combined" & (!is.na(new_ocp$operation_within_5years_Q121) | !is.na(new_ocp$operation_within_5years_Q122))] <- NA
new_ocp$common_medications_column[new_ocp$common_medications_column == "progestogens" & (!is.na(new_ocp$operation_within_5years_Q121) | !is.na(new_ocp$operation_within_5years_Q122))] <- NA

#this line is only here!
new_ocp$common_medications_column[new_ocp$common_medications_column == "progestogens1"] <-"local_progestogens"


OCP_medication_only$common_medications_column[OCP_medication_only$current_OCP == 0 & is.na(OCP_medication_only$common_medications_column)] <- "none"

for (i in 1:nrow(new_ocp)) {
  m <- which(OCP_medication_only$FID == new_ocp$FID[i])
  OCP_medication_only[m,"common_medications_column"] <- new_ocp[i,"common_medications_column"]
}
####################factor them by numbers now and for those reporting taking OCP without the medication name I have to put them as NA
#OCP_medication_only$common_medications_column[is.na(OCP_medication_only$current_OCP)] <- NA
OCP_medication_only$common_medications_column[OCP_medication_only$current_OCP == 1 & (is.na(OCP_medication_only$common_medications_column) | OCP_medication_only$common_medications_column == "none")] <- NA


med <- OCP_medication_only[,c("FID","current_OCP","common_medications_column")]

OCP_medication_only$common_medications_column[OCP_medication_only$common_medications_column == "progestogens"] <-"oral_progestogens"

OCP_medication_only$common_medications_column[OCP_medication_only$common_medications_column == "injection_progestin" | OCP_medication_only$common_medications_column == "implant_progestin"] <-"impinjection_progestogens"
OCP_medication_only$common_medications_column[OCP_medication_only$common_medications_column == "Combined"] <-"oral_combined"

#############################################################################################
confounders_data <- OCP_medication_only[c("Age","BMI","Smoking","common_medications_column")]
confounders_data$common_medications_column <- factor(confounders_data$common_medications_column)


omics <- OCP_medication_only[,296:3218]

###change the omics names 
library(openxlsx)
library(tidyverse)
coding_pro <- vroom::vroom("/proj/sens2017538/nobackup/Tijana/proteins_name.tsv")

coding_pro <- coding_pro %>%
  separate('meaning', into = c("Protein","meaning"), sep = ";" )



head(coding_pro[1756,])

nf = OCP_medication_only %>% 
  select(1, 296:3218) %>% 
  pivot_longer(., cols = -1, names_to = "code", values_to = "npx")

nf
#add prot names

nf2 = nf %>% 
  rename(coding = code) %>% 
  mutate(across(coding, as.numeric)) %>% 
  left_join(., coding_pro)
head(nf2)

omics = nf2 %>% 
  select(1, npx, Protein) %>% 
  pivot_wider(., names_from = "Protein", values_from = "npx")

head(omics)[,1:5]
omics <- omics[,-1]

#############batch and processing time confounders
time_to_process <- OCP_medication_only[,3219:6141]

#batch <- OCP_medication_only[,6142:9064]


#linear regression analysis
models2 <- list()

for (i in 1:ncol(omics)) {
  protein <- colnames(omics)[i]
  print(paste("working on protein ", protein))
  pt_col <- paste0(protein, "_pt")
  data <- data.frame(
    omics_value = omics[[protein]],
    processing_time = time_to_process[[pt_col]],
    Age = confounders_data$Age,
    BMI = confounders_data$BMI,
    smoking = confounders_data$Smoking,
    common_medications_column = confounders_data$common_medications_column)
  
  data$common_medications_column <- relevel(data$common_medications_column, ref = "none")
  
  # Fit the linear model
  model = lm(data[,1] ~ common_medications_column + Age + I(Age^2) + BMI + smoking + processing_time , data = data)
  
  
  models2[[protein]] <- model
  
} 

# Access and analyze the models
for (i in seq_along(models2)) {
  omics_variable <- colnames(omics)[i]
  print(paste("Model for", omics_variable))
  print(summary(models2[[i]]))
  
}




#extract the p-values
# Initialize an empty list to store p-values
p_values_list <- list()
# Loop through the list of models
for (i in seq_along(models2)) {
  model2 <- models2[[i]]
  p_values <- summary(model2)$coefficients[,4]
  p_values_list[[i]] <- p_values
}

# Convert the list of p-values to a data frame
p_values_df <- data.frame(do.call(rbind, p_values_list))

# Save to an Excel file
#library(openxlsx)

# Specify the file path
#file_path <- "/proj/sens2017538/nobackup/Tijana/P_values.xlsx"

# Save the data frame with row and column names
#write.xlsx(p_values_df, file = file_path, colNames = TRUE)




### 
#library(openxlsx)
#p_values_df1 <- read.xlsx("/proj/sens2017538/nobackup/Tijana/P_values.xlsx")

par(mfrow=c(1,2))
library(qqman)
qq(p_values_df$common_medications_columnimpinjection_progestogens , main='implant/injection_progestogens contraceptives')
qq(p_values_df$common_medications_columnlocal_progestogens , main='local_progestogens contraceptives')
qq(p_values_df$common_medications_columnoral_combined , main='oral combined contraceptives')
qq(p_values_df$common_medications_columnoral_progestogens , main='oral progestogens contraceptives')

#multiple testing
# Assuming your data frame is named my_data
# Assuming your data frame is named my_data
threshold <- 0.05/(length(omics))  # Set your desired threshold here
p_value_omic_oralcombined <- p_values_df["common_medications_columnoral_combined"]
p_value_omic_oralcombined$features <- as.numeric(rownames(p_value_omic_oralcombined))

p_value_omic_implant.injection_progestin <- p_values_df["common_medications_columnimpinjection_progestogens"]
p_value_omic_implant.injection_progestin$features <- as.numeric(rownames(p_value_omic_implant.injection_progestin))

p_value_omic_local_progestin <- p_values_df["common_medications_columnlocal_progestogens"]
p_value_omic_local_progestin$features <- as.numeric(rownames(p_value_omic_local_progestin))

p_value_omic_oral_progestogens <- p_values_df["common_medications_columnoral_progestogens"]
p_value_omic_oral_progestogens$features <- as.numeric(rownames(p_value_omic_oral_progestogens))

# Extract rows higher than the threshold
#Higher_significant_combined <- p_value_omic_combined[p_value_omic_combined$common_medications_columnCombined > threshold , ]
#Higher_significant_features_combined <- Higher_significant_combined[, c("features", "common_medications_columnCombined")]
#Higher_significant_progestogens <- p_value_omic_progestogens[p_value_omic_progestogens$common_medications_columnprogestogens > threshold , ]
#Higher_significant_features_progestogens <- Higher_significant_progestogens[, c("features", "common_medications_columnprogestogens")]

# Extract rows lower than or equal to the threshold
lower_significant_oralcombined <- p_value_omic_oralcombined[p_value_omic_oralcombined$common_medications_columnoral_combined <= threshold, ]
lower_significant_features_oralcombined <- lower_significant_oralcombined[, c("features", "common_medications_columnoral_combined")]

lower__significant_implant.injection_progestin <- p_value_omic_implant.injection_progestin[p_value_omic_implant.injection_progestin$common_medications_columnimpinjection_progestogens <= threshold, ]
lower_significant_features_implant.injection_progestin <- lower__significant_implant.injection_progestin[, c("features", "common_medications_columnimpinjection_progestogens")]

lower__significant_local_progestin <- p_value_omic_local_progestin[p_value_omic_local_progestin$common_medications_columnlocal_progestogens <= threshold, ]
lower_significant_features_local_progestin <- lower__significant_local_progestin[, c("features", "common_medications_columnlocal_progestogens")]

lower__significant_oral_progestogens <- p_value_omic_oral_progestogens[p_value_omic_oral_progestogens$common_medications_columnoral_progestogens <= threshold, ]
lower_significant_features_oral_progestogens <- lower__significant_oral_progestogens[, c("features", "common_medications_columnoral_progestogens")]


library(dplyr)

colnames_omics <- colnames(omics[1:length(omics)])
# Create a new column 'new_column' by extracting names based on 'features'
lower_significant_oralcombined <- lower_significant_oralcombined %>%
  mutate(feature_name = colnames_omics[features])

lower__significant_implant.injection_progestin <- lower__significant_implant.injection_progestin %>%
  mutate(feature_name = colnames_omics[features])

lower__significant_local_progestin <- lower__significant_local_progestin %>%
  mutate(feature_name = colnames_omics[features])

lower__significant_oral_progestogens <- lower__significant_oral_progestogens %>%
  mutate(feature_name = colnames_omics[features])

#extract the coeficient
# Initialize an empty list to store coeficient
coeficient_list <- list()
# Loop through the list of models
for (i in seq_along(models2)) {
  model2 <- models2[[i]]
  coeficient <- summary(model2)$coefficients[, 1]
  coeficient_list[[i]] <- coeficient
}

# Convert the list of coeficient to a data frame
coeficient_df <- data.frame(do.call(rbind, coeficient_list))


# Initialize an empty list to store std
std_list <- list()
# Loop through the list of models
for (i in seq_along(models2)) {
  model2 <- models2[[i]]
  std <- summary(model2)$coefficients[, 2]
  std_list[[i]] <- std
}

# Convert the list of std to a data frame
std_df <- data.frame(do.call(rbind, std_list))

# Specify the file path
#file_path <- "/proj/sens2017538/nobackup/tayebe/New_test/3000protein/OCPanalysis_without_confounders/OCPanalysis_without_excluding_medications/estimates.xlsx"

# Save the data frame with row and column names
#write.xlsx(coeficient_df, file = file_path, colNames = TRUE)



#file_path <- "/proj/sens2017538/nobackup/tayebe/New_test/3000protein/OCPanalysis_without_confounders/OCPanalysis_without_excluding_medications/lower_significant_features_progestogen.xlsx"

# Save the data frame with row and column names
#write.xlsx(lower_significant_features_progestogens, file = file_path, colNames = TRUE)

#file_path <- "/proj/sens2017538/nobackup/tayebe/New_test/3000protein/OCPanalysis_without_confounders/OCPanalysis_without_excluding_medications/lower_significant_features_combined.xlsx"

# Save the data frame with row and column names
#write.xlsx(lower_significant_features_combined, file = file_path, colNames = TRUE)

######################################

omic_col <- data.frame(colnames(omics))
whole_result <- cbind(omic_col, p_values_df,coeficient_df, std_df)


whole_result$common_medications_columnoral_combined <- ifelse(whole_result$common_medications_columnoral_combined == 0.000000e+00 , 5e-324 , whole_result$common_medications_columnoral_combined)
whole_result$common_medications_columnimpinjection_progestogens <- ifelse(whole_result$common_medications_columnimpinjection_progestogens == 0.000000e+00 , 5e-324 , whole_result$common_medications_columnimpinjection_progestogens)
whole_result$common_medications_columnlocal_progestogens <- ifelse(whole_result$common_medications_columnlocal_progestogens == 0.000000e+00 , 5e-324 , whole_result$common_medications_columnlocal_progestogens)
whole_result$common_medications_columnoral_progestogens <- ifelse(whole_result$common_medications_columnoral_progestogens == 0.000000e+00 , 5e-324 , whole_result$common_medications_columnoral_progestogens)


file_path <- "/proj/sens2017538/nobackup/Tijana/HC_Routesofadmin/whole_results_routeadmin_HC.xlsx"
write.xlsx(whole_result, file = file_path, colNames = TRUE)
############################
low_localprogestogens <- whole_result[lower_significant_features_local_progestin$features,]
low_implantinjectionprogestin <- whole_result[lower_significant_features_implant.injection_progestin$features,]
low_oralprogestin <- whole_result[lower_significant_features_oral_progestogens$features,]
low_oralcombined <- whole_result[lower_significant_features_oralcombined$features,]

file_path <- "/proj/sens2017538/nobackup/Tijana/HC_Routesofadmin/top_significant_implant_injection_progestogen.xlsx"
write.xlsx(low_implantinjectionprogestin, file = file_path, colNames = TRUE)

file_path <- "/proj/sens2017538/nobackup/Tijana/HC_Routesofadmin/top_significant_oralcombined.xlsx"
write.xlsx(low_oralcombined, file = file_path, colNames = TRUE)

file_path <- "/proj/sens2017538/nobackup/Tijana/HC_Routesofadmin/top_significant_oralprogestin.xlsx"
write.xlsx(low_oralprogestin, file = file_path, colNames = TRUE)

file_path <- "/proj/sens2017538/nobackup/Tijana/HC_Routesofadmin/top_significant_localprogestin.xlsx"
write.xlsx(low_localprogestogens, file = file_path, colNames = TRUE)



########plot  volcanoplot

oralcombinedhc <- whole_result[,c(1,5,15)]
colnames(oralcombinedhc) <- c("a","pval","b")



oralcombined_th <- 0.05 / 2923

ordered_pvalues <- data.frame(oralcombinedhc[order(oralcombinedhc$pval), "a"])
oralcombinedhc$delabel <- ifelse(oralcombinedhc$a %in% ordered_pvalues$oralcombinedhc.order.oralcombinedhc.pval....a..[1:6] , oralcombinedhc$a, NA)


# Modify the dataframe to create a new categorical variable
oralcombinedhc <- oralcombinedhc %>%
  mutate(expression_status = case_when(
    pval >= oralcombined_th ~ "Not associated with oral combined",
    pval < oralcombined_th ~ "Associated with oral combined"
  ))

oralcombinedhc <- oralcombinedhc %>%
  mutate(label_text = ifelse(expression_status %in% c("Associated with oral combined"), 
                             delabel, NA))  # Assign NA for non-significant proteins

library(ggrepel)

plot1 <- ggplot(data = oralcombinedhc, 
                aes(x = b, y = -log10(pval), col = expression_status, label = delabel)) + 
  geom_point(size = 2, alpha = 0.7, na.rm = TRUE) +  # Ignore NA values
  scale_color_manual(values = c("Not associated with oral combined" = "grey", 
                                "Associated with oral combined" = "#00AFBB")) + 
  coord_cartesian(ylim = c(0, 330), xlim = c(-2.4, 3)) + 
  labs(color = 'Association status', 
       size = '-log10(p-value)',  # Label the size legend
       x = "Estimates", 
       y = "-log10(p-value)") + 
  theme_minimal(base_size = 18) +  # White background (minimal theme)
  theme(
    legend.position = "bottom",  # Position legends at the bottom
    legend.box = "vertical",  # Arrange legends vertically
    axis.title = element_text(size = 24),  # Increase axis title font size
    axis.text = element_text(size = 18),   # Increase axis text font size
    legend.text = element_text(size = 18),  # Increase legend text size
    legend.title = element_text(size = 20),  # Increase legend title size
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # Black border
    panel.background = element_rect(fill = "white", color = NA),  # White background
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  ) +
  geom_hline(yintercept = -log10(oralcombined_th), linetype = "dashed", color = "black", linewidth = 1) +  # Threshold line
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +  # X-axis center line
  geom_label_repel(aes(label = label_text), 
                   max.overlaps = 10, 
                   size = 5, 
                   box.padding = 0.5, 
                   label.padding = 0.3, 
                   fill = "white", 
                   color = "black", 
                   segment.color = "grey50")  

# Print the plot
print(plot1)


#########progestin
implant.injechc <- whole_result[,c(1,3,13)]
colnames(implant.injechc) <- c("protein","pval","b")



implant.injechc_th <- 0.05 / 2923

ordered_pvalues <- data.frame(implant.injechc[order(implant.injechc$pval), "protein"])
implant.injechc$delabel <- ifelse(implant.injechc$protein %in% ordered_pvalues$implant.injechc.order.implant.injechc.pval....protein..[1:6] , implant.injechc$protein, NA)


# Modify the dataframe to create a new categorical variable
implant.injechc <- implant.injechc %>%
  mutate(expression_status = case_when(
    pval >= implant.injechc_th ~ "Not associated with implant/injection HC",
    pval < implant.injechc_th ~ "Associated with implant/injection HC"
  ))

implant.injechc <- implant.injechc %>%
  mutate(label_text = ifelse(expression_status %in% c("Associated with implant/injection HC"), 
                             delabel, NA))  # Assign NA for non-significant proteins


plot2 <- ggplot(data = implant.injechc, 
                aes(x = b, y = -log10(pval), col = expression_status, label = delabel)) + 
  geom_point(size = 2, alpha = 0.7, na.rm = TRUE) +  # Ignore NA values
  scale_color_manual(values = c("Not associated with implant/injection HC" = "grey", 
                                "Associated with implant/injection HC" = "#FF69B4")) + 
  coord_cartesian(ylim = c(0, 20), xlim = c(-2.5, 1)) + 
  labs(color = 'Association status', 
       size = '-log10(p-value)',  # Label the size legend
       x = "Estimates", 
       y = "-log10(p-value)") + 
  theme_minimal(base_size = 18) +  # White background (minimal theme)
  theme(
    legend.position = "bottom",  # Position legends at the bottom
    legend.box = "vertical",  # Arrange legends vertically
    axis.title = element_text(size = 24),  # Increase axis title font size
    axis.text = element_text(size = 18),   # Increase axis text font size
    legend.text = element_text(size = 18),  # Increase legend text size
    legend.title = element_text(size = 20),  # Increase legend title size
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # Black border
    panel.background = element_rect(fill = "white", color = NA),  # White background
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  ) +
  geom_hline(yintercept = -log10(implant.injechc_th), linetype = "dashed", color = "black", linewidth = 1) +  # Threshold line
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +  # X-axis center line
  geom_label_repel(aes(label = label_text), 
                   max.overlaps = 10, 
                   size = 5, 
                   box.padding = 0.5, 
                   label.padding = 0.3, 
                   fill = "white", 
                   color = "black", 
                   segment.color = "grey50")  

# Print the plot
print(plot2)

##########oralprogestin
oralcombinedmht <- whole_result[,c(1,6,16)]
colnames(oralcombinedmht) <- c("protein","pval","b")




oralcombinedmht_th <- 0.05 / 2923

ordered_pvalues <- data.frame(oralcombinedmht[order(oralcombinedmht$pval), "protein"])
oralcombinedmht$delabel <- ifelse(oralcombinedmht$protein %in% ordered_pvalues$oralcombinedmht.order.oralcombinedmht.pval....protein..[1:6] , oralcombinedmht$protein, NA)


oralcombinedmht <- oralcombinedmht %>%
  mutate(expression_status = case_when(
    pval >= oralcombinedmht_th ~ "Not associated with oral progestin",
    pval < oralcombinedmht_th ~ "Associated with oral progestin"
  ))

oralcombinedmht <- oralcombinedmht %>%
  mutate(label_text = ifelse(expression_status %in% c("Associated with oral progestin"), 
                             delabel, NA))  # Assign NA for non-significant proteins


plot3 <- ggplot(data = oralcombinedmht, 
                aes(x = b, y = -log10(pval), col = expression_status, label = delabel)) + 
  geom_point(size = 2, alpha = 0.7, na.rm = TRUE) +  # Ignore NA values
  scale_color_manual(values = c("Not associated with oral progestin" = "grey", 
                                "Associated with oral progestin" = "purple")) + 
  coord_cartesian(ylim = c(0, 37), xlim = c(-2, 0.5)) + 
  labs(color = 'Association status', 
       size = '-log10(p-value)',  # Label the size legend
       x = "Estimates", 
       y = "-log10(p-value)") + 
  theme_minimal(base_size = 18) +  # White background (minimal theme)
  theme(
    legend.position = "bottom",  # Position legends at the bottom
    legend.box = "vertical",  # Arrange legends vertically
    axis.title = element_text(size = 24),  # Increase axis title font size
    axis.text = element_text(size = 18),   # Increase axis text font size
    legend.text = element_text(size = 18),  # Increase legend text size
    legend.title = element_text(size = 20),  # Increase legend title size
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # Black border
    panel.background = element_rect(fill = "white", color = NA),  # White background
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  ) +
  geom_hline(yintercept = -log10(oralcombinedmht_th), linetype = "dashed", color = "black", linewidth = 1) +  # Threshold line
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +  # X-axis center line
  geom_label_repel(aes(label = label_text), 
                   max.overlaps = 10, 
                   size = 5, 
                   box.padding = 0.5, 
                   label.padding = 0.3, 
                   fill = "white", 
                   color = "black", 
                   segment.color = "grey50")  


# Print the plot
print(plot3)
#########local progestin
oralestrogenmht <- whole_result[,c(1,4,14)]
colnames(oralestrogenmht) <- c("protein","pval","b")


oralestrogenmht_th <- 0.05 / 2923

ordered_pvalues <- data.frame(oralestrogenmht[order(oralestrogenmht$pval), "protein"])
oralestrogenmht$delabel <- ifelse(oralestrogenmht$protein %in% ordered_pvalues$oralestrogenmht.order.oralestrogenmht.pval....protein..[1:6] , oralestrogenmht$protein, NA)


# Modify the dataframe to create a new categorical variable
oralestrogenmht <- oralestrogenmht %>%
  mutate(expression_status = case_when(
    pval >= oralestrogenmht_th ~ "Not associated with local progestin",
    pval < oralestrogenmht_th ~ "Associated with local progestin"
  ))

oralestrogenmht <- oralestrogenmht %>%
  mutate(label_text = ifelse(expression_status %in% c("Associated with local progestin"), 
                             delabel, NA))  # Assign NA for non-significant proteins


plot4 <- ggplot(data = oralestrogenmht, 
                aes(x = b, y = -log10(pval), col = expression_status, label = delabel)) + 
  geom_point(size = 2, alpha = 0.7, na.rm = TRUE) +  # Ignore NA values
  scale_color_manual(values = c("Not associated with local progestin" = "grey", 
                                "Associated with local progestin" = "#bb0c00")) + 
  coord_cartesian(ylim = c(0, 7), xlim = c(-0.5, 0.5)) + 
  labs(color = 'Association status', 
       size = '-log10(p-value)',  # Label the size legend
       x = "Estimates", 
       y = "-log10(p-value)") + 
  theme_minimal(base_size = 18) +  # White background (minimal theme)
  theme(
    legend.position = "bottom",  # Position legends at the bottom
    legend.box = "vertical",  # Arrange legends vertically
    axis.title = element_text(size = 24),  # Increase axis title font size
    axis.text = element_text(size = 18),   # Increase axis text font size
    legend.text = element_text(size = 18),  # Increase legend text size
    legend.title = element_text(size = 20),  # Increase legend title size
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # Black border
    panel.background = element_rect(fill = "white", color = NA),  # White background
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  ) +
  geom_hline(yintercept = -log10(oralestrogenmht_th), linetype = "dashed", color = "black", linewidth = 1) +  # Threshold line
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +  # X-axis center line
  geom_label_repel(aes(label = label_text), 
                   max.overlaps = 10, 
                   size = 5, 
                   box.padding = 0.5, 
                   label.padding = 0.3, 
                   fill = "white", 
                   color = "black", 
                   segment.color = "grey50")  


# Print the plot
print(plot4)

########make them one plot
library(patchwork)  # Ensure you have this package for arranging plots
library(ggpubr)

# Modify the first plot with annotation "a"
plot1 <- plot1 + 
  annotate("text", x = -Inf, y = Inf, label = "a", hjust = -0.5, vjust = 1.5, size = 6, fontface = "bold")

# Modify the second plot with annotation "b"
plot2 <- plot2 + 
  annotate("text", x = -Inf, y = Inf, label = "b", hjust = -0.5, vjust = 1.5, size = 6, fontface = "bold")

plot3 <- plot3 + 
  annotate("text", x = -Inf, y = Inf, label = "c", hjust = -0.5, vjust = 1.5, size = 6, fontface = "bold")

plot4 <- plot4 + 
  annotate("text", x = -Inf, y = Inf, label = "d", hjust = -0.5, vjust = 1.5, size = 6, fontface = "bold")

# Remove legends from all plots
plot11 <- plot1 + theme(axis.title.x = element_blank(), legend.position = "none")
plot21 <- plot2 + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), legend.position = "none")
plot31 <- plot3 + theme(axis.title.y = element_blank(), legend.position = "none")  
plot41 <- plot4 + theme(legend.position = "none")  # Keep y-axis title

# Adjust axis labels
plot11 <- plot11 + labs(y = "-log10(p-value)")  
plot31 <- plot31 + labs(x = "Estimates", y = "-log10(p-value)")  # ✅ Keep y-axis label
plot41 <- plot41 + labs(x = "Estimates")  

# Create an empty plot for alignment
empty_plot <- ggplot() + theme_void()

# Arrange plots in a 3x2 layout with equal width and height across rows
combined_plot <- (plot11 | plot21) / 
  (plot31 | plot41) +
  plot_layout(widths = c(1, 1, 1), heights = c(1, 1))  # Equal widths and heights


# Display the combined plot
print(combined_plot)






