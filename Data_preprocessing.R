load('/proj/sens2017538/nobackup/UKBB_41143_Data/All.Phenos.20210212.Females_WiteCau.RData')

wh_w <- temp[,c("f.eid","f.21000.0.0","f.21022.0.0", "f.53.0.0","f.21001.0.0","f.3581.0.0","f.48.0.0", "f.49.0.0", "f.2784.0.0", "f.2724.0.0", "f.20116.0.0", "f.3591.0.0", "f.2834.0.0")]

library(dplyr)
wh_w <- wh_w %>% mutate(WHR = f.48.0.0 / f.49.0.0)
white_woman <- cbind(wh_w[,1:6],wh_w[,14],wh_w[9:13])


colnames(white_woman)<-c("FID","ethinicity","Age", "Date_of_attending","BMI","WHR","Age_at_menopause","oral contraceptive pill","Menopause","Smoking", "Hysterectomy","Bilateral_oophorectomy")
colnames(white_woman)


med_col <- select(temp, starts_with("f.20003.0"))
name_col <- paste0("medication_", 1:ncol(med_col))
names(med_col) <- name_col

operation_col <- select(temp, starts_with("f.41272.0"))
name_col <- paste0("operation_" , 1:ncol(operation_col))
names(operation_col) <- name_col

operation_dates_col <- select(temp, starts_with("f.41282.0"))
name_col <- paste0("operation_date_" , 1:ncol(operation_dates_col))
names(operation_dates_col) <- name_col

white_woman <- cbind(white_woman, med_col, operation_col, operation_dates_col)


###########################################################################################
#load proteomics data
Data <- read.table("/proj/sens2017538/nobackup/UKBB_41143_Data/olink_data_41143_20231122.txt",  header = TRUE)

library(tidyverse)
Data_wide <- pivot_wider(Data, names_from = protein_id, values_from = result)

print(unique(Data_wide$ins_index))

data_filtered <- subset(Data_wide, ins_index != 2)
data_filtered <- subset(data_filtered, ins_index != 3)

print(unique(data_filtered$ins_index))


proteomics <- data_filtered[, -which(names(data_filtered) == "ins_index")]
head(proteomics)

colnames(proteomics)[colnames(proteomics) == "eid"] <- "FID"


merged_dataset <- merge(white_woman,proteomics, by = "FID")

##########################################################################################
# (Yes: 1, No: 0, Not sure: -5, Prefer not to answer: -3) put -3 and -5 as NA for Bilateral_oophorectomy
merged_dataset$Bilateral_oophorectomy <- ifelse(merged_dataset$Bilateral_oophorectomy == "-3", NA, merged_dataset$Bilateral_oophorectomy)
merged_dataset$Bilateral_oophorectomy <- ifelse(merged_dataset$Bilateral_oophorectomy == "-5", NA, merged_dataset$Bilateral_oophorectomy)

# (Yes: 1, No: 0, Not sure - had a hysterectomy: 2, Not sure - other reasons: 3, Prefer not to answer: -3) put -3 and 3 and 2 as NA for Menopause
merged_dataset$Menopause <- ifelse(merged_dataset$Menopause == "2", NA, merged_dataset$Menopause)
merged_dataset$Menopause <- ifelse(merged_dataset$Menopause == "3", NA, merged_dataset$Menopause)
merged_dataset$Menopause <- ifelse(merged_dataset$Menopause == "-3", NA, merged_dataset$Menopause)

# (Prefer not to answer: -3, Never: 0, Previous: 1, current: 2) put previous users as non-user
merged_dataset$Smoking <- ifelse(merged_dataset$Smoking == "-3", NA, merged_dataset$Smoking)
merged_dataset$Smoking <- ifelse(merged_dataset$Smoking == "1", 0, merged_dataset$Smoking)
merged_dataset$Smoking <- ifelse(merged_dataset$Smoking == "2", 1, merged_dataset$Smoking)

# (Yes: 1, No: 0, Not sure: -5, Prefer not to answer: -3) put -3 and -5 as NA for Hysterectomy
merged_dataset$Hysterectomy <- ifelse(merged_dataset$Hysterectomy == "-3", NA, merged_dataset$Hysterectomy)
merged_dataset$Hysterectomy <- ifelse(merged_dataset$Hysterectomy == "-5", NA, merged_dataset$Hysterectomy)
unique(merged_dataset$Menopause)

##Contraceptive use
OC<-temp$f.2784.0.0 ## this crates a new variable with OC ever use info
OC[which(temp$f.2784.0.0 == -1)]<-NA
#To be continued, just happened to press enter...
OC[which(temp$f.2784.0.0 == -3)]<-NA
OC[which(temp$f.2784.0.0 == 1)]<-"previous"
OC[which(temp$f.2784.0.0 == 0)]<-"never"
#Now you have basic info on ever/never
#TO also add current use you need the second variable:
OC[which(temp$f.2804.0.0 == -11)]<-"current"
OC[which(temp$f.2804.0.0 == -3)]<-NA
ocp <- data.frame(OC)

ocp$OC <- ifelse(ocp$OC == "current", 1, ocp$OC)
ocp$OC <- ifelse(ocp$OC == "never", 0, ocp$OC)
ocp$OC <- ifelse(ocp$OC == "previous", 0, ocp$OC)
final_oc <- cbind(temp[,1],ocp)
names(final_oc)[1] <- "FID"

########################################################
#final_data <- merge(merged_dataset,final_mht, by="FID")
final_data <- merge(merged_dataset,final_oc, by="FID")
final_data <- final_data[,c(1:7,3218,8:3217)]
View(final_data)
names(final_data)[9] <- "pre_OCP"
names(final_data)[8] <- "current_OCP"
#names(final_data)[8] <- "pre_MHT"
#names(final_data)[9] <- "current_MHT"

#table(final_data$current_MHT, useNA = "ifany")
table(final_data$current_OCP, useNA = "ifany")

#######################################################

file_path <- "/proj/sens2017538/nobackup/Tijana/New_proteins_and_participants.tsv"

# Save as a TSV file
write.table(final_data, file = file_path, sep = "\t", quote = FALSE, row.names = FALSE)

