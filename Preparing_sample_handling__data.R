fffinal_data  <- vroom::vroom("/proj/sens2017538/nobackup/Tijana/New_proteins_and_participants.tsv")

#this data has plateID as a column that we need
PlatID <- vroom::vroom("/proj/sens2017538/proj_41143/ukb674499.tab", col_select = c(f.eid, "f.30900.0.0","f.30901.0.0","f.30902.0.0"))

#extract the column of PlatID with the individual ID
PlatID_data <- PlatID[,c("f.eid", "f.30900.0.0","f.30901.0.0","f.30902.0.0")]
c <- c("FID", "Number_of_proteins_measured", "PlateID", "well_used")
colnames(PlatID_data) <- c
PlatID_data$PlateID <- as.numeric(PlatID_data$PlateID)




library(readxl)
#load batch number data
excel_file_path <- "/proj/sens2017538/nobackup/Tijana/Olink_batchnumber.xlsx"

# Read the Excel file into a data frame
Olink_batchnumber <- read_excel(excel_file_path, col_types = "text")


#load processing time data
excel_file_path <- "/proj/sens2017538/nobackup/Tijana/ProteomicsProcessing_date.xlsx"

# Read the Excel file into a data frame
Proteomics_proceesing_data <- read_excel(excel_file_path)

Proteomics_proceesing_data$PlateID <- as.numeric(Proteomics_proceesing_data$PlateID)




merged_data1 <- merge(PlatID_data, Olink_batchnumber, by = "PlateID", all.x = TRUE)
merged_data1 <- merge(merged_data1, PlatID_data, by = "PlateID", all.x = TRUE)


##################################assay,assay version and lot number
excel_file_path <- "/proj/sens2017538/nobackup/Tijana/assay.xlsx"

# Read the Excel file into a data frame
assay <- vroom::vroom(excel_file_path)




