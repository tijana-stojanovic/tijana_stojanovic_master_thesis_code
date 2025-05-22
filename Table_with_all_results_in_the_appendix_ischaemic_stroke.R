HC= rio::import("/proj/sens2017538/nobackup/Tijana/HC_Routesofadmin/whole_results_routeadmin_HC.xlsx")


colnames(HC)[1] <- "protein"

library(tidyr)
library(dplyr)
HC= HC %>%
  separate('protein', into = c("protein", "ProteinID"), sep = ";")

#colnames(HC) <- HC[2,] 
#HC <- HC[-c(1,2),]
#colnames(HC)[1] <- "Proteins"
#colnames(HC)[2] <- "ProteinID"
proteins_to_use <- c("PTPRK", "PTPRM", "CST6")




new <- HC[HC$protein %in% proteins_to_use, ]

new1 <- new[,c(1,4:7,15:18)]


newmodes= rio::import("/proj/sens2017538/nobackup/Tijana/mr_results_stroke_proteins.xlsx")

MRresults <- newmods[which(newmods$method == "Inverse variance weighted"),]

MRresults$OR <- exp(MRresults$b)
MRresults$CI_lower <- exp(MRresults$b - 1.96 * MRresults$se)
MRresults$CI_upper <- exp(MRresults$b + 1.96 * MRresults$se)

MRresults <- MRresults[,c("protein","OR","CI_lower","CI_upper")]

maindata <- merge(new1,MRresults, by = "protein")


library(writexl)

# Example: Save dataframe `mydata` to Excel file
write_xlsx(maindata, "Stroketable.xlsx")