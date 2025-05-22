#stroke
protein1 <- c("IL18BP", "SELE", "FASLG")
protein2 <- data.frame(protein1)

#the same as above for HC
OCP_results <- rio::import("C:/Users/tikas/Documents/Master's programme/Degree project/Coding/whole_results_routeadmin_HC.xlsx")
#check you have the protein column name as "proteins"

colnames(OCP_results)[1] <- 'proteins'

library(dplyr)
sig_proteins_OCP <- OCP_results %>%
  filter(proteins %in% protein1)

# original names of the columns in the code: "common_medications_columnoral_Combined...16","common_medications_columnimplant.injection_progestin...14"; are my changes okay?
imp_OCP_col_beta <- sig_proteins_OCP %>% 
  select("proteins","common_medications_columnimpinjection_progestogens...13","common_medications_columnlocal_progestogens...14","common_medications_columnoral_combined...15","common_medications_columnoral_progestogens...16") %>% 
  rename(b_implant.injection.progestin = common_medications_columnimpinjection_progestogens...13, 
         b_local.progestin = common_medications_columnlocal_progestogens...14,
         b_oral.combined = common_medications_columnoral_combined...15,
         b_oral.progestin = common_medications_columnoral_progestogens...16)


# original names of the columns in the code:"common_medications_columnoral_Combined...27","common_medications_columnimplant.injection_progestin...25"; are my changes okay?
imp_OCP_col_SE <- sig_proteins_OCP %>% 
  select("proteins","common_medications_columnimpinjection_progestogens...23","common_medications_columnlocal_progestogens...24", "common_medications_columnoral_combined...25", "common_medications_columnoral_progestogens...26") %>% 
  rename(se_implant.injection.progestin = common_medications_columnimpinjection_progestogens...23, 
         se_local.progestin = common_medications_columnlocal_progestogens...24,
         se_oral.combined = common_medications_columnoral_combined...25,
         se_oral.progestin = common_medications_columnoral_progestogens...26)



oralcombined_HC <- left_join(imp_OCP_col_beta[,c(1,4)], imp_OCP_col_SE[,c(1,4)])
#oralcombined_HC$OR <- exp(oralcombined_HC$b_oral.combined)
oralcombined_HC$CI_lower <- oralcombined_HC$b_oral.combined - 1.96 * oralcombined_HC$se_oral.combined
oralcombined_HC$CI_upper <- oralcombined_HC$b_oral.combined + 1.96 * oralcombined_HC$se_oral.combined
oralcombined_HC$category <- "oralcombined"
names(oralcombined_HC)[2] <- "Beta"

implant.injrction_progestinonly_HC <- left_join(imp_OCP_col_beta[,1:2], imp_OCP_col_SE[,1:2])
#implant.injrction_progestinonly_HC$OR <- exp(implant.injrction_progestinonly_HC$b_implant.injection.progestin)
implant.injrction_progestinonly_HC$CI_lower <- implant.injrction_progestinonly_HC$b_implant.injection.progestin - 1.96 * implant.injrction_progestinonly_HC$se_implant.injection.progestin
implant.injrction_progestinonly_HC$CI_upper <- implant.injrction_progestinonly_HC$b_implant.injection.progestin + 1.96 * implant.injrction_progestinonly_HC$se_implant.injection.progestin
implant.injrction_progestinonly_HC$category <- "implant.injection_progestinonly"
names(implant.injrction_progestinonly_HC)[2] <- "Beta"

oral_progestinonly_HC <- left_join(imp_OCP_col_beta[,c(1,5)], imp_OCP_col_SE[,c(1,5)])
#oral_progestinonly_HC$OR <- exp(oral_progestinonly_HC$b_oral.progestin)
oral_progestinonly_HC$CI_lower <- oral_progestinonly_HC$b_oral.progestin - 1.96 * oral_progestinonly_HC$se_oral.progestin
oral_progestinonly_HC$CI_upper <- oral_progestinonly_HC$b_oral.progestin + 1.96 * oral_progestinonly_HC$se_oral.progestin
oral_progestinonly_HC$category <- "oral_progestinonly"
names(oral_progestinonly_HC)[2] <- "Beta"

local_progestinonly_HC <- left_join(imp_OCP_col_beta[,c(1,3)], imp_OCP_col_SE[,c(1,3)])
#local_progestinonly_HC$OR <- exp(local_progestinonly_HC$b_local.progestin)
local_progestinonly_HC$CI_lower <- local_progestinonly_HC$b_local.progestin - 1.96 * local_progestinonly_HC$se_local.progestin 
local_progestinonly_HC$CI_upper <- local_progestinonly_HC$b_local.progestin + 1.96 * local_progestinonly_HC$se_local.progestin
local_progestinonly_HC$category <- "local_progestinonly"
names(local_progestinonly_HC)[2] <- "Beta"

HC_results <- rbind(oralcombined_HC[,c(1,2,4:6)],oral_progestinonly_HC[,c(1,2,4:6)],local_progestinonly_HC[,c(1,2,4:6)] ,implant.injrction_progestinonly_HC[,c(1,2,4:6)])

#all_data <- rbind(MHT_results,HC_results)

# I changed all_data to just HC_results
# Convert factors
HC_results$category <- factor(HC_results$category, levels = c("implant.injection_progestinonly", "local_progestinonly", "oral_progestinonly", "oralcombined"))
HC_results$proteins <- factor(HC_results$proteins, levels = protein2$protein1)
HC_results$category <- factor(HC_results$category, levels = rev(unique(HC_results$category)))
HC_results$proteins <- factor(HC_results$proteins, levels = protein2$protein1)

#all_data <- all_data %>%
# mutate(CT_upper = pmin(CT_upper, 3))  # Any value above 3 is set to 3

# Plotting
library(ggplot2)
library(dplyr)

# Ensure row_id is created before plotting
HC_results <- HC_results %>%
  mutate(row_id = as.numeric(factor(proteins)),  # Assign numeric values to proteins
         bg_color = ifelse(row_id %% 2 == 0, "gray90", "white"))  # Alternate row colors

my_forest_plot1 <- ggplot(HC_results, aes(x = Beta, y = factor(proteins, levels = sort(unique(proteins))), xmin = CI_lower, xmax = CI_upper, color = category)) +
  # Add alternating background colors
  geom_rect(aes(ymin = row_id - 0.5, ymax = row_id + 0.5, xmin = -Inf, xmax = Inf, fill = bg_color), 
            inherit.aes = FALSE, alpha = 0.5) +  
  geom_point(size = 3, position = position_dodge(width = 0.5), pch = 18) +
  geom_errorbarh(height = 0.2, position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0, linetype = 5) +
  scale_x_continuous(
    limits = c(min(HC_results$CI_lower, na.rm = TRUE), max(HC_results$CI_upper, na.rm = TRUE)),
    breaks = seq(-1, 0.5, by = 0.1), 
    labels = seq(-1, 0.5, by = 0.1),
    expand = c(0, 0)
  ) +
  labs(x = "Estimate (95% CI)", y = "", title = "") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 12, hjust = 0.5, vjust = 1),
    legend.position = "bottom",
    legend.box = "rect",  
    legend.background = element_rect(color = "black", fill = "white"),  
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 12, margin = margin(r = 10)),
    axis.text.x = element_text(size = 12),
    plot.title = element_text(size = 16, hjust = 0.5),  
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    # Add black border around the plot
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) +
  scale_color_manual(
    name = "Category",  # Capitalized legend title
    values = c(
      "implant.injection_progestinonly" = "#00429d",
      "local_progestinonly" = "#ff9500",
      "oral_progestinonly" = "#007f5f",
      "oralcombined" = "#d90368"
    ),
    labels = c(
      "implant.injection_progestinonly" = "Progestin-Only Implant/Injection",
      "local_progestinonly" = "Local progestin-only",
      "oral_progestinonly" = "Oral progestin-only",
      "oralcombined" = "Oral Combined"
    )
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 4), nrow = 3),
    size = guide_legend(nrow = 1)
  ) +
  scale_fill_identity()  # Ensures the fill color is used directly

# Error: Discrete value supplied to continuous scale
# Print the plot
print(my_forest_plot1)






data1 <- HC_results
##############plot a forest plot for MR results
strok <- rio::import("/Users/tikas/Documents/Master's programme/Degree project/Results/MR VTE, stroke and CAD results/mr_results_stroke_proteins.xlsx")


strok[, 7:13] <- data.frame(lapply(strok[, 7:13], as.numeric))

# Error in exp(strok$b) : non-numeric argument to mathematical function
strok$OR <- exp(strok$b)
strok$CI_lower <- exp(strok$b - 1.96 * strok$se)
strok$CI_upper <- exp(strok$b + 1.96 * strok$se)



sig_proteins <- strok %>%
  filter(protein %in% protein1)


Inverse_variance_weighted <- sig_proteins[which(sig_proteins$method == "Inverse variance weighted"),]
Inverse_variance_weighted <- Inverse_variance_weighted[match(oralcombined_HC$proteins, Inverse_variance_weighted$protein), ]

MR_Egger <- sig_proteins[which(sig_proteins$method == "MR Egger"),]
MR_Egger <- MR_Egger[match(oralcombined_HC$proteins, MR_Egger$protein), ]

Weighted_median <- sig_proteins[which(sig_proteins$method == "Weighted median"),]
Weighted_median <- Weighted_median[match(oralcombined_HC$proteins, Weighted_median$protein), ]

# Combine data
HC_results <- bind_rows(Inverse_variance_weighted, MR_Egger, Weighted_median)
# Convert factors
HC_results$method <- factor(HC_results$method, levels = c("Inverse variance weighted", "MR Egger", "Weighted median"))
HC_results$protein <- factor(HC_results$protein, levels = protein2$protein1)
HC_results$method <- factor(HC_results$method, levels = rev(unique(HC_results$method)))
HC_results$protein <- factor(HC_results$protein, levels = protein2$protein1)



# Ensure row_id is created before plotting
HC_results <- HC_results %>%
  mutate(row_id = as.numeric(factor(protein)),  # Assign numeric values to proteins
         bg_color = ifelse(row_id %% 2 == 0, "gray90", "white"))  # Alternate row colors

# Create forest plot
my_forest_plot <- ggplot(HC_results, aes(x = OR, y =factor(protein, levels = sort(unique(protein))), xmin = CI_lower, xmax = CI_upper, color = method)) +
  # Add alternating background colors
  geom_rect(aes(ymin = row_id - 0.5, ymax = row_id + 0.5, xmin = -Inf, xmax = Inf, fill = bg_color), 
            inherit.aes = FALSE, alpha = 0.5) +  
  geom_point(size = 3, position = position_dodge(width = 0.5), pch = 18) +
  geom_errorbarh(height = 0.2, position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 1, linetype = 5) +
  scale_x_continuous(
    limits = c(min(HC_results$CI_lower, na.rm = TRUE), max(HC_results$CI_upper, na.rm = TRUE)),
    breaks = seq(0, 5.5, by = 0.5), 
    labels = seq(0, 5.5, by = 0.5),
    expand = c(0, 0)
  ) +
  labs(x = "Odds ratio (95% CI)", y = "", title = "") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 12, hjust = 0.5, vjust = 1),
    legend.position = "bottom",
    legend.box = "rect",
    legend.background = element_rect(color = "black", fill = "white"),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 12, margin = margin(r = 10)),
    axis.text.x = element_text(size = 12),
    plot.title = element_text(size = 16, hjust = 0.5),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    # Add black border around the plot
    panel.border = element_rect(color = "black", fill = NA, size = 1)  
  ) +
  scale_color_manual(
    name = "Method",  # Capitalized legend title
    values = c(
      "Inverse variance weighted" = "purple",
      "MR Egger" = "orange",
      "Weighted median" = "#00AFBB"
    )
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 4), nrow = 2),
    size = guide_legend(nrow = 1)
  ) +
  scale_fill_identity()  # Ensure the fill color is used directly

# Print the plot
print(my_forest_plot)



#############put 2 plots together
library(ggplot2)
library(patchwork)
library(cowplot)  # Needed for get_legend()

# Modify the second plot to remove study names
my_forest_plot2 <- my_forest_plot + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# Define a list of your seven plots
plots <- list(my_forest_plot1, my_forest_plot2)


# Empty plot to fill grid space
empty_plot <- ggplot() + theme_void()




library(ggpubr)

# Arrange the plots side by side with legends aligned at the bottom
final_plot <- ggarrange(
  my_forest_plot1, my_forest_plot2,
  ncol = 2,  # Arrange in one row
  common.legend = FALSE,  # Keep individual legends
  legend = "bottom",  # Place legends at the bottom
  align = "hv"  # Align both horizontally and vertically
)

print(final_plot)


library(ggplot2)
library(ggpubr)

# # Modify the first plot with annotation "a"
# my_forest_plot1 <- my_forest_plot1 + 
#   annotate("text", x = -Inf, y = Inf, label = "a", hjust = -0.5, vjust = 1.5, size = 6, fontface = "bold")
# 
# # Modify the second plot with annotation "b"
# my_forest_plot2 <- my_forest_plot2 + 
#   annotate("text", x = -Inf, y = Inf, label = "b", hjust = -0.5, vjust = 1.5, size = 6, fontface = "bold")

# Arrange the plots side by side with legends aligned at the bottom
final_plot <- ggarrange(
  my_forest_plot1, my_forest_plot2,
  ncol = 2,  # Arrange in one row
  common.legend = FALSE,  # Keep individual legends
  legend = "bottom",  # Place legends at the bottom
  align = "hv"  # Align both horizontally and vertically
)

# Print the final plot
print(final_plot)

