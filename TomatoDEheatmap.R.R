#Pseudo:
#Create a heatmap to track logFC over four time points for each genotype
#Since excel file isn't in .csv format, and it contains data in separate spreadsheets, make a new file combining only the relevant columns
#Keep row and column orders
#Disperse colors with mid = val 0 = white, follow heat indicators
#Rename logFC colum to ID value per time point
#Columns to be carried over to the new sheet: gene, logFC 15m, logFC 30m, logFC 45m, logFC 75m
#One excel file per genotype to avoid join errors
#Figures should contain headers, title, log and heat legend
#According to Anthony's analysis relevant DE genes are distributed as: VF36 = 54 genes, OE3 = 55 genes, are = 129 genes
#Column Q = adjusted p-value, filtered by <=0.049
#Files names: are_CvT_DEgenes.xlsx , OE3_CvT_DEgenes.xlsx , VF36_CvT_DEgenes.xlsx
#Session path for test: ~Desktop/HotTot
#R project: HotTot
#Bitbucket repository: https://bitbucket.org/hotpollen/flavonoid-rnaseq/src/main/72_F3H_PollenTube/




#Install and load libraries 
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("ComplexHeatmap")
install.packages("openxlsx")

library(ComplexHeatmap)
library(openxlsx)
library(circlize)
library(grid)

#are
#Read excel file
are_data <- read.xlsx("are_CvT_DEgenes.xlsx", sheet = "are_Rel_Data")

#Set gene column as row name
rownames(are_data) <- are_data$gene
are_data_matrix <- as.matrix(are_data[, c("logFC_15min", "logFC_30min", "logFC_45min", "logFC_75min")])

#Rename the columns of the matrix
colnames(are_data_matrix) <- c("logFC 15min", "logFC 30min", "logFC 45min", "logFC 75min")

#QC data correctly loaded and transformed
print(str(are_data_matrix))

#Replace NA values with 0
are_data_matrix[is.na(are_data_matrix)] <- 0

#Define two gradients
blue_to_white <- colorRampPalette(c("steelblue3", "white"))(50)
white_to_red <- colorRampPalette(c("white", "tomato2"))(50)

#Combine the gradients
excel_colors <- c(blue_to_white, white_to_red[-1])

#Define breaks
break_interval1 <- (-3 - 0) / 49
break_interval2 <- (8 - 0) / 49

blue_white_breaks <- seq(-3, 0, length.out = 50)
white_red_breaks <- seq(0, 8, length.out = 50)[2:50]

excel_breaks <- c(blue_white_breaks, white_red_breaks)

#Create color mapping
color_fun <- circlize::colorRamp2(breaks = excel_breaks, colors = excel_colors)

#Draw and save heatmap
png("are_heatmap.png", width = 3000, height = 2000)

#Create heatmap with custom features
heatmap_plot <- Heatmap(are_data_matrix, 
                        col = color_fun,
                        name = "logFC ",
                        show_row_names = FALSE,
                        row_names_gp = gpar(fontsize = 40, fontface = "bold"),
                        column_names_gp = gpar(fontsize = 40, fontface = "bold"),
                        column_names_rot = 0,
                        cluster_rows = FALSE,
                        cluster_columns = FALSE,  
                        height = nrow(are_data_matrix) * 34,
                        width = ncol(are_data_matrix) * 0.1,
                        top_annotation = NULL,
                        heatmap_legend_param = list(title_gp = gpar(fontsize = 46, fontface = "bold"),
                                                    labels_gp = gpar(fontsize = 46),
                                                    title = NULL,
                                                    at = c(-3, 0, 8))
)

# Print heatmap 
draw(heatmap_plot, annotation_legend_side = "bot")

dev.off()




#OE3
#Read excel file
OE3_data <- read.xlsx("OE3_CvT_DEgenes.xlsx", sheet = "OE3_Rel_Data")

#Set gene column as row name
rownames(OE3_data) <- OE3_data$gene
OE3_data_matrix <- as.matrix(OE3_data[, c("logFC_15min", "logFC_30min", "logFC_45min", "logFC_75min")])

#Rename the columns of the matrix
colnames(OE3_data_matrix) <- c("logFC 15min", "logFC 30min", "logFC 45min", "logFC 75min")

#QC data correctly loaded and transformed
print(str(OE3_data_matrix))

#Replace NA values with 0
OE3_data_matrix[is.na(OE3_data_matrix)] <- 0

#Create color mapping
breaks = c(-8, 0, 8)
colors = c("steelblue3", "white", "tomato2")

#Draw and save heatmap with adjusted width and height
png("OE3_heatmap.png", width = 3000, height = 2000)

#Create a heatmap with custom features
heatmap_plot <- Heatmap(OE3_data_matrix, 
                        col = colors,
                        name = "logFC ",
                        show_row_names = TRUE,
                        row_names_gp = gpar(fontsize = 20),
                        column_names_gp = gpar(fontsize = 20),
                        column_names_rot = 0,
                        height = nrow(OE3_data_matrix) * 15.5,
                        width = ncol(OE3_data_matrix) * 0.9  
)

#Print heatmap
draw(heatmap_plot, annotation_legend_side = "bot")

#Add title back with adjusted positioning
grid.text("Differential Gene Expression in OE3 Tomato", x = 0.5, y = 1.1, gp = gpar(fontface = "bold", fontsize = 22))

dev.off()



#VF36
#Read excel file
VF36_data <- read.xlsx("VF36_CvT_DEgenes.xlsx", sheet = "VF36_Rel_Data")

#Set gene column as row name
rownames(VF36_data) <- VF36_data$gene
VF36_data_matrix <- as.matrix(VF36_data[, c("logFC_15min", "logFC_30min", "logFC_45min", "logFC_75min")])

#Rename the columns of the matrix
colnames(VF36_data_matrix) <- c("logFC 15min", "logFC 30min", "logFC 45min", "logFC 75min")

#QC data correctly loaded and transformed
print(str(VF36_data_matrix))

#Replace NA values with 0
VF36_data_matrix[is.na(VF36_data_matrix)] <- 0

#Create color mapping
breaks = c(-8, 0, 8)
colors = c("steelblue3", "white", "tomato2")

#Draw and save heatmap with adjusted width and height
png("VF36_heatmap.png", width = 3000, height = 2000)

#Create a heatmap with custom features
heatmap_plot <- Heatmap(VF36_data_matrix, 
                        col = colors,
                        name = "logFC ",
                        show_row_names = TRUE,
                        row_names_gp = gpar(fontsize = 20),
                        column_names_gp = gpar(fontsize = 20),
                        column_names_rot = 0,
                        height = nrow(VF36_data_matrix) * 15.5,
                        width = ncol(VF36_data_matrix) * 0.9 
) 

#Print heatmap
draw(heatmap_plot, annotation_legend_side = "bot")

#Add title back with adjusted positioning
grid.text("Differential Gene Expression in VF36 Tomato", x = 0.5, y = 1.1, gp = gpar(fontface = "bold", fontsize = 22))

dev.off()

