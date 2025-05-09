# load libraries
library(devtools)
library(ggplot2)
library(ggbiplot)

# Change workspace to directory containing files
setwd("C:/Users/Darcy/PCA Analysis/PCA")

# load log FC data set containing standardized array values (i.e. normalized array value - average normalized array value per gene across all sample groups)
PCAdata = read.csv(file="PCA_IPF_COPD_CTRL_norm.csv", header = TRUE, sep=",")

# display first three rows of the data set
head(PCAdata, 3)

# cut the disease.state and label column, leaving the genes only (PCAdata[, 1:n]) where n is the number of genes
PCAdata.genes <- PCAdata[, 1:15]

# create a vector with the disease.state column in PCAdata[, n+1] where n is the number of genes
disease.state <- PCAdata[, 16]

# apply PCA on log transformed data (PCAdata) 
# scale. = TRUE is highly advisable, but default is FALSE. 
matrix.pca <- prcomp(PCAdata.genes, center = TRUE, scale. = TRUE)

# print method returns the standard deviation of each of the four PCs
# and their rotation (or loadings), which are the coefficients of the 
# linear combinations of the continuous variables
print(matrix.pca)

# plot method returns a plot of the variances (y-axis) associated with the PCs (x-axis)
plot(matrix.pca, type = "l")

# summary method describes the importance of the PCs
# The first row describe again the standard deviation associated with each PC
# The second row shows the proportion of the variance in the data explained by each component 
# the third row describe the cumulative proportion of explained variance
summary(matrix.pca) 

# create a biplot
# biplot generated by the function ggbiplot of the ggbiplot package
PCAplot <- ggbiplot::ggbiplot(matrix.pca, obs.scale = 1, var.scale = 2, groups = PCAdata$Disease.state, labels = NULL, point.size = 1.25, ellipse = TRUE, ellipse.fill = FALSE, ellipse.linewidth = 0.5,
                              var.axes=FALSE, loadings= FALSE)
print(PCAplot)
final_plot<- PCAplot+ theme_bw()
print(final_plot)
theme <-
  theme(
    panel.grid = element_line(),
    axis.line.x = element_line(),
    axis.line.y = element_line())
PCA_white <- PCAplot + theme_bw() + theme(legend.direction = 'horizontal', legend.position = "top", legend.text = element_text(size=12), legend.key.size = unit(1, 'cm'), legend.key.height = unit(1, 'cm'), legend.key.width = unit(0.5, 'cm')) + theme +  geom_hline(yintercept=0, alpha=0.4) + geom_vline(xintercept=0, alpha=0.4) 
PCA_axis <- PCA_white + theme(axis.title = element_text(size = 14, color = "black", face = "bold"), 
                              axis.text = element_text(color = "black", size=12, face = NULL))
# display PCA plot
print(PCA_axis)

# use a palette modified for persons with colorblindness (https://davidmathlogic.com/colorblind/#%23648FFF-%23785EF0-%23DC267F-%23FE6100-%23FFB000) 
PCA_colors <- PCA_axis + scale_color_manual(values = c("COPD" = "#FFB000","CTRL" = "#648FFF", "IPF" = "#DC267F")) 

# display updated PCA plot
print(PCA_colors)

# save file (PCA_colors)
ggsave("1D.tiff", plot = last_plot(), width = 5, height = 4, dpi = 400)
