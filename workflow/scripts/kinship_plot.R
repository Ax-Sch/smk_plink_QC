# Load necessary libraries
library(ggplot2)
library(plotly)
library(dplyr)
library(rmarkdown)

#kinship_data_1_path<- "results/kinship/Geno05_CR_sex_snp_qc_snpqc2_EUR.kin0"
#kinship_data_2_path<- "results/kinship/last_check_kinship.kin0"
# Load the kinship files (paths provided by Snakemake)
kinship_data_1_path= "VCR_CR_sex_snp_qc_snpqc2_EUR.kin0"
kinship_data_2_path= "last_check_kinship.kin0"
output_graphs_path= "kinship_scatter_plots.html"
kinship_data_1 <- read.table(file=kinship_data_1_path, header = TRUE)
kinship_data_2 <- read.table(file=kinship_data_2_path, header = TRUE)

# Plot IBS0 vs Kinship for the first kinship file (interactive with plotly)
p1 <- ggplot(kinship_data_1, aes(x = IBS0, y = Kinship, text = paste("FID1:", FID1, "ID1:", ID1, "FID2:", FID2, "ID2:", ID2))) +
    geom_point(alpha = 0.7, color = "blue") +
    labs(title = "Kinship File 1: IBS0 vs Kinship", x = "IBS0", y = "Kinship") +
    theme_minimal()

# Convert the ggplot object to a plotly object for interactivity
p1_interactive <- ggplotly(p1, tooltip = "text")

# Plot IBS0 vs Kinship for the second kinship file (interactive with plotly)
p2 <- ggplot(kinship_data_2, aes(x = IBS0, y = Kinship, text = paste("FID1:", FID1, "ID1:", ID1, "FID2:", FID2, "ID2:", ID2))) +
    geom_point(alpha = 0.7, color = "green") +
    labs(title = "Kinship File 2: IBS0 vs Kinship", x = "IBS0", y = "Kinship") +
    theme_minimal()

# Convert the ggplot object to a plotly object for interactivity
p2_interactive <- ggplotly(p2, tooltip = "text")

# Display both interactive plots
subplot(p1_interactive, p2_interactive, nrows = 1, titleX = TRUE, titleY = TRUE)

ggsave("kinship_plot1.pdf", plot = p1, width = 6, height = 5)
ggsave("kinship_plot2.pdf", plot = p2, width = 6, height = 5)
