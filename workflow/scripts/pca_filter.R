library(optparse)
library(data.table)

# Define command line arguments
option_list <- list(
  make_option(c("--pc1_min"), type = "numeric", default = -0.01,
              help = "Minimum value for PC1", metavar = "numeric"),
  make_option(c("--pc1_max"), type = "numeric", default = -0.005,
              help = "Maximum value for PC1", metavar = "numeric"),
  make_option(c("--pc2_min"), type = "numeric", default = -0.017,
              help = "Minimum value for PC2", metavar = "numeric"),
  make_option(c("--pc2_max"), type = "numeric", default = -0.0128,
              help = "Maximum value for PC2", metavar = "numeric"),
  make_option(c("--eigenvec_file"), type = "character", default = "results/PCA/pca.eigenvec",
              help = "Path to the eigenvec file", metavar = "character"),
  make_option(c("--eur_outfile"), type = "character", default = "eur_ids.txt",
              help = "Output file for EUR sample IDs", metavar = "character"),
  make_option(c("--non_eur_outfile"), type = "character", default = "non_eur_ids.txt",
              help = "Output file for non-EUR sample IDs", metavar = "character")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Assign arguments to variables
pc1_min <- opt$pc1_min
pc1_max <- opt$pc1_max
pc2_min <- opt$pc2_min
pc2_max <- opt$pc2_max
eigenvec_file <- opt$eigenvec_file
eur_outfile <- opt$eur_outfile
non_eur_outfile <- opt$non_eur_outfile

# Read in the eigenvec data
eigenvec_data <- fread(eigenvec_file, header= FALSE)

# Subset the data based on the provided PC1 and PC2 ranges
subset_data <- eigenvec_data[eigenvec_data$V3 >= pc1_min & eigenvec_data$V3 <= pc1_max &
                               eigenvec_data$V4 >= pc2_min & eigenvec_data$V4 <= pc2_max, ]

# Extract sample IDs
sample_ids <- subset_data$V2
all_samples <- eigenvec_data$V2 
non_eur_samples <- setdiff(all_samples, sample_ids)

# Write the output files
write.table(sample_ids, eur_outfile, row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(non_eur_samples, non_eur_outfile, row.names = FALSE, col.names = FALSE, quote = FALSE)

# Print the number of samples
cat("Number of EUR samples in the subset: ", length(sample_ids), "\n")
cat("Number of non-EUR samples in the subset: ", length(non_eur_samples), "\n")

