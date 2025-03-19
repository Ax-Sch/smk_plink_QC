library(optparse)
library(ggplot2)
library(plotly)
library(dplyr)
library(rmarkdown)

#option_list = list(
#    make_option(c("-a", "--sexcheckfile"), type="character"),
#    make_option(c("-b", "--problem_sex_info"), type="character"),
#    make_option(c("-c", "--problem_sex_ids"), type="character"),
#    make_option(c("-g", "--f_coeff_graph"), type="character")
#)


# Parse command-line arguments
#pt_parser <- OptionParser(option_list = option_list)
#opt <- parse_args(pt_parser)

# Read the sex check file
sexcheck_data <- read.table("VCR_SCR.sexcheck",, header = TRUE)

# Extract rows with "PROBLEM" in the STATUS column
problem_sex_info <- subset(sexcheck_data, STATUS == "PROBLEM")

# Write the problem sex info to a file
write.table(problem_sex_info, file = "sex_check_problem_sex.txt", row.names = FALSE, quote = FALSE, sep = "\t")

# Extract FID and IID for problem samples and write to a file
problem_sex_ids <- problem_sex_info[, c("FID", "IID")]
write.table(problem_sex_ids, file = "sex_check_problem_ids.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

########graph
# Define custom colors
sexcheck_data <- sexcheck_data %>%
    mutate(color = case_when(
        STATUS == "PROBLEM" ~ "violet",
        PEDSEX == 1 ~ "lightblue",  # Male
        PEDSEX == 2 ~ "lightpink"   # Female
   ))
# Create a ggplot object
p <- ggplot(sexcheck_data, aes(x = FID, y = F, color = color, text = paste("ID:", IID))) +
	geom_point(size = 3, alpha = 0.8) +
	scale_color_identity() +
	labs(title = "F-Coefficient by FID",
	    x = "FID",
	    y = "F Coefficient") +
	theme_minimal()

# Convert the ggplot to a plotly object for interactivity
p_interactive <- ggplotly(p, tooltip = "text")

#Display the interactive plot
p_interactive
#htmlwidgets::saveWidget(p_interactive, opt$f_coeff_graph)


