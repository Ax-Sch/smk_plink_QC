library(tidyverse)
library(optparse)
library(ggplot2)

option_list = list(
  make_option(c("-a", "--cases"), type="character"),
  make_option(c("-b", "--conts"), type="character"),
  make_option(c("-l", "--lower_cutoff"), type="numeric", default= -0.2),
  make_option(c("-u", "--upper_cutoff"), type="numeric", default= 0.2),
  make_option(c("-o", "--outfile"), type="character"),
  make_option(c("-c", "--outgraph_cases"), type="character"),
  make_option(c("-t", "--outgraph_controls"), type="character") 
)

pt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(pt_parser)

het_file_cases=read.table(file = opt$cases, sep="", header=T)
het_file_controls=read.table(file = opt$conts, sep="", header=T)

F_lco<-opt$lower_cutoff #lower cutoff
F_uco<-opt$upper_cutoff #upper cutoff

cases_not_passing<-het_file_cases %>%
  filter(F > F_uco | F < F_lco)

controls_not_passing<-het_file_controls %>%
  filter(F > F_uco | F < F_lco)

all_not_passing<-rbind(cases_not_passing,
                       controls_not_passing)

print("not passing samples:")
print(all_not_passing)

write.table(x=all_not_passing,
            file=opt$outfile)

##############graphs

#extract relevant columns (IID and F coefficient)
cases_data <- data.frame(IID = het_file_cases$IID, F = het_file_cases$F)
controls_data <- data.frame(IID = het_file_controls$IID, F = het_file_controls$F)

# Function to create and save a histogram
plot_histogram <- function(data, title, lower_cutoff, upper_cutoff, output_file) {
    p <- ggplot(data, aes(x = F)) +
    geom_histogram(binwidth = 0.02, fill = "lightblue", color = "black") +
    geom_vline(xintercept = lower_cutoff, color = "red", linetype = "dashed", size = 1) +
    geom_vline(xintercept = upper_cutoff, color = "red", linetype = "dashed", size = 1) +
    labs(title = title, x = "F coefficient", y = "Count") +
    theme_minimal() +
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14))
				    
    # Save the plot to a PDF file
    ggsave(output_file, plot = p, width = 8, height = 6)
}

# Create and save the histogram for cases
plot_histogram(cases_data, "Cases - F Coefficient Histogram", F_lco, F_uco, opt$outgraph_cases)

# Create and save the histogram for controls
plot_histogram(controls_data, "Controls - F Coefficient Histogram", F_lco, F_uco, opt$outgraph_controls)

