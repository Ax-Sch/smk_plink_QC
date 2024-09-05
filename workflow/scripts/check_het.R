library(tidyverse)
library(optparse)


option_list = list(
  make_option(c("-a", "--cases"), type="character"),
  make_option(c("-b", "--conts"), type="character"),
  make_option(c("-l", "--lower_cutoff"), type="numeric", default= -0.2),
  make_option(c("-u", "--upper_cutoff"), type="numeric", default= 0.2),
  make_option(c("-o", "--outfile"), type="character")
  
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
