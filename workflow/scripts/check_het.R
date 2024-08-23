library(tidyverse)
library(optparse)


option_list = list(
  make_option(c("-a", "--cases"), type="character"),
  make_option(c("-b", "--conts"), type="character"),
  make_option(c("-o", "--outfile"), type="character"),
  make_option(c("-l", "--lower_cutoff"), type="numeric"),
  make_option(c("-u", "--upper_cutoff"), type="numeric")
  
)
opt = parse_args(OptionParser(option_list=option_list))


het_file_cases=read.table(file = opt$cases, sep="", header=T)
het_file_controls=read.table(file = opt$conts, sep="", header=T)

F_cutoff=as.numeric(opt$cutoff)
####!!!!!!





cases_not_passing<-het_file_cases %>%
  filter(F > 0.2 | F < -0.2)

controls_not_passing<-het_file_controls %>%
  filter(F > 0.2 | F < -0.2)

all_not_passing<-rbind(cases_not_passing,
                       controls_not_passing)

print("not passing samples:")
print(all_not_passing)

write.table(x=all_not_passing,
            file=opt$outfile)
