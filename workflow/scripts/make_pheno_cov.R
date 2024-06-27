library(tidyverse)
library(optparse)

option_list = list(
  make_option(c("-b", "--pheno_buffy"), type="character", default="resources/pheno_age_buffy_coat.tsv"),
  make_option(c("-l", "--pheno_longCovid"), type="character", default="resources/pheno_age_longCovids.tsv"),
  make_option(c("-p", "--pcs"), type="character", default="results/covar_PCA/pcs.txt"),
  make_option(c("-c", "--cov"), type="character", default="results/pheno/cov.cov"),
  make_option(c("-t", "--pheno"), type="character", default="results/pheno/pheno.pheno")
)

opt = parse_args(OptionParser(option_list=option_list))

pheno_longcovid<-read_tsv(opt$pheno_longCovid)
pheno_buffy<-read_tsv(opt$pheno_buffy)
pcs<-read_tsv(opt$pcs)

# model:
# Phenotype ~ variant + age + age2 + sex + age*sex + PCs + study_specific_covariates

# create phenotypes
pheno_buffy_mod<-pheno_buffy %>%
  mutate(NQ1_3=NA,
         NQ2_3=0)

pheno_longcovid_mod<-pheno_longcovid%>%
  mutate(NQ1_3=ifelse(long_covid=="ja",1,
                     ifelse(long_covid=="nein", 0, NA) )
  )%>%
  mutate(NQ2_3=NQ1_3) %>% 
  select(-long_covid)

pheno_all<-rbind(pheno_buffy_mod,
                 pheno_longcovid_mod )

# create covariates
pheno_all_cov<-pheno_all %>%
  mutate(age2=age*age,
         age_sex=age*sex) %>%
  distinct(USI, .keep_all=T)

# add to pcs
pcs_w_pheno_cov<-pcs %>%
  left_join(pheno_all_cov, by=c("IID"="USI"))


# splice out pheno and cov files:
final_pheno_file<-pcs_w_pheno_cov %>% 
  select(FID, IID, NQ1_3, NQ2_3)

write_tsv(x=final_pheno_file,
          file=opt$pheno)

final_cov_file<-pcs_w_pheno_cov %>% 
  select(-NQ1_3, -NQ2_3)

write_tsv(x=final_cov_file,
          file=opt$cov)
