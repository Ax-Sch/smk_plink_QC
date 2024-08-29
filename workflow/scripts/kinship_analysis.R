library(tidyverse)
library(optparse)

option_list = list(
  make_option(c("-f", "--fam_file"), type="character", default="results/PCA/cohort_eur.fam"),
  make_option(c("-k", "--kinship_file1"), type="character", default="results/kinship/Geno05_CR_sex_snp_qc_snpqc2_EUR.kin0"),
  make_option(c("-i", "--kinship_file2"), type="character", default="results/kinship/Geno05_CR_sex_snp_qc_snpqc2_EUR.kin0"),
  make_option(c("-c", "--kinship_cutoff"), type="numeric", default= 0.04),
  make_option(c("-o", "--outfile"), type="character", default="results/kinship/exclude_ids.txt")
)

opt = parse_args(OptionParser(option_list=option_list))

#parent-offspring (PO) pairs, full sibling (FS) 

fam_file<-read_tsv(file=opt$fam_file, col_names = FALSE)

cases<-fam_file[fam_file$X6==2,]$X2

kinship_matrix1<-read_tsv(opt$kinship_file1)
kinship_matrix2<-read_tsv(opt$kinship_file2)
kinship_cutoff<- opt$kinship_cutoff

common_colnames_ind <-colnames(kinship_matrix1) %in% colnames(kinship_matrix2)
common_colnames<-colnames(kinship_matrix1)[common_colnames_ind]

kinship_matrix<-rbind(kinship_matrix1 %>% select(any_of(common_colnames)),
                      kinship_matrix2 %>% select(any_of(common_colnames)))

relateds<-kinship_matrix %>% 
  filter(Kinship>kinship_cutoff) %>%
  distinct()

# make sorted related list
ids_relateds<-tibble(
  id=c(relateds$ID1, relateds$ID2)   )

ids_relateds<-ids_relateds %>% 
  mutate(case=id %in% cases) %>% 
  group_by(id, case) %>%
  summarise(n_relateds=n()) %>%
  arrange(case, -n_relateds)

#remove the controls with many relateds first
exclude=ids_relateds[0,]
related_list_track=relateds
for (exclude_candidate in ids_relateds$id){
  if (exclude_candidate %in% c(related_list_track$ID1,related_list_track$ID2) ){
    related_list_track<-related_list_track %>%
      filter(ID1 != exclude_candidate) %>%
      filter(ID2 != exclude_candidate)
    
    exclude=rbind(exclude, 
                  ids_relateds[ids_relateds$id==exclude_candidate, ])
  }
}

print("excluding:")
print(exclude)

# write to file
fam_for_exclude<-fam_file %>% filter(!X2 %in% exclude$id)
write_delim(file=opt$outfile, fam_for_exclude, delim=" ", col_names = FALSE)

#### run Plink
#### run KING
#### check if there are still relateds:

#nrow(read_tsv("kinship/Geno05_CR_sex_snp_qc_snpqc2_EUR_remaining.kin0") %>% filter(Kinship>0.04))
#nrow(read_tsv("kinship/Geno05_CR_sex_snp_qc_snpqc2_EUR_remaining.kin") %>% filter(Kinship>0.04))


