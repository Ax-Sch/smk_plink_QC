library(data.table)

pc1_min= -0.01
pc1_max= -0.005
pc2_min= -0.017
pc2_max= -0.0128

eigenvec_file <- "/mnt/int2/lc_covid_imputed_jan_24/with_HW/pca/pca.eigenvec"
eigenvec_data <- fread(eigenvec_file, header= FALSE)
subset_data <- eigenvec_data[eigenvec_data$V3 >= pc1_min & eigenvec_data$V3 <= pc1_max &
                               eigenvec_data$V4 >= pc2_min & eigenvec_data$V4 <= pc2_max, ]

sample_ids <- subset_data$V2
all_samples <- eigenvec_data$V2 
non_eur_samples<- setdiff(all_samples, sample_ids)



write.table(sample_ids, "eur_ids.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(non_eur_samples, "non_eur_ids.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
cat("Number of eur samples in the subset: ", length(sample_ids), "\n")
cat("Number of non_eur samples in the subset: ", length(non_eur_samples), "\n")



### double check
library(tidyverse)
ggplot(eigenvec_data)+
  geom_point(aes(x=V3, y=V4))+
  geom_point(data=subset_data, aes(x=V3, y=V4), color="red")

ggplot(eigenvec_data)+
  geom_point(aes(x=V5, y=V6))+
  geom_point(data=subset_data, aes(x=V5, y=V6), color="red")

ggplot(eigenvec_data)+
  geom_point(aes(x=V7, y=V8))+
  geom_point(data=subset_data, aes(x=V7, y=V8), color="red")
