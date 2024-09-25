# Load necessary libraries
library(randomForest)
library(tidyverse)
library(plotly)
library(rmarkdown)

# Get command line arguments
args=c('pca.eigenvec', '../../../resources/H_1000G/H_ped/20130606_g1k.ped', 'population_PCA.html')

# Define input and output files
eigenvec_file <- args[1]
ped_file <- args[2]
output_file <- args[3]

# Read in the eigenvectors, produced in PLINK
eigenvec <- read.table(eigenvec_file, header = FALSE, skip=0, sep = ' ')
eigenvec <- eigenvec[,2:ncol(eigenvec)]
colnames(eigenvec) <- c("Individual.ID",paste('PC', c(1:20), sep = ''))

# Read in the PED data
PED <- read.table(ped_file, header = TRUE, skip = 0, sep = '\t')

# Build data frame for random forest classifier
dataRF <- merge(eigenvec, PED[, c("Individual.ID", "Population")], all.x=TRUE)

# Build plot
dataRF$Population <- factor(dataRF$Population, levels=c(
  "ACB","ASW","ESN","GWD","LWK","MSL","YRI",
  "CLM","MXL","PEL","PUR",
  "CDX","CHB","CHS","JPT","KHV",
  "CEU","FIN","GBR","IBS","TSI",
  "BEB","GIH","ITU","PJL","STU"))

dataRF$Continental <- rep(NA_character_, nrow(dataRF))
dataRF$Continental[which(dataRF$Population %in% c("ACB","ASW","ESN","GWD","LWK","MSL","YRI"))]<-"AFR"
dataRF$Continental[which(dataRF$Population %in% c("CLM","MXL","PEL","PUR"))]<-"AMR"
dataRF$Continental[which(dataRF$Population %in% c("CDX","CHB","CHS","JPT","KHV"))]<-"EAS"
dataRF$Continental[which(dataRF$Population %in% c("CEU","FIN","GBR","IBS","TSI"))]<-"EUR"
dataRF$Continental[which(dataRF$Population %in% c("BEB","GIH","ITU","PJL","STU"))]<-"SAS"
dataRF$Continental<-as.factor(dataRF$Continental)

rf_classifier = randomForest(Continental ~ ., 
                             data=dataRF[which(!is.na(dataRF$Continental)), c("PC1","PC2","PC3","PC4","PC5","PC6", "Continental")], 
                             ntree=3000, importance=TRUE)

# Predict population in your cohort
dataPred<-dataRF[which(is.na(dataRF$Continental)),]
dataRF<-dataRF[which(!is.na(dataRF$Continental)),]
dataPred$Prediction<-rep(NA, nrow(dataPred))
dataPred$Prediction<-predict(rf_classifier,dataPred[,c("PC1","PC2","PC3","PC4","PC5","PC6")])

PC1_2_plot <- ggplot() +
  geom_point(data=dataRF, aes(x=PC1, y=PC2, shape=Continental), alpha=0.5, size=3) + 
  geom_point(data=dataPred, aes(x=PC1, y=PC2, color=Prediction, text=Individual.ID), shape=1, alpha=0.8, size=3) +
  theme_bw() +
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10))

ggsave(filename="PC1_2_plot.pdf", plot=PC1_2_plot, width=4, height=3.5)
ggplotly(PC1_2_plot)

ggplotly(ggplot()+
	            geom_point(data=dataRF, aes(x=PC3,y=PC4,shape=Continental), alpha=0.5,size=3)+ 
		               geom_point(data=dataPred, aes(x=PC3,y=PC4,color=Prediction, text=Individual.ID),shape=1, alpha=0.8,size=2 ) )

ggplotly(ggplot()+
	            geom_point(data=dataRF, aes(x=PC5,y=PC6,shape=Continental), alpha=0.5,size=3)+ 
		               geom_point(data=dataPred, aes(x=PC5,y=PC6,color=Prediction, text=Individual.ID),shape=1, alpha=0.8,size=2 ) )

ggplotly(ggplot()+
	            geom_point(data=dataRF, aes(x=PC7,y=PC8,shape=Continental), alpha=0.5,size=3)+ 
		               geom_point(data=dataPred, aes(x=PC7,y=PC8,color=Prediction, text=Individual.ID),shape=1, alpha=0.8,size=2 ) )

# Save the results
write_tsv(dataPred %>% select(Individual.ID, Prediction), file="populations.txt")

# Render the R Markdown file to HTML
