# install required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install("affy")
BiocManager::install("affyPLM")
BiocManager::install("sva")
BiocManager::install("AnnotationDbi")
BiocManager::install("hgu133plus2.db")
install.packages("ggfortify")
install.packages("tidyverse")
install.packages("plotly")
#--------------------------------------------
library("affy") ##load the affy package
library("affyPLM")
library("sva")
library("AnnotationDbi")
library("hgu133plus2.db")
library("ggfortify")
library("tidyverse")
library("plotly")

# set the working directory and read data
setwd("/home/eikthedragonslayer/Desktop/BF528/project01/CEL_files/")
Data <- ReadAffy()
eset <- rma(Data) #create eSet object

# Fit into the model
Pset <- fitPLM(Data,normalize=TRUE,background=TRUE)

# shorten sample names
nuse_plot <- boxplot(Pset)
names <- nuse_plot$names
split_names <- strsplit(names,"_")
basic_names <- list()
for (i in 1:length(split_names)) {
  temp <- substr(split_names[[i]][1],4,9)
  basic_names[[i]] <- temp
}

# compute NUSE and NLE
NUSE(Pset, names=basic_names,ylab="Intensities",main="NUSE Plot",outline=TRUE,notch=TRUE,las = 2) #TODO: CHANGE THE PLOT FORMATTING
RLE(Pset,names=basic_names,ylab="Log Expression",main="RLE Plot",outline=TRUE,notch=TRUE,las = 2)

# correct batch effects
annotations <- read.csv(file="../proj_metadata.csv")
edata <- exprs(eset) #retrieve expression from eSet object
batch <- annotations$normalizationcombatbatch
combatmod <- model.matrix(~annotations$normalizationcombatmod)
combat_edata = ComBat(dat=edata, batch=batch, mod=combatmod[,2], par.prior=TRUE, prior.plots=FALSE) #mod has 2 indicators, discard the intercept
# write to the output
write.csv(combat_edata,'normalized_expression.csv')

# perform PCA, transpose the matrix for gene scaling
transpose = t(combat_edata)
scaled_transpose = scale(transpose)
scaled_edata = t(scaled_transpose)
pca <- prcomp(scaled_edata,center=FALSE,scale=FALSE)
summary(pca) # view pca in a tabular format
names(pca)
sub_pca <- pca$x[1:134,1:2] # this frame will be ploted

# plot PC1 VS. PC2
percent_variance <- pca$sdev^2/sum(pca$sdev^2)
plot <- pca$x %>% 
  as.data.frame %>%
  ggplot(aes(x=PC1,y=PC2,color='Genes')) + geom_point(aes(text=row.names(pca$x)),size=0.5) +
  theme_bw(base_size=20) + 
  labs(x=paste0("PC1: ",round(percent_variance[1]*100,1),"%"),
       y=paste0("PC2: ",round(percent_variance[2]*100,1),"%")) +
  theme(legend.position="none")
ggplotly(plot)

# section 4 Noise filtering & dimension reduction
# filter by row (probe), for each row, at lease 20% cols should be greater than log2(15)
expres_data <- data.frame(combat_edata) # make a deep copy for future analysis meanwhile preserve the original data
test <- expres_data[which(rowSums(expres_data > log2(15)) > 0.2*134),]
# for each probe calculate the variance and get the median variance
row_vars <- rowVars(test)
mean_var <- mean(row_vars)
required_var <- qchisq(0.01,df=length(row_vars)-1,ncp=mean_var ) # find the 99th percentile of the chi-squared distribution with n degrees of freedom
# filter by row, for each row, the coefficient of variation is greater than 0.186
test <- test[which(cv(test) > 0.186),]

