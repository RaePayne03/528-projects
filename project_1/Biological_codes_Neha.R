#installing required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("hgu133plus2.db")
install.packages('tidyverse')
BiocManager::install("affy")
BiocManager::install("GSEABase")

#calling the required packages
library(BiocManager)
library(hgu133plus2.db)
library(tidyverse)
library(GSEABase)


#setting the directory
setwd("C:/Users/Neha/Desktop/BU/bf528/project1")

#reading csv file and assign appropriate column names
diff_exp_results <- read.csv('welch_5_6.csv', col.names = c('num','PROBEID', 't', 'p', 'p_adj'))
#remove column 1 as it had row numbers
diff_exp_results <- diff_exp_results[, -1]
head(diff_exp_results)

## 1.Using the select() function of the bioconductor package hgu133plus2.db, map the probeset IDs to gene symbols by specifying the appropriate key and column arguments. 
gene_symbols <- AnnotationDbi::select(hgu133plus2.db, diff_exp_results$PROBEID, c('SYMBOL'))
View(gene_symbols)

#function to merge 2 gene symbols with same probe id
merger <- function(x) {x %>% unique %>% sort %>% paste(collapse = '|') }

#running merger on the genesymbols of our data
gene_symbols <- gene_symbols %>% group_by(PROBEID) %>% summarise_each(funs(merger)) %>% ungroup
View(gene_symbols)


#joining diff_exp_results with this gene symbols
merged_results <- merge(diff_exp_results, gene_symbols, on ='PROBEID')
View(merged_results)

#dropping the null rows
merged_results <- merged_results[!(is.na(merged_results$SYMBOL) | merged_results$SYMBOL == ""), ]

#finding and counting most significant probe id for each gene symbol that has multiple probe ids
significant <- merged_results %>% group_by(SYMBOL) %>% count() %>% filter(n>1)
significant[order(significant$n), ]
View(significant)

#for each repaeated symbol, taking min p_adjusted value to consider for most significant among the repeated ones
diff <- data.frame(PROBEID = character(), t = numeric(), p = numeric(), p_adj = numeric(), SYMBOL = character() )
for (i in significant$SYMBOL) 
  {
  x <- merged_results[merged_results$SYMBOL == i, ]
  x <- x[x$p_adj == min(x$p_adj),  ]
  merged_results <- merged_results[!merged_results$SYMBOL == i, ]
  diff <- rbind(diff, x)
}

merged_results <- rbind(merged_results, diff)


#importing gene data sets
go <- getGmt('c5.go.v7.2.symbols.gmt')
kegg <- getGmt('c2.cp.kegg.v7.2.symbols.gmt')
hallmark <- getGmt('h.all.v7.2.symbols.gmt')

#number of genes in each gene sets
length(go) #10271
length(kegg) #186
length(hallmark) #50


#for up and down regulated, arrangge all the t-stats values in decreaing order
merged_results <- merged_results[order(merged_results$t, decreasing = T), ]
head(merged_results)
tail(merged_results)


#extract top and bottom 1000 for up and down regualted genes
top_1000_up <- head(merged_results, n = 1000)
top_1000_down <- tail(merged_results, n=1000)

#for report, selecting top 10 out of those 1000 for each up and down
top_10_up <- head(top_1000_up, n=10)
top_10_down <- tail(top_1000_down, n =10)

#writing the results for top 10 up and down in cv files
write.csv(top_10_up, 'Top_10_upregulated.csv')
write.csv(top_10_down, 'Top_10_downregulated.csv')


#genes not expressed up and down
up_not_diff_exp <- subset(merged_results, !merged_results$SYMBOL %in% top_1000_up$SYMBOL)
down_not_diff_exp <- subset(merged_results, !merged_results$SYMBOL %in% top_1000_down$SYMBOL)

#constrcuting function for fisher test and contigency table
fisher_table <- function(gene_list, gene_set, not_diff_exp)
{
  diff_exp_in_gene_set <- length(intersect(gene_list, gene_set))
  diff_exp_notin_gene_set <- length(gene_list) - diff_exp_in_gene_set
  not_diff_exp_ingene_set <-  length(intersect(not_diff_exp, gene_set))
  not_diff_exp_not_ingene_set <- length(not_diff_exp) - not_diff_exp_ingene_set
  
  return(c(diff_exp_in_gene_set, diff_exp_notin_gene_set, not_diff_exp_ingene_set, not_diff_exp_not_ingene_set))
}

#creating dataframe to store table results
go_result <- data.frame(name = character(), pvalue = numeric(), estimate = numeric(), expression = character(), stringsAsFactors = F)
kegg_result <- data.frame(name = character(), pvalue = numeric(), estimate = numeric(), expression = character(), stringsAsFactors = F)
hallmark_result <- data.frame(name = character(), pvalue = numeric(), estimate = numeric(), expression = character(), stringsAsFactors = F)

#computing fisher test and contigency table for each gene set
#GO
for (i in 1:length(go))
{
  gene_id <- geneIds(go[i])
  fisher_up <- fisher_table(top_1000_up$SYMBOL, gene_id[[names(gene_id)]], up_not_diff_exp$SYMBOL)
  fisher_down <- fisher_table(top_1000_down$SYMBOL, gene_id[[names(gene_id)]], down_not_diff_exp$SYMBOL)
  up <- fisher.test(matrix(fisher_up,nrow=2))
  down <- fisher.test(matrix(fisher_down, nrow=2))
  go_result[nrow(go_result) +1, ] <- c(names(gene_id), up$p.value, up$estimate, 'UP')
  go_result[nrow(go_result) +1, ] <- c(names(gene_id), down$p.value, down$estimate, 'Down')
}

go_result <- go_result %>% mutate(pvalue = as.numeric(pvalue), estimate = as.numeric(estimate))
View(go_result)


#KEGG
for (i in 1:length(kegg))
{
  gene_id <- geneIds(kegg[i])
  fisher_up <- fisher_table(top_1000_up$SYMBOL, gene_id[[names(gene_id)]], up_not_diff_exp$SYMBOL)
  fisher_down <- fisher_table(top_1000_down$SYMBOL, gene_id[[names(gene_id)]], down_not_diff_exp$SYMBOL)
  up <- fisher.test(matrix(fisher_up,nrow=2))
  down <- fisher.test(matrix(fisher_down, nrow=2))
  kegg_result[nrow(kegg_result) +1, ] <- c(names(gene_id), up$p.value, up$estimate, 'UP')
  kegg_result[nrow(kegg_result) +1, ] <- c(names(gene_id), down$p.value, down$estimate, 'Down')
}

kegg_result <- kegg_result %>% mutate(pvalue = as.numeric(pvalue), estimate = as.numeric(estimate))
View(kegg_result)

##Hallmarks
for (i in 1:length(hallmark))
{
  gene_id <- geneIds(hallmark[i])
  fisher_up <- fisher_table(top_1000_up$SYMBOL, gene_id[[names(gene_id)]], up_not_diff_exp$SYMBOL)
  fisher_down <- fisher_table(top_1000_down$SYMBOL, gene_id[[names(gene_id)]], down_not_diff_exp$SYMBOL)
  up <- fisher.test(matrix(fisher_up,nrow=2))
  down <- fisher.test(matrix(fisher_down, nrow=2))
  hallmark_result[nrow(hallmark_result) +1, ] <- c(names(gene_id), up$p.value, up$estimate, 'UP')
  hallmark_result[nrow(hallmark_result) +1, ] <- c(names(gene_id), down$p.value, down$estimate, 'Down')
}

hallmark_result <- hallmark_result %>% mutate(pvalue = as.numeric(pvalue), estimate = as.numeric(estimate))
View(hallmark_result)



#performing Benjamini_Hochberg method
go_result$BH <- p.adjust(go_result$pvalue, method = "BH", n = length(go_result$pvalue))
kegg_result$BH <- p.adjust(kegg_result$pvalue, method = "BH", n = length(kegg_result$pvalue))
hallmark_result$BH <- p.adjust(hallmark_result$pvalue, method = "BH", n = length(hallmark_result$pvalue))


#wriitng final files
write.csv(go_result, 'Final_go_results.csv')
write.csv(kegg_result, 'Final_kegg_results.csv')
write.csv(hallmark_result, 'Final_hallmark_results.csv')

