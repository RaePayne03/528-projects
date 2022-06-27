#7.3 Heatmap


#load differential expressed genes 
de_genes <- read.delim("/projectnb/bf528/users/swiss_cheese2/project_2/cuffdiff_out/gene_exp.diff")

#find significant genes, create new top 100 deferentially expressed genes subset
de_genes  <- de_genes %>% arrange(q_value)  %>% slice_head(n=1000)
de_genes  <- de_genes[de_genes$significant =='yes',]
top  <- de_genes %>% arrange(q_value)  %>% slice_head(n=50)
deg_top <- top$gene

#subset to top 50 genes 
fpkm_combined_top <- fpkm_combined[fpkm_combined$gene_short_name %in% deg_top, ]

#create final_fpkm_combined matrix 
final_fpkm_combined <- fpkm_combined_top[-1:-2]
#make genes row names
rownames(final_fpkm_combined) <- make.names(fpkm_combined[,2], unique = TRUE)
rownames(final_fpkm_combined) <- make.name(fpkm_combined$gene_short_name)
final_fpkm_combined <- final_fpkm_combined[apply(final_fpkm_combined[,-1], 1, function(x) !all(x==0)),]


#convert to matrix for heatmap
data <- as.matrix(final_fpkm_combined[,])

#create heatmap
colors = brewer.pal(n = 11, name = "GnBu")
colors = colorRampPalette(colors)(50)
colors = rev(colors)

heatmap <- pheatmap(data, scale = "row", color = colors,fontsize_row = 4,border_color = NA, clustering_distance_rows="euclidean",
                    clustering_distance_cols="euclidean", main = "Top 1000 DE Genes ")
