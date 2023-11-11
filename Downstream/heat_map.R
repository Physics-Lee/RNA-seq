library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(org.Mm.eg.db)
library(RColorBrewer)

# import data
my_gene_expression_matrix <- read.csv("./data_concatenated/concatenated_data.csv")
my_count_data <- my_gene_expression_matrix[,-1]
rownames(my_count_data) <- my_gene_expression_matrix[,1]

# to DGE-list
group_name <- c("wt","wt","wt","an3","an3","an3","bcl7ab","bcl7ab","bcl7ab",
                "brd1213","brd1213","brd1213","brip12","brip12","brip12","brm","brm","brm",
                "minu12","minu12","minu12","pms2ab","pms2ab","pms2ab","swi3d","swi3d","swi3d",
                "swp73b","swp73b","swp73b","syd","syd","syd","sys123","sys123","sys123")
my_dge_list <- DGEList(counts=my_count_data,group=factor(group_name))

# screen genes using counts per million.
keep <- rowSums(cpm(my_dge_list)>0.5) >= 6 # 6 of 36 columns
my_dge_list <- my_dge_list[keep,]
my_dge_list$samples$lib.size <- colSums(my_dge_list$counts) # drop un-passed genes

# normalize the data (I don't know why normalizing step is done after screening step)
my_dge_list <- calcNormFactors(my_dge_list)

# Different Expressed Genes
my_dge_list_2 <- estimateCommonDisp(my_dge_list, verbose=T)
names(my_dge_list_2)
my_dge_list_2 <- estimateTagwiseDisp(my_dge_list_2)
names(my_dge_list_2)
plotBCV(my_dge_list_2)

# create a list to store tables for all comparisons
result_list <- list()
FDR_and_logFC_screened_list <- list()

# Perform exact tests for each comparison and count up- and down-regulated genes
for (i in 1:11){
  et <- exactTest(my_dge_list_2, pair=c(12, i)) # log(group i/ group 12)
  result_of_negative_binomial_test <- topTags(et, n = dim(my_dge_list$counts)[1],
               adjust.method = "BH", sort.by = "none", p.value = 1)
  result_table <- result_of_negative_binomial_test$table
  
  FDR_screened <- result_table[result_table[, 4] < 0.05, ]
  FDR_and_logFC_screened <- result_table[, 4] < 0.05 & 
    abs(result_table[, 1]) > 1
  result_list[[i]] <- result_table
  FDR_and_logFC_screened_list[[i]] <- FDR_and_logFC_screened
}

# Combine all FDR and logFC masks to create a final mask
result_mask_2 <- Reduce(`|`, FDR_and_logFC_screened_list)
sum(result_mask_2)

# Extract logFC values for all genes passing the FDR and logFC filter
all_logFC = result_list[[1]]
for (i in 2:11){
  all_logFC[,i] = result_list[[i]][,1]
}
colnames(all_logFC) <- c("an3","bcl7ab","brd1213","brip12",
                         "brm","minu12","pms2ab","swi3d",
                         "swp73b","syd","sys123")

all_logFC_screened <- all_logFC[result_mask,c(3,4,5,6,7,8,10,11)]

#######################################################
library(pheatmap)
all_logFC_screened <- all_logFC_screened[order(-all_logFC_screened[,1]),]

# Define custom color palette
my_palette <- colorRampPalette(c("#0000FF", "#FFFFFF", "#FF0000"), space = "Lab")(100)

# Set breaks and color palette for heatmap
breaks <- seq(from = -4, to = 4, length.out = 101)
col_colors <- my_palette

# Create heatmap with custom colors
pheatmap(all_logFC_screened,
         show_rownames = FALSE,
         color = col_colors,
         breaks = breaks)
#######################################################

# Create a vector specifying the desired column order
new_order <- c(3,1,2,7,6,8,4,5)
all_logFC_screened <- all_logFC_screened[, new_order]
cor_matrix <- cor(all_logFC_screened)
library(gplots)
heatmap.2(cor_matrix, trace = "none", dendrogram = "none",
          col = colorRampPalette(c("blue", "white", "red"))(100),
          Rowv = FALSE, Colv = FALSE)
######################################################
gene_name <- c("AT5G67060","AT2G45210","AT3G45230","AT5G64170","AT2G28550","AT3G15354",
               "AT1G10870","AT3G13960","AT1G01470","AT3G15510","AT3G06490","AT3G57880",
               "AT2G28610","AT1G68480","AT5G55820","AT3G07780","AT3G06430",
               "AT5G65080","AT5G10140","AT3G22886","AT3G54340","AT2G46830",
               "AT1G65480","AT1G69490","AT4G15248","AT4G24540")
gene_name_popular <- c("HEC1","SAUR36","APAP1","LNK1","RAP2.7","SPA3",
                       "VAL3","GRF5","LEA14","NAC056","MYB108","MCTP3",
                       "WOX3","JAG","WYRD","OBE1","PPR2",
                       "MAF5","FLC","MIR167A","AP3","CCA1",
                       "FT","NAC029","MIP1A","AGL24")
old_to_new_names <- setNames(gene_name_popular, gene_name)
matrix_of_fig_5 <- all_logFC_screened[rownames(all_logFC_screened) %in% gene_name,]
matrix_of_fig_5 <- matrix_of_fig_5[order(match(rownames(matrix_of_fig_5), gene_name)),]
rownames(matrix_of_fig_5) <- old_to_new_names[rownames(matrix_of_fig_5)]
matrix_of_fig_5_transpose <- t(matrix_of_fig_5)

# Define the color palette
my_palette <- colorRampPalette(c("#0000FF", "#FFFFFF", "#FF0000"), space = "Lab")(100)

# Create the heatmap with row and column clustering disabled
pheatmap(matrix_of_fig_5_transpose,
         color = my_palette,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_colnames = TRUE,
         show_rownames = TRUE,
         )

pheatmap(matrix_of_fig_5_transpose)
heatmap.2(matrix_of_fig_5_transpose)