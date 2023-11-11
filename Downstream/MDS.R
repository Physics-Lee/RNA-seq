rm(list = ls())

library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(org.Mm.eg.db)
library(RColorBrewer)

# just for fun
print("haha")

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
my_dge_list
dim(my_dge_list)

# screen genes using counts per million.
keep <- rowSums(cpm(my_dge_list)>0.5) >= 6 # 6 of 36 columns
my_dge_list <- my_dge_list[keep,]
dim(my_dge_list)
my_dge_list$samples$lib.size <- colSums(my_dge_list$counts)

barplot(my_dge_list$samples$lib.size, las=2)
title("Barplot of library sizes")
par(las=2)  # rotate the x-axis labels
text(x = barplot(dgeObj$samples$lib.size, add=TRUE), y = -5, 
     labels = colnames(dgeObj), xpd=TRUE, srt=90, adj=1, cex=0.6)

# normalize the data (I don't know why normalizing step is done after screening step)
my_dge_list <- calcNormFactors(my_dge_list)
my_dge_list

# MDS plot
color_vector <- c("red", "green", "blue", "orange", "purple", "gray",
                  "pink", "brown", "black", "cyan", "magenta", "yellow")
cell_types <- factor(rep(c("type1", "type2", "type3", "type4", "type5", "type6",
                           "type7", "type8", "type9", "type10", "type11", "type12"
                           ), each = 3))
col_vector <- color_vector[as.numeric(cell_types)]
plotMDS(my_dge_list, col=col_vector)
title("MDS plot")