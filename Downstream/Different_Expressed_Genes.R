rm(list = ls())

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

# Create empty vectors to store the results
num_up_genes <- numeric(11)
num_down_genes <- numeric(11)

# Perform exact tests for each comparison and count up- and down-regulated genes
for (i in 1:11){
  et <- exactTest(my_dge_list_2, pair=c(12, i)) # log(group i/ group 12)
  x <- topTags(et, n = 10000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
  x_2 <- x$table
  x_3 <- x_2[x_2[, 4] < 0.05, ]
  num_up_genes[i] <- sum(x_3[, 1] > 1)
  num_down_genes[i] <- sum(x_3[, 1] < -1)
}

# Print the results
print(num_up_genes)
print(num_down_genes)

# Create a matrix with the number of up- and down-regulated genes for each comparison
gene_counts <- cbind(num_up_genes, num_down_genes)

x_ticks <- c("an3","bcl7ab","brd1213","brip12",
             "brm","minu12","pms2ab","swi3d",
             "swp73b","syd","sys123")

# Create a color matrix with two colors for each comparison
color_matrix <- matrix(
  c(rep(c("red", "blue"), each = 1)),
  ncol = 2,
  byrow = TRUE
)

# Create the bar plot
bp <- barplot(
  t(gene_counts),
  col = color_matrix,
  xlab = "Group comparison",
  ylab = "Number of genes",
  names.arg = x_ticks,
  ylim = c(0, max(abs(gene_counts))),
  beside = TRUE,
  main = "Differential Expressed Genes",
  cex.axis = 0.8,
  las = 2,  # <-- Set las = 2 to orient x_ticks vertically
  cex.names = 0.8
)

text(x = bp, y = -1, labels = x_ticks, srt = 90, adj = c(1, 0.5), cex = 0.8)

# Add a legend with larger font size
legend(
  "topright",
  legend = c("Up-regulated", "Down-regulated"),
  fill = c("red", "blue"),
  border = NA,
  cex = 0.6,
)

#######################################################
library(pheatmap)
gene_expression_matrix <- my_dge_list$counts
pheatmap(gene_expression_matrix)
