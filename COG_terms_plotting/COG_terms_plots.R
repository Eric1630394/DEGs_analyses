library(ggplot2)
library("tidyr")
library("dplyr")

COG_terms <- read.table("Comparison_whole_DEGs.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)

data_long <- COG_terms %>%
  tidyr::pivot_longer(
    cols = c(Whole.genome, DEGs),
    names_to = "Count_Type",
    values_to = "Count"
  )

#Barplot:
barplot <- ggplot(data_long, aes(x = X, y = Count, fill = Count_Type, alpha = Count_Type)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7,color="black",linewidth=0.5) +
  labs(
    title = "COG terms (Whole genome vs DEGs)",
    x = "COG_term",
    y = "Proportion to the whole genome"
  ) +
  theme_linedraw() +
  theme(
    legend.position = "top",
    legend.box = "horizontal",
    legend.spacing.x = unit(0.5, "cm"),
    plot.title = element_text(size=15,hjust=0.5),
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12),
  ) +
  guides(
    fill = guide_legend(title.position = "top"),
    alpha = guide_legend(title.position = "top")
  )

ggsave(filename = "Barplot_comparison.png", plot = barplot, width = 10, height = 10, bg = "white")

#Enrichment analysis:
cog_data <- data.frame(
  COG = c("C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Z"),
  Whole_genome = c(74,35,79,56,59,27,34,187,326,98,107,35,289,469,13,566,96,238,21,8,1,25,1),
  DEGs = c(48,19,40,28,33,13,17,125,144,39,51,15,147,235,8,290,46,108,14,4,1,3,1),
  Total_genes = 2844,
  Total_DEGs = 1429
)

Category <- c()
p_values <- c()

for (i in 1:nrow(cog_data)) {
  Category <- c(Category,cog_data$COG[i])
  contingency_table <- matrix(c(cog_data$DEGs[i],cog_data$Total_DEGs[i]-cog_data$DEGs[i],cog_data$Whole_genome[i],cog_data$Total_genes[i]-cog_data$Whole_genome[i]),nrow=2)
  print(contingency_table)
  test <- fisher.test(contingency_table, alternative = "greater")
  p_values <- c(p_values,test$p.value)
}

results <- data.frame(Category,p_values)
results$significant <- results$p_values < 0.05
results$significance_level <- cut(results$p_values,
                                 breaks = c(0, 0.001, 0.01, 0.05, 1),
                                 labels = c("***", "**", "*", "ns"))

write.csv(results, "COG_terms_statistics.csv", row.names = FALSE)
