library("ggplot2")
library("tidyverse")
library("ggrepel")

### VOLCANO PLOT ### 
table <- read.table("Mutant_vs_WildType_results.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)

table$"log2(FoldChange)" <- NA
table$"-log10(padj)" <- NA
table$DEGs <- NA
table$Significance <- NA
table$Type <- NA

for (i in 1:nrow(table)) {
  value1 <- log2(table$FoldChange[i])
  value2 <- -log10(table$padj[i])
  table$`log2(FoldChange)`[i] <- value1
  table$`-log10(padj)`[i] <- value2
  if (table$`log2(FoldChange)`[i]>=1){
    table$DEGs[i] <- "Overexpressed"
  } else if (table$`log2(FoldChange)`[i]<=-1){
    table$DEGs[i] <- "Underexpressed"
  } else {
    table$DEGs[i] <- "Not_DEG"
  }
}    

for (i in 1:nrow(table)) {
  if (table$`-log10(padj)`[i]>=1.3){
    table$Significance[i] <- "Significant"
  } else {
    table$Significance[i] <- "Non_Significant"
  }
}   

for (i in 1:nrow(table)){
  if (table$DEGs[i] == "Overexpressed" && table$Significance[i] == "Significant"){
    table$Type[i] <- "Up-regulated"
  } else if (table$DEGs[i] == "Underexpressed" && table$Significance[i] == "Significant"){
    table$Type[i] <- "Down-regulated"
  } else {
    table$Type[i] <- "Not_DEG"
  }
}

volcano_plot <- ggplot(data = table, mapping = aes(x = log2(FoldChange), y = -log10(padj), colour = Type)) +
  geom_point(size=0.7) +
  geom_vline(xintercept = -1, linetype = "dashed", color = "black") +  # Vertical line at x=3
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "black") +
  labs(title = "Volcano_plot",
       x = "Log2(fold_change)",
       y = "-log10(padj)") +
  geom_text_repel(
    data = subset(table, abs(log2(FoldChange)) > 1 & -log10(padj) > 90),  # Optional filter
    aes(label = Locus_tag),  # Column with labels
    size = 3,
    color = "black",
    box.padding = 1,  # Adjust spacing
    max.overlaps = 15  # Increase if some labels are missing
  ) +
  scale_color_manual(values = c("green", "blue", "red")) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 15,hjust = 0.5),
    axis.title.y = element_text(size = 15,hjust = 0.5),
    plot.title = element_text(size = 20,hjust = 0.5),
    legend.text = element_text(size=10),
    legend.title = element_text(size=12),
    legend.position = c(0.85, 0.85),
    legend.background = element_rect(fill = "white", color = "gray")
  )
print(volcano_plot)
ggsave(filename = "volcano_plot.png", plot = volcano_plot, width = 10, height = 10, bg = "white")
