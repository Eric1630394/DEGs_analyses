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

### BAR PLOT ###
subset_up_regulated <- subset(table, Type == "Up-regulated")
barplot_up <- ggplot(subset_up_regulated, aes(x = Locus_tag, y = log2(FoldChange), fill = Type)) +
  geom_bar(stat = "identity",color="black",linewidth=0.2) +
  theme_minimal() +
  labs(title = "Logarithmic fold change for Up-regulated genes", x = "Locus tags", y = "Values") +
  theme(
    plot.title = element_text(size = 10,hjust = 0.5),
    axis.text.x = element_text(angle = 75, hjust = 1, size = 5, margin = margin(t = 10, b = 40)),
  )

print(barplot_up)
ggsave(filename = "Barplot_up.png", plot = barplot_up, width = 10, height = 10, bg = "white",limitsize=FALSE)

subset_down_regulated <- subset(table, Type == "Down-regulated")
barplot_down <- ggplot(subset_down_regulated, aes(x = Locus_tag, y = log2(FoldChange), fill = Type)) +
  geom_bar(stat = "identity",color="black",linewidth=0.2) +
  theme_minimal() +
  labs(title = "Logarithmic fold change for Down-regulated genes", x = "Locus tags", y = "Values") +
  theme(
    plot.title = element_text(size = 10,hjust = 0.5),
    axis.text.x = element_text(angle = 75, hjust = 1, size = 5, margin = margin(t = 10, b = 40)),
  )

print(barplot_down)
ggsave(filename = "Barplot_down.png", plot = barplot_down, width = 49, height = 49, bg = "white",limitsize=FALSE)

write.csv(table, "Mutant_vs_WildType_results_v2.txt", row.names = FALSE)

### DISPERSION DIAGRAM ###
table_fpkm <- read.table("gene_count_fpkm.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
table_fpkm$Mutant_fpkm_average <- NA
table_fpkm$WildType_fpkm_average <- NA
table_fpkm$"log10(Mutant1_fpkm)" <- NA
table_fpkm$"log10(Mutant2_fpkm)" <- NA
table_fpkm$"log10(Mutant3_fpkm)" <- NA
table_fpkm$"log10(Mutant_fpkm_average)" <- NA
table_fpkm$"log10(WildType1_fpkm)" <- NA
table_fpkm$"log10(WildType2_fpkm)" <- NA
table_fpkm$"log10(WildType3_fpkm)" <- NA
table_fpkm$"log10(WildType_fpkm_average)" <- NA

for (i in 1:nrow(table_fpkm)) {
  # Mutant
  average_mutant <- (table_fpkm$Mutant_fpkm_replica1[i]+table_fpkm$Mutant_fpkm_replica2[i]+table_fpkm$Mutant_fpkm_replica3[i])/3
  table_fpkm$Mutant_fpkm_average[i] <- average_mutant
  log_average_mutant <- log10(average_mutant)
  table_fpkm$"log10(Mutant1_fpkm)"[i] <- log10(table_fpkm$Mutant_fpkm_replica1[i])
  table_fpkm$"log10(Mutant2_fpkm)"[i] <- log10(table_fpkm$Mutant_fpkm_replica2[i])
  table_fpkm$"log10(Mutant3_fpkm)"[i] <- log10(table_fpkm$Mutant_fpkm_replica3[i])
  table_fpkm$"log10(Mutant_fpkm_average)"[i] <- log_average_mutant
  # Wild_type
  average_wild_type <- (table_fpkm$WildType_fpkm_replica1[i]+table_fpkm$WildType_fpkm_replica2[i]+table_fpkm$WildType_fpkm_replica3[i])/3
  table_fpkm$WildType_fpkm_average[i] <- average_wild_type
  log_average_wt <- log10(average_wild_type)
  table_fpkm$"log10(WildType1_fpkm)"[i] <- log10(table_fpkm$WildType_fpkm_replica1[i])
  table_fpkm$"log10(WildType2_fpkm)"[i] <- log10(table_fpkm$WildType_fpkm_replica2[i])
  table_fpkm$"log10(WildType3_fpkm)"[i] <- log10(table_fpkm$WildType_fpkm_replica3[i])
  table_fpkm$"log10(WildType_fpkm_average)"[i] <- log_average_wt
}

dispersion_diagram <- ggplot(data = table_fpkm, mapping = aes(x = `log10(WildType_fpkm_average)`, y = `log10(Mutant_fpkm_average)`)) +
  geom_abline(intercept=0, slope=1, color="black", linetype="dashed",linewidth=0.5) +
  geom_point(size=1) + 
  labs(title = "Dispersion diagram (Mutant vs Wild Type)",
       x = "log10(WildType_fpkm_average)",
       y = "log10(Mutant_fpkm_average)") +
  theme_linedraw() +
  theme(
    plot.title = element_text(size = 15,hjust = 0.5),
    axis.title.x = element_text(size = 12),  # Add space below x-axis title
    axis.title.y = element_text(size = 12),   # Add space to the left of y-axis title
  )

ggsave(filename = "Dispersion_diagram.png", plot = dispersion_diagram, width = 8, height = 8, bg = "white",limitsize=FALSE)

write.csv(table_fpkm, "fpkm_log.txt", row.names = FALSE)

### BOX PLOT ###
long_data <- table_fpkm %>%
  select(starts_with("log10")) %>%
  pivot_longer(cols = everything(), names_to = "Sample", values_to = "log10_FPKM") %>%
  mutate(
    Type = ifelse(str_detect(Sample, "Mutant"), "Mutant", "WildType"),
    Data_Type = ifelse(str_detect(Sample, "average"), "Average", "Replicate"),
    Replicate = case_when(
      str_detect(Sample, "1_fpkm") ~ "Rep1",
      str_detect(Sample, "2_fpkm") ~ "Rep2",
      str_detect(Sample, "3_fpkm") ~ "Rep3",
      TRUE ~ "Average"
    ),
    Group = paste(Type, Replicate, sep = "_"),
    Group = factor(Group, levels = c(
      "WildType_Average", "WildType_Rep1",
      "WildType_Rep2", "WildType_Rep3",
      "Mutant_Average", "Mutant_Rep1",
      "Mutant_Rep2", "Mutant_Rep3"
    ))
  )

# Define custom fill colors: Different for Replicate vs. Average
fill_palette <- c(
  "Mutant.Replicate" = "slateblue",
  "WildType.Replicate" = "lightgreen",
  "Mutant.Average" = "navy",    
  "WildType.Average" = "forestgreen"       
)

box_plot <- ggplot(long_data, aes(x = Group, y = log10_FPKM)) +
  geom_boxplot(
    aes(fill = interaction(Type, Data_Type)),
    color = "black",                          
    alpha = 0.8                               
  ) +
  scale_fill_manual(
    name = "Sample Type",
    values = fill_palette,
    labels = c(
      "Mutant_Replicate" = "Mutant (Replicate)",
      "WildType_Replicate" = "WildType (Replicate)",
      "Mutant_Average" = "Mutant (Average)",
      "WildType_Average" = "WildType (Average)"
    )
  ) +
  labs(
    title = "Box plot on log10(fpkm) values",
    x = "Sample",
    y = "log10(FPKM)"
  ) +
  theme_linedraw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
  )

ggsave(filename = "Boxplot.png", plot = box_plot, width = 10, height = 10, bg = "white",limitsize=FALSE)

up <- subset(table, `log2(FoldChange)` > 1 & `-log10(padj)` > 1.3)
down <- subset(table, `log2(FoldChange)` < -1 & `-log10(padj)` > 1.3)
