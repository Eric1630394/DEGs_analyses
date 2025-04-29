library("ggplot2")
library("tidyverse")
library("ggrepel")

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
