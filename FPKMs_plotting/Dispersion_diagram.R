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
