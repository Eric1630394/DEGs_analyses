# Problem 2: RNA-seq analysis 🧬
For this problem, we had several data files, such as BAM files📄, txt files with read counts and FPKM correction📄 and a Differential Expression analysis file📄. Some of the steps needed to solve the diverse questions and objectives we were asked to adress were solved with the following:

## 📁 FPKMS_plotting
- 📄 `Mutant_vs_WildType_results.txt`: file with normalized FPKMs for each replica, with calculated fold change (FC) and p-value.
- 💻 `Box_plot.R`: script that plots the box plot📈 of FPKMs in all the replicas and the average.
- 💻 `Dispersion_diagram.R`: script that plots the dispersion diagram📈 of FPKMs in all the replicas.
- 🎨 `Boxplot.png` and `Dispersion_diagram.png`: images of the Box plot and the Dispersions Diagram. 


## 📁 Volcano_plot
- 📄 `Mutant_vs_WildType_results.txt`: file with normalized FPKMs for each replica, with calculated fold change (FC) and p-value.
- 💻 `Volcano_plot.R`: script that creates the volcano plot📈 of the fold changes and the corresponding p-values.
- 🎨 `volcano_plot.png`: volcano plot for the Fold Change and p-adj values.  


## 📁 COG_terms_plotting
- 📄 `DEGs_woth_COGs.csv`: COG term for every differentially expressed gene associated to a locus_tag. 
- 💻 `COG_terms_plots.R`: script that mathches the COG terms to our DEGs, creates a barplot and computes Fisher's exact test on data. 
- 🎨 `Barplot_comparison_significance.png`: barplot for the different COG terms distinguished by DEGs and original count. 
