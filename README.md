#circosplot
conda create -n syny
conda activate syny
run_syny.pl -a *.gbff.gz -o finaloutput directory
run_syny.pl -a *.gbff -g 0 1 5 -r A_solani -o $tika --circos all --labels names --label_size 20 --label_font semibold
run_syny.pl -a *.gbff -g 0 1 5 -e 1e-10 -r Alinariea25chrfinal -o tikafinal --circos all --no_skews --labels names --label_size 20 --label_font semibold --image_size 3000,3000 --output_format png --dpi 300
conda activate bio-env
nano ./fasta_to_gbff.py
./fasta_to_gbff.py
run_syny.pl -a *.gbff -o tika directory





# Alternaria-linaries
Rajan
quality control 
fastqc *.fastq -o qc_reports
# Load necessary libraries
library(ggplot2)
# Load necessary libraries
library(ggplot2)

# Load necessary libraries
library(ggplot2)
library(viridis)

# Create a data frame with your data
pathway_data <- data.frame(
  Pathway = c(
    "Protein families: Genetic information processing",
    "Genetic information processing",
    "Carbohydrate metabolism",
    "Protein families: Signaling and cellular processes",
    "Cellular processes",
    "Unclassified: Metabolism",
    "Protein families: Metabolism",
    "Amino acid metabolism",
    "Environmental information processing",
    "Lipid metabolism",
    "Energy metabolism",
    "Metabolism of cofactors and vitamins",
    "Glycan biosynthesis and metabolism",
    "Nucleotide metabolism",
    "Organismal systems",
    "Human diseases",
    "Metabolism of other amino acids",
    "Xenobiotics biodegradation and metabolism",
    "Metabolism of terpenoids and polypeptides",
    "Unclassified: Signaling and cellular processing",
    "Biosynthesis of other Secondary Metabolites",
    "Unclassified: Genetic information processing",
    "Unclassified"
  ),
  Number_of_genes = c(
    761, 715, 360, 255, 229, 199, 155, 170, 146, 
    123, 104, 102, 70, 62, 61, 50, 28, 25, 25, 
    21, 14, 8, 31
  )
)

# Create the bar graph
ggplot(pathway_data, aes(x = reorder(Pathway, -Number_of_genes), y = Number_of_Genes, fill = Number_of_Genes)) +
  geom_bar(stat = "identity") +
  labs(title = "Number of Genes by Biological Pathway",
       x = "Pathway",
       y = "Number of genes") +
  theme(
        axis.text.x = element_text(angle = 0, hjust = 1, family = "Times New Roman", size = 14, color = "black"),
        axis.text.y = element_text(colour = "black", size = 14),
        axis.title.y = element_text(family = "Times New Roman", size = 14, color = "black"),
        axis.title.x = element_text(family = "Times New Roman", size = 14, color = "black"),
        plot.title = element_text(family = "Times New Roman", size = 14, hjust = 0.5)) +
   scale_fill_gradientn(colours = pal) + 
  coord_flip() +  # Flip the coordinates to have pathways on the y-axis
  guides(fill = "none") +  # Remove the legend
  theme_minimal()



# Create the bar graph
ggplot(pathway_data, aes(x = reorder(Pathway, -Number_of_Genes), y = Number_of_Genes, fill = Number_of_Genes)) +
  geom_bar(stat = "identity") +
  labs(title = "Number of Genes by Biological Pathway",
       x = "Pathway",
       y = "Number of Genes") +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, family = "Times New Roman", size = 14, color = "black"),
        axis.text.x = element_text(angle = 0, hjust = 1, family = "Times New Roman", size = 14, color = "black"),
        axis.title.y = element_text(family = "Times New Roman", size = 14, color = "black"),
        axis.title.x = element_text(family = "Times New Roman", size = 14, color = "black"),
        plot.title = element_text(family = "Times New Roman", size = 14, hjust = 0.5, color = "black")) +
  scale_fill_gradientn(colours = pal) + 
  coord_flip() +  # Flip the coordinates to have pathways on the y-axis
  guides(fill = "none") +  # Remove the legend
  theme_minimal()

# Create the Cleveland dot plot
ggplot(pathway_data, aes(x = Number_of_genes, y = reorder(Pathway, -Number_of_genes), color = Number_of_genes, add = "segments", add.params = list(color = "lightgray", size = 2)), +
  geom_point(size = 3) + # Dots colored based on Number_of_Genes
  scale_color_gradient(low = "blue", high = "red")+
  labs(
    title = "Cleveland Dot Plot: Pathways vs. Number of genes",
    x = "Number of genes",
    y = "Pathway",
    color = "Number of genes"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(family = "Times New Roman", size = 12, color = "black"),
    axis.text.x = element_text(family = "Times New Roman", size = 12, color = "black"),
    axis.title = element_text(family = "Times New Roman", size = 12, color = "black"),
    plot.title = element_text(family = "Times New Roman", size = 12, color = "black", hjust = 0.5),
    legend.text = element_text(family = "Times New Roman", size = 12, color = "black"),
    legend.title = element_text(family = "Times New Roman", size = 12, color = "black"))
)

 

# Sort data in ascending order
pathway_data <- pathway_data[order(pathway_data$Number_of_genes), ]

# Create the Cleveland dot plot
ggplot(pathway_data, aes(x = Number_of_genes, y = reorder(Pathway, Number_of_genes))) +
  geom_point(size = 3, color = "black") + # Black dots
  labs(
    title = "Cleveland Dot Plot: Pathways vs. Number of Genes",
    x = "Number of genes",
    y = "Pathway"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(family = "Times New Roman", size = 12, color = "black"),
    axis.text.x = element_text(family = "Times New Roman", size = 12, color = "black"),
    axis.title = element_text(family = "Times New Roman", size = 12, color = "black"),
    plot.title = element_text(family = "Times New Roman", size = 12, color = "black", hjust = 0.5),
    legend.text = element_text(family = "Times New Roman", size = 12, color = "black"),
    legend.title = element_text(family = "Times New Roman", size = 12, color = "black")
  ) 



library(wesanderson)
install.packages("wesanderson")
pal <- wes_palette("Zissou1", 100, type = "continuous")


# Create a data frame with your data
go_data <- data.frame(
  GO = c(
    "rRNA metabolic process", "rRNA processing", "Ribosome biogenesis",
    "Ribonucleoprotein complex biogenesis", "RNA processing",
    "Cellular component biogenesis", "RNA metabolic process",
    "Cellular component organization or biogenesis", "Nucleic acid metabolic process",
    "Gene expression", "Nucleobase-containing compound metabolic process",
    "Heterocycle metabolic process", "Cellular component organization",
    "Cellular aromatic compound metabolic process", "Cellular nitrogen compound metabolic process",
    "Organic cyclic compound metabolic process", "Macromolecule metabolic process",
    "Nitrogen compound metabolic process", "Primary metabolic process",
    "Cellular metabolic process", "Nuclear chromosome", "Nucleolus", 
    "Nuclear lumen", "Nucleoplasm", "Ribonucleoprotein complex", 
    "Protein-containing complex", "Catalytic complex", "Membrane-enclosed lumen", 
    "Organelle lumen", "Intracellular organelle lumen", "Non-membrane-bounded organelle", 
    "Intracellular non-membrane-bounded organelle", "Organelle membrane", 
    "Endomembrane system", "Nucleus", "Intracellular membrane-bounded organelle", 
    "Membrane-bounded organelle", "Intracellular organelle", "Organelle", 
    "Cytoplasm", "Translation initiation factor activity", "GTPase binding", 
    "Translation regulator activity", "Translation regulator activity, nucleic acid binding", 
    "Translation factor activity, rna binding", "Lipid binding", "RNA binding", 
    "ATPase activity", "Protein binding", "Nucleoside-triphosphatase activity", 
    "Ribonucleotide binding", "Purine nucleotide binding", 
    "Purine ribonucleotide binding", "Purine ribonucleoside triphosphate binding", 
    "Adenyl nucleotide binding", "Nucleic acid binding", "Anion binding", 
    "Organic cyclic compound binding", "Heterocyclic compound binding", 
    "Binding"
  ),
  Gene_number = c(
    183, 180, 258, 310, 451, 599, 696, 1262, 932, 
    895, 1125, 1245, 1032, 1258, 1479, 1343, 2139, 
    2508, 3005, 3069, 190, 208, 557, 203, 365, 
    1435, 500, 773, 773, 773, 861, 861, 866, 
    783, 1760, 3352, 3552, 3679, 3808, 3717, 40, 
    68, 60, 57, 54, 86, 450, 222, 571, 311, 
    743, 718, 715, 712, 617, 886, 1064, 1936, 
    1923, 3193
  ),
  Category = c(
    rep("Biological process", 20),
    rep("Cellular components", 20),
    rep("Molecular function", 20)
  )
)

# Create a factor for GO terms within categories
go_data$GO <- factor(go_data$GO, levels = go_data$GO)

# Create the bar graph
ggplot(go_data, aes(x = Category, y = Gene_number, fill = Category)) +
  geom_bar(stat = "identity", position = position_dodge()) +  # Group by category
  geom_bar(stat = "identity", aes(x = GO, group = Category), position = position_dodge()) +  # Display GO terms
  scale_fill_manual(values = c("Biological process" = "skyblue", 
                               "Cellular components" = "lightgreen", 
                               "Molecular function" = "salmon")) +  # Manual color mapping
  labs(title = "Gene Numbers by GO Terms",
       x = "GO terms",
       y = "Number of genes") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, family = "Times New Roman", size = 12, color = "black"),
    axis.text.y = element_text(family = "Times New Roman", size = 9, color = "black"),
    axis.title.y = element_text(family = "Times New Roman", size = 9, color = "black"),
    axis.title.x = element_text(family = "Times New Roman", size = 9, color = "black"),
    plot.title = element_text(family = "Times New Roman", size = 9, hjust = 0.5, color = "black")
  ) +
  coord_flip() +  # Flip the coordinates to have categories on the y-axis
  guides(fill = guide_legend(title = "Category"))  # Add legend if needed



library(ggplot2)
library(dbplyr)
library(reshape2)  
# Data setup
data <- data.frame(
  Species = c("A. alternata SRC1k2F", "A. aborecenes FERA675", 
              "A. tomatophila BMP32002", "A. linariae 25", 
              "A. solani NL30003", "A. brassicicola Abra43", 
              "A. dauci A2006"),
  Proteins = c(13466, 12924, 12811, 11768, 11026, 10688, 10026),
  Clusters = c(11729, 11630, 13001, 11152, 10730, 8812, 9545),
  Singletons = c(1305, 831, 2421, 347, 110, 1540, 146))

# Reshape data for ggplot
data_melt <- melt(data, id.vars = "Species")

# Add position for text labels
data_melt <- data_melt %>%
  group_by(Species) %>%
  mutate(position = cumsum(value) - 0.5 * value)

# Define custom labels with italicized species names
species_labels <- expression(
  italic("A. alternata") ~ "SRC1k2F",
  italic("A. aborecenes") ~ "FERA675",
  italic("A. tomatophila") ~ "BMP32002",
  italic("A. linariae") ~ "25",
  italic("A. solani") ~ "NL30003",
  italic("A. brassicicola") ~ "Abra43",
  italic("A. dauci") ~ "A2006"
)



# Data setup
data <- data.frame(
  Species = c("A. alternata SRC1k2F", "A. aborecenes FERA675", 
              "A. tomatophila BMP32002", "A. linariea 25", 
              "A. solani NL30003", "A. brassicicola Abra43", 
              "A. dauci A2006"),
  Proteins = c(12048, 12895, 13001, 11768, 11026, 12456, 10026),
  Clusters = c(11729, 11630, 12811, 11152, 10730, 8812, 9545),
  Singletons = c(1305, 831, 2421, 347, 110, 1540, 146)
)

# Reshape data for ggplot
data_melt <- melt(data, id.vars = "Species")

# Reorder levels to place Singletons at the end of each stack
data_melt$variable <- factor(data_melt$variable, levels = c("clusters", "Proteins", "Singletons"))

# Add position for text labels
data_melt <- data_melt %>%
  group_by(Species) %>%
  mutate(position = cumsum(value) - 0.5 * value)

# Define custom labels with italicized species names
species_labels <- expression(
  italic("A. alternata") ~ "SRC1k2F",
  italic("A. aborecenes") ~ "FERA675",
  italic("A. tomatophila") ~ "BMP32002",
  italic("A. linariae") ~ "25",
  italic("A. solani") ~ "NL30003",
  italic("A. brassicicola") ~ "Abra43",
  italic("A. dauci") ~ "A2006"
)

# Plot with values in each stack and species names on the y-axis
ggplot(data_melt, aes(y = Species, x = value, fill = variable)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = value, x = position), color = "black", size = 4) +
  scale_y_discrete(labels = species_labels) +
  labs(y = "Species", x = "Counts", fill = "Category") +
  theme(
    axis.text.x = element_text(angle = 1, hjust = 1, family = "Times New Roman", size = 12, color = "black"),
    axis.text.y = element_text(family = "Times New Roman", size = 12, color = "black"),
    axis.title.y = element_text(family = "Times New Roman", size = 12, color = "black"),
    axis.title.x = element_text(family = "Times New Roman", size = 12, color = "black"),
    plot.title = element_text(family = "Times New Roman", size = 12, hjust = 0.5, color = "black")
  ) +
  ggtitle("Proteins, Clusters, and Singletons for Different Species")


# Data setup
data <- data.frame(
  Species = c("A. tomatophila BMP32002", "A. solani NL30003", "A. linariae 25", "A. dauci A2006", "A. brassicicola Abra43", "A. alternata SRC1k2F", "A. aborecenes FERA675"),
  Singletons = c(1305, 831, 2421, 347, 110, 1540, 146),
  Clusters = c(12048, 12895, 13001, 11768, 11026, 12456, 10026),
  Proteins = c(11729, 11630, 12811, 11152, 10730, 8812, 9545)
)

# Reshape data for ggplot
data_melt <- melt(data, id.vars = "Species")

# Reorder levels to place Singletons at the end of each stack
data_melt$variable <- factor(data_melt$variable, levels = c("Clusters", "Proteins", "Singletons"))

# Add position for text labels
data_melt <- data_melt %>%
  group_by(Species) %>%
  mutate(position = cumsum(value) - 0.5 * value)

# Modify Species labels to be italicized using ggplot parsing
species_labels <- c(
  "italic('A. alternata')~'SRC1k2F'",
  "italic('A. aborecenes')~'FERA675'",
  "italic('A. tomatophila')~'BMP32002'",
  "italic('A. linariae')~'25'",
  "italic('A. solani')~'NL30003'",
  "italic('A. brassicicola')~'Abra43'",
  "italic('A. dauci')~'A2006'"
)

# Plot with values in each stack and italicized species names on the y-axis
ggplot(data_melt, aes(y = Species, x = value, fill = variable)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = value, x = position), color = "black", size = 4) +
  scale_y_discrete(labels = parse(text = species_labels)) +
  labs(y = "Species", x = "Counts", fill = "Category") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 1, hjust = 1, family = "Times New Roman", size = 12, color = "black"),
    axis.text.y = element_text(family = "Times New Roman", size = 12, color = "black"),
    axis.title.y = element_text(family = "Times New Roman", size = 12, color = "black"),
    axis.title.x = element_text(family = "Times New Roman", size = 12, color = "black"),
    plot.title = element_text(family = "Times New Roman", size = 12, hjust = 0.5, color = "black")
  ) +
  ggtitle("Proteins, Clusters, and Singletons for Different Species")

# Load ggplot2 library
library(ggplot2)

# Create the data frame
data <- data.frame(
  Field1 = c("Apoplastic", "Apoplastic/Cytoplasmic", "Cytoplasmic"),
  Count = c(131, 47, 25)
)

# Plot the bar graph
ggplot(data, aes(x = Field1, y = Count, fill = Field1)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Counts of Field1 Categories", x = "Effector type", y = "Number of effectors") +
  scale_fill_manual(values = c("brown", "cyan", "orange")) +
theme(
  axis.text.x = element_text(angle = 45, hjust = 1, family = "Times New Roman", size = 12, color = "black"),
  axis.text.y = element_text(family = "Times New Roman", size = 12, color = "black"),
  axis.title.y = element_text(family = "Times New Roman", size = 12, color = "black"),
  axis.title.x = element_text(family = "Times New Roman", size = 12, color = "black")
)

# Define the similarity matrix
similarity_matrix <- matrix(c(
  100, 93, 54, 65, 65, 69, 69,
  91, 100, 55, 66, 66, 70, 71,
  49, 52, 100, 44, 45, 47, 47,
  68, 69, 50, 100, 87, 91, 92,
  67, 70, 50, 87, 100, 93, 99,
  66, 69, 50, 86, 88, 100, 94,
  66, 69, 49, 85, 92, 92, 100
), nrow = 7, byrow = TRUE)

# Set row and column names
rownames(similarity_matrix) <- c("A. arborescens FERA675", "A. alternata SRC1lrK2f", "A. brasicicola Abra43", 
                                 "A. dauci A2016", "A. linariae 25", "A. solani NL03003", "A. tomatophila BMP2032")
colnames(similarity_matrix) <- c("A. arborescens FERA675", "A. alternata SRC1lrK2f", "A. brasicicola Abra43", 
                                 "A. dauci A2016", "A. linariae 25", "A. solani NL03003", "A. tomatophila BMP2032")
)

# Load the required libraries
library(pheatmap)
library(RColorBrewer)

# Create the heatmap with clustering
pheatmap(similarity_matrix,
         clustering_distance_rows = "euclidean", # Euclidean distance for rows
         clustering_distance_cols = "euclidean", # Euclidean distance for columns
         clustering_method = "complete",  # Method of clustering (complete linkage)
         color = colorRampPalette(brewer.pal(9, "Blues"))(100),  # Color scheme
         scale = "none",  # No scaling of the data
         main = "Similarity Heatmap of Alternaria Strains",  # Title of the plot
         fontsize = 10)  # Font size for the labels



# Load the required libraries
library(pheatmap)
library(RColorBrewer)

# Define the similarity matrix
similarity_matrix <- matrix(c(
  100, 97, 86, 88, 87, 88, 88,
  97, 100, 86, 88, 88, 88, 88,
  86, 86, 100, 86, 86, 86, 86,
  88, 88, 86, 100, 94, 94, 94,
  88, 88, 86, 94, 100, 97, 99,
  88, 88, 86, 94, 98, 100, 99,
  88, 88, 86, 94, 99, 97, 100
), nrow = 7, byrow = TRUE)

# Set row and column names (species names)
species_names <- c("A. arborescens FERA675", "A. alternata SRC1lrK2f", "A. brasicicola Abra43", 
                   "A. dauci A2016", "A. solani NL03003", "A. tomatophila BMP2032", "A. linariae 25")

rownames(similarity_matrix) <- species_names
colnames(similarity_matrix) <- species_names

italic_colnames <- c(
  expression(italic("A. arborescens") ~ FERA675),
  expression(italic("A. alternata") ~ SRC1lrK2f),
  expression(italic("A. brassicicola") ~ Abra43),
  expression(italic("A. dauci") ~ A2016),
  expression(italic("A. linariae") ~ 25),
  expression(italic("A. solani") ~ NL03003),
  expression(italic("A. tomatophila") ~ BMP2032)
)

italic_rownames <- c(
  expression(italic("A. arborescens") ~ FERA675),
  expression(italic("A. alternata") ~ SRC1lrK2f),
  expression(italic("A. brassicicola") ~ Abra43),
  expression(italic("A. dauci") ~ A2016),
  expression(italic("A. linariae") ~ 25),
  expression(italic("A. solani") ~ NL03003),
  expression(italic("A. tomatophila") ~ BMP2032)
)



# Round the values in the matrix to remove the decimals (.00)
similarity_matrix <- round(similarity_matrix)

# Italicize species names
italic_rownames <- sapply(species_names, function(x) paste0("*", x, "*"))
italic_colnames <- sapply(species_names, function(x) paste0("*", x, "*"))
custom_red_palette <- colorRampPalette(c("white", "pink"))(200)  # 200 shades from white to dark red
# Create the heatmap with customization
pheatmap(similarity_matrix,
         color = custom_red_palette,
         scale = "none",  # No scaling of the data
         main = "Similarity Heatmap of Alternaria Strains",  # Title of the plot
         fontsize = 12,  # Font size for the labels
         angle_col = 45,  # Rotate column labels for better readability
         angle_row = 0,  # Keep row labels horizontal
         show_rownames = TRUE,  # Show row names (species)
         show_colnames = TRUE,  # Show column names (species)
         legend = TRUE,  # Show the legend
         annotation_legend = FALSE,  # Disable annotation legend
         cluster_rows = F,  # Do not cluster rows
         cluster_cols = F,  # Do not cluster columns
         fontfamily = "Times", # Times New Roman font
         display_numbers = T, # Show numbers in cells
         labels_col = italic_colnames,  # Italicize column names
         labels_row = italic_rownames,  # Italicize row names
         cellwidth = 30, # Adjust cell width
         cellheight = 30, # Adjust cell height
         number_format = "%.0f"
)

# Load libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr) # Ensure tidyr is loaded for pivot_longer

# Create data frame
data <- data.frame(
  Species = c("A. solani NL03003", "A. linariae 25", "A. tomatophila BMP2032", "A. alternata SRC1lrK2f"),
  Apoplastic = c(157, 154, 102, 163),
  Cytoplasmic = c(32, 27, 30, 33),
  Apoplastic_cytoplasmic = c(42, 40, 39, 45)
)


# Add totals row
data <- data %>% 
  mutate(Total = Apoplastic + Cytoplasmic + Apoplastic_cytoplasmic)

# Reshape data for ggplot
data_long <- data %>%
  select(Species, Apoplastic, Cytoplasmic, Apoplastic_cytoplasmic) %>%
  pivot_longer(cols = c(Apoplastic, Cytoplasmic, Apoplastic_cytoplasmic),
               names_to = "EffectorType", values_to = "Count")

# Define custom labels with italicized species names
species_labels <- expression(
  italic("A. alternata") ~ "SRC1k2F",
  italic("A. linariae") ~ "25",
  italic("A. solani") ~ "NL30003",
  italic("A. tomatophila") ~ "BMP32002"
)

# Plot
ggplot(data_long, aes(x = Species, y = Count, fill = EffectorType)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(labels = parse(text = species_labels))+
  scale_fill_manual(values = c("Apoplastic" = "brown", "Cytoplasmic" = "orange", "Apoplastic_cytoplasmic" = "cyan")) +
  labs(y = "Number of effectors", x = NULL) +
  theme_minimal(base_size = 12, base_family = "Times New Roman") +
    theme(
      axis.text.x = element_text(angle = 1, hjust = 1, family = "Times New Roman", size = 12, color = "black"),
      axis.text.y = element_text(family = "Times New Roman", size = 12, color = "black"),
      axis.title.y = element_text(family = "Times New Roman", size = 12, color = "black"),
      axis.title.x = element_text(family = "Times New Roman", size = 12, color = "black"),
    legend.title = element_blank(),
    legend.position = "right" # Position legend to the right
  ) +
  coord_flip()
# Load necessary libraries
library(ggplot2)
library(tidyr)  # For data reshaping

# Load necessary libraries
library(ggplot2)
library(tidyr)  # For data reshaping

# Create the data frame
data <- data.frame(
  Enzyme = c("AA6", "AA8", "CBM42", "GH45", "GH53", "GH55", "GT25", "GT50", "GT59", "GT76", "PL11", "PL26", "PL42", "CBM36", "CBM37", "CBM38", "CBM39", "CBM63", "GH33"),
  A_linariae_25 = c(0, 0, 1, 1, 1, 0, 4, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0),
  A_tomatophila_BMP2032 = c(0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0),
  A_solani_NL03003 = c(1, 1, 1, 2, 1, 4, 4, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0),
  A_alternata_SRC1lrK2f = c(1, 1, 1, 2, 1, 4, 4, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1),
  A_daucii_CBS107_38 = c(0, 1, 1, 2, 1, 4, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0),
  A_brassicicola_Abra43 = c(1, 1, 1, 2, 1, 4, 3, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0),
  A_arborescens_FERA675 = c(1, 1, 1, 2, 1, 4, 3, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0)
)

# Reshape data to long format for ggplot
data_long <- data %>%
  pivot_longer(cols = -Enzyme, names_to = "Species", values_to = "Count")

# Define labels with italicized species names
species_labels <- c(
  A_linariae_25 = expression(italic("A. linariae") ~ "25"),
  A_tomatophila_BMP2032 = expression(italic("A. tomatophila") ~ "BMP2032"),
  A_solani_NL03003 = expression(italic("A. solani") ~ "NL03003"),
  A_alternata_SRC1lrK2f = expression(italic("A. alternata") ~ "SRC1lrK2f"),
  A_daucii_CBS107_38 = expression(italic("A. daucii") ~ "CBS107.38"),
  A_brassicicola_Abra43 = expression(italic("A. brassicicola") ~ "Abra43"),
  A_arborescens_FERA675 = expression(italic("A. arborescens") ~ "FERA675")
)

# Plot the bubble chart
ggplot(data_long, aes(x = Species, y = Enzyme, size = Count)) +
  geom_point(color = "skyblue", alpha = 0.7) +
  scale_size(range = c(1, 10)) +  # Adjust bubble size
  scale_x_discrete(labels = species_labels) +  # Apply italicized labels
  labs(x = "Species", y = "Enzyme", size = "Count") +
  theme_minimal(base_size = 12, base_family = "Times New Roman") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, family = "Times New Roman", color = "black"),
    axis.text.y = element_text(family = "Times New Roman", color = "black"),
    axis.title = element_text(family = "Times New Roman", color = "black")
  ) +
  coord_flip()

# Load necessary libraries
library(ggplot2)
library(tidyr)  # For data reshaping

# Create the data frame
data <- data.frame(
  Enzyme = c("AA6", "AA8", "CBM42", "GH45", "GH53", "GH55", "GT25", "GT50", "GT59", "GT76", "PL11", "PL26", "PL42", "CBM36", "CBM37", "CBM38", "CBM39", "CBM63", "GH33"),
  A_linariae_25 = c(0, 0, 1, 1, 1, 0, 4, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0),
  A_tomatophila_BMP2032 = c(0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0),
  A_solani_NL03003 = c(1, 1, 1, 2, 1, 4, 4, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0),
  A_alternata_SRC1lrK2f = c(1, 1, 1, 2, 1, 4, 4, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1),
  A_daucii_CBS107_38 = c(0, 1, 1, 2, 1, 4, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0),
  A_brassicicola_Abra43 = c(1, 1, 1, 2, 1, 4, 3, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0),
  A_arborescens_FERA675 = c(1, 1, 1, 2, 1, 4, 3, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0)
)

# Reshape data to long format for ggplot
data_long <- data %>%
  pivot_longer(cols = -Enzyme, names_to = "Species", values_to = "Count")

# Define labels with italicized species names
species_labels <- c(
  A_linariae_25 = expression(italic("A. linariae") ~ "25"),
  A_tomatophila_BMP2032 = expression(italic("A. tomatophila") ~ "BMP2032"),
  A_solani_NL03003 = expression(italic("A. solani") ~ "NL03003"),
  A_alternata_SRC1lrK2f = expression(italic("A. alternata") ~ "SRC1lrK2f"),
  A_daucii_CBS107_38 = expression(italic("A. daucii") ~ "CBS107.38"),
  A_brassicicola_Abra43 = expression(italic("A. brassicicola") ~ "Abra43"),
  A_arborescens_FERA675 = expression(italic("A. arborescens") ~ "FERA675")
)

# Plot the bubble chart
ggplot(data_long, aes(x = Species, y = Enzyme, size = Count)) +
  geom_point(color = "orange", alpha = 0.7) +
  scale_size(range = c(1, 10)) +  # Adjust bubble size
  scale_x_discrete(labels = species_labels) +  # Apply italicized labels for species
  labs(x = "Species", y = "CAZyme", size = "Count") +
  theme_minimal(base_size = 12, base_family = "Times New Roman") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, family = "Times New Roman", size = 12, color = "black"),
    axis.text.y = element_text(family = "Times New Roman", size = 12,color = "black"),
    axis.title = element_text(family = "Times New Roman",size = 12, color = "black")
  )
# Load necessary libraries
library(ggplot2)
library(dplyr)
# Sample data frame based on your provided data
data <- data.frame(
  GO = c("rRNA metabolic process", "rRNA processing", "Ribosome biogenesis", 
         "Ribonucleoprotein complex biogenesis", "RNA processing", 
         "Cellular component biogenesis", "RNA metabolic process", 
         "Cellular component organization or biogenesis", 
         "Nucleic acid metabolic process", "Gene expression", 
         "Nucleobase-containing compound metabolic process", 
         "Heterocycle metabolic process", "Cellular component organization", 
         "Cellular aromatic compound metabolic process", 
         "Cellular nitrogen compound metabolic process", 
         "Organic cyclic compound metabolic process", 
         "Macromolecule metabolic process", "Nitrogen compound metabolic process", 
         "Primary metabolic process", "Cellular metabolic process", 
         "Nuclear chromosome", "Nucleolus", "Nuclear lumen", "Nucleoplasm", 
         "Ribonucleoprotein complex", "Protein-containing complex", 
         "Catalytic complex", "Membrane-enclosed lumen", "Organelle lumen", 
         "Intracellular organelle lumen", "Non-membrane-bounded organelle", 
         "Intracellular non-membrane-bounded organelle", "Organelle membrane", 
         "Endomembrane system", "Nucleus", "Intracellular membrane-bounded organelle", 
         "Membrane-bounded organelle", "Intracellular organelle", "Organelle", 
         "Cytoplasm", "Translation initiation factor activity", "GTPase binding", 
         "Translation regulator activity", "Translation regulator activity, nucleic acid binding", 
         "Translation factor activity, rna binding", "Lipid binding", "RNA binding", 
         "ATPase activity", "Protein binding", "Nucleoside-triphosphatase activity", 
         "Ribonucleotide binding", "Purine nucleotide binding", 
         "Purine ribonucleotide binding", "Purine ribonucleoside triphosphate binding", 
         "Adenyl nucleotide binding", "Nucleic acid binding", "Anion binding", 
         "Organic cyclic compound binding", "Heterocyclic compound binding", 
         "Binding"),
  Gene_number = c(183, 180, 258, 310, 451, 599, 696, 1262, 932, 895, 1125, 1245, 
                  1032, 1258, 1479, 1343, 2139, 2508, 3005, 3069, 190, 208, 557, 203, 
                  365, 1435, 500, 773, 773, 773, 861, 861, 866, 783, 1760, 3352, 3552, 
                  3679, 3808, 3717, 40, 68, 60, 57, 54, 86, 450, 222, 571, 311, 743, 
                  718, 715, 712, 617, 886, 1064, 1936, 1923, 3193),
  Color = c(rep("Biological process", 20), rep("Cellular component", 20), rep("Molecular function", 20))
)

# Convert the Color column to a factor and order it explicitly
data$Color <- factor(data$Color, levels = c("Biological process", "Cellular component", "Molecular function"))

# Sort the data by Color and Gene_number (optional: descending order of Gene_number)
data <- data[order(data$Color, data$Gene_number, decreasing = TRUE), ]
# Plot the bar graph, sorted by Color
ggplot(data, aes(x = Gene_number, y = GO, fill = Color)) + 
  geom_bar(stat = "identity") + 
  coord_flip() +  # Flip the axes to make GO terms readable on Y-axis
  theme_minimal() + 
  theme(
    text = element_text(family = "Times New Roman", size = 12, color = "black"), 
    axis.text.y = element_text(size = 10, color = "black"),  # Adjust size of GO term labels
    axis.title.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 12, color = "black")
  ) +
  labs(x = "GO terms", y = "Gene number") +
  scale_fill_manual(values = c("Biological process" = "blue", 
                               "Cellular components" = "green", 
                               "Molecular function" = "red")) +
  theme(axis.text.x = element_text(angle = 50))  # Ensure text is horizontal for better readability


 order within each group)
data <- data[order(data$Color, data$Gene_number), ]

# Create the bar plot
ggbarplot(
  data, 
  x = "GO", 
  y = "Gene_number", 
  fill = "Color", 
  color = "black",  # Add border around bars
  add = "segments",                             # Add segments from y = 0 to dots
  sort.val = "none",  # Disable automatic sorting, as data is already sorted
  sort.by.groups = FALSE,  # Retain our custom sorting
  x.text.angle = 0,  # Keep X-axis text horizontal
  orientation = "horiz",# Flip axes to make GO terms on Y-axis
  ggtheme = theme_pubr()
) +
  theme_minimal() +
  theme(
    text = element_text(family = "Times New Roman", size = 12, color = "black"), 
    axis.text.y = element_text(size = 12, color = "black"), 
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
  ) +
  labs(x = "GO terms", y = "Gene number") +
  scale_fill_manual(values = c("Biological process" = "blue", 
                               "Cellular component" = "green", 
                               "Molecular function" = "red")) +
  guides(fill = guide_legend(title = "GO category"))

