library(tidyverse)
library(tidyr)

library(ggtree)
library(treeio)
library(ape)

library(gridExtra)
library(cowplot)

### Plots settings
# theme code adapted from https://rpubs.com/koundy/71792 
base_size=14
base_family="helvetica"

# Created this vector based on tree order
genome_order <- rev(c("Pleery1", "Pletuo1", "Pleflo1", "Plecor1", "Plepla1", "PleostPC15", "PleostPC9", "PlepulLGAM",  "PlepulSS2",
                  "PlepulSS13", "PlepulSS5", "Plesap1", "Pledja1", "Plesal1","Plecit1","Plegig1","Pletub1"))

#### 1. Tree ####
# Code adapted from following
# https://guangchuangyu.github.io/ggtree-book
# https://yulab-smu.top/treedata-book/index.html
# I flip some clades order following Li et al. (2020) reference tree

tree <- read.iqtree("data/genomic/concat_onlyUFBOOT.treefile")

# Turn to phylo object and root
tree_phylo <- as.phylo(tree)
outgroup_tips <- c("H_portegna_HKAS74040", "H_unguicularis_HKAS90443")
ougroup_node <- getMRCA(tree_phylo, outgroup_tips)
tree_rooted <- root(tree_phylo, node = ougroup_node, resolve.root = TRUE)

# Turn it to ggtree object
p <- ggtree(tree_rooted)
# bin UFBOOT (>95 is confident based on IQTREE documentation)
p$data <- p$data %>%
  mutate(ufboot = as.numeric(label),
         support_cat = case_when(ufboot >= 95 ~ "95–100", ufboot >= 80 ~ "80–94", TRUE ~ "<80"))

# Fix names for display
p$data$label <- gsub("_", " ", p$data$label ) %>% gsub("P ", "P. ", .) %>% gsub("H ", "H. ", .) %>% 
  gsub(" sp", " sp. ", .) %>% gsub("var fer", "var ferulae", .) %>% gsub("var ery", "var eryngii", .) %>% 
  gsub("var ela", "var elaeoselini", .) %>%  gsub("P. sp.", "Pleurotus sp.", .)

tree_plot <- p %>% # if you want to color the branches , aes(color = as.numeric(label))
  # here I flip the pulmonarius and ostreatus clades for comparison with Li et al. (), no real reason
  flip(73, 75) %>%
  # here flip the LGAM and italian pulmonarius clades so that the three Italian pulmonarius are consecutive
  flip(117, 118) +
  # if you want to display values: geom_text2(aes(subset = !isTip & as.numeric(label) >= 80, label = label), size = 2.5, hjust = -0.3, color = "black") +
  geom_point2(aes(subset = !isTip, color = support_cat), size = 1.5) +
  # fontface makes them all italics
  geom_tiplab(size = 3, color = "black", family = base_family, fontface = 3) +
  geom_treescale() + 
  # if you don't want to bin the UFBOOT could use a continuous theme 
  #scale_color_viridis_c( name = "UFBOOT", limits = c(18, 100), oob = scales::squish) +
  scale_color_manual(name = "UFBoot", values = c("95–100" = "#0072B2","80–94"  = "#E69F00","<80" = "grey70")) +
  theme(legend.position = "right", legend.text = element_text(size = base_size), 
        legend.key.size= unit(0.6, "cm"), legend.title = element_text(size = base_size))

tree_plot


#### 2. BUSCO plot ####
# Proteome source indicates where the predicted proteome is downloaded from
# In case there was no available proteome it was predicted with Genemark-ES optimized for fungal genomes 

df <- read_csv("data/genomic/pleurotus_genomes.csv", show_col_types = FALSE)
df <- df %>% 
  filter(!is.na(`Proteome source`)) %>%
  mutate(Complete = as.numeric(sub("^C:([0-9\\.]+)%.*", "\\1",  `BUSCO proteome`)),
         Single = as.numeric(sub(".*S:([0-9\\.]+)%.*", "\\1",  `BUSCO proteome`)),
         Duplicate = as.numeric(sub(".*D:([0-9\\.]+)%.*", "\\1",  `BUSCO proteome`)),
         Fragmented = as.numeric(sub(".*F:([0-9\\.]+)%.*", "\\1",  `BUSCO proteome`)),
         Missing = as.numeric(sub(".*M:([0-9\\.]+)%.*", "\\1", `BUSCO proteome`)))
busco_df <- df %>%
  select(ShortName, Single, Duplicate, Fragmented, Missing) %>%
  pivot_longer(cols = c(Single, Duplicate, Fragmented, Missing), names_to = "Type", values_to = "Value") %>%
  mutate(Type = factor(Type, levels = c("Missing" ,"Fragmented", "Duplicate",  "Single" )))  # Set order

busco_df$ShortName <- factor(busco_df$ShortName, levels = genome_order)  

busco_plot <- ggplot(busco_df, aes(x = ShortName, y = Value, fill = Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Single" = "#0072B2",  "Duplicate" = "#56B4E9", "Fragmented" = "#E69F00" , "Missing" = "#D55E00")) +
  coord_flip() +
  theme_minimal() +
  labs(title = "BUSCO", y = "Percentage", x = "") + 
  theme(
    text = element_text(size = base_size, family = base_family ),
    plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
    axis.title = element_text(face = "bold",size = rel(1)),
    axis.title.y = element_blank(), #axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    axis.title.x = element_text(vjust = -0.2),
    legend.key = element_rect(colour = NA),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.size= unit(0.6, "cm"),
    legend.spacing  = unit(0, "cm"),
    legend.title = element_blank(),
    plot.margin=unit(c(10,5,5,5),"mm")
    )
busco_plot



#### 3. OrthoFinder ####
# Get count of proteins in Core, Accessory and Signletons
# Core will be defined as present in orthogroups with n-1 genomes 

### 3.1 Combine orthogroups and unassigned groups
gene_counts <- read.csv("data/genomic/Orthogroups.GeneCount.tsv.gz", sep = "\t")
# Read the unassigned subpart and make it numeric, assign 0 to empty strings and 1 to non-empty strings
unassigned <- read.csv("data/genomic/Orthogroups_UnassignedGenes.tsv.gz", sep = "\t")
unassigned <- unassigned %>%
  mutate(across(-1, ~ ifelse(. == "", 0, 1))) %>%
  mutate(Total = 1)
#combine them
gene_counts_all<- rbind(gene_counts, unassigned)
n_orthogroups_all <- nrow(gene_counts_all)


### 3.2 Some numbers
# how many orthogroups?
n_orthogroups <- nrow(gene_counts)
n_orthogroups

#how many strains?
cols_to_exclude<- c("Orthogroup", "Total", "number_genomes")
strains_only <- gene_counts_all[,!names(gene_counts_all) %in% cols_to_exclude]
strain_names <- names(strains_only)
ngenomes <- length(unique(strain_names)) 

# Get genome count for each group
gene_counts_all$number_genomes <- rowSums( strains_only > 0)

# Assign to core, accesory and singletons
cutoff <- ngenomes - 1
gene_counts_all <- gene_counts_all %>% mutate(Status = ifelse(number_genomes >= cutoff,  "Core", "Accessory"),
                                              Status = ifelse(number_genomes == 1,  "Singleton", Status))

paste0(round(100 * nrow(gene_counts_all %>% filter(Status == "Core")) / n_orthogroups,2), " % are core orthogroups")
paste0(round(100 * nrow(gene_counts_all %>% filter(Status == "Accessory")) / n_orthogroups,2), " % are accessory orthogroups")

### Check changes if you remove the low part busco proteomes
low_buscos <- c("Pledja1","Plesap1", "Plecor1" )
gene_counts_all_noLowBuscos <- gene_counts_all[,!names(gene_counts_all) %in% c(cols_to_exclude, low_buscos, "Status")]
gene_counts_all_noLowBuscos$number_genomes <- rowSums( gene_counts_all_noLowBuscos > 0)
cutoff_noLowBuscos <- cutoff - length(low_buscos)
gene_counts_all_noLowBuscos <- gene_counts_all_noLowBuscos %>% mutate(Status = ifelse(number_genomes >= cutoff_noLowBuscos,  "Core", "Accessory"),
                                              Status = ifelse(number_genomes == 1,  "Singleton", Status))
paste0(round(100 * nrow(gene_counts_all_noLowBuscos %>% filter(Status == "Core")) / n_orthogroups,2), " % are core orthogroups")
paste0(round(100 * nrow(gene_counts_all_noLowBuscos %>% filter(Status == "Accessory")) / n_orthogroups,2), " % are accessory orthogroups")

### 3.3 Transform dataframe and plot

ortho_df_wide <- gene_counts_all[,!names(gene_counts_all) %in% cols_to_exclude] %>%
  group_by(Status) %>%
  dplyr::summarise(across(everything(), sum))

ortho_df_long <- gene_counts_all[,!names(gene_counts_all) %in% cols_to_exclude] %>%
  group_by(Status) %>%
  summarise(across(everything(), sum)) %>%
  pivot_longer(cols = -Status, names_to = "ShortName", values_to = "Value") %>%
  mutate(Status = factor(Status, levels = c("Singleton" ,"Accessory", "Core")))

ortho_df_long$ShortName <- factor(ortho_df_long$ShortName, levels = genome_order)

# Plot horizontal barplot
gene_plot <- ggplot(ortho_df_long, aes(x = ShortName, y = Value, fill = Status)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("Core" = "#0072B2", "Accessory" = "#56B4E9","Singleton"= "#E69F00")) +
  coord_flip() +
  theme_minimal() +
  labs(title = "OrthoFinder", y = "Gene Count", χ = "") + 
  theme(plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
        text = element_text(size = base_size, family = base_family ),
        axis.title = element_text(face = "bold",size = rel(1)),
        axis.title.y = element_blank(), #axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.x = element_text(vjust = -0.2), axis.ticks.x = element_line(),
        panel.grid.major = element_line(colour="gray80"), panel.grid.minor = element_blank(),
        legend.key = element_rect(colour = NA),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.size= unit(0.6, "cm"),
        legend.spacing = unit(0, "cm"),
        legend.title = element_blank(),
        plot.margin=unit(c(10,5,5,5),"mm"),
        strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
        strip.text = element_text(face="bold")
  )
  
gene_plot



#### 4. Combine plots ####
plots <- plot_grid(tree_plot,busco_plot, gene_plot,
                   rel_widths = c(2,1,1), labels = c("a)", "b)", "c)"), ncol = 3)
plots
#ggsave(file="genomePlots.svg", plot=plots, device = "svg", width=16, height=10)
