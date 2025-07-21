library(tidyverse)
library(stringr)
library(ggplot2)
library(gridExtra)
library(MSnbase)
library(ggVennDiagram)
library(patchwork)
library(ggpattern)

base_size=12
base_family="helvetica"

#### 1. Variable modifications ####
# I won't count contaminants

### 1.1 Read mzTab and categorize
varMods <- c("M\\(Oxidation\\)", "W\\(Oxidation\\)", "\\(Acetyl\\)", "N\\(Deamidated\\)", "\\(Gln->pyro-Glu\\)")
extra_varMods <- c("W\\(Oxidation\\)", "N\\(Deamidated\\)", "\\(Gln->pyro-Glu\\)")

mzt <- "data/MS/sdrf_openms_design_openms.mzTab.gz"
psms <- data.frame(psms(MzTab(mzt))) %>% 
  filter(!grepl("CONTAMINANT", accession)) %>%
  mutate(peptidoform = opt_global_cv_MS.1000889_peptidoform_sequence) %>% 
  select(accession, sequence, peptidoform)

psms_no_var_mods <- psms %>% filter(!grepl(paste(varMods, collapse = "|"), peptidoform))
psms_standard_var_mods <- psms %>% filter(grepl(paste(varMods, collapse = "|"), peptidoform)) %>% filter(!grepl(paste(extra_varMods, collapse = "|"), peptidoform))
psms_extra_var_mods <- psms %>% filter(grepl(paste(extra_varMods, collapse = "|"), peptidoform))

### 1.2 Compute counts and percentages

# Compute total counts for normalization
total_counts <- c(nrow(psms), 
                  length(unique(psms$peptidoform)), 
                  length(unique(psms$sequence)))

nr_peptides_only_standard_varMods <- length(unique(setdiff(psms_standard_var_mods$sequence, psms_no_var_mods$sequence)))
nr_peptides_only_extra_varMods <- length(unique(setdiff(psms_extra_var_mods$sequence, union(psms_no_var_mods$sequence, psms_standard_var_mods$sequence))))

# Compute percentages
counts <- data.frame(
  Level = rep(c("PSMs", "Peptidoforms", "Peptides"), each = 3),
  Category = rep(c("no", "standard", "extra"), times = 3),
  Percentage = c(
    
    # PSM Level
    nrow(psms_no_var_mods) / total_counts[1],
    nrow(psms_standard_var_mods) / total_counts[1],
    nrow(psms_extra_var_mods) / total_counts[1],
    
    # Peptidoform Level
    length(unique(psms_no_var_mods$peptidoform)) / total_counts[2],
    length(unique(psms_standard_var_mods$peptidoform)) / total_counts[2],
    length(unique(psms_extra_var_mods$peptidoform)) / total_counts[2] ,
    
    # Peptide Level
    length(unique(psms_no_var_mods$sequence)) / total_counts[3],
    nr_peptides_only_standard_varMods / total_counts[3],
    nr_peptides_only_extra_varMods / total_counts[3] 
  )
)
counts$Level <- factor(counts$Level, levels = c("PSMs", "Peptidoforms","Peptides" ))
counts$Category <- factor(counts$Category, levels = c( "extra", "standard","no"))


### 1.3 Generate the percentage stacked bar plot

category_colors <- c("no" = "#E69F00",
                     "standard" = "#56B4E9", 
                     "extra" = "#0072B2")

mods_plot <- ggplot(counts, aes(x = Level, y = Percentage, fill = Category)) +
  geom_bar(stat = "identity", position = "fill") +  # 'fill' normalizes to percentages
  scale_fill_manual(values = category_colors,  name = "Variable PTMs") +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(scale = 100), position = "right") +  
  theme_minimal(base_size=base_size, base_family=base_family) +
  labs(title = "", x = "", y = "Percentage") +
  theme(
    plot.title = element_text(face = "bold",  size = rel(1.2), hjust = 0.5),
    text = element_text(size = base_size, family = base_family ),
    axis.title = element_text(face = "bold",size = rel(1)),
    axis.title.y = element_text(angle=90,vjust =2),
    axis.title.x = element_text(vjust = -0.2),
    axis.text = element_text(), 
    panel.grid.major = element_line(colour="gray80"),
    panel.grid.minor = element_blank(),
    legend.key = element_rect(colour = NA),
    legend.position = "bottom", legend.direction = "horizontal",
    legend.key.size= unit(0.6, "cm"),
    legend.spacing = unit(0, "cm"),
    #legend.title = element_blank(),
    plot.margin=unit(c(10,5,5,5),"mm"),
    strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
    strip.text = element_text(face="bold")
  )
mods_plot

#### 2. Venn diagram ####
df <- read.csv("results/MS/MSresults.csv.gz") %>% filter(MS_status == "quant")

venn_data <- list(
  xyl = df %>% filter(!is.na(avgAbd_xyl)) %>% pull(protein_id),
  cs = df %>% filter(!is.na(avgAbd_cs)) %>% pull(protein_id),
  bw = df %>% filter(!is.na(avgAbd_bw)) %>% pull(protein_id)
)

venn_plot <- ggVennDiagram(venn_data, set_size = 6, label_alpha = 0, label_size = 4, label_percent_digit = 1) +  
  scale_fill_gradient(low = "lightsteelblue1", high = "#4981BF") +
  theme(legend.title  = element_blank(), legend.text = element_text(size = 12))
venn_plot


#### 3. DEPs barplot ####
### 3.1 Read data

df <- read.csv("results/MS/prolfqua_unfiltered_moderated_v1_mc.csv.gz")
df$contrast <- gsub("-", "vs", df$contrast) 
contrast_order <- c("bw vs xyl", "cs vs xyl", "bw vs cs")
# mark CAZymes
annot <- read.csv("results/MS/MSresults.csv.gz")
cazymes <- annot %>% filter(!is.na(CAZy) & CAZy != "" & CAZy != "EXPN") %>% pull(protein_id)
# up and down-regulated
df <- df %>%
  mutate(
    CAZy = ifelse(protein_Id %in% cazymes, "CAZyme", "Other"),
    DE = case_when(
      diff > 0 & FDR < 0.05 ~ "Up",
      diff < 0 & FDR < 0.05 ~ "Down")) %>%
    filter(DE %in% c("Up", "Down"))

### 3.2 Count DEPs and CAZymes
count_DEPs <- df %>% count(contrast, DE, name = "Count")
count_cazy <- df %>% filter(CAZy == "CAZyme") %>% count(contrast, DE, name = "CAZymes")

count_DEPs <- left_join(count_DEPs, count_cazy, by = c("contrast", "DE")) %>%
  mutate(CAZymes = replace_na(CAZymes, 0),
         'non CAZymes' = Count - CAZymes) 

count_DEPs_long <- count_DEPs %>%
  mutate(contrast = factor(contrast, levels = contrast_order)) %>% 
  pivot_longer(cols = c(CAZymes,  'non CAZymes'), names_to = "Type", values_to = "Value") %>%
  mutate(Type = factor(Type, levels = c("non CAZymes", "CAZymes")))


### 3.3 Plot
deps_plot <- ggplot(count_DEPs_long, aes(x = interaction(contrast, DE), y = Value, fill = DE, pattern = Type)) +
  geom_bar_pattern(stat = "identity", position = "stack", color = "black", pattern_fill = "black", pattern_density = 0.2, pattern_spacing = 0.02) +
  scale_fill_manual(values = c("Up" = "#E69F00", "Down" = "#56B4E9")) +
  scale_pattern_manual(values = c("CAZymes" = "stripe", "non CAZymes" = "none")) +
  scale_x_discrete(labels = function(x) sub("\\..*", "", x)) +   # Clean up x-axis labels to show only the contrast name
  theme_minimal() +
  labs(title = "" , x = "Contrast", y = "Protein Count") +
  theme(plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
        text = element_text(size = base_size, family = base_family ),
        axis.title = element_text(face = "bold",size = rel(1)),
        axis.title.y = element_text(angle=90,vjust =2),
        axis.title.x = element_text(vjust = -0.2),
        axis.text = element_text(), 
        panel.grid.major.y = element_line(colour="gray80"), panel.grid.major.x = element_blank(),panel.grid.minor = element_blank(),
        legend.key = element_rect(colour = NA),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.size= unit(0.6, "cm"),
        legend.spacing = unit(0, "cm"),
        legend.title = element_blank(),
        plot.margin=unit(c(10,5,5,5),"mm"),
        strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
        strip.text = element_text(face="bold"))


#### 4. Combine plots ####

library(cowplot)
down_plots <- plot_grid(venn_plot, deps_plot, labels = c("b)", "c)"), ncol = 2, nrow = 1)
plot_grid(mods_plot, down_plots, 
          labels = c("a)", ""),
          ncol = 1, nrow = 2)
  
#ggsave(file="../results/figures/proteome/proteome_plots.jpg", width=7, height=8)
#ggsave(file="../results/figures/proteome/proteome_plots.svg", width=7, height=8)
