source("scripts/generic_functions.R")
source("scripts/parse_mzTab.R")
source("scripts/prolfqua_functions.R")
library(PCAtools)
library(DEP)
library(ggVennDiagram)

########################################### 1. Read and prepare data #######################################################
quantms_figdir = "secretome/quantms/figures/"
mzt <-"secretome/quantms/proteomicslfq/sdrf_openms_design_openms.mzTab.gz"

### Take mzTab info
prot_info <- extract_protein_stats(mzt, output_file = paste(quantms_figdir, "../proteomicslfq/unfiltered_MS_data.csv.gz", sep = ''))
proteomicslfq_ids_to_keep <- get_proteinIds_from_proteomicslfq(mzt)
plot_ptms <- count_peptide_modifications(mzt,  plot_type = "separate")
plot_ptms
ggsave(paste(quantms_figdir, "PTMs_separate.jpg", sep = ''), width=10, height=8)

### prepeare for prolfqua
data <- read.csv(paste(quantms_figdir, "../proteomicslfq/sdrf_openms_design_msstats_in.csv.gz", sep = ""), header = TRUE, sep = ',')
data$Condition <- gsub("corn stover", "cornstover", data$Condition)
replacements <- c(
  "E26698_2p_50uPAC12_trap10_PRC-5442_2.mzML" = 	"xylose_1",
  "E26704_2p_50uPAC12_trap10_PRC-5442_5.mzML" = 	"xylose_2",
  "B28551_Ap_IonOpt_PRC-6063__newprep_3_10ul.mzML" = 	"xylose_3",
  "E28115_1p_50uPAC13__trap9_PRC-5590_1.mzML" = 	"cornstover_1",
  "B28547_Ap_IonOpt_PRC-6063_newprep_1_10ul.mzML" = 	"cornstover_2",
  "B28549_Ap_IonOpt_PRC-6063__newprep_2_10ul.mzML" = 	"cornstover_3",
  "B28553_Ap_IonOpt_PRC-6063__newprep_4_10ul.mzML" = 	"beachwood_1",
  "B28543_Ap_IonOpt_PRC-6063__newprep_5_10ul.mzML" = 	"beachwood_2",
  "B28545_Ap_IonOpt_PRC-6063__newprep_6_10ul.mzML" = 	"beachwood_3"
)
data$Reference <- replacements[data$Reference]

########################################### 2. prolfqua peptide level #######################################################
lfqdata <- create_lfqdata(data)
lfqdata$hierarchy_counts()

density_raw_pep <- lfqdata$get_Plotter()$intensity_distribution_density()  + labs(tag = "A) Raw")

p1 <- lfqdata$get_Summariser()$plot_missingness_per_group() + labs(tag = "A)")
p2 <-  lfqdata$get_Plotter()$missigness_histogram() + labs(tag = "B)")
png(paste(quantms_figdir, "MVs_per_Culture.png", sep = ''), width = 700, height = 900)
gridExtra::grid.arrange(p1, p2)
dev.off()

stats <- lfqdata$get_Stats()
prolfqua::table_facade( stats$stats_quantiles()$wide, paste0("quantile of ",stats$stat ))
p1 <- stats$density_median() + labs(tag = "A)")
p2 <-  stats$violin() + labs(tag = "B)")
stdm_raw <- stats$stdv_vs_mean(size = 10000) +
  ggplot2::scale_x_log10() +
  ggplot2::scale_y_log10()
png(paste(quantms_figdir, "CV.png", sep = ''), width = 700, height = 900)
gridExtra::grid.arrange(p1, p2)
dev.off()

### Normalize peptide intensities
lt <- lfqdata$get_Transformer()
lfqdataPeptideNorm <- lt$log2()$robscale()$lfq
density_norm_pep <- lfqdataPeptideNorm$get_Plotter()$intensity_distribution_density() + labs(tag = "B) robscale")
png(paste(quantms_figdir, "prolfqua/peptide_density.png", sep = ''), width = 900, height = 700)
gridExtra::grid.arrange(density_raw_pep, density_norm_pep)
dev.off()
pl <-  lfqdataPeptideNorm$get_Plotter()
png(paste(quantms_figdir, "prolfqua/PCA_peptide_robscaled.png", sep = ''), width = 900, height = 700)
pl$pca() +  labs(tag = "PCA of robscaled peptide intensities")
dev.off()
pep_data <- data.frame(lfqdataPeptideNorm$to_wide()$data)
save_gzipped_csv(data.frame(pep_data),paste(quantms_figdir, "../prolfqua/robscaled_peptide_intensities.csv.gz", sep = ''))

########################################### 3. Protein aggregation and EDA #######################################################

### 3.1 Protein aggregation
# - Use Tukey's median polish protein aggregation
# - Remove artifacts with same value in all samples 
# - Replace the artifact intensities with top3 method
# - Merge the two methods and import in prolfqua for DE analysis

mp_df <- protein_aggregation(lfqdataPeptideNorm)
mp_df <- filter_for_leading_protein(mp_df, proteomicslfq_ids_to_keep)
save_gzipped_csv(data.frame(mp_df), "secretome/quantms/prolfqua/uncorected_protein_intensities.csv.gz")
rownames(mp_df) <- mp_df$protein_Id
mp_df <- mp_df[,c(2:(length(data$Reference %>% unique()) + 1))]
metadata <- data.frame(row.names = colnames(mp_df))
metadata$Condition <- gsub("^(.*)_.*$",  "\\1", colnames(mp_df))
metadata$Batch <- c("e2_m2","e2_m2","e2_m2","e1_m1", "e1_m2","e1_m2","e1_m1","e1_m1","e2_m2")
p <- PCAtools::pca(na.omit(mp_df), metadata = metadata, removeVar = 0.1)
PCAtools::biplot(p, colby = "Condition", shape = "Batch", title = "PCA of prolfqua protein intensities",  legendPosition = "left",
                 subtitle = "Without filtering for proteins in >=2 replicates", labSize = 5, pointSize = 5)

### 3.2 EDA plots after filtering for reliably quantified proteins
mp_df_quant <- filter_proteins_by_replicates(mp_df, c("beachwood", "cornstover", "xylose"), 2)
quantified_proteins <- rownames(mp_df_quant)
p <- PCAtools::pca(na.omit(mp_df_quant), metadata = metadata, removeVar = 0.1)
PCAtools::biplot(p, colby = "Condition", shape = "Batch", title = "PCA of prolfqua robscaled - TMP or top3 summarized protein intensities", 
                 subtitle = "Filtered for proteins in >=2 replicates\nSame value proteins are summarized with top3 instead of TMP",
                 labSize = 5,legendPosition = "left", pointSize = 5)#, xlim = c(-40, 40), ylim = c(-30, 30))
ggsave( paste(quantms_figdir, "prolfqua/PCA_protein.png", sep = ''), width=10, height=8)

horn <- parallelPCA(na.omit(mp_df_quant))
elbow <- findElbowPoint(p$variance)
PCAtools::screeplot(p,
                    components = getComponents(p, 1:20),
                    vline = c(horn$n, elbow)) +
  
  geom_label(aes(x = horn$n + 0.5, y = 50,
                 label = 'Horn\'s', vjust = -1, size = 8)) +
  geom_label(aes(x = elbow + 1, y = 30,
                 label = 'Elbow method', vjust = -1, size = 8))
ggsave( paste(quantms_figdir, "prolfqua/screeplot.png", sep = ''), width=10, height=8)

png(paste(quantms_figdir, "prolfqua/pearson.png", sep = ''), width = 900, height = 700)
GGally::ggcorr(mp_df_quant, method = c("pairwise", "pearson"), label = TRUE, label_round = 2, label_size = 3)
dev.off()

transposed_data <- t(mp_df_quant)
dist_matrix <- dist(transposed_data, method = "euclidean")
hclust_result <- hclust(dist_matrix, method = "complete")
png(paste(quantms_figdir, "prolfqua/dendrogram.png", sep = ''), width = 600, height = 600)
plot(hclust_result, hang = -1)
dev.off()

### 3.3 DEP
### 3.3.1 Sample_level
protein_intensities <- tibble::rownames_to_column(mp_df_quant, "proteinId")
condition_cols <- c(2:10)
data_se <- DEP:::make_unique(protein_intensities,"proteinId","proteinId",delim=";")
exp_design <- data.frame(
  label =rownames(metadata),
  condition = metadata$Condition,
  replicate = rep(c(1:3),3)
)
data_se<-DEP:::make_se(data_se, condition_cols ,exp_design)

plot_frequency(data_se) + labs(title = "Protein identification overlap") + xlab("Number of samples")
ggsave(file=paste(quantms_figdir, "DEP_sample_overlap.png", sep = ''), width=10, height=8)

plot_numbers(data_se) + labs(title = "Proteins per sample") + ylab("Number of proteins") + scale_x_discrete(labels = exp_design$label) +
  theme(legend.position = "none")
ggsave(file=paste(quantms_figdir, "DEP_proteins_per_sample.png", sep = ''), width=10, height=8)

plot_missval(data_se)
plot_detect(data_se)


### 3.3.2 Substrate level
protein_intensities <- tibble::rownames_to_column(mp_df_quant, "proteinId")

# Prepare substrate dataframe ####
print("A protein is detected if present one replicate")
protein_intensities$beachwood_s1 <- rowMeans(protein_intensities[, 2:4], na.rm = TRUE)
protein_intensities$cornstover_s1 <- rowMeans(protein_intensities[, 5:7], na.rm = TRUE)
protein_intensities$xylose_s1 <- rowMeans(protein_intensities[, 8:10], na.rm = TRUE)

condition_cols <- c(11:13)
conditions <- c("beachwood", "cornstover", "xylose")
data_se <- DEP:::make_unique(protein_intensities,"proteinId","proteinId",delim=";")
exp_design <- data.frame(
  label =paste(unique(conditions), "_s1", sep = ""),
  condition = unique(conditions),
  replicate = rep(c(1),3)
)
data_se<-DEP:::make_se(data_se, condition_cols ,exp_design)

plot_frequency(data_se) + labs(title = "Protein identification overlap") + xlab("Number of substrates")
ggsave(file=paste(quantms_figdir, "DEP_substrate_overlap.png", sep = ''), width=10, height=8)

plot_numbers(data_se) + labs(title = "Proteins per substrate") + ylab("Number of proteins") + scale_x_discrete(labels = conditions) +
  theme(legend.position = "none")
ggsave(file=paste(quantms_figdir, "DEP_proteins_per_substrate.png", sep = ''), width=10, height=8)


### 3.4 ggVenn diagram 

# Define quantified_proteins list
quantified_proteins <- list(
  beachwood = NULL,
  cornstover = NULL,
  xylose = NULL
)

# Create a list of conditions and corresponding column names
conditions <- list(
  beachwood = "beachwood_s1",
  cornstover = "cornstover_s1",
  xylose = "xylose_s1"
)

# Loop through each element in quantified_proteins
for (protein in names(quantified_proteins)) {
  # Get the condition column name for the protein
  condition_column <- conditions[[protein]]
  
  # Check if the condition column is not all NA
  if (!all(is.na(protein_intensities[[condition_column]]))) {
    # Populate quantified_proteins with proteinId values
    quantified_proteins[[protein]] <- protein_intensities$proteinId[!is.na(protein_intensities[[condition_column]])]
  }
}

ggVennDiagram(quantified_proteins) +
  theme(legend.position = 'none',
        legend.text = element_text(size = 24)) 
ggsave(file=paste(quantms_figdir, "ggvenn_substrate_overlap.jpg", sep = ''), width=5, height=4)


########################################### 5. prolfqua DE #######################################################
# Will run with and without filtering/covariate

### 5.1 Unfiltered proteins

long_df <- data.frame(mp_df) %>%
  rownames_to_column("protein_Id") %>%  # Convert rownames to a new 'proteinId' column
  pivot_longer(cols = -protein_Id, names_to = "Reference", 
               values_to = "Intensity", values_drop_na = TRUE)
head(long_df)

annot <- data.frame(colnames(mp_df)) %>%
  mutate(
    Reference = colnames.mp_df., 
    Run = c(1:9),
    Condition = gsub("^(.*)_.*$", "\\1", colnames.mp_df.),
    replicate = gsub("^.*_r([1-3])$", "\\1", colnames.mp_df.),
    Batch = c("e2_m2","e2_m2","e2_m2","e1_m1", "e1_m2","e1_m2","e1_m1","e1_m1","e2_m2")
  ) %>% select(Reference, Run, Condition, replicate, Batch)

startdata <- dplyr::inner_join(long_df, annot, by = "Reference")
lfqdata <- create_lfqdata(startdata, response_level = "protein", proteinId_column = "protein_Id", extra_factor = "Batch" )

# EDA plots
p1 <- lfqdata$get_Summariser()$plot_missingness_per_group() + labs(tag = "A)") 
p2 <- pl$missigness_histogram() + labs(tag = "B)")
png(paste(quantms_figdir, "prolfqua/MVs_protein_level.png", sep = ''), width = 900, height = 700)
gridExtra::grid.arrange(p1, p2)
dev.off()

### 5.1.1 Without covariate

formula_Condition <-  strategy_lm("Intensity ~ condition_")

Contrasts <- c("cs - xyl" = "condition_cornstover - condition_xylose",
               "bw - xyl" = "condition_beachwood - condition_xylose",
               "bw - cs" = "condition_beachwood - condition_cornstover")

mod <- prolfqua::build_model(
  lfqdata$data,
  formula_Condition,
  subject_Id = lfqdata$config$table$hierarchy_keys())
contr <- prolfqua::Contrasts$new(mod, Contrasts)
v1 <- contr$get_Plotter()$volcano()
mC <- ContrastsMissing$new(lfqdata = lfqdata, contrasts = Contrasts)
merged  <- prolfqua::merge_contrasts_results(prefer = contr,add = mC)$merged
merged <- prolfqua::ContrastsModerated$new(merged)
unfiltered_nobatch <- merged$get_Plotter()$volcano()
unfiltered_nobatch$FDR


### 5.1.2 With covariate
formula_Condition <-  strategy_lm("Intensity ~ condition_ + extrafactor_", model_name = "lm")


mod <- prolfqua::build_model(
  lfqdata$data,
  formula_Condition,
  subject_Id = lfqdata$config$table$hierarchy_keys())
contr <- prolfqua::Contrasts$new(mod, Contrasts)
v1 <- contr$get_Plotter()$volcano()
mC <- ContrastsMissing$new(lfqdata = lfqdata, contrasts = Contrasts)
merged  <- prolfqua::merge_contrasts_results(prefer = contr,add = mC)$merged
merged <- prolfqua::ContrastsModerated$new(merged)
unfiltered_batch <- merged$get_Plotter()$volcano()


### 5.2 Filtered for quantified proteins

long_df <- data.frame(mp_df_quant) %>%
  rownames_to_column("protein_Id") %>%  # Convert rownames to a new 'proteinId' column
  pivot_longer(cols = -protein_Id, names_to = "Reference", 
               values_to = "Intensity", values_drop_na = TRUE)

startdata <- dplyr::inner_join(long_df, annot, by = "Reference")
lfqdata <- create_lfqdata(startdata, response_level = "protein", proteinId_column = "protein_Id", extra_factor = "Batch" )

# EDA plots
p1 <- lfqdata$get_Summariser()$plot_missingness_per_group() + labs(tag = "A)") 
p2 <- pl$missigness_histogram() + labs(tag = "B)")
png(paste(quantms_figdir, "prolfqua/MVs_quant-protein_level.png", sep = ''), width = 900, height = 700)
gridExtra::grid.arrange(p1, p2)
dev.off()

### 5.2.1 Without covariate

formula_Condition <-  strategy_lm("Intensity ~ condition_")

Contrasts <- c("cs - xyl" = "condition_cornstover - condition_xylose",
               "bw - xyl" = "condition_beachwood - condition_xylose",
               "bw - cs" = "condition_beachwood - condition_cornstover")

mod <- prolfqua::build_model(
  lfqdata$data,
  formula_Condition,
  subject_Id = lfqdata$config$table$hierarchy_keys())
contr <- prolfqua::Contrasts$new(mod, Contrasts)
v1 <- contr$get_Plotter()$volcano()
mC <- ContrastsMissing$new(lfqdata = lfqdata, contrasts = Contrasts)
merged  <- prolfqua::merge_contrasts_results(prefer = contr,add = mC)$merged
merged <- prolfqua::ContrastsModerated$new(merged)
filtered_nobatch <- merged$get_Plotter()$volcano()

### 5.2.2 With covariate
formula_Condition <-  strategy_lm("Intensity ~ condition_ + extrafactor_", model_name = "lm")

mod <- prolfqua::build_model(
  lfqdata$data,
  formula_Condition,
  subject_Id = lfqdata$config$table$hierarchy_keys())
contr <- prolfqua::Contrasts$new(mod, Contrasts)
v1 <- contr$get_Plotter()$volcano()
mC <- ContrastsMissing$new(lfqdata = lfqdata, contrasts = Contrasts)
merged  <- prolfqua::merge_contrasts_results(prefer = contr,add = mC)$merged
merged <- prolfqua::ContrastsModerated$new(merged)
filtered_batch <- merged$get_Plotter()$volcano()

### 5.3 Gather results and compare

comparisons <- c("cs - xyl", "bw - xyl", "bw - cs")

dataframes <- list(
  unfiltered_nobatch$FDR$data,
  unfiltered_batch$FDR$data,
  filtered_nobatch$FDR$data,
  filtered_batch$FDR$data
)

dataframe_names <- c(
  "unfiltered_nobatch",
  "unfiltered_batch",
  "filtered_nobatch",
  "filtered_batch"
)

results <- data.frame(
  DataFrame = character(),
  Comparison = character(),
  DEPs = integer()
)

count_diff_proteins <- function(df, comparison) {
  df %>%
    filter(contrast == comparison & abs(diff) > 1 & FDR < 0.05) %>%
    nrow()
}

for (i in seq_along(dataframes)) {
  for (comp in comparisons) {
    diff_proteins <- count_diff_proteins(dataframes[[i]], comp)
    results <- rbind(results, data.frame(
      DataFrame = dataframe_names[i],
      Comparison = comp,
      DEPs = diff_proteins
    ))
  }
}

print(results)

# Generate a plot from the results
ggplot(results, aes(x = Comparison, y = DEPs, fill = DataFrame)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() + scale_fill_viridis_d() +
  labs(title = "Count of DEPs per design",
       x = "Comparison",
       y = "Number of DEPs",
       fill = "Design")
ggsave(paste(quantms_figdir, "prolfqua/DEPs_per_design.jpg", sep = ''), width=10, height=8)


de_test <- filtered_batch$FDR$data
write.csv(data.frame(de_test),gzfile(paste(quantms_figdir, "../prolfqua/prolfqua_filtered_batch_moderated_v1_mc.csv.gz", sep = '')), row.names = FALSE)






