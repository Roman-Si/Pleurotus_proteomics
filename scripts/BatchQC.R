library(BatchQC)
library(imputeLCMD)
library(gridExtra)
library(tidyverse)
library(limma)
library(PCAtools)


#### 1. Functions ####
impute_se <- function(data, impute_method) {
  
  if (impute_method == "MinProb") {
    data <- impute.MinProb(data)
  }
  if (impute_method == "MinDet") {
    data <- impute.MinDet(data)
  }
  
  se_object <- BatchQC::summarized_experiment(data, batches)
  SummarizedExperiment::assayNames(se_object) <- "log_intensity"
  
  se_object <- BatchQC::batch_correct(se = se_object, method = "ComBat", assay_to_normalize = "log_intensity", batch = "comb_batch",
                                      covar = c("condition"), output_assay_name = "Combat_comb")
  se_object <- BatchQC::batch_correct(se = se_object, method = "ComBat", assay_to_normalize = "log_intensity", batch = "ms_batch",
                                      covar = c("condition"), output_assay_name = "Combat_ms")
  se_object <- BatchQC::batch_correct(se = se_object, method = "ComBat", assay_to_normalize = "log_intensity", batch = "exp_batch",
                                      covar = c("condition"), output_assay_name = "Combat_exp")
  
  return(se_object)
}

plot_variation <- function(se_object) {
  p1 <- BatchQC::EV_plotter(batchqc_ev = batchqc_explained_variation(se = se_object,batch = "comb_batch", 
                                                                     condition = "condition", assay_name = "log_intensity")[[1]]) +  labs(title = "Combined batch uncorrected")
  p3 <- BatchQC::EV_plotter(batchqc_ev = batchqc_explained_variation(se = se_object,batch = "ms_batch", 
                                                                     condition = "condition", assay_name = "log_intensity")[[1]])  + labs(title = "MS batch uncorrected")
  p5 <- BatchQC::EV_plotter(batchqc_ev = batchqc_explained_variation(se = se_object,batch = "exp_batch", 
                                                                     condition = "condition", assay_name = "log_intensity")[[1]]) +   labs(title = "Experiment batch uncorrected")
  p2 <- BatchQC::EV_plotter(batchqc_ev = batchqc_explained_variation(se = se_object,batch = "comb_batch", 
                                                                     condition = "condition", assay_name = "Combat_comb")[[1]]) +  labs(title = "Combined batch corrected")
  p4 <- BatchQC::EV_plotter(batchqc_ev = batchqc_explained_variation(se = se_object,batch = "ms_batch", 
                                                                     condition = "condition", assay_name = "Combat_ms")[[1]])  + labs(title = "MS batch corrected")
  p6 <- BatchQC::EV_plotter(batchqc_ev = batchqc_explained_variation(se = se_object,batch = "exp_batch", 
                                                                     condition = "condition", assay_name = "Combat_exp")[[1]]) +   labs(title = "Experiment batch corrected")
  
  gridExtra::grid.arrange(p1, p2, p3, p4, p5, p6)
}

#### 2. Read and prepare data ####
# First filter for peptides present in at least 2 samples
batch_figure_dir = "secretome/quantms/figures/prolfqua/batches/"
pep_data <- read.csv("secretome/quantms/prolfqua/uncorected_robscaled_peptide_intensities.csv.gz")
pep_df <- pep_data[,c(4:12)]
rownames(pep_df) <- pep_data$peptide_Id
# filter for 2 NaN
pep_df <- pep_df[rowSums(is.na(pep_df)) <= 7, ]

# batches
batches <- data.frame(row.names = colnames(pep_df))
batches$condition  <- as.factor(rep(c("bw", "cs", "xyl"),each=3)) 
batches$exp_batch <- as.factor(c(2,2,2,1, 1,1,1,1,2))
batches$ms_batch <- as.factor(c(2,2,2,1, 2,2,1,1,2))
batches$comb_batch <- as.factor(c("e2_m2","e2_m2","e2_m2","e1_m1", "e1_m2","e1_m2","e1_m1","e1_m1","e2_m2"))

se_object <- BatchQC::summarized_experiment(pep_df, batches)
SummarizedExperiment::assayNames(se_object) <- "log_intensity"
BatchQC::confound_metrics(se = se_object, batch = "condition")


p <- PCAtools::pca(na.omit(pep_df), metadata = batches, removeVar = 0.1)
PCAtools::biplot(p, colby = "condition", shape = 'comb_batch', title = "PCA of robscaled peptide intensities",
                 labSize = 5, pointSize = 5)


#### 3. Imputations ####
# There are a lot of missing values of course but ComBat can not handle them so will impute with 2 different methods. MinDet which replaces missing values with the minimum value of the sample and MinProb which generates a (normal?) distribution around the minimum value. First threw peptides present in less than 2 samples. This will affect samples with many missing values and completely change their distribution, let's see number of NaN values in each sample

cat("There are", nrow(pep_df), "peptides in total.\nLet's count NaN in each sample")
colSums(is.na(pep_df))

cols =   c("darkorange1", "coral2", "coral4", "yellow", "greenyellow", "green4", "purple1", "slateblue", "slategrey")
limma::plotDensities(pep_df, col =cols)
limma::plotDensities(impute.MinDet(pep_df), col = cols, legend = FALSE)
limma::plotDensities(impute.MinProb(pep_df), col = cols, legend = FALSE)

se_object_md <- impute_se(pep_df, "MinDet")
se_object_mp <- impute_se(pep_df, "MinProb")

#### 4. Plots of BatchQC ####
plot_variation(se_object_md)
plot_variation(se_object_mp)

pca_plot_md <- BatchQC::PCA_plotter(se = se_object_md, nfeature = 20, color = "condition",
                                    shape = "comb_batch", assays = c("log_intensity", "Combat_comb"),
                                    xaxisPC = 1, yaxisPC = 2, log_option = FALSE)
pca_plot_mp <- BatchQC::PCA_plotter(se = se_object_mp, nfeature = 20, color = "condition",
                                    shape = "comb_batch", assays = c("log_intensity", "Combat_comb"),
                                    xaxisPC = 1, yaxisPC = 2, log_option = FALSE)

pca_plot_md$plot
pca_plot_mp$plot

heatmaps <- BatchQC::heatmap_plotter(se = se_object_md, assay = "log_intensity",
                                     nfeature = 38, annotation_column = c("comb_batch", "condition"),
                                     log_option = "FALSE")
heatmap_cor_md <- heatmaps$correlation_heatmap
dev.off()
heatmap_cor_md

heatmaps <- BatchQC::heatmap_plotter(se = se_object_mp, assay = "log_intensity",
                                     nfeature = 38, annotation_column = c("comb_batch", "condition"),
                                     log_option = "FALSE")
heatmap_cor_mp <- heatmaps$correlation_heatmap
dev.off()
heatmap_cor_mp


#### 5. Check TMP protein intensities ####

prot_df <- read_csv("secretome/quantms/prolfqua/protein_intensities.csv.gz")

p <- PCAtools::pca(na.omit(prot_df[,c(2:10)]), metadata = batches, removeVar = 0.1)
PCAtools::biplot(p, colby = "condition", shape = 'comb_batch', title = "PCA of protein intensities",
                             labSize = 5, pointSize = 5)

horn <- parallelPCA(na.omit(prot_df[,c(2:10)]))
elbow <- findElbowPoint(p$variance)
PCAtools::screeplot(p,
                                  components = getComponents(p, 1:20),
                                  vline = c(horn$n, elbow)) +
  
  geom_label(aes(x = horn$n + 0.5, y = 50,
                 label = 'Horn\'s', vjust = -1, size = 8)) +
  geom_label(aes(x = elbow + 1, y = 30,
                 label = 'Elbow method', vjust = -1, size = 8))




###### OLD STUFF with shiny app ################

nbatch <- 3
ncond <- 3
npercond <- 3
ms_batch=c(2,2,2,1, 2,2,1,1,2)
exp_batch = c(2,2,2,1, 1,1,1,1,2)
comb_batch = c(3,3,3,1, 2,2,1,1,3)

condition <- rep(c(3, 2, 1),each=3)
batchQC(pep_df , batch=comb_batch, condition=condition, 
        report_file="quantms/figures/prolfqua/batch_correction/peptide_level/BatchQC/batchqc_CombinedBatch/batchqc_CombinedBatch.html", 
        report_dir="/home/roman/Documents/Bioinfo/ntua_fungi/Zerba/Proteomics/Citrinopileatus/quantms/figures/prolfqua/batch_correction/peptide_level/BatchQC/batchqc_CombinedBatch/", 
        report_option_binary="111111111",
        view_report=FALSE, interactive=TRUE, batchqc_output=TRUE)
batchQC(pep_df , batch=ms_batch, condition=condition, 
        report_file="quantms/figures/prolfqua/batch_correction/peptide_level/BatchQC/batchqc_MSbatch/batchqc_MSbatch.html", 
        report_dir="/home/roman/Documents/Bioinfo/ntua_fungi/Zerba/Proteomics/Citrinopileatus/quantms/figures/prolfqua/batch_correction/peptide_level/BatchQC/batchqc_MSbatch/", 
        report_option_binary="111111111",
        view_report=FALSE, interactive=TRUE, batchqc_output=TRUE)
batchQC(pep_df , batch=exp_batch, condition=condition, 
        report_file="quantms/figures/prolfqua/batch_correction/peptide_level/BatchQC/batchqc_EXPbatch/batchqc_EXPbatch.html", 
        report_dir="/home/roman/Documents/Bioinfo/ntua_fungi/Zerba/Proteomics/Citrinopileatus/quantms/figures/prolfqua/batch_correction/peptide_level/BatchQC/batchqc_EXPbatch/", 
        report_option_binary="111111111",
        view_report=FALSE, interactive=TRUE, batchqc_output=TRUE)