library(SKAT)

# 3 per minute
## Running SKAT, no annotations
args <- commandArgs(trailingOnly = TRUE)
CHRO_NB <- as.numeric(1) # chromosome number
TYPE_ <- as.character("OLIG2") # cell type: OLIG2, LHX2, NEUN, peripheralPU1nuclei, coding
ITERATION <- as.character("1") # iteration that this run is or folder that this should be saved as
OUTPUT <- as.character("hi")
RACE <- as.character("ALL")
path = "/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/Alzheimer-RV/data/SKAT/"

# Load phenotype and covariate information
PHENOTYPE_PATH = "/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/Alzheimer-RV/data/phenotypes/36K_QC_filtered_final.csv"
XY <- read.table(file = PHENOTYPE_PATH, sep = ',', header = TRUE)
if (RACE != "ALL") {
  keep = which(XY$predicted_ancestry == RACE)
  XY = XY[XY$predicted_ancestry == RACE,]
}

Y <- XY[c("Diagnosis")]
y = c(Y)$Diagnosis
X <- XY[c('Sex', 'Age', 'apoe_e4', 'apoe_e2', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10', 'Illumina_HiSeq_2000', 'Illumina_HiSeqX', 'Illumina_NovaSeq', 'Illumina', 'USUHS', 'USUHS.Miami', 'NYGC', 'MEDGENOME', 'Baylor', 'Broad', 'WashU', 'PRS_5')] # Choose covariates
# Without PCs
# X <- XY[c('Sex','Age','apoe_e4','apoe_e2', 'Illumina_HiSeq_2000', 'Illumina_HiSeqX','Illumina_NovaSeq', 'Illumina', 'USUHS', 'USUHS.Miami','NYGC', 'MEDGENOME', 'Baylor', 'Broad', 'WashU', 'PRS_5')] # Choose covariates

X$age_age = X$Age * X$Age
X$age_sex = X$Age * X$Sex
X$age_sex2 = X$Age * X$Sex * X$Sex
X <- scale(X) # normalize our covariates, X
X <- X[, colSums(is.na(X)) == 0] # no NANs
X <- cbind(X, 1) # Intercept term

# Run the null model (y ~ X) from covariates
obj.b <- SKAT_Null_Model(y ~ X, out_type = "D")

# Get gene matrices and prep output
Gs <- list.files(path = paste0('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/Alzheimer-RV/data/gene_matrices_maf/', TYPE_, '/chr', CHRO_NB), pattern = "*geno_imputed.csv", full.names = TRUE, recursive = FALSE)

# Process the first gene only
if (length(Gs) > 0) {
  G <- read.csv(Gs[1])
  if (RACE != "ALL") {
    G <- G[keep, ]
  }
  
  # Calculate column sums
  column_sums <- colSums(G)
  # Filter columns where the sum is greater than or equal to 2
  G_filtered <- G[, column_sums >= 10]
  # Check the dimensions of filtered_data
  print(dim(G_filtered))
  
  gene_name = strsplit(basename(Gs[1]), "_")[[1]][1]
  
  # Ensure the filtered data has correct dimensions
  if (ncol(G_filtered) > 0 && nrow(G_filtered) == length(y)) {
    burden <- SKAT(as.matrix(G_filtered), obj.b, method = "Burden")$p.value
    skato <- SKAT(as.matrix(G_filtered), obj.b, method = "SKATO")$p.value
    skat <- SKAT(as.matrix(G_filtered), obj.b, method = "SKAT")$p.value
    
    # Print the p-values
    print(paste("Gene:", gene_name))
    print(paste("Burden p-value:", burden))
    print(paste("SKATO p-value:", skato))
    print(paste("SKAT p-value:", skat))
  } else {
    print("Filtered data dimensions are inconsistent")
  }
} else {
  print("No gene matrices found.")
}
