packages <- c("Rcpp", "RcppArmadillo", "GMMAT", "Matrix", "devtools")
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}
invisible(sapply(packages, install_if_missing))

library(Rcpp)
library(RcppArmadillo)
library(GMMAT)
library(Matrix)
library(GENESIS)
library(devtools)

Sys.setenv(PATH = paste("/gpfs/commons/home/adas/R/cmake-3.29.0-linux-x86_64/bin", Sys.getenv("PATH"), sep=":"))
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GENESIS")library(GENESIS)
devtools::install_github("xihaoli/STAAR")

library(STAAR)


args <- commandArgs(trailingOnly = TRUE)
CHRO_NB <-as.numeric(args[1]) # chromosome number
TYPE_ <- as.character(args[2]) # cell type: OLIG2, LHX2, NEUN, peripheralPU1nuclei, coding
ITERATION <- as.character(args[3]) # iteration that this run is or folder that this should be saved as
OUTPUT <- as.character(args[4])

print(paste("---- STARTING CHROMOSOME ",CHRO_NB, "FOR ", TYPE_, "VARIANTS", ITERATION, "NUMBER ITERATION", OUTPUT, "OUTPUT"))

PHENOTYPE_PATH = "/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/Alzheimer-RV/data/phenotypes/36K_QC_filtered_final.csv"
XY <- read.table(file = PHENOTYPE_PATH, sep = ',', header = TRUE)
XY <- XY[c('Diagnosis','Sex','Age','apoe_e4','apoe_e2','PC1', 'PC2','PC3','PC4','PC5','PC6', 'PC7','PC8','PC9','PC10', 'Illumina_HiSeq_2000', 'Illumina_HiSeqX','Illumina_NovaSeq', 'Illumina', 'USUHS', 'USUHS.Miami','NYGC', 'MEDGENOME', 'Baylor', 'Broad', 'WashU', 'PRS_5')] # Choose covariates
XY$age_age = XY$Age * XY$Age
XY$age_sex =  XY$Age * XY$Sex
XY$age_sex2 = XY$Age * XY$Sex * XY$Sex
#XY <- scale(XY) # normalize our covariates, X
XY <- XY[ , colSums(is.na(XY)) == 0] # no NANs
XY <- cbind(XY, 1) # Intercept term
XY <- data.frame(XY)

null_model = fit_null_glm(Diagnosis~Sex + Age + apoe_e4 + apoe_e2	+ PC1	+PC2+	PC3+	PC4+PC5+	PC6+PC7+	PC8+PC9+PC10+	Illumina_HiSeq_2000 +	Illumina_HiSeqX	+Illumina_NovaSeq	+Illumina+USUHS+USUHS.Miami+NYGC+MEDGENOME+	Baylor+Broad+WashU+PRS_5 , data = XY, "binomial")

Zs <- list.files(path=paste0('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/Alzheimer-RV/data/gene_matrices_maf/', TYPE_, '/chr', CHRO_NB), pattern="*_anno_nmf.csv", full.names=TRUE, recursive=FALSE)
Gs <- list.files(path=paste0('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/Alzheimer-RV/data/gene_matrices_maf/', TYPE_, '/chr', CHRO_NB), pattern="*geno_imputed.csv", full.names=TRUE, recursive=FALSE)
Z <- read.csv(Zs[1])
features_selected = c(colnames(Z)[3:length(colnames(Z))]) # basically all annotations -- this can be NMF too


has_fewer_than_3_unique <- function(col) {
  length(unique(col)) < 3
}

staar_output <- data.frame(matrix(ncol = 0, nrow = 0))

count <-0
for(gene in 1:length(Zs)){
  count<-count+1
  print(paste("Gene",count,"over",length(Zs)))
  tryCatch({
    Z <- read.csv(Zs[gene])
    Z <- Z[features_selected] 
    #Z <- as.matrix(Z[ , colSums(is.na(scale(Z)))==0])
    Z <- Z[, !apply(Z, 2, has_fewer_than_3_unique)]
    G <- as.matrix(read.csv(Gs[gene]))
    Z <- cbind(intercept = c(rep(1,length(dim(Z)[1]))), Z) # intercept
    
    result <- STAAR(genotype = G, obj_nullmodel = null_model, annotation_phred = Z, rare_maf_cutoff = 0.05, rv_num_cutoff = 2) # maybe increase rv num cutoff
    gene_name <- substr(basename(Zs[gene]), 1, nchar(basename(Zs[gene]))-13)
    staar_output[nrow(staar_output) + 1, c("Gene", "STAAR-O", "ACAT-O")] <- c(gene_name, result$results_STAAR_O, result$results_ACAT_O)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

path = "/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/Alzheimer-RV/data/STAAR"
dir.create(paste0(path,OUTPUT))
dir.create(paste0(path,OUTPUT,"/", TYPE_))
dir.create(paste0(path,OUTPUT,"/", TYPE_, "/", ITERATION))
write.table(staar_output, file=paste0(path,OUTPUT,"/", TYPE_, "/", ITERATION,'/by_region_chr',CHRO_NB,'.csv'), quote=FALSE, sep='\t')





