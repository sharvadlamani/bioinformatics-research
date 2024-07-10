library(SKAT)
#Sys.setenv(PATH = paste("/gpfs/commons/home/adas/R/cmake-3.29.0-linux-x86_64/bin", Sys.getenv("PATH"), sep=":"))

# 3 per minute
## Running SKAT, no annotations
args <- commandArgs(trailingOnly = TRUE)
CHRO_NB <-as.numeric(args[1]) # chromosome number
TYPE_ <- as.character(args[2]) # cell type: OLIG2, LHX2, NEUN, peripheralPU1nuclei, coding
ITERATION <- as.character(args[3]) # iteration that this run is or folder that this should be saved as
MAF <- numeric(args[4]) # Desired MAF
OUTPUT <- as.character(args[5]) # Output path
RACE <- as.character(args[6])

path = "/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/Alzheimer-RV/data/SKAT/"

# Load phenotype and covariate information

PHENOTYPE_PATH = "/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/Alzheimer-RV/data/phenotypes/36K_QC_filtered_final.csv"
XY <- read.table(file = PHENOTYPE_PATH, sep = ',', header = TRUE)
if(RACE != "ALL"){
  keep = which(XY$predicted_ancestry == RACE)
  XY = XY[XY$predicted_ancestry == RACE,]
}

Y <- XY[c("Diagnosis")]
y = c(Y)$Diagnosis
X <- XY[c('Sex','Age','apoe_e4','apoe_e2','PC1', 'PC2','PC3','PC4','PC5','PC6', 'PC7','PC8','PC9','PC10', 'Illumina_HiSeq_2000', 'Illumina_HiSeqX','Illumina_NovaSeq', 'Illumina', 'USUHS', 'USUHS.Miami','NYGC', 'MEDGENOME', 'Baylor', 'Broad', 'WashU', 'PRS_5')] # Choose covariates


X$age_age = X$Age * X$Age
X$age_sex =  X$Age * X$Sex
X$age_sex2 = X$Age * X$Sex * X$Sex
X <- scale(X) # normalize our covariates, X
X <- X[ , colSums(is.na(X)) == 0] # no NANs
X <- cbind(X, 1) # Intercept term


# Run the null model (y ~ X) from covariates
obj.b<-SKAT_Null_Model(y ~ X, out_type = "D")

# Get gene matrices and prep output
Gs <- list.files(path=paste0('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/Alzheimer-RV/data/gene_matrices_maf/', TYPE_, '/chr', CHRO_NB), pattern="*geno_imputed.csv", full.names=TRUE, recursive=FALSE)
df_p_value <- data.frame(matrix(ncol = 0, nrow = 0))

# Run SKAT, SKAT0, and Burden tests for all genes in given chromosome and cell type
count <-0
for(gene in 1:length(Gs)){
  #for(gene in 1:2){
  count <- count + 1
  print(paste("Gene", count, "over", length(Gs), "MAF:", MAF))
  tryCatch({
    G <- read.csv(Gs[gene])
    S <- read.csv(Gs[gene])
    if (RACE != "ALL") {
      G <- G[keep, ]
    }
    
    # # Calculate column sums
    column_sums <- colSums(G)
    # # Filter columns where the maf is less than or equal to the argument
    max_instances <- floor(nrow(G)*2*MAF)
    filtered_data_G <- G[, column_sums <= max_instances]
    
    #filtered_data_G <-G
    print(dim(filtered_data_G))
    
    gene_name = strsplit(basename(Gs[gene]), "_")[[1]][1]
    
    # Ensure the filtered data has correct dimensions
    if (ncol(filtered_data_G) > 0 && nrow(filtered_data_G) == length(y)) {
      burden = SKAT(as.matrix(filtered_data_G), obj.b, method = "Burden")$p.value
      skato = SKAT(as.matrix(filtered_data_G), obj.b, method = "SKATO")$p.value
      skat = SKAT(as.matrix(filtered_data_G), obj.b, method = "SKAT")$p.value
      
      df_p_value[nrow(df_p_value) + 1, c("Gene", "Burden", "SKATO", "SKAT")] <- c(gene_name, burden, skato, skat)
    } else {
      print("Filtered data dimensions are inconsistent")
    }
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

my_output_path = "/gpfs/commons/home/svadlamani/" #saving output to a folder in Shar's home directory

# Save output
dir.create(paste0(my_output_path,OUTPUT))
dir.create(paste0(my_output_path,OUTPUT,"/", MAF))
dir.create(paste0(my_output_path,OUTPUT,"/", MAF, "/",TYPE_))
dir.create(paste0(my_output_path,OUTPUT,"/", MAF, "/", TYPE_, "/", ITERATION))
write.table(df_p_value, file=paste0(my_output_path,OUTPUT,"/", MAF,"/", TYPE_, "/", ITERATION,'/chr',CHRO_NB,'.csv'), quote=FALSE, sep='\t', row.names = FALSE)