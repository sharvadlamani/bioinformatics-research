library(dplyr)
library(FSTpackage)
library(data.table)
library(sqldf)
start.time <- Sys.time()

# the FST test requires the following tables:
# Y: outcomes, n by 1 matrix where n is the total number of observations aka individuals √
# X: covariates, n by d matrix √
# G: genotype matrix, n by p matrix where n is the total number of subjects, p the number of variants in a set
# Z: functional annotation matrix, p by q matrix, p is number of variants in set, q is their annotations. First column in Z should be 1 to use SKAT/burden test

# If we want to run the gene set test
# GeneSetID: p*2 matrix indicating the genes in which the variants are located, column 1 = gene name, column 2 = variant names

args <- commandArgs(trailingOnly = TRUE)
CHRO_NB <-as.numeric(args[1]) # chromosome number
TYPE_ <- as.character(args[2]) # cell type: OLIG2, LHX2, NEUN, peripheralPU1nuclei
ITERATION <- as.character(args[3]) # iteration that this run is or folder that this should be saved as
OUTPUT <- as.character(args[4])
RACE <- as.character(args[5])

# FOR TESTING
#CHRO_NB = 22
#TYPE_ = "promoter"
#ITERATION = "1"
print(paste("---- STARTING CHROMOSOME ",CHRO_NB, "FOR ", TYPE_, "VARIANTS", ITERATION, "NUMBER ITERATION", OUTPUT, "OUTPUT"))

PHENOTYPE_PATH = "/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/Alzheimer-RV/data/phenotypes/36K_QC_filtered_final.csv"
XY <- read.table(file = PHENOTYPE_PATH, sep = ',', header = TRUE)
if(RACE != "ALL"){
  keep = which(XY$predicted_ancestry == RACE)
  XY = XY[XY$predicted_ancestry == RACE,]
}

# Permutations ## COMMENT OUT IF NOT RUNNING PERMUTATIONS
#perm_order <- sample(nrow(XY))
#XY$Diagnosis <- XY$Diagnosis[perm_order]
#XY$PC1 <- XY$PC1[perm_order]
#XY$PC2 <- XY$PC2[perm_order]
#XY$PC3 <- XY$PC3[perm_order]
#XY$PC4 <- XY$PC4[perm_order]
#XY$PC5 <- XY$PC5[perm_order]
#XY$PC6 <- XY$PC6[perm_order]
#XY$PC7 <- XY$PC7[perm_order]
#XY$PC8 <- XY$PC8[perm_order]
#XY$PC9 <- XY$PC9[perm_order]
#XY$PC10 <- XY$PC10[perm_order]

Y <- XY[c("Diagnosis")]
rownames(Y) <- XY$SampleID
X <- XY[c('Sex','Age','apoe_e4','apoe_e2','PC1', 'PC2','PC3','PC4','PC5','PC6', 'PC7','PC8','PC9','PC10', 'Illumina_HiSeq_2000', 'Illumina_HiSeqX','Illumina_NovaSeq', 'Illumina', 'USUHS', 'USUHS.Miami','NYGC', 'MEDGENOME', 'Baylor', 'Broad', 'WashU', 'PRS_5')] # Choose covariates
X$age_age = X$Age * X$Age
X$age_sex =  X$Age * X$Sex
X$age_sex2 = X$Age * X$Sex * X$Sex
X <- scale(X) # normalize our covariates, X
X <- X[ , colSums(is.na(X)) == 0] # no NANs
X <- cbind(X, 1) # Intercept term
rownames(X) <- XY$SampleID

# Preliminary data management
result.prelim<-FST.prelim(Y,X=X,out_type='D')


Zs <- list.files(path=paste0('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/Alzheimer-RV/data/gene_matrices_maf/', TYPE_, '/chr', CHRO_NB), pattern="*_anno_nmf.csv", full.names=TRUE, recursive=FALSE)
Gs <- list.files(path=paste0('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/Alzheimer-RV/data/gene_matrices_maf/', TYPE_, '/chr', CHRO_NB), pattern="*geno_imputed.csv", full.names=TRUE, recursive=FALSE)
Z <- read.csv(Zs[1])
features_selected = c(colnames(Z)[3:length(colnames(Z))]) # basically all annotations -- this can be NMF too

# prep dataframes to be filled
df_dispersion <- data.frame(matrix(ncol = length(features_selected)+3, nrow = 0))
df_burden <- data.frame(matrix(ncol = length(features_selected)+3, nrow = 0))
df_d_b <- data.frame(matrix(ncol = length(features_selected)+3, nrow = 0))
df_p_value_individal_variant <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(df_dispersion) <- c("intercept", c(features_selected), "minP",'fisher')
colnames(df_burden) <- c("intercept", c(features_selected), "minP",'fisher')
colnames(df_d_b) <- c("intercept", c(features_selected), "minP",'fisher')
new_list <- c()
nb_variants <- 0 

#loop through genes
count <-0
for(gene in 1:length(Zs)){
  count<-count+1
  print(paste("Gene",count,"over",length(Zs)))
  tryCatch({
    Z <- read.csv(Zs[gene])
    G <- read.csv(Gs[gene])
    rownames(G) <- XY$SampleID
    rownames(Z) <- colnames(G) 
    #weight_maf <- 1 / as.numeric(Z[,'gnomAD_genomes_POPMAX_AF']) ## using GNOMAD
    #weight_maf <- 1 / as.numeric(colSums(G, na.rm = TRUE) / dim(G)[1]) ## using MAF from data
    Z <- Z[features_selected] # USE THIS WHEN INCLUDING ANNOTATIONS
    ## NEED TO include correctly scaled Z matrix, no case or adsp maf
    #Z$control_maf = -log(Z$control_maf)
    #Z$control_maf[is.infinite(Z$control_maf)] = 14
    #Z$control_maf = (Z$control_maf - 0.6931472) / (14 - 0.6931472)
    Z <- as.matrix(Z[ , colSums(is.na(scale(Z)))==0]) # removes columns w no variation # INCLUDE WHEN USING ANNOTATIONS
    #Z <- subset(Z, select = -c(case_maf, adsp_maf))
    Z <- cbind(intercept = c(rep(1,length(rownames(Z)))), Z) # intercept
    result <- FST.test(result.prelim,G, Z, B=5000) # could use weights = weight_maf, if nothing specificed uses beta
    df_dispersion[nrow(df_dispersion) + 1, c(colnames(Z), "minP","fisher")] <- c(result$p.value[1,])
    df_burden[nrow(df_burden) + 1, c(colnames(Z), "minP","fisher")] <- c(result$p.value[2,])
    df_d_b[nrow(df_d_b) + 1, c(colnames(Z), "minP","fisher")] <- c(result$p.value[3,])
    gene_name <- substr(basename(Zs[gene]), 1, nchar(basename(Zs[gene]))-13) # OR subtract 9 if not NMF
    df_p_single <- as.data.frame(result$p.single)
    df_p_single$gene <- rep(gene_name, length(rownames(df_p_single)))
    df_p_value_individal_variant <- rbind(df_p_value_individal_variant, df_p_single)
    new_list <-c(new_list, gene_name)
    nb_variants <- nb_variants + length(result$p.single)
    print(paste("number of variant in this gene: ",result$n.marker))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

print(paste("number genes skipped bc does not include any variant",length(Zs)-length(new_list)))
colnames(df_p_value_individal_variant) <-c("p-value",	"gene_name")
rownames(df_dispersion) <- new_list
rownames(df_burden) <- new_list
rownames(df_d_b) <- new_list


# group df
df_burden$type_test<-rep("Burden", length(rownames(df_burden)))
df_burden$gene_index<-1:length(rownames(df_burden))

df_dispersion$type_test<-rep("Dispersion", length(rownames(df_dispersion)))
df_dispersion$gene_index<-1:length(rownames(df_dispersion))

df_d_b$type_test<-rep("D+B", length(rownames(df_d_b)))
df_d_b$gene_index<-1:length(rownames(df_d_b))

df_p_value<- rbind(df_burden, df_dispersion, df_d_b)
df_p_value_raw <- data.frame(df_p_value)


# Ajust p_value with total number of genes

end.time <- Sys.time()
time.taken <- end.time - start.time
print('SUMMARY')
print(paste("number of genes tested for chromosome", CHRO_NB,":", length(new_list)))
print(paste("number of variants tested over all genes for chromosome", CHRO_NB,":", nb_variants))
print(time.taken)



path = "/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/Alzheimer-RV/data/FST/"
dir.create(paste0(path,OUTPUT))
dir.create(paste0(path,OUTPUT,"/", TYPE_))
dir.create(paste0(path,OUTPUT,"/", TYPE_, "/", ITERATION))
write.table(df_p_value, file=paste0(path,OUTPUT,"/", TYPE_, "/", ITERATION,'/by_region_chr',CHRO_NB,'.csv'), quote=FALSE, sep='\t')
write.table(df_p_value_individal_variant, file=paste0(path,OUTPUT,"/", TYPE_, "/", ITERATION,'/by_variant_chr',CHRO_NB,'.csv'), quote=FALSE, sep='\t')



#coding_region <- read.csv("/gpfs/commons/home/adas/Rare-Variant-AD/data/FST/results/all_coding_region.csv")
#noncoding_region <- read.csv("/gpfs/commons/home/adas/Rare-Variant-AD/data/FST/results/all_noncoding_region.csv")
#promoter_region <- read.csv("/gpfs/commons/home/adas/Rare-Variant-AD/data/FST/results/all_promoter_region.csv")
#all_regions <- rbind(coding_region, noncoding_region, promoter_region)

#df = all_regions[all_regions['type_test']=="D+B",]
#df = df[,c('gene','CHR','BP','fisher')]
#colnames(df) <- c('SNP','CHR','BP','P')
#manhattan(df, annotatePval = 1e-8)


