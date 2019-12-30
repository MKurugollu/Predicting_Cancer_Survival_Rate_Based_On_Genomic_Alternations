# Tutorial: SA of TCGA patients integrating gene expression
# https://www.biostars.org/p/153013/

install.packages(c("survival", "survminer", "BiocManager"))

#set directory to where the RNA and Clinical folders are
setwd("D:/Documents/lv4_Project/Predicting_Cancer_Survival_Rate_Based_On_Genomic_Alternations/Survival Analysis/Tutorials/Data/KIRC")

getwd() #check to see if we are in correct directory
library(survival)

rna <- read.table('RNA/KIRC.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt',nrows=20533, header=T,row.names=1,sep='\t')
rna<-rna[-1,] #don't need 1st row

clinical <- t(read.table('Clinical/KIRC.merged_only_clinical_clin_format.txt', header = T, row.names = 1, sep = '\t'))


# first remove genes whose expression is ==0 in more than 50% of the samples
rem <- function(x){
  x <- as.matrix(x)
  x <- t(apply(x,1,as.numeric))
  r <- as.numeric(apply(x,1,function(i) sum(i == 0)))
  remove <- which(r > dim(x)[2]*0.5)
  return(remove)
}

#apply 
remove <- rem(rna)
rna <- rna[-remove,]


table(substr(colnames(rna),14,14))

#get index of the normal/control samples
n_index <-which(substr(colnames(rna),14,14) == '1')
t_index <-which(substr(colnames(rna),14,14) == '0')


BiocManager::install("limma")
library(limma)
# apply voom function from limma package to normalize the data
vm <- function(x){
  cond <- factor(ifelse(seq(1,dim(x)[2],1) %in% t_index, 1,  0))
  d <- model.matrix(~1+cond)
  x <- t(apply(x,1,as.numeric))
  ex <- voom(x,d,plot=F)
  return(ex$E)
}

rna_vm <- vm(rna)
colnames(rna_vm) <- gsub('\\.', '-', substr(colnames(rna),1,12))


hist(rna_vm)

rm(rna)

#z = [(value gene X in tumor Y)-(mean gene X in normal)]/(standard deviation X in normal)

# calculate z-scores
scal <- function(x,y){
  mean_n <- rowMeans(y)  # mean of normal
  sd_n <- apply(y,1,sd)  # SD of normal
  # z score as (value - mean normal)/SD normal
  res <- matrix(nrow=nrow(x), ncol=ncol(x))
  colnames(res) <- colnames(x)
  rownames(res) <- rownames(x)
  for(i in 1:dim(x)[1]){
    for(j in 1:dim(x)[2]){
      res[i,j] <- (x[i,j]-mean_n[i])/sd_n[i]
    }
  }
  return(res)
}
z_rna <- scal(rna_vm[,t_index],rna_vm[,n_index])
# set the rownames keeping only gene name
rownames(z_rna) <- sapply(rownames(z_rna), function(x) unlist(strsplit(x,'\\|'))[[1]])

rm(rna_vm) #we don't need it anymore



# match the patient ID in clinical data with the colnames of z_rna
barcodes <- clinical[,'patient.bcr_patient_barcode']


#clinical$IDs <- toupper(clinical$patient.bcr_patient_barcode)
clinical$IDs <- toupper(clinical[,'patient.bcr_patient_barcode'])
sum(clinical$IDs %in% colnames(z_rna)) # we have 529 patients that we could use

# get the columns that contain data we can use: days to death, new tumor event, last day contact to....
ind_keep <- grep("days_to_new_tumor_event_after_initial_treatment",colnames(clinical))

ind_keep

new_tum <- as.matrix(clinical[,ind_keep])
new_tum_collapsed <- c()
for (i in 1:dim(new_tum)[1]){
  if ( sum ( is.na(new_tum[i,])) < dim(new_tum)[2]){
    m <- min(new_tum[i,],na.rm=T)
    new_tum_collapsed <- c(new_tum_collapsed,m)
  } else {
    new_tum_collapsed <- c(new_tum_collapsed,'NA')
  }
}