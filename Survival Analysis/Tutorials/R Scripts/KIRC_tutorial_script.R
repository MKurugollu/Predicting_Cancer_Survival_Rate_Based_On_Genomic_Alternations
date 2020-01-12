# Tutorial: SA of TCGA patients integrating gene expression
# https://www.biostars.org/p/153013/

#install.packages(c("survival", "survminer", "BiocManager"))
#BiocManager::install("limma")

library(survival)
library(limma)


#set directory to where the RNA and Clinical folders are
setwd("D:/Documents/lv4_Project/Predicting_Cancer_Survival_Rate_Based_On_Genomic_Alternations/Survival Analysis/Tutorials/Data/KIRC")



rna <- read.table('RNA/KIRC.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt',nrows=20533, header=T,row.names=1,sep='\t')
rna<-rna[-1,] #don't need 1st row

clinical <- as.data.frame(t(read.table('Clinical/KIRC.merged_only_clinical_clin_format.txt', header = T, row.names = 1, sep = '\t')))
#warning: incomplete final line found by readTableHeader on 'Clinical/KIRC.merged_only_clinical_clin_format.txt'


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


# do the same to death
ind_keep <- grep('days_to_death',colnames(clinical))
death <- as.matrix(clinical[,ind_keep])
death_collapsed <- c()
for (i in 1:dim(death)[1]){
  if ( sum ( is.na(death[i,])) < dim(death)[2]){
    m <- max(death[i,],na.rm=T)
    death_collapsed <- c(death_collapsed,m)
  } else {
    death_collapsed <- c(death_collapsed,'NA')
  }
}

# and days last follow up here we take the most recent which is the max number
ind_keep <- grep('days_to_last_followup',colnames(clinical))
fl <- as.matrix(clinical[,ind_keep])
fl_collapsed <- c()
for (i in 1:dim(fl)[1]){
  if ( sum (is.na(fl[i,])) < dim(fl)[2]){
    m <- max(fl[i,],na.rm=T)
    fl_collapsed <- c(fl_collapsed,m)
  } else {
    fl_collapsed <- c(fl_collapsed,'NA')
  }
}

# and put everything together
all_clin <- data.frame(new_tum_collapsed,death_collapsed,fl_collapsed)
colnames(all_clin) <- c('new_tumor_days', 'death_days', 'followUp_days')
                       
# create vector with time to new tumor containing data to censor for new_tumor
all_clin$new_time <- c()
for (i in 1:length(as.numeric(as.character(all_clin$new_tumor_days)))){
  all_clin$new_time[i] <- ifelse ( is.na(as.numeric(as.character(all_clin$new_tumor_days))[i]),
                                   as.numeric(as.character(all_clin$followUp_days))[i],as.numeric(as.character(all_clin$new_tumor_days))[i])
}


# create vector time to death containing values to censor for death
all_clin$new_death <- c()
for (i in 1:length(as.numeric(as.character(all_clin$death_days)))){
  all_clin$new_death[i] <- ifelse ( is.na(as.numeric(as.character(all_clin$death_days))[i]),
                                    as.numeric(as.character(all_clin$followUp_days))[i],as.numeric(as.character(all_clin$death_days))[i])
}


# create vector for death censoring
table(clinical$patient.vital_status)


all_clin$death_event <- ifelse(clinical$patient.vital_status == 'alive', 0,1)
rownames(all_clin) <- clinical$IDs

# create event vector for RNASeq data
event_rna <- t(apply(z_rna, 1, function(x) ifelse(abs(x) > 1.96,1,0)))


# since we need the same number of patients in both clinical and RNASeq data take the indices for the matching samples
ind_tum <- which(unique(colnames(z_rna)) %in% rownames(all_clin))
ind_clin <- which(rownames(all_clin) %in% colnames(z_rna))

# pick your gene of interest
ind_gene <- which(rownames(z_rna) == 'TP53')

table(event_rna[ind_gene,])


# run survival analysis
s <- survfit(Surv(as.numeric(as.character(all_clin$new_death))[ind_clin],all_clin$death_event[ind_clin])~event_rna[ind_gene,ind_tum])
s1 <- tryCatch(survdiff(Surv(as.numeric(as.character(all_clin$new_death))[ind_clin],all_clin$death_event[ind_clin])~event_rna[ind_gene,ind_tum]), error = function(e) return(NA))

pv <- ifelse ( is.na(s1),next,(round(1 - pchisq(s1$chisq, length(s1$n) - 1),3)))[[1]]

p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)

print(s1)
print(s)

plot(survfit(Surv(as.numeric(as.character(all_clin$new_death))[ind_clin],all_clin$death_event[ind_clin])~event_rna[ind_gene,ind_tum]),
     col=c(1:3), frame=F, lwd=2,main=paste('KIRK',rownames(z_rna)[ind_gene],sep='\n'))

# add lines for the median survival
x1 <- ifelse ( is.na(as.numeric(summary(s)$table[,'median'][1])),'NA',as.numeric(summary(s)$table[,'median'][1]))
x2 <- as.numeric(summary(s)$table[,'median'][2])
if(x1 != 'NA' & x2 != 'NA'){
  lines(c(0,x1),c(0.5,0.5),col='blue')
  lines(c(x1,x1),c(0,0.5),col='black')
  lines(c(x2,x2),c(0,0.5),col='red')
}


# add legend
legend(1800,0.995,legend=paste('p.value = ',pv[[1]],sep=''),bty='n',cex=1.4)
legend(max(as.numeric(as.character(all_clin$death_days)[ind_clin]),na.rm = T)*0.7,0.94,
       legend=c(paste('NotAltered=',x1),paste('Altered=',x2)),bty='n',cex=1.3,lwd=3,col=c('black','red'))

for(i in rownames(z_rna)){
  if(i != "?"){
    print(i)
    ind_gene <- which(rownames(z_rna) == i)
    s <- survfit(Surv(as.numeric(as.character(all_clin$new_death))[ind_clin],all_clin$death_event[ind_clin])~event_rna[ind_gene,ind_tum])
    s1 <- tryCatch(survdiff(Surv(as.numeric(as.character(all_clin$new_death))[ind_clin],all_clin$death_event[ind_clin])~event_rna[ind_gene,ind_tum]), error = function(e) return(NA))
    
    pv <- ifelse ( is.na(s1),next,(round(1 - pchisq(s1$chisq, length(s1$n) - 1),3)))[[1]]
    print(s1)
    print(pv)
  }
}
