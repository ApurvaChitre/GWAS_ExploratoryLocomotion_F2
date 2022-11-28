# file constants ----

HEADER_FILE <- "/projects/ps-palmer/apurva/round8/gcta_assoc_header.txt"
QTL_DIR <- "QTL_locuszoom/"
COV_FILE <- "temp.txt"
ASSOC_FILE <- "temp_gcta_cond.assoc.txt.mlma"
FAM_FILE <- paste0("/projects/ps-palmer/apurva/u01_huda_akil/genotypes/",
                   "u01_huda_akil_genotypes")
PHENO_DIR <- "QTL_locuszoom/pheno/"
ASSOC_DIR <- "QTL_locuszoom/assoc/"
TOTAL_N_FILE <- "/projects/ps-palmer/apurva/u01_huda_akil/total_N.csv"
DOSAGES_FILE <- paste0("/projects/ps-palmer/apurva/u01_huda_akil/genotypes",
                       "/u01_huda_akil_genotypes_dosages_n534.txt")
GCTA_LOC <- paste0("/projects/ps-palmer/software/local/src/gcta_1.93.2beta",
                   "/gcta64")
GRM_LOC <- "/projects/ps-palmer/apurva/u01_huda_akil/grm/"
CHR_DIR <- "QTL_locuszoom/conditional_analysis/chr"

# script ----

args <- commandArgs(trailingOnly = TRUE)
exp_dir <- args[1]
expname <- args[2]
trait <- args[3]
chr <- args[4]

exp_dir<-paste0(exp_dir,"/")

header<-read.table(HEADER_FILE,header=T,stringsAsFactors=F)
load(paste0(exp_dir,"residuals/",expname,"_traits.RData"))

library(data.table)

chr <- as.numeric(gsub("chr","",chr))
cov_file=paste0(exp_dir,QTL_DIR,chr,"_",trait,"_",COV_FILE)
assoc_gcta=paste0(exp_dir,QTL_DIR,chr,"_",trait,"_",ASSOC_FILE)
fam=read.table(paste0(FAM_FILE, ".fam"),header=F,stringsAsFactors = F)
fam=fam[,c("V1","V2")]
if (file.exists(cov_file)) {
  #Delete file if it exists
  file.remove(cov_file)
}

if (file.exists(assoc_gcta)) {
  #Delete file if it exists
  file.remove(assoc_gcta)
}

pheno<-paste0(exp_dir,PHENO_DIR,trait,".txt")  
assoc<-paste0(exp_dir,ASSOC_DIR,trait,".loco.mlma")  

total_N=read.csv(TOTAL_N_FILE,header=T,stringsAsFactors=F)

expsub<- paste0(expname,"_")
trait_name_n<-gsub(expsub,"",trait)

sample_size<-total_N[which(total_N$trait==trait_name_n),"N"]

ifelse(sample_size>300,threshold<-6.664,threshold<-6.254 )

assoc_for_pval<-list()

topSNP_part2_list<- list()
topSNP_part2<-data.frame()
cov_counter<- 0
loop_run <- TRUE
assoc_counter<-1

while(loop_run){
  if(cov_counter==0){
    CHR<-fread(paste0('cat ',assoc,' | cut -f 1,2,3,7,9 '),
               header=T,stringsAsFactors =F)
    column_names<-c('Chr','SNP','bp','b','p')
    setnames(CHR, column_names)
    CHR$log10P = -log10(CHR$p)
    
    CHR<-CHR[Chr==chr] 
    assoc_for_pval[[assoc_counter]]<-CHR
    
  }else{
    CHR=fread(paste0('cat ',exp_dir,QTL_DIR,chr,'_',trait,'_',ASSOC_FILE,
                     ' | cut -f 1,2,3,7,9 '),header=T,stringsAsFactors = F)  
    
    column_names<-c('Chr','SNP','bp','b','p')
    setnames(CHR, column_names)
    CHR$log10P = -log10(CHR$p)
    assoc_for_pval[[assoc_counter]]<-CHR
  }
  
  if(!cov_counter==0){
    file.remove(assoc_gcta) 
  }
  
  topsnp<-CHR$SNP[which(CHR$log10P==max(CHR$log10P,na.rm=T))]
  
  if(length(topsnp)>1){
    dups<-CHR[which(CHR$SNP %in% topsnp),]
    dups<-dups[order(dups$bp),]
    topsnp<-dups$SNP[which.max((abs(dups$b)))]  
  }
  
  topsnp_log10P<-CHR$log10P[which(CHR$SNP== topsnp)]
  
  #update topSNP_part2 table
  trait=trait
  
  ifelse(cov_counter==0,covariate <- NA,covariate <- cond_topsnp)
  if(cov_counter>0){
    p_val_covar <- vector(mode="character", length=cov_counter)
    for(m in 1:length(p_val_covar)){
      p_val_covar[m]<-assoc_for_pval[[m]][
        which(assoc_for_pval[[m]]$SNP==covariate),"log10P"]
      
    }  
    covar_pval<-paste(p_val_covar, sep="",collapse="_")  
  }else{
    covar_pval<-NA
  }
  
  topsnp=topsnp
  topsnp_log10P=topsnp_log10P
  
  df<-data.frame()
  df <- data.frame(trait,topsnp,topsnp_log10P,covariate,
                   covar_pval,stringsAsFactors = F)
  
  topSNP_part2  <- rbind (topSNP_part2,df)
  
  covariate <- topsnp
  
  ds_command <- paste0("grep -w '",topsnp,"' ", DOSAGES_FILE)
  
  if(cov_counter == 0){
    cov=fread(cmd=ds_command,header=F,stringsAsFactors = F)
    #first 3 columns should be removed
    colsToDelete <-  colnames(cov)[c(1,2,length(colnames(cov)))]
    cov[, (colsToDelete) := NULL]  
    cov<-t(cov)  
    cov <- as.data.frame(cov)  
    cov<-cbind(fam,cov)
    cov_file=paste0(exp_dir,QTL_DIR,chr,"_",trait,"_","temp.txt")
    write.table(cov,file=cov_file,row.names=F,col.names=F,quote=F)
  }else{
    c = fread(cmd=ds_command,header=F,stringsAsFactors = F)
    colsToDelete <-  colnames(c)[c(1,2,length(colnames(c)))]
    c[, (colsToDelete) := NULL] 
    c <- t(c)
    c <- as.data.frame(c) 
    cov<-cbind(cov,c)
    cov_file=paste0(exp_dir,QTL_DIR,chr,"_",trait,"_","temp.txt")
    write.table(cov,file=cov_file,row.names=F,col.names=F,quote=F)
  }
  
  if(cov_counter == 0){
    f<-paste0(exp_dir,"QTL_locuszoom/code/",trait,"_",chr,".sh")
    if (file.exists(f)) file.remove(f)
    d1<-paste0("#!/bin/bash")
    ##This needs to be changed for each pheno file name
    d2<-paste0("#PBS -N ",trait,"_",chr)  
    d3<-paste0("#PBS -S /bin/bash")
    d4<-paste0("#PBS -l walltime=6:00:00")
    d5<-paste0("#PBS -l nodes=1:ppn=5")
    d6<-paste0("#PBS -j oe")
    d7<-paste0("#PBS -o ",exp_dir,"QTL_locuszoom/pbs_log/",trait,"_",chr,".out") 
    d8<-paste0("#PBS -q condo")  
    part1 <- paste0(GCTA_LOC, " --mlma --grm ", 
                    GRM_LOC, "u01_huda_akil_genotypes")
    part2 <- paste0("--mlma-subtract-grm ",GRM_LOC,"chr",chr,
                    ".u01_huda_akil_genotypes --bfile ",
                    FAM_FILE," --chr ",chr," ")
    part3<- paste0("--pheno ",pheno," --qcovar ",cov_file," ")
    part4<- paste0("--out ",exp_dir,"QTL_locuszoom/",chr,
                   "_",trait,"_","temp_gcta_cond.assoc.txt ")
    part5<- paste0("--thread-num 5")
    command <- paste0(part1,part2,part3,part4,part5)
    line=c(d1,d2,d3,d4,d5,d6,d7,d8,command)
    write(line,file=paste0(exp_dir,"QTL_locuszoom/code/",trait,"_",chr,".sh"),
          append=TRUE) 
    
    system(paste0("qsub ",exp_dir,"QTL_locuszoom/code/",trait,"_",chr,".sh"),
           ignore.stdout = T,ignore.stderr = T,wait = T)
  }else{
    system(paste0("qsub ",exp_dir,"QTL_locuszoom/code/",trait,"_",chr,".sh"),
           ignore.stdout = T,ignore.stderr = T,wait = T)    
  } 

  #check if pbs_log file exists

  while (!file.exists(paste0(exp_dir,QTL_DIR,chr,'_',trait,'_',ASSOC_FILE))) {
    Sys.sleep(100)
  }
  
  CHR=fread(paste0('cat ',exp_dir,QTL_DIR,chr,'_',trait,'_',ASSOC_FILE,
                   ' | cut -f 1,2,3,7,9 '),header=T,stringsAsFactors = F)  
  column_names<-c('Chr','SNP','bp','b','p')
  setnames(CHR, column_names)
  CHR$log10P = -log10(CHR$p)
  
  nrow_next<-nrow(CHR[which(CHR$log10P>threshold),])
  
  if(nrow_next==0){
    loop_run<-FALSE
    rm(cov,ds_command)
    
  }else{
    cond_topsnp<-topsnp
    cov_counter=cov_counter+1  
    assoc_counter=assoc_counter+1
  }
}

save(topSNP_part2,
     file=paste0(exp_dir,CHR_DIR,chr,"_",trait,"_QTL_conditional_analysis.RData"))
write.csv(topSNP_part2,
          paste0(exp_dir,CHR_DIR,chr,"_",trait,"_QTL_conditional_analysis.csv"),
          row.names=F,quote=F)

topSNP_part2_list[[1]]<- topSNP_part2

no_covar<-do.call("rbind",topSNP_part2_list)
write.csv(no_covar,
          paste0(exp_dir,CHR_DIR,chr,"_",trait,"_QTL_conditional_analysis.csv"),
          row.names=F,quote=F)

rm(list=ls())