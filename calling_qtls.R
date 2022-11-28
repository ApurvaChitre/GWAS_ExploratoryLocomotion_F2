# file constants ----

TRAITS_DATA <- paste0("/projects/ps-palmer/apurva/u01_huda_akil/residuals",
                      "/u01_huda_akil_traits.RData")
ASSOC_DIR <- "QTL_locuszoom/assoc/"
PHENO_DIR <- "QTL_locuszoom/pheno/"
CODE_DIR <- "QTL_locuszoom/code/"
PBS_LOG_DIR <- "QTL_locuszoom/pbs_log/"
TOTAL_N_CSV <- "/projects/ps-palmer/apurva/u01_huda_akil/total_N.csv"
WD <- "/projects/ps-palmer/apurva/u01_huda_akil/r2/"
OASIS_LOC <- "/oasis/tscc/scratch/aschitre/software/plink-1.90/plink"
LOCAL_LOC <- paste0("/projects/ps-palmer/apurva/u01_huda_akil/genotypes/",
                    "u01_huda_akil_genotypes")
MLMA_FILE <- "QTL_locuszoom/temp_gcta_cond.loco.mlma"
TEMP_FILE <- "QTL_locuszoom/temp.txt"
RECODE_LOC <- "/projects/ps-palmer/apurva/u01_huda_akil/r2/recodeA"
OASIS_GCTA <- "/oasis/tscc/scratch/aschitre/gcta/gcta64"
LOCAL_GRM <- paste0("projects/ps-palmer/apurva/u01_huda_akil/grm",
                    "/u01_huda_akil_genotypes")

# end ----

expname <- "u01_huda_akil"

base_dir <- paste0("/projects/ps-palmer/apurva/",expname)

load(TRAITS_DATA)

traitlist<-traits

library(data.table)

setwd(paste0(exp_dir))

system(paste0("mkdir ",exp_dir,ASSOC_DIR))
system(paste0("mkdir ",exp_dir,PHENO_DIR))
system(paste0("mkdir ",exp_dir,CODE_DIR))
system(paste0("mkdir ",exp_dir,PBS_LOG_DIR))
setwd(paste0(exp_dir,ASSOC_DIR))

all_traits_chr<-list.files(path=paste0(exp_dir,"/results"),full.names = T)

filenames<-basename(all_traits_chr)
filenames<-gsub("allChr_","",filenames)
filenames<-gsub(".loco.mlma","",filenames)

column_names<-c('Chr','SNP','bp','b','p')

readdata<-function(x){
  data<-read.table(x,header=T,stringsAsFactors =F)
  data<-data[,c(1,2,3,7,9)]
  colnames(data)<-c('Chr','SNP','bp','b','p')
  data=data[which(!is.nan(data$p)),]
  data$log10P = -log10(data$p)
  return(data)
}

all_data<-lapply(all_traits_chr,readdata)

names(all_data)<-filenames

all_traits_chr<-filenames

topSNP_part1<-data.frame()

total_N=read.csv(TOTAL_N_CSV,header=T,stringsAsFactors=F)

for(i in seq_along(all_traits_chr)){
  trait_name<-gsub(paste0(expname,"_"),"",all_traits_chr[i])
  
  sample_size<-total_N$N[which(total_N$trait==trait_name)]
  
  ifelse(sample_size>300,threshold<-6.664,threshold<-6.254 )
  
  if(length(which(all_data[[all_traits_chr[i]]]$log10P > threshold  )) > 0){
    trait=all_traits_chr[i]
    above_threshold="yes"
    
    df <- data.frame(trait,above_threshold,stringsAsFactors = F)
  }else{
    trait=all_traits_chr[i]
    above_threshold="no"
    df <- data.frame(trait,above_threshold,stringsAsFactors = F)
  }
  
  topSNP_part1  <- rbind (topSNP_part1,df)
}

signi_traits_n<-length(which(topSNP_part1$above_threshold=="yes"))
signi_traits<-topSNP_part1$trait[which(topSNP_part1$above_threshold=="yes")]

#this will contain chrs with signi results
signi_chrs <- setNames(replicate(signi_traits_n,data.frame()),signi_traits)

for(i in seq_along(signi_traits)){
  trait_name<-gsub(paste0(expname,"_"),"",signi_traits[i])
  
  sample_size<-total_N$N[which(total_N$trait==trait_name)]
  
  ifelse(sample_size>300,threshold<-6.664,threshold<-6.254 )
  
  chrs<-unlist(unique(all_data[[signi_traits[i]]][
    which(all_data[[signi_traits[i]]]$log10P > threshold) ,"Chr" ]))
  chrs<-unname(chrs)
  
  signi_chrs[[signi_traits[i]]]<-chrs
}

summary_QTL <- list()
topSNP_part2_list <- list ()

system(paste("mkdir", WD))
setwd(WD)

#define parameters that can be soft coded
#This is log10P - 1.5 (to see if there are any supporting SNPs)
#sup_P_term<-1.5
sup_P_term<-2
#This is 0.5 MB flanking dist- to see if there are any supporting snps for 
#top snp or if it is a rogue snp
sup_dist<-500000
#LD command parameters for plink
ld_window_kb=11000
ld_r2=0.4
#QTL boundary distance 
qtl_dist<-1000000

for(i in seq_along(signi_traits)){
  topSNP_part2_list <- list ()  
  chrs <-signi_chrs[[signi_traits[i]]]
  
  trait_name<-gsub(paste0("_",expname),"",signi_traits[i])
  
  sample_size<-total_N$N[which(total_N$trait==trait_name)]
  
  ifelse(sample_size>300,threshold<-6.664,threshold<-6.254 )
  
  for(j in seq_along(chrs)){
    CHR <- all_data[[signi_traits[i]]][
      which(all_data[[signi_traits[i]]]$Chr == chrs[j]),]
    
    topSNP_part2<-data.frame()
    
    loop_run <- TRUE
    while(loop_run){
      #run until there are no topsnps
      #find top snp
      
      all_topsnps<-CHR[which((CHR$log10P > threshold)),]
      topsnp<-
        all_topsnps$SNP[which(all_topsnps$log10P==max(all_topsnps$log10P))]
      
      if(length(topsnp)>1){
        dups<-CHR[which(CHR$SNP %in% topsnp),]
        dups<-dups[order(dups$bp),]
        topsnp<-dups$SNP[which.max((abs(dups$b)))]  
      }
      
      topsnp_log10P<-CHR$log10P[which(CHR$SNP== topsnp)]
      sup_P<-(topsnp_log10P - sup_P_term)
      
      chr_topsnp<-gsub("chr","",strsplit(topsnp,split=":")[[1]][1])
      ps_topsnp<-strsplit(topsnp,split=":")[[1]][2]
      ps_topsnp<-as.numeric(ps_topsnp)
      
      start<- ps_topsnp -sup_dist
      stop<- ps_topsnp +sup_dist
      
      ##This evaluates if it is not a rogue SNP
      if(nrow(CHR[which((CHR$bp %in% seq(start,stop)) & 
                        (CHR$log10P > sup_P)),])>2){
        trait=signi_traits[i]  
        topsnp=topsnp
        QTL="yes"
        topsnp_log10P=topsnp_log10P
        df<-data.frame()
        df <- data.frame(trait,topsnp,QTL,topsnp_log10P,stringsAsFactors = F)
      }else{
        trait=signi_traits[i]  
        topsnp=topsnp
        QTL="no"
        topsnp_log10P=topsnp_log10P
        df<-data.frame()
        df <- data.frame(trait,topsnp,QTL,topsnp_log10P,stringsAsFactors = F)
      }
      
      topSNP_part2  <- rbind (topSNP_part2,df)
      
      system(paste0(OASIS_LOC," --bfile ",LOCAL_LOC," --chr ",chr_topsnp,
                    " --nonfounders --r2  --ld-snp ",topsnp,
                    " --ld-window 1000000 -ld-window-kb ",ld_window_kb,
                    " --ld-window-r2 ",ld_r2," --out ",WD,"temp_qtl_n"),
             ignore.stdout = T,ignore.stderr = T,wait = T)
      
      #no need to sleep because system command prints output 
      #and then moves onto next line
      
      #find next topsnp on the same chr which is not in the following list
      
      f<-paste0(WD,"temp_qtl_n.ld")
      if (file.exists(f)){
        qtl_snps<-fread(paste(f,"| grep -v R2 | awk '{ print $6 }'"),sep=" ",
                        header = F,stringsAsFactors = F)
        
        #0.5 
        #1MB
        sub_start<-ps_topsnp-qtl_dist
        sub_stop<-ps_topsnp+qtl_dist
        
        nrow_next<-nrow(CHR[which((!CHR$SNP %in% c(qtl_snps$V1,topsnp)) & 
                                    (CHR$log10P>threshold) & 
                                    (!CHR$bp %in% seq(sub_start,sub_stop))),])
        
        if(!(nrow_next ==0)){
          CHR<-CHR[which((!CHR$SNP %in% c(qtl_snps$V1,topsnp)) & 
                           (!CHR$bp %in% seq(sub_start,sub_stop))),]
        }else{
          loop_run<-FALSE
        }
      }
      else{
        nrow_next<-nrow(CHR[which((!CHR$SNP %in% topsnp) & 
                                    (CHR$log10P>threshold) ),])
        
        if(!(nrow_next ==0)){
          CHR<-CHR[which((!CHR$SNP %in% topsnp)),] 
        }else{
          loop_run<-FALSE 
        }
      }
    }
    
    topSNP_part2_list[[j]]<-topSNP_part2
  }
  
  summary_QTL[[i]] <- topSNP_part2_list
}

names(summary_QTL) <- signi_traits

temp <- unlist(summary_QTL, recursive = FALSE)
final <- do.call("rbind", temp)
write.csv(final,paste0(exp_dir,"QTL_locuszoom/",expname,"_QTL.csv"),
          row.names=F,quote=F)
