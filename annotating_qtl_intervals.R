# file constants ----

TOPSPNS_CSV <- paste0("/projects/ps-palmer/apurva/u01_huda_akil/QTL_locuszoom",
                      "/u01_huda_akil_final_QTL_list.csv")
WD <- "/projects/ps-palmer/apurva/u01_huda_akil/r2/"
OASIS_LOC <- "/oasis/tscc/scratch/aschitre/software/plink-1.90/plink"
LOCAL_LOC <- paste0("/projects/ps-palmer/apurva/u01_huda_akil/genotypes",
                    "/u01_huda_akil_genotypes")
OUTPUT_CSV <- paste0("/projects/ps-palmer/apurva/u01_huda_akil/QTL_locuszoom",
                     "/u01_huda_akil_table3.csv")
LIB_LOC <- "/home/aschitre/R_libs/"
URL <- "https://rest.rgd.mcw.edu/rgdws/genes/"

# end ----

library(data.table)

topsnps=read.csv(TOPSNPS_CSV,header=T,stringsAsFactors = F)

topsnps$chr<-sapply(strsplit(topsnps$topsnp,split=":"),`[`,1)
topsnps$pos<-sapply(strsplit(topsnps$topsnp,split=":"),`[`,2)

column_names<-c("snp1","snp2","rsquare")
LD_start_stop<-function(x,trait,chr,pos){
  data<-fread(x,header=T,stringsAsFactors =F,select=c(3,6,7))
  setnames(data, column_names)
  data$dprime<-rep(NA,nrow(data))
  data<-data[which(data$rsquare>=0.6),]
  
  #split snp1
  data[, c("chr", "pos") := tstrsplit(snp2, ":", fixed=TRUE)]
  
  data$pos <- as.numeric(data$pos)
  data<-data[order(pos),]    
  LD_interval_start<-data$pos[1]
  LD_interval_stop<-data$pos[nrow(data)]    
  output<-data.frame(LD_interval_start,LD_interval_stop)
  return(output)  
}

setwd(WD)
out<-data.frame()

for(i in 1:nrow(topsnps)){
  trait<-topsnps$trait[i]  
  topsnp<-topsnps$topsnp[i]
  chr<-topsnps$chr[i]
  pos<-topsnps$pos[i]
  
  system(paste0(OASIS_LOC," --bfile ",LOCAL_LOC," --chr ",gsub("chr","",chr),
                " --nonfounders --r2  --ld-snp ",topsnp," ",
                "--ld-window 100000 -ld-window-kb 6000 --ld-window-r2 0 --out ",
                WD,"temp_qtl_n"),ignore.stdout = T,ignore.stderr = T,wait = T)
  
  test<-LD_start_stop("temp_qtl_n.ld",trait = trait,chr=chr,pos=pos)
  out=rbind(out,test)
}

final<-cbind(topsnps,out)

final$LD_interval_size_bp <- final$LD_interval_stop - final$LD_interval_start

setwd(paste0(exp_dir,"QTL_locuszoom/assoc"))

column_names<-c("Chr","SNP","A1","A2","Freq","b","se","p")

getassocAfBetaSe<-function(trait,topsnp,chr){
  assoc_name<-paste0(trait,".loco.mlma")  
  data  =read.table(assoc_name,header=T,stringsAsFactors=F)
  
  data=data[,c(1,2,4,5,6,7,8,9)]
  
  colnames(data)<-column_names
  
  data<-data[which(data$Chr==chr),] 
  rs<-data$SNP[which(data$SNP == topsnp)]  
  af<-data$Freq[which(data$SNP == topsnp)]  
  beta<-data$b[which(data$SNP == topsnp)]  
  se<-data$se[which(data$SNP == topsnp)]  
  allele1<-data$A1[which(data$SNP == topsnp)]
  allele2<-data$A2[which(data$SNP == topsnp)]
  p_score<-data$p[which(data$SNP == topsnp)]
  output<-data.frame(trait,rs,af,beta,se,allele1,allele2,p_score)
  return(output)  
}

out<-NULL

for(i in 1:nrow(final)){
  
  trait<-final$trait[i]
  topsnp<-final$topsnp[i]
  chr<-gsub("chr","",final$chr[i])
  test<-getassocAfBetaSe(trait=trait,topsnp=topsnp,chr=chr)
  out=rbind(out,test)  
}

out$trait_snp<-NA
out$trait_snp = paste0(out$trait,"_",out$rs)

final$trait_snp <- paste0(final$trait,"_",final$topsnp)

merged=merge(final,out,by="trait_snp",all=F)

write.csv(merged,file=OUTPUT_CSV,row.names=F,quote=F)

library(httr)
library(RJSONIO,lib.loc=LIB_LOC)
library(openxlsx,lib.loc=LIB_LOC)

library(purrr,lib.loc=LIB_LOC)

raw=merged

raw$genes<-paste0(URL,gsub("chr","",raw$chr),'/',raw$LD_interval_start,'/',
                  raw$LD_interval_stop,"/360")
raw$num_of_genes<-NA

for(i in 1:nrow(raw)){
  raw_genes<-fromJSON(raw$genes[i])
  raw$num_of_genes[i]<-length(raw_genes)
  temp<-paste(map_chr(raw_genes,"symbol"),collapse = ",")
  raw$gene_symbol[i]<-temp
}

raw$genes<-NULL

write.xlsx(raw, paste0(exp_dir,"QTL_locuszoom/",expname,"_table3_final.xlsx"))
write.csv(raw, paste0(exp_dir,"QTL_locuszoom/",expname,"_table3_final.csv"),
          row.names=F,quote=T)