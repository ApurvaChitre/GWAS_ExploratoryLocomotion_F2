# file constants ----

TRAITS_DATA <- paste0("/projects/ps-palmer/apurva/u01_huda_akil/residuals",
                      "/u01_huda_akil_traits.RData")
RESIDUALS_DATA <- paste0("/projects/ps-palmer/apurva/u01_huda_akil/residuals",
                         "/residuals_u01_huda_akil.RData")
LIB_LOC <- "/home/aschitre/R_libs/"
TOPSNPS_CSV <- paste0("/projects/ps-palmer/apurva/u01_huda_akil/QTL_locuszoom",
                      "/u01_huda_akil_final_QTL_list.csv")
INPUT_DIR <- paste0("/oasis/tscc/scratch/aschitre/round8/unpruned",
                    "/snp_annotations/snpeff_input/")
OASIS_LOC <- "/oasis/tscc/scratch/aschitre/software/plink-1.90/plink"
LOCAL_LOC <- paste0("/projects/ps-palmer/apurva/u01_huda_akil/genotypes",
                    "/u01_huda_akil_genotypes")
TEMP_DIR <- "/projects/ps-palmer/apurva/u01_huda_akil/r2/"
RESULTS_WD <- "/projects/ps-palmer/apurva/u01_huda_akil/results"
OUTPUT_DIR <- paste0("/oasis/tscc/scratch/aschitre/round8/unpruned",
                     "/snp_annotations/snpeff_output/")
SNPEFF_DIR <- "/projects/ps-palmer/software/local/src/snpEff/"
JAVA_LOC <- "/usr/lib/jvm/jre-1.8.0/bin/java"

# end ----

args <- commandArgs(trailingOnly = TRUE)
exp_dir <- args[1]
expname <- args[2]

exp_dir<-paste0(exp_dir,"/")


load(TRAITS_DATA)
load(RESIDUALS_DATA)
load(paste0(exp_dir,"residuals/",expname,"_traits.RData"))

.libPaths(c(LIB_LOC, .libPaths()))

library("data.table", lib.loc=LIB_LOC)

library("BSgenome.Rnorvegicus.UCSC.rn6",lib.loc=LIB_LOC)

#VEP's region REST endoint requires variants are described as 
#[chr]:[start]-[end]:[strand]/[allele]. 
#This follows the same conventions as the default input format described above, 
#with the key difference being that this format does not require the reference 
#(REF) allele to be included; VEP will look up the reference allele using either 
#a provided FASTA file (preferred) or Ensembl core database. Strand is optional 
#and defaults to 1 (forward strand).

topsnps=read.csv(TOPSNPS_CSV,header=T,stringsAsFactors=F)

topsnps$chr<-sapply(strsplit(topsnps$topsnp,split=":"),`[`,1)
topsnps$pos<-sapply(strsplit(topsnps$topsnp,split=":"),`[`,2)

qtls<-topsnps

#create directory
dir.create(paste0(INPUT_DIR,expname), showWarnings = FALSE)

#LD file

LD_ref_snps<-function(x,trait,chr,pos){
  data<-fread(x,header=T,stringsAsFactors =F,select=c(6,7))
  column_names_ld<-c("snp","rsquare")
  setnames(data, column_names_ld)
  data<-data[which(data$rsquare>=0.6),]
  #get reference allele
  data$reference_allele<-NA
  for(j in 1:nrow(data)){
    chr_ref<-strsplit(data$snp[j],split=":")[[1]][1]
    pos_ref<-as.numeric(strsplit(data$snp[j],split=":")[[1]][2])
    data$reference_allele[j]<-as.character(
      getSeq(BSgenome.Rnorvegicus.UCSC.rn6, chr_ref,start=pos_ref,end=pos_ref)) 
  }
  
  return(data)  
}

results_dir = paste0(exp_dir,"results/")
readAssoc<-function(x){
  data<-fread(paste0('cat ',x,' | cut -f 2,4,5 '),header=T,stringsAsFactors =F)
  column_names<-c("snp","allele1","allele2")
  setnames(data, column_names)
  return(data)
}

for(i in 1:nrow(qtls)){
  trait<-qtls$trait[i]  
  topsnp<-qtls$topsnp[i]
  chr<-qtls$chr[i]
  pos<-qtls$pos[i]
  assoc_name = paste0(trait,".loco.mlma")
  
  system(paste0(OASIS_LOC," --bfile ",LOCAL_LOC," --chr ",gsub("chr","",chr),
                " --nonfounders --r2  --ld-snp ",topsnp," ",
                "--ld-window 100000 -ld-window-kb 6000 --ld-window-r2 0 --out ",
                TEMP_DIR,expname,"_temp_qtl_n"),
         ignore.stdout = T,ignore.stderr = T,wait = T)
  
  LD_df<-LD_ref_snps(paste0(TEMP_DIR,expname,"_temp_qtl_n.ld"),
                     trait = trait,chr=chr,pos=pos)
  
  setwd(RESULTS_WD)
  assoc_df<- readAssoc(assoc_name)
  
  merged<-merge(LD_df,assoc_df,by="snp",all=F)
  
  merged$vep_allele<-NA
  
  for(k in 1:nrow(merged)){
    chr_ref<-strsplit(merged$snp[k],split=":")[[1]][1]
    chr_ref=paste0("chr",chr_ref)
    pos_ref<-as.numeric(strsplit(merged$snp[k],split=":")[[1]][2])
    
    if(merged$reference_allele[k] %in% merged$allele1[k]){
      merged$vep_allele[k]<-merged$allele2[k]
    }else{
      merged$vep_allele[k]<-merged$allele1[k]  
    }
  }
  
  ##CHROM  POS ID      REF ALT QUAL    FILTER  INFO
  setwd(paste0(INPUT_DIR,expname))
  colnames(merged)[grepl("vep_allele",colnames(merged))]<-"ALT"
  
  merged$CHROM<-NA
  merged$POS<-NA
  merged$CHROM<-sapply(strsplit(merged$snp,split=":"),`[`,1)
  merged$POS<-sapply(strsplit(merged$snp,split=":"),`[`,2)
  merged$ID<-paste0(merged$CHR,"_",merged$POS)
  
  colnames(merged)[grepl("reference_allele",colnames(merged))]<-"REF"
  
  colnames(merged)[grepl("CHROM",colnames(merged))]<-"##CHROM"
  
  cat("##fileformat=VCFv4.0 \n",file=paste0(trait,"_",chr,"_",pos,".vcf"))
  
  merged$QUAL<-"."
  merged$FILTER<-"."
  merged$INFO<-"."
  
  write.table(
    merged[,c("##CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")],
    paste0(trait,"_",chr,"_",pos,".vcf"),
    col.names=T,append=T,row.names=F,quote=F,sep="\t",na = ".")
}

#run snpeff for one file
dir.create(paste0(OUTPUT_DIR,expname), showWarnings = FALSE)

qtls$vcf_snpeff<-paste0(INPUT_DIR,expname,"/",qtls$trait,
                        "_",qtls$chr,"_",qtls$pos,".vcf")

if (file.exists(paste0(SNPEFF_DIR,expname,".sh"))) {
  file.remove(paste0(SNPEFF_DIR,expname,".sh"))
}

for(i in 1:nrow(qtls)){
  input<-paste0(INPUT_DIR,expname,"/",qtls$trait[i],"_",
                qtls$chr[i],"_",qtls$pos[i],".vcf")
  part1<-paste0(JAVA_LOC," -Xms40g -Xmx40g -jar ",SNPEFF_DIR,
                "snpEff.jar Rnor_6.0.99  -no-intergenic -no-intron -noStats ")
  part2<-" > "
  output<-paste0(OUTPUT_DIR,expname,"/",qtls$trait[i],"_",
                 qtls$chr[i],"_",qtls$pos[i],"_anno.vcf")
  
  command<- paste0(part1,input,part2,output)
  
  write(command,file=paste0(SNPEFF_DIR,expname,".sh"),append=TRUE)  
}

system(paste0("chmod +x ",SNPEFF_DIR,expname,".sh"))
system(paste0(SNPEFF_DIR,expname,".sh"),wait=T)

#run .sh command on an interactive node

if (file.exists(paste0(SNPEFF_DIR,expname,"_filter.sh"))) {
  file.remove(paste0(SNPEFF_DIR,expname,"_filter.sh"))
}

for(i in 1:nrow(qtls)){
  part1<-paste0("cat ",OUTPUT_DIR,expname,"/",qtls$trait[i],"_",
                qtls$chr[i],"_",qtls$pos[i],"_anno.vcf \\")
  part2<-paste0("    | ",SNPEFF_DIR,"scripts/vcfEffOnePerLine.pl \\")
  part3<-paste0("| ",JAVA_LOC," -Xms512m -Xmx512m -jar ",SNPEFF_DIR,
                "SnpSift.jar extractFields - CHROM POS REF ALT ",
                "\"ANN[*].ALLELE\" \"ANN[*].ANNOTATION\" \"ANN[*].IMPACT\" ",
                "\"ANN[*].GENE\" \"ANN[*].GENEID\" \"ANN[*].FEATURE\" \"ANN[*]",
                ".FEATUREID\" \"ANN[*].BIOTYPE\" \"ANN[*].RANK\" \"ANN[*].",
                "HGVS_C\" \"ANN[*].HGVS_P\" \"ANN[*].CDNA_POS\" \"ANN[*].",
                "CDNA_LEN\" \"ANN[*].CDS_POS\" \"ANN[*].CDS_LEN\" \"ANN[*].",
                "AA_POS\" \"ANN[*].AA_LEN\" \"ANN[*].DISTANCE\" \"ANN[*].",
                "ERRORS\" >")
  output<- paste0(OUTPUT_DIR,expname,"/",qtls$trait[i],"_",
                  qtls$chr[i],"_",qtls$pos[i],"_fields.txt")
  output<-paste0(part3,output)
  command<-paste(part1,part2,output,sep = "\n")
  write(command,file=paste0(SNPEFF_DIR,expname,"_filter.sh"),append=TRUE)  
}

system(paste0("chmod +x ",SNPEFF_DIR,expname,"_filter.sh"))
system(paste0(SNPEFF_DIR,expname,"_filter.sh"),wait=T)

.libPaths(c(LIB_LOC, .libPaths()))
library("data.table", lib.loc=LIB_LOC)

setwd(paste0(OUTPUT_DIR,expname))

anno_files<-list.files(path=".",full.names = F,pattern=paste0("*_fields.txt_*"))

readdata<-function(x){
  data<-fread(file=x,header=T,stringsAsFactors =F,na.strings = "")
  column_names=c('chr','pos','ref','alt','snp_allele','snp_effect','snp_impact',
                 'snp_gene','snp_geneid','snp_feature','snp_featureid',
                 'snp_biotype','snp_rank','snp_hgvs_c','snp_hgvs_p',
                 'snp_cdna_pos','snp_cdna_length','snp_cds_pos','snp_cds_len',
                 'snp_aa_pos','snp_aa_length','snp_dist_to_feature',
                 'snp_errors')
  setnames(data, column_names)
  data=data[snp_impact %in% c("HIGH","MODERATE")]
  trait_snp<-gsub("_fields.txt","",x)
  data$trait_snp<-trait_snp
  return(data)
}

anno_snpeff<-lapply(anno_files,readdata)
anno_files<-gsub("_fields.txt","",anno_files)
names(anno_snpeff)<-anno_files

annotations<-do.call("rbind",anno_snpeff)

write.csv(annotations,
          paste0(exp_dir,"QTL_locuszoom/",expname,
                 "_snp_annotations_final_withWarnings.csv"),
          row.names=F,quote=T)



# file constants ----

LIB_LOC <- "/home/aschitre/R_libs/"
URL <- "https://rgd.mcw.edu/rgdweb/report/gene/main.html?id="
OASIS_LOC <- "/oasis/tscc/scratch/aschitre/software/plink-1.90/plink"
LOCAL_LOC <- paste0("/projects/ps-palmer/apurva/u01_huda_akil/genotypes",
                    "/u01_huda_akil_genotypes")
TEMP_DIR <- "/projects/ps-palmer/apurva/u01_huda_akil/r2/"

# end ----

args <- commandArgs(trailingOnly = TRUE)
exp_dir <- args[1]
expname <- args[2]

exp_dir<-paste0(exp_dir,"/")

library("data.table", lib.loc=LIB_LOC)

load(paste0(exp_dir,"residuals/",expname,"_traits.RData"))

annotations<-read.csv(paste0(exp_dir,"QTL_locuszoom/",expname,
                             "_snp_annotations_final_withWarnings.csv"),
                      header=T,stringsAsFactors = F)
annotations<-annotations[is.na(annotations$snp_errors),]

annotations$gene_name<-NA
annotations$RGD_link<-NA
annotations$gene_long_name<-NA

ensembl_genes<-unique(annotations$snp_geneid)
ensembl_genes<-ensembl_genes[!is.na(ensembl_genes)]

.libPaths(c(LIB_LOC, .libPaths()))

library(vctrs, lib.loc=LIB_LOC)
library(annotables, lib.loc=LIB_LOC)

if(nrow(rnor6[which(rnor6$ensgene %in% ensembl_genes),])>0){
  
  for(i in 1:length(ensembl_genes))  {
    if(ensembl_genes[i] %in% rnor6$ensgene){
      annotations$gene_name[
        which(annotations$snp_geneid %in% ensembl_genes[i])]<-
        as.character(rnor6[which(rnor6$ensgene %in% ensembl_genes[i]),"symbol"])
      acc<-gsub("]","",strsplit(as.character(rnor6[
        which(rnor6$ensgene %in% ensembl_genes[i]),"description"]),
        split="Acc:")[[1]][2])
      if(is.na(acc)){
        annotations$RGD_link[
          which(annotations$snp_geneid %in% ensembl_genes[i])]<-NA
        annotations$gene_long_name[
          which(annotations$snp_geneid %in% ensembl_genes[i])]<-NA
      }else{
        annotations$RGD_link[
          which(annotations$snp_geneid %in% ensembl_genes[i])]<-paste0(URL,acc)  
        annotations$gene_long_name[
          which(annotations$snp_geneid %in% ensembl_genes[i])]<-
          strsplit(as.character(rnor6[
            which(rnor6$ensgene %in% ensembl_genes[i]),"description"]),
            split=" \\[")[[1]][1]
      }
    }
  }
}

write.csv(annotations,
          paste0(exp_dir,"QTL_locuszoom/",expname,"_snp_annotations_final.csv"),
          row.names=F,quote=T)
annotations$r2_with_trait_topsnp<-NA
annotations$dprime_with_trait_topsnp<-NA

for(i in 1:nrow(annotations)){
  topsnp<-gsub(".*_chr","",annotations$trait_snp[i]) 
  chr<-strsplit(topsnp,split="_")[[1]][1]
  trait_topsnp<-paste0("chr",chr,":",strsplit(topsnp,split="_")[[1]][2])
  snp<-paste0(annotations$chr[i],":",annotations$pos[i])
  command<-paste0(OASIS_LOC," --bfile ",LOCAL_LOC," --chr ",chr,
                  " --nonfounders --r2  --ld ",trait_topsnp," ",snp," --out ",
                  TEMP_DIR,expname,"_temp_snpeff")
  system(command,wait=T)
  r2_log=fread(paste0("cat ",TEMP_DIR,expname,"_temp_snpeff.log | grep R-sq"),
               header=F,stringsAsFactors = F)
  if(nrow(r2_log)<2){
    annotations$r2_with_trait_topsnp[i]<-r2_log$V3
    annotations$dprime_with_trait_topsnp[i]<-r2_log$V6
  }else{
    annotations$r2_with_trait_topsnp[i]<-
      as.numeric(r2_log[which.max(r2_log$V3),"V3"])
    annotations$dprime_with_trait_topsnp[i]<-
      as.numeric(r2_log[which.max(r2_log$V3),"V6"])
  }
}

annotations$r2_with_trait_topsnp<-unlist(annotations$r2_with_trait_topsnp)
annotations$dprime_with_trait_topsnp<-
  unlist(annotations$dprime_with_trait_topsnp)
write.csv(annotations,
          paste0(exp_dir,"QTL_locuszoom/",expname,"_snp_annotations_final.csv"),
          row.names=F,quote=T)