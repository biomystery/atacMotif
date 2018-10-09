source("./bin/aux_funs.R")
dat <- readRDS("./test/appDat.rds")


# by supfamily  ------------------------------------------------------
pd.sub <- dat$motif%>%
  #filter(log2FE>0)%>%
  #filter(adj_p.value > 10)%>%
  dplyr::select(sample,subfamily,log2FE)%>%
  group_by(sample,subfamily)%>%
  summarise(mlog2FE=max(log2FE))%>%
  spread(key=sample,value = mlog2FE)%>%
  as.data.frame()%>%
  column_to_rownames("subfamily")

# order by level 
rord <- pd.sub %>% 
  rownames_to_column("tf")%>%
  arrange(desc(alpha_1),desc(alpha_2),
          desc(beta_1),desc(beta_2),
          desc(delta_1),desc(delta_2),
          desc(endothelial_1),desc(endothelial_2),
          desc(exocrine),desc(gamma),
          desc(glial),desc(immune),desc(stellate))

pheatmap(pd.sub[rord$tf,],scale = "none",cluster_rows = F,cluster_cols = F,na_col = "grey",
         color = colorRampPalette(brewer.pal(9,"Blues"))(21)[7:21],
         show_rownames = F,fontsize_row = 6)


# sanity check list -------------------------------------------------------
dic<- readRDS("./db/dic_jaspar_tfclass.rds")
tf.check <- read.csv("./test/snATAC_check_list.csv",header = T,stringsAsFactors = F)
tf.check <- tf.check[,-c(4,7,9)]
tf.check<-tf.check%>%
  gather(key = "cellType",value = "TF")%>%
  filter(TF!="")%>%
  group_by(TF)%>%
  summarise(cellType=paste0(cellType,collapse = ";"))

# tf list 
tf.all <- as.vector(as.matrix(tf.check$TF,nrow=1))


# find ensemble id
require(biomaRt)
human<-  useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
name.dic <- getBM(attributes =c("ensembl_gene_id", "external_gene_name"),
                  filters = "external_gene_name",
                  tf.all,
                  mart = human)
tf.all[!tf.all %in% name.dic$external_gene_name]

# add ensemble id 
tf.check <- tf.check%>%
  left_join(name.dic,by = c("TF"="external_gene_name"))

# find annotation 
tf.check%>% 
  left_join(dat$tfclass[,c(1:3,)],by=c("ensembl_gene_id"="ensembl.id"))

tf.all.anno <- tf.all.anno[!is.na(tf.all.anno$external_gene_name),]
sum(res$tfclass$ensembl.id%in% name.dic$ensembl_gene_id )
sum(name.dic$ensembl_gene_id%in%  res$tfclass$ensembl.id)
sum(name.dic$external_gene_name%in%  res$tfclass$genus.name)
dup.id <- duplicated(tf.all.anno$ensembl_gene_id)
dup.tf <- tf.all.anno[dup.id,]
tf.all.anno <- tf.all.anno[!dup.id,]
rownames(tf.all.anno) <- tf.all.anno$external_gene_name

# for each celltype 
res <-lapply(colnames(tf.check),function(nm){
  cell.tfs <- tf.check[,nm]
  cell.tfs <- cell.tfs[cell.tfs!=""]
  cell.sfs <- tf.all.anno[cell.tfs,]
  cell.sfs <- cell.sfs[!is.na(cell.sfs$ensembl_gene_id),]
  cell.sfs <- cell.sfs[(paste(cell.sfs$subfamily.id,cell.sfs$subfamily.name,sep = "_") %in% rownames(pd.sub)),]
  data.frame(round(pd.sub[paste(cell.sfs$subfamily.id,cell.sfs$subfamily.name,sep = "_"),],2),type=nm)
})
res <- do.call(rbind,res)
write.csv(res,"check.csv",quote = F)
