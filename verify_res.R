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


# sanity check list (clean) -------------------------------------------------------
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
tf.check<-tf.check%>% 
  left_join(dic$merged,by=c("ensembl_gene_id"="ensembl.id"))

write.csv(tf.check,"./test/snATAC_check_list_anno.csv",quote = T)
system("open ./test/snATAC_check_list_anno.csv")

# filter TF with no sf annotation 
table(is.na(tf.check$subfamily.id)) #22 vs. 2 
tf.check <- tf.check %>% filter(!is.na(subfamily.id))


# check results  ----------------------------------------------------------

tf.check.res <- left_join(tf.check,dat$motif%>%
                            filter(genus.id %in% tf.check$genus.id[!is.na(tf.check$genus.id)])%>%
                            dplyr::select(sample,genus.id,adj_p.value)%>%
                            mutate(adj_p.value=round(adj_p.value,2))%>%
                            group_by(sample,genus.id)%>%
                            spread(key=sample,value = adj_p.value))%>%
  arrange(genus.id)

write.csv(tf.check.res,"./test/snATAC_check_list_res_padj.csv",quote = T)
system("open ./test/snATAC_check_list_res_padj.csv")

