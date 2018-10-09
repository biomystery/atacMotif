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
  if(T){
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
  system("open ./test/snATAC_check_list_anno.csv")}

tf.check <- read.csv("./test/snATAC_check_list_anno.csv",header = T,stringsAsFactors = F)

# filter TF with no sf annotation 
table(is.na(tf.check$subfamily.id)) #22 vs. 2 
#tf.check <- tf.check %>% filter(!is.na(subfamily.id))


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

# extra check 
dat$motif%>%
  filter(genus.id=="4.1.1.6.2")%>%
  dplyr::select(sample,genus.id,adj_p.value)%>%
  mutate(adj_p.value=round(adj_p.value,2))%>%
  group_by(sample,genus.id)%>%
  spread(key=sample,value = adj_p.value)



# check SF rank  ----------------------------------------------------------
if(T){pd.sub.padj <- dat$motif%>%
  #filter(log2FE>0)%>%
  #filter(adj_p.value > 10)%>%
  dplyr::select(sample,subfamily,adj_p.value )%>%
  group_by(sample,subfamily)%>%
  summarise(madj_p.value =max(adj_p.value ))%>%
  spread(key=sample,value = madj_p.value)%>%
  as.data.frame()%>%
  column_to_rownames("subfamily")


rord <- pd.sub.padj %>% 
  rownames_to_column("tf")%>%
  arrange(desc(alpha_1),desc(alpha_2),
          desc(beta_1),desc(beta_2),
          desc(delta_1),desc(delta_2),
          desc(endothelial_1),desc(endothelial_2),
          desc(exocrine),desc(gamma),
          desc(glial),desc(immune),desc(stellate))
}

# sf check 
sf.check <- tf.check %>% 
  filter(!is.na(subfamily.id))%>%
  group_by(subfamily.id,subfamily.name,family.name)%>%
  summarise(cellType=paste(cellType,collapse = ";"))%>%
  data.frame(stringsAsFactors = F)
sf.check[17,2] <- as.character(sf.check[17,3])



# check at sf level --------------------------------------------------------

pd.sub.padj <- pd.sub.padj[rord$tf,]

tar.sf <- "1.1.3.1"
anno_row <- data.frame(sf= factor(rep(0,nrow(pd.sub)),levels = c(0,1)))
rownames(anno_row)<-rownames(pd.sub)

anno_row$sf[grepl(tar.sf,rownames(anno_row))] <- 1
ann_colors <- list(sf=c("0"="white", "1"="firebrick"))
if(T){
  setEPS()
  i<-1
  postscript(paste0("sf_hm_",i,".eps"),onefile = F,width = 2,height = 6)
  pheatmap(pd.sub,scale = "none",cluster_rows = F,
           cluster_cols = F,na_col = "grey",cellwidth = 8,
           color = colorRampPalette(brewer.pal(9,"Blues"))(21)[7:21],
           show_rownames = F,show_colnames = T,fontsize   = 8,annotation_row = anno_row,
           annotation_colors = ann_colors,annotation_legend = F,annotation_names_row = F)
  dev.off()
}

# check sanity list for cell type & get ranks 
sf.check.celltype<-(sf.check%>%
  filter(subfamily.id==tar.sf))$cellType
sf.check.celltype<- unique(unlist(strsplit(sf.check.celltype,split=";")))
sf.check.motif <- pd.sub[grepl(tar.sf,rownames(anno_row)),]

pd.sub.rank <- round(apply(-pd.sub.padj,2,rank))
pd.sub.rank.tot <- apply(pd.sub.padj,2,function(x) sum(!is.na(x)))
pd.sub.rank.anno <- pd.sub.rank
for(nm in colnames(pd.sub.rank))   pd.sub.rank.anno[,nm] <- sapply(pd.sub.rank[,nm],function(x) paste(x,pd.sub.rank.tot[nm],sep = "/"))
pd.sub.rank.anno[is.na(pd.sub)] <- NA


# save sf.check rank 
pd.sub.rank.anno.check  <- pd.sub.rank.anno[unlist(sapply(sf.check$subfamily.id, function(x) grep(x,rownames(pd.sub.padj)))),]%>%
  as.data.frame()%>%
  rownames_to_column("sf")%>%
  separate(sf,into = c("sf.id","sf.nm"),sep = "_")%>%
  left_join(sf.check,by = c("sf.id"="subfamily.id"))%>%
  arrange(sf.id)
write.csv(pd.sub.rank.anno.check,file = "./test/pd.sub.rank.anno.check.csv")
system( "open ./test/pd.sub.rank.anno.check.csv")


pd.sub.rank.anno.check  <- round(pd.sub.padj[unlist(sapply(sf.check$subfamily.id, function(x) grep(x,rownames(pd.sub.padj)))),],2)%>%
  as.data.frame()%>%
  rownames_to_column("sf")%>%
  separate(sf,into = c("sf.id","sf.nm"),sep = "_")%>%
  left_join(sf.check,by = c("sf.id"="subfamily.id"))%>%
  arrange(sf.id)
write.csv(pd.sub.rank.anno.check,file = "./test/pd.sub.padj.anno.check.csv")
system( "open ./test/pd.sub.padj.anno.check.csv")
#  save sf.check padj
pd.sub.padj


# expression 
sf.tfs <- dic$merged%>% 
  filter(subfamily.id==tar.sf)%>%
  filter(grepl("ENSG",ensembl.id ))
pd.tfs <- dat$promoter.cpm%>%
  right_join(sf.tfs[,c("ensembl.id","genus.id","tf.symbol")],by = c("ensembl_gene_id"="ensembl.id"))%>%
  arrange(genus.id)%>%
  column_to_rownames("tf.symbol")
if(T){
  setEPS()
  postscript(paste0("sf_exp_hm_",i,".eps"),onefile = F,width = 4,height = 4)
  pheatmap(pd.tfs[,3:16],scale = "row",cluster_rows = F,cluster_cols = F,
           cellwidth = 8,fontsize   = 8,
           main = paste(unique(sf.tfs$subfamily.id),
                        "SF:",unique(sf.tfs$subfamily.name),"F:",unique(sf.tfs$family.name)))
  dev.off()
}
