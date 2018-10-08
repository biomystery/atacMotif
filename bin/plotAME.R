
source("./bin/aux_funs.R")
tfclass <- readRDS("./db/tfclass.rds")


# input parameters --------------------------------------------------------
ame_path <- "./test/ame_2kbg_all.res.txt"
th <- .00001

# ame’s output ------------------------------------------------------------
ame_res <- read.table(ame_path,
                      header = T,
                      stringsAsFactors = F)


# load jaspar_to_TFclass dic---------------------------------------------------------

# filter adj_p.value <1e-5
ame_res.anno <- ame_res%>% 
  filter(adj_p.value<=th)

# jaspa dic 
dic$jasparTOensembl <- unique(ame_res.anno$motif_ID)
names(dic$jasparTOensembl) <- unique(ame_res.anno$motif_ID)

if(F){
  # grab jaspar info from web 
  for(i in 1:length(dic.jasparTOensembl)) {
    print(i);
    dic.jasparTOensembl[i]=getEnsemblFromJaspar(dic.jasparTOensembl[i])}
  
  # tfclass dic 
  dic.ensemblTOsubfamily <- tfclass$merge%>%
    separate_rows(tf.id,sep = ";") %>%
    distinct(tf.id,.keep_all = T)%>%
    column_to_rownames("tf.id")
  
  # jaspa to tfclass 
  dic.jasparTOtfclass <-data.frame(ensembl.id = dic.jasparTOensembl)%>%
    rownames_to_column("motif_ID")%>%
    separate_rows(ensembl.id,sep = ";") %>%
    right_join(dic.ensemblTOsubfamily%>%rownames_to_column("ensembl.id"))%>%
    filter(!is.na(motif_ID))
  
  saveRDS(list(jasparTOensembl=dic.jasparTOensembl,
               ensemblTOtfclass=dic.ensemblTOsubfamily,
               jasparTOtfclass=dic.jasparTOtfclass),
          file = "./db/dic_jaspar_tfclass.rds")
}

dic<- readRDS("./db/dic_jaspar_tfclass.rds")
names(dic)

# add annotation -----------------------------------------------------------------

tmp <- sapply(unique(ame_res[ame_res$adj_p.value<=th,"motif_ID"]),function(a)
  (dic$jasparTOtfclass%>%
     filter(row_number()==grep(a,dic$jasparTOtfclass$jaspar.id)[1])),simplify = F)
tmp<-do.call(rbind,tmp)

# 17 out of 430 are not in the dic. 
all(unique(ame_res$motif_ID) %in% tmp$jaspar.id)
sum(!unique(ame_res$motif_ID) %in% tmp$jaspar.id)
length(unique(ame_res$motif_ID))
unique(ame_res$motif_ID) [!unique(ame_res$motif_ID) %in% tmp$jaspar.id]

pd.ame_res <- ame_res%>%
  filter(adj_p.value<=th)
pd.ame_res<-cbind(pd.ame_res,tmp[pd.ame_res$motif_ID,-1])%>%
  #unite(tf,1:2)%>%
  mutate(tf=motif_alt_ID)%>%
  mutate(subfamily.name=ifelse(is.na(subfamily.name),as.character(family.name),subfamily.name))%>%
  unite(family,family.id,family.name)%>%
  unite(subfamily,subfamily.id,subfamily.name)%>%
  mutate(p.value=-log10(p.value),
         adj_p.value=-log10(adj_p.value),
         E.value=-log10(E.value),
         log2FE = log2((X.TP+0.1)/(X.FP+0.1)))


res.sum <- data.frame(pd.ame_res[,c("sample","FASTA_max")]%>%distinct(),
                      nmotifs=sapply(unique(pd.ame_res$sample), 
                                     function(x) nrow(pd.ame_res%>% filter(sample==x))),
                      nmotifs.FE=sapply(unique(pd.ame_res$sample), 
                                     function(x) nrow(pd.ame_res%>%
                                                        filter(log2FE>1)%>%
                                                        filter(sample==x)))
                      )
res.sum%>% mutate(frac=round(nmotifs/FASTA_max,2))

write.table(data.frame(pd.ame_res[,c("sample","FASTA_max")]%>%distinct(),
                         nmotifs=sapply(unique(pd.ame_res$sample), function(x) nrow(pd.ame_res%>% filter(sample==x)))),
            file = "n.csv",sep = ',',quote = F,row.names = F)

ggplot(pd.ame_res%>%filter(log2FE>1),aes(sample,tf))+
  geom_tile(aes(fill=log2FE))+
  theme(axis.text.y = element_blank())
  #+scale_fill_gradientn(colours =  wes_palette("Zissou1", 21, type = "continuous"))


# thresholding ------------------------------------------------------------
pd <- pd.ame_res%>%
  filter(log2FE>1)%>%
  filter(adj_p.value > 10)%>%
  dplyr::select(sample,tf,log2FE)%>%
  spread(key=sample,value = log2FE)%>%
  column_to_rownames("tf")

pd.2 <- pd; pd.2[is.na(pd)]<-0

# hcluster ----------------------------------------------------------------
require(cluster)
distfunc <- function(x) daisy(x,metric="gower")
d <- distfunc(pd.2)
dend <- as.dendrogram(hclust(d))
dend <- as.dendrogram(hclust(dist(pd.2)))
plot(dend)
pheatmap(pd[order.dendrogram(dend),],scale = "none",cluster_rows = F,cluster_cols = F,na_col = "grey",
         color = colorRampPalette(brewer.pal(9,"Blues"))(21),
         show_rownames = F)



# ranking -----------------------------------------------------------------
# order by column
rord <- pd %>% 
  rownames_to_column("tf")%>%
  arrange(desc(alpha_1),desc(alpha_2),
          desc(beta_1),desc(beta_2),
          desc(delta_1),desc(delta_2),
          desc(endothelial_1),desc(endothelial_2),
          desc(exocrine),desc(gamma),
          desc(glial),desc(immune),desc(stellate))

pheatmap(pd[rord$tf,],scale = "none",cluster_rows = F,cluster_cols = F,na_col = "grey",
         color = colorRampPalette(brewer.pal(9,"Blues"))(21)[7:21],
         show_rownames = F,fontsize_row = 6)

require(heatmaply)
heatmaply(pd[rord$tf,],scale="none",
             Rowv = NULL,Colv = NULL,
             colors = colorRampPalette(brewer.pal(9,"Blues"))(21)[7:21])

# customized ranking ------------------------------------------------------------------
custRank_all <- function(input.pd=pd.2,na.val=0){
  custRank <- function(pd.2,v="alpha_1"){
    y <- pd.2[,v] 
    idx.nna <- which(y!=na.val); idx.na <- which(y==na.val)
    idx.idx <- order(apply(pd.2[idx.nna,],1,mean),decreasing = T)  
    res.list <- list()
    res.list$ordered <- pd.2[idx.nna[idx.idx],]
    res.list$un_ordered <- pd.2[idx.na,]
    res.list
  }
  
  tmp.pd <- input.pd;final.pd <- input.pd[-(1:nrow(input.pd)),]
  for(x in colnames(input.pd)){
    tmp.res <- custRank(pd.2 = tmp.pd,v=x)  
    tmp.pd <- tmp.res$un_ordered
    final.pd<- rbind(final.pd,tmp.res$ordered)
  }
  final.pd[-1,]
}

final.pd <- custRank_all()
pheatmap(pd[rownames(final.pd),],scale = "none",cluster_rows = F,cluster_cols = F,na_col = "grey",
         color = colorRampPalette(brewer.pal(9,"Blues"))(21)[7:21],
         show_rownames = F,fontsize_row = 6)

# digital customized ranking ----------------------------------------------
bks <- seq(min(pd, na.rm = T), max(pd,na.rm = T), length.out = 15 + 1)
mat = as.matrix(pd)
pd.3 <- matrix(as.numeric(cut(as.vector(mat), breaks = bks, include.lowest = T)),
               nrow(mat), ncol(mat), dimnames = list(rownames(mat), colnames(mat)))
pd.3[is.na(pd.3)]<- 0
final.pd <- custRank_all(input.pd = pd.3)

pheatmap(pd[rownames(final.pd),],scale = "none",cluster_rows = F,cluster_cols = F,na_col = "grey",
         color = colorRampPalette(brewer.pal(9,"Blues"))(21)[7:21],
         show_rownames = F,fontsize_row = 6)

# kmeans ----------------------------------------------
getClusts <- function(pd.2,...){
  require(NbClust)
  nb <- NbClust(pd.2,method="complete",...)
  require(factoextra)
  fviz_nbclust(nb) + theme_minimal()
  ords <-(sapply(1:max(nb$Best.partition), function(x) which(nb$Best.partition==x)))
  brks <- cumsum(lapply(ords, length))
  ords <- unlist(ords)
  list(ords=ords,brks=brks)  
}

pheatmap(pd[ords,],scale = "none",cluster_rows = F,cluster_cols = F,na_col = "grey",
         color = colorRampPalette(brewer.pal(9,"Blues"))(21)[7:21],gaps_row = brks,
         show_rownames = F,fontsize_row = 6)

heatmaply(pd[ords,],scale="none",
          Rowv = NULL,Colv = NULL,gaps_row=brks,
          colors = colorRampPalette(brewer.pal(9,"Blues"))(21)[7:21])

# row-wise scale  ---------------------------------------------------------
pd.3 <- t(scale(t(pd.2)))
pd.3[pd.3>1.96] <- 1.96; pd.3[pd.3 < -1.96] <- -1.96
d <- distfunc(pd.3)
dend <- as.dendrogram(hclust(d))

pd.3[is.na(pd)]<- NA
pheatmap(pd.3[order.dendrogram(dend),],scale = "none",cluster_rows = F,cluster_cols = T,na_col = "grey",
         show_rownames = F,color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                     "RdYlBu")))(21))

pd.3[is.na(pd)]<- -1.96


# plot by supfamily  ------------------------------------------------------
pd.sub <- pd.ame_res%>%
  filter(log2FE>0)%>%
  filter(adj_p.value > 10)%>%
  dplyr::select(sample,subfamily,log2FE)%>%
  group_by(sample,subfamily)%>%
  summarise(mlog2FE=max(log2FE))%>%
  spread(key=sample,value = mlog2FE)%>%
  as.data.frame()%>%
  column_to_rownames("subfamily")

# hclust
pd.sub <- pd.sub[apply(pd.sub, 1, function(x) sum(is.na(x))<ncol(pd.sub)),]
pd.sub.2 <- pd.sub;pd.sub.2[is.na(pd.sub.2)] <- 0   
dend <- as.dendrogram(hclust(dist(pd.sub.2)))
pheatmap(pd.sub[order.dendrogram(dend),],scale = "none",cluster_rows = F,cluster_cols = T,na_col = "grey",
         color = colorRampPalette(brewer.pal(9,"GnBu"))(11),
         show_rownames = T,fontsize_row = 6)


# nbClust
ords <- getClusts(pd.sub.2,min.nc=5)
pheatmap(pd.sub[ords$ords,],scale = "none",cluster_rows = F,cluster_cols = F,na_col = "grey",
         color = colorRampPalette(brewer.pal(9,"Blues"))(21)[7:21],gaps_row = ords$brks,
         show_rownames = F,fontsize_row = 6)

# cust
final.pd <- custRank_all(input.pd = pd.sub.2)

heatmaply(pd.sub[rownames(final.pd),],scale="none",
          Rowv = NULL,Colv = NULL,gaps_row=brks,
          colors = colorRampPalette(brewer.pal(9,"Blues"))(21)[7:21])

# by family  -----------------------------------------------------------
pd.fa <- pd.ame_res%>%
  filter(log2FE>1)%>%
  filter(adj_p.value > 10)%>%
  dplyr::select(sample,family,log2FE)%>%
  group_by(sample,family)%>%
  summarise(mlog2FE=max(log2FE))%>%
  spread(key=sample,value = mlog2FE)%>%
  as.data.frame()%>%
  column_to_rownames("family")

pd.fa <- pd.fa[apply(pd.fa, 1, function(x) sum(is.na(x))<ncol(pd.fa)),]
pd.fa.2 <- pd.fa;pd.fa.2[is.na(pd.fa)] <- 0   
dend <- as.dendrogram(hclust(dist(pd.fa.2)))
pheatmap(pd.fa[order.dendrogram(dend),],scale = "none",cluster_rows = F,cluster_cols = T,na_col = "grey",
         color = colorRampPalette(brewer.pal(9,"GnBu"))(11),
         show_rownames = T,fontsize_row = 8)
na.value <- -ncol(pd.fa)
pd.fa.2 <- pd.fa;pd.fa.2[is.na(pd.fa)] <- na.value
final.pd <- custRank_all(input.pd = pd.fa.2,na.val = na.value)
final.pd <- final.pd[-grep("NA",rownames(final.pd)),]
heatmaply(pd.fa[rownames(final.pd),],scale="none",
          Rowv = NULL,Colv = NULL,gaps_row=brks,
          colors = colorRampPalette(brewer.pal(9,"Blues"))(21)[7:21])

# choose a quantity to plot  -------------------------------------------------------------------
col.max <- apply(pd.ame_res[,grepl("sample",colnames(pd.ame_res))],
                 2,max)
col.max.col <- cut(col.max,seq(min(col.max)-.01,max(col.max)+.01,length.out = 20))
col.map <- colorRampPalette(brewer.pal(9,"Blues"))(19)
names(col.map)<- levels(col.max.col)
col.map.col <- col.map[col.max.col]
names(col.map.col) <- colnames(pd.ame_res)[grepl("sample",colnames(pd.ame_res))]

normalise <- function(x) (x-min(x))/(max(x)-min(x))
pd <- apply(pd.ame_res[,grepl("sample",colnames(pd.ame_res))],2,
            normalise)

pd <- pd.ame_res[grepl("ENSG",pd.ame_res$ensembl.id),
                   grepl("sample",colnames(pd.ame_res))]
pd.ame_res.hsap <- pd.ame_res[grepl("ENSG",pd.ame_res$ensembl.id),]

pd.scale <- scale(pd)
require(heatmaply)
require(RColorBrewer)
pd.2 <-apply(pd,2,normalise) %>% 
  as.data.frame %>%
  arrange(sample_1,sample_2,sample_3,sample_4,sample_5,sample_6,sample_7,sample_8,sample_9,sample_10,sample_11,sample_12,sample_13,sample_14,sample_15)
rownames(pd.2)<- rownames(pd)
heatmaply(apply(pd,2,normalise),
          seriate="mean",
          scale = "none",  
          colors=rev(brewer.pal(11,"RdYlBu")),
          k_col = 5,k_row = 15,
          RowSideColors = pd.ame_res[grepl("ENSG",pd.ame_res$ensembl.id),
                                       c("subfamily","family")],
          row_side_colors = brewer.pal(9,"Set1"),
          #ColSideColors = col.max.col,
          #col_side_palette = (col.map.col),
          column_text_angle = 90,
          key.title = "enrichment",margins = c(80,0,0,0)
          ) %>% layout(width=500,height=800)
pd.dist <- as.matrix(dist(t(apply(pd,2,normalise)),diag=T,upper=T))
pd.dist <-1- (pd.dist-min(pd.dist))/(max(pd.dist)-min(pd.dist))
heatmaply(as.matrix(pd.2),Rowv = NULL,scale = 'none',margins = c(80,0,NA,0),
          colors=rev(brewer.pal(11,"RdYlBu")))







saveRDS(list(pd.chu_ame=pd.chu_ame,
             pd.chu_ame.hsap=pd.chu_ame.hsap,
             pd.chu_ame.mmus=pd.chu_ame.mmus,
             pd.ame_res=pd.ame_res),file = './data/ame.res.rds')





# TOP5 TFs  ----------------------------------------------
pd.ame_res.top5<- pd.ame_res.hsap%>% 
  as.tibble()%>% 
  gather(key="sample",value = "log10_adj_p.value",1:15)%>%
  group_by(sample)%>%
  top_n(5,log10_adj_p.value)%>%
  spread(sample,log10_adj_p.value)

# print
pd.ame_res.top5<-pd.ame_res.top5[,c(colnames(pd.ame_res.top5)[1:7],paste0("sample_",1:15))]%>% 
  arrange(sample_1,sample_2,sample_3,sample_4,sample_5,sample_6,sample_7,sample_8,sample_9,sample_10,sample_11,sample_12,sample_13,sample_14,sample_15)%>%
  print(n=28)


# ensemble -> gene symbol 
require(biomaRt)
human<-  useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
name.dic <- getBM(attributes =c("ensembl_gene_id", "external_gene_name"),
      filters = "ensembl_gene_id",
      pd.ame_res.top5$ensembl.id,
      mart = human)
pd.ame_res.top5<- right_join(name.dic,pd.ame_res.top5,by=c("ensembl_gene_id"="ensembl.id"))


# Gene expression connection ----------------------------------------------
require(RSQLite)
sqlitePath <- "./data/zhang-rnaseq-db.sqlite"
mydb <- dbConnect(SQLite(), sqlitePath)
query <- sprintf("SELECT * FROM %s", "samples")
metadata.samples <- dbGetQuery(mydb, query)
query <- sprintf("SELECT gene_id,description FROM %s", "genes_no_mt_tpm_rescale")
genes <- dbGetQuery(mydb, query)
rownames(genes) <- genes$gene_id

getExp <-  function(gg){
  data.frame(TPM= as.numeric(dbGetQuery(mydb,'SELECT * FROM genes_no_mt_tpm_rescale WHERE "gene_id" = :g',
                                        params=list(g=gg))[metadata.samples$id]),
             metadata.samples,Gene=gg)%>%dplyr::select(id,TPM)%>%spread(key = id,value = TPM)
}

exp.dat <- sapply(pd.ame_res.top5$external_gene_name,
                  getExp,simplify = F)
all.equal(names(exp.dat),
          pd.ame_res.top5$external_gene_name)
exp.dat<-do.call(rbind,exp.dat)
exp.dat <- exp.dat[,c(9:10,1:8,11:18)]
pd.ame_res.top5<- right_join(pd.ame_res.top5,exp.dat%>%
                                 rownames_to_column("external_gene_name"))

# save 
write.csv(pd.ame_res.top5,file = "./data/zhang_top5_enrich_motif.csv")

# heatmap 

p1 <- ggplot(pd.ame_res.top5%>% 
  unite(gene,external_gene_name,subfamily)%>%
  dplyr::select("gene",starts_with("sample"))%>%
    mutate(gene=factor(gene,levels = gene)) %>%
  gather("sample","log10_adj_p.value",2:16)%>%
  mutate(sample=factor(sample,levels = paste0("sample_",1:15))),
  aes(sample,gene,fill=log10_adj_p.value))%+% geom_tile()%+%
  scale_fill_gradientn(colours =rev(brewer.pal(7,"YlGnBu")),na.value="NA")+
  theme(plot.background=element_blank(),
        axis.text=element_text(face="bold"),
        panel.grid = element_line(colour = "black"))+
  theme_bw()
  ggplot(pd.ame_res.top5%>% 
                 unite(gene,external_gene_name,subfamily)%>%
                 dplyr::select("gene",
                               starts_with("Me"),
                               starts_with("AEC"),
                               starts_with("VEC"))%>%
                 mutate(gene=factor(gene,levels = gene)) %>%
                 gather("sample","logTPM",2:18),
               aes(sample,gene,fill=logTPM))+ geom_tile()+
    geom_tile(colour="white",size=0.25)+
  scale_fill_gradientn(colours =rev(brewer.pal(7,"YlGnBu")),na.value="NA")
  
  theme(plot.background=element_blank(),
        axis.text=element_text(face="bold"),
        panel.grid = element_line(colour = "black"))+
  theme_bw()
  
# heatmap expression  

pd <- pd.ame_res.top5%>% 
  unite(gene,external_gene_name,subfamily)%>%
  dplyr::select("gene",starts_with("sample"))%>%
  column_to_rownames("gene")
b <- simplot(pd)


pd.2.dup <-right_join(pd.ame_res.top5[,1:2],
           pd.ame_res.hsap%>% 
             filter(ensembl.id %in% pd.ame_res.top5$ensembl_gene_id),
           by=c('ensembl_gene_id'="ensembl.id"))%>%
  unite(gene,external_gene_name,subfamily)%>%
  filter(duplicated(gene))
pd.2.dup <- pd.2.dup[c()]
pd.2 <- right_join(pd.ame_res.top5[,1:2],
                   pd.ame_res.hsap%>% 
  filter(ensembl.id %in% pd.ame_res.top5$ensembl_gene_id),
  by=c('ensembl_gene_id'="ensembl.id"))%>%
  unite(gene,external_gene_name,subfamily)%>%
  filter(!duplicated(gene))%>%
  dplyr::select("gene",starts_with("sample"))%>%
  column_to_rownames("gene")

pd.2[pd.2.dup$gene[5],1:15] <-pd.2.dup[5,3:17]
pd.2 <- pd.2[rownames(pd),paste0("sample_",1:15)]
b.2 <- simplot(pd.2)
b.3 <- simplot(apply(pd.2[rownames(pd),paste0("sample_",1:15)],2,normalise))

pd <- pd.ame_res.top5%>% 
    unite(gene,external_gene_name,subfamily)%>%
    dplyr::select("gene",
                  starts_with("Me"),
                  starts_with("AEC"),
                  starts_with("VEC"))%>%
    column_to_rownames("gene")

a<- simplot(pd)

  
a2<- simplot(t(apply(pd, 1,normalise)))


pd.2 <- pd %>%rownames_to_column("gene")%>%
  gather(key=sample,value = TPM,2:(ncol(pd)+1))%>%
  separate(sample,into =c("Cell","Day","Rep")) %>% 
  group_by(gene,Cell,Day)%>%
  summarise(avgTPM=mean(TPM))%>%
  unite(sample,Cell,Day)%>%
  spread(key = sample,value = avgTPM)%>%
  ungroup()%>%
  mutate(gene=factor(gene,levels = rownames(pd)))%>%
  arrange(gene)%>%as.data.frame()%>%
  column_to_rownames("gene")
pd.2 <- pd.2[,c(5,1:4,6:9)]

c<- simplot(pd.2)
c2<- simplot(t(apply(pd.2, 1, normalise)))


require(gridExtra)
grid.arrange(b+theme(axis.text.x = element_blank()),
             a+theme(axis.text = element_blank()),ncol=2)
  
grid.arrange(b+theme(axis.text.x = element_blank()),
             a2+theme(axis.text = element_blank()),ncol=2)

grid.arrange(b+theme(axis.text.x = element_blank()),
             c+theme(axis.text = element_blank()),ncol=2)

grid.arrange(b+theme(axis.text.x = element_blank()),
             b.3+theme(axis.text = element_blank()),
             c2+theme(axis.text = element_blank()),ncol=3)

grid.arrange(b+theme(axis.text.x = element_blank()),
             b.2+theme(axis.text = element_blank()),
             c+theme(axis.text = element_blank()),ncol=3)


# chu ---------------------------------------------------------------------

pd.chu_ame <- chu_ame%>% unite(tf,1:3)%>%
  filter(adj_p.value<=0.00001)%>%
  mutate(log10.adj_p.value=-log10(adj_p.value))%>%
  dplyr::select(-adj_p.value)%>%
  spread(key=sample,value = log10.adj_p.value,fill = 1)%>%
  column_to_rownames("tf")

pd.chu_ame.hsap <- chu_ame%>% unite(tf,1:3)%>%
  filter(adj_p.value<=0.00001)%>%
  filter(grepl("hsap",sample))%>%
  mutate(log10.adj_p.value=-log10(adj_p.value))%>%
  dplyr::select(-adj_p.value)%>%
  spread(key=sample,value = log10.adj_p.value,fill = 1)%>%
  column_to_rownames("tf")

pd.chu_ame.mmus <- chu_ame%>% unite(tf,1:3)%>%
  filter(adj_p.value<=0.00001)%>%
  filter(!grepl("hsap",sample))%>%
  mutate(log10.adj_p.value=-log10(adj_p.value))%>%
  dplyr::select(-adj_p.value)%>%
  spread(key=sample,value = log10.adj_p.value,fill = 1)%>%
  column_to_rownames("tf")


require(d3heatmap)
require(heatmaply)

heatmaply(apply(pd.chu_ame,2,normalise),
          seriate="mean",
          scale = "none",  
          colors=rev(brewer.pal(11,"RdYlBu")),
          k_col = 5,k_row = 15,column_text_angle = 90,
          key.title = "enrichment",margins = c(80,0,NA,0)
) %>% layout(width=600,height=1200)

d3heatmap(apply(pd.chu_ame,2,normalise),
          colors=rev(brewer.pal(11,"RdYlBu")),k_col = 15,k_row = 5)

d3heatmap(log2(pd.chu_ame.hsap),
          colors=rev(brewer.pal(11,"RdYlBu")),k_col = 15,k_row = 5)
d3heatmap(log2(pd.chu_ame.mmus),
          colors=rev(brewer.pal(11,"RdYlBu")),k_col = 15,k_row = 5)
d3heatmap(log2(pd.ame_res),
          colors=rev(brewer.pal(11,"RdYlBu")),k_col = 15,k_row = 5)
