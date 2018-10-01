
source("./aux_funs.R")
tfclass <- readRDS("./data/tfclass.rds")

# ame’s output ------------------------------------------------------------
zhang_ame <- read.table("./data/ame.all.res.txt",
                      col.names = c('jaspar.id',"tf.name","tf.seq","padj","clust"),
                      stringsAsFactors = F)


chu_ame <- read.table("./data/chu_ame_res.txt",
                      col.names = c('jaspar.id',"tf.name","tf.seq","padj","clust"),
                      stringsAsFactors = F)

#jaspar_to_TFclass dic---------------------------------------------------------

# filter padj <1e-5
th <- .01
zhang_ame.anno <- zhang_ame%>% 
  filter(padj<=th)

# jaspa dic 
dic.jasparTOensembl <- unique(zhang_ame.anno$jaspar.id)
names(dic.jasparTOensembl) <- unique(zhang_ame.anno$jaspar.id)
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
  rownames_to_column("jaspar.id")%>%
  separate_rows(ensembl.id,sep = ";") %>%
  right_join(dic.ensemblTOsubfamily%>%rownames_to_column("ensembl.id"))%>%
  filter(!is.na(jaspar.id))
  

saveRDS(list(jasparTOensembl=dic.jasparTOensembl,
             ensemblTOtfclass=dic.ensemblTOsubfamily,
             jasparTOtfclass=dic.jasparTOtfclass),
        file = "./data/dic_jaspar_tfclass.rds")
dic<- readRDS("./data/dic_jaspar_tfclass.rds")
attach(dic)
# heatmap -----------------------------------------------------------------
tmp <- sapply(unique(zhang_ame[zhang_ame$padj<=th,"jaspar.id"]),function(a)
  (dic.jasparTOtfclass%>%
     filter(row_number()==grep(a,dic.jasparTOtfclass$jaspar.id)[1])),simplify = F)
tmp<-do.call(rbind,tmp)


pd.zhang_ame <- zhang_ame%>%
  filter(padj<=th)%>%
  mutate(log10.padj=-log10(padj))%>%
  dplyr::select(-padj)%>%
  spread(key=clust,value = log10.padj,fill = 0)
#  semi_join(tmp[pd.zhang_ame$jaspar.id,],copy = T)

pd.zhang_ame<-cbind(pd.zhang_ame,tmp[pd.zhang_ame$jaspar.id,])%>%
  unite(tf,1:3)%>%
  remove_rownames()%>%
  column_to_rownames("tf")%>%
  unite(family,family.id,family.name)%>%
  unite(subfamily,subfamily.id,subfamily.name)

col.max <- apply(pd.zhang_ame[,grepl("clust",colnames(pd.zhang_ame))],
                 2,max)
col.max.col <- cut(col.max,seq(min(col.max)-.01,max(col.max)+.01,length.out = 20))
col.map <- colorRampPalette(brewer.pal(9,"Blues"))(19)
names(col.map)<- levels(col.max.col)
col.map.col <- col.map[col.max.col]
names(col.map.col) <- colnames(pd.zhang_ame)[grepl("clust",colnames(pd.zhang_ame))]

normalise <- function(x) (x-min(x))/(max(x)-min(x))
pd <- apply(pd.zhang_ame[,grepl("clust",colnames(pd.zhang_ame))],2,
            normalise)

pd <- pd.zhang_ame[grepl("ENSG",pd.zhang_ame$ensembl.id),
                   grepl("clust",colnames(pd.zhang_ame))]
pd.zhang_ame.hsap <- pd.zhang_ame[grepl("ENSG",pd.zhang_ame$ensembl.id),]

pd.scale <- scale(pd)
require(heatmaply)
require(RColorBrewer)
pd.2 <-apply(pd,2,normalise) %>% 
  as.data.frame %>%
  arrange(clust_1,clust_2,clust_3,clust_4,clust_5,clust_6,clust_7,clust_8,clust_9,clust_10,clust_11,clust_12,clust_13,clust_14,clust_15)
rownames(pd.2)<- rownames(pd)
heatmaply(apply(pd,2,normalise),
          seriate="mean",
          scale = "none",  
          colors=rev(brewer.pal(11,"RdYlBu")),
          k_col = 5,k_row = 15,
          RowSideColors = pd.zhang_ame[grepl("ENSG",pd.zhang_ame$ensembl.id),
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
             pd.zhang_ame=pd.zhang_ame),file = './data/ame.res.rds')





# TOP5 TFs  ----------------------------------------------
pd.zhang_ame.top5<- pd.zhang_ame.hsap%>% 
  as.tibble()%>% 
  gather(key="clust",value = "log10_padj",1:15)%>%
  group_by(clust)%>%
  top_n(5,log10_padj)%>%
  spread(clust,log10_padj)

# print
pd.zhang_ame.top5<-pd.zhang_ame.top5[,c(colnames(pd.zhang_ame.top5)[1:7],paste0("clust_",1:15))]%>% 
  arrange(clust_1,clust_2,clust_3,clust_4,clust_5,clust_6,clust_7,clust_8,clust_9,clust_10,clust_11,clust_12,clust_13,clust_14,clust_15)%>%
  print(n=28)


# ensemble -> gene symbol 
require(biomaRt)
human<-  useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
name.dic <- getBM(attributes =c("ensembl_gene_id", "external_gene_name"),
      filters = "ensembl_gene_id",
      pd.zhang_ame.top5$ensembl.id,
      mart = human)
pd.zhang_ame.top5<- right_join(name.dic,pd.zhang_ame.top5,by=c("ensembl_gene_id"="ensembl.id"))


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

exp.dat <- sapply(pd.zhang_ame.top5$external_gene_name,
                  getExp,simplify = F)
all.equal(names(exp.dat),
          pd.zhang_ame.top5$external_gene_name)
exp.dat<-do.call(rbind,exp.dat)
exp.dat <- exp.dat[,c(9:10,1:8,11:18)]
pd.zhang_ame.top5<- right_join(pd.zhang_ame.top5,exp.dat%>%
                                 rownames_to_column("external_gene_name"))

# save 
write.csv(pd.zhang_ame.top5,file = "./data/zhang_top5_enrich_motif.csv")

# heatmap 

p1 <- ggplot(pd.zhang_ame.top5%>% 
  unite(gene,external_gene_name,subfamily)%>%
  dplyr::select("gene",starts_with("clust"))%>%
    mutate(gene=factor(gene,levels = gene)) %>%
  gather("clust","log10_padj",2:16)%>%
  mutate(clust=factor(clust,levels = paste0("clust_",1:15))),
  aes(clust,gene,fill=log10_padj))%+% geom_tile()%+%
  scale_fill_gradientn(colours =rev(brewer.pal(7,"YlGnBu")),na.value="NA")+
  theme(plot.background=element_blank(),
        axis.text=element_text(face="bold"),
        panel.grid = element_line(colour = "black"))+
  theme_bw()
  ggplot(pd.zhang_ame.top5%>% 
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

pd <- pd.zhang_ame.top5%>% 
  unite(gene,external_gene_name,subfamily)%>%
  dplyr::select("gene",starts_with("clust"))%>%
  column_to_rownames("gene")
b <- simplot(pd)


pd.2.dup <-right_join(pd.zhang_ame.top5[,1:2],
           pd.zhang_ame.hsap%>% 
             filter(ensembl.id %in% pd.zhang_ame.top5$ensembl_gene_id),
           by=c('ensembl_gene_id'="ensembl.id"))%>%
  unite(gene,external_gene_name,subfamily)%>%
  filter(duplicated(gene))
pd.2.dup <- pd.2.dup[c()]
pd.2 <- right_join(pd.zhang_ame.top5[,1:2],
                   pd.zhang_ame.hsap%>% 
  filter(ensembl.id %in% pd.zhang_ame.top5$ensembl_gene_id),
  by=c('ensembl_gene_id'="ensembl.id"))%>%
  unite(gene,external_gene_name,subfamily)%>%
  filter(!duplicated(gene))%>%
  dplyr::select("gene",starts_with("clust"))%>%
  column_to_rownames("gene")

pd.2[pd.2.dup$gene[5],1:15] <-pd.2.dup[5,3:17]
pd.2 <- pd.2[rownames(pd),paste0("clust_",1:15)]
b.2 <- simplot(pd.2)
b.3 <- simplot(apply(pd.2[rownames(pd),paste0("clust_",1:15)],2,normalise))

pd <- pd.zhang_ame.top5%>% 
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
  filter(padj<=0.00001)%>%
  mutate(log10.padj=-log10(padj))%>%
  dplyr::select(-padj)%>%
  spread(key=clust,value = log10.padj,fill = 1)%>%
  column_to_rownames("tf")

pd.chu_ame.hsap <- chu_ame%>% unite(tf,1:3)%>%
  filter(padj<=0.00001)%>%
  filter(grepl("hsap",clust))%>%
  mutate(log10.padj=-log10(padj))%>%
  dplyr::select(-padj)%>%
  spread(key=clust,value = log10.padj,fill = 1)%>%
  column_to_rownames("tf")

pd.chu_ame.mmus <- chu_ame%>% unite(tf,1:3)%>%
  filter(padj<=0.00001)%>%
  filter(!grepl("hsap",clust))%>%
  mutate(log10.padj=-log10(padj))%>%
  dplyr::select(-padj)%>%
  spread(key=clust,value = log10.padj,fill = 1)%>%
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
d3heatmap(log2(pd.zhang_ame),
          colors=rev(brewer.pal(11,"RdYlBu")),k_col = 15,k_row = 5)
