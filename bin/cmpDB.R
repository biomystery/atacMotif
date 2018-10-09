
source("./aux_funs.R")

jaspar <- read.table("./db/jaspar2018_vertebrates.txt",skip = 1,stringsAsFactors = F)[,2:3] #579
tfclass <- readRDS("./db/tfclass.rds") 
dic<- readRDS("./db/dic_jaspar_tfclass.rds")


# update dic 1---------------------------------------------------------

# jaspa dic 
jasparTOensembl <- setdiff(jaspar$V2,names(dic$jasparTOensembl))# 138
names(jasparTOensembl) <- jasparTOensembl

for(i in 1:length(jasparTOensembl)) {
  print(i);
  jasparTOensembl[i]=getEnsemblFromJaspar(jasparTOensembl[i])}

# tfclass dic 
ensemblTOsubfamily <- tfclass$merge%>%
  separate_rows(tf.id,sep = ";") %>%
  distinct(tf.id,.keep_all = T)%>%
  column_to_rownames("tf.id")

# jaspa to tfclass 
jasparTOtfclass <-data.frame(ensembl.id = jasparTOensembl)%>%
  rownames_to_column("jaspar.id")%>%
  separate_rows(ensembl.id,sep = ";") %>%
  right_join(ensemblTOsubfamily%>%rownames_to_column("ensembl.id"))%>%
  filter(!is.na(jaspar.id))
  

# added 
dic$jasparTOensembl <- c(dic$jasparTOensembl,jasparTOensembl)
dic$ensemblTOtfclass  <- rbind(dic$ensemblTOtfclass,ensemblTOsubfamily)
dic$jasparTOtfclass <- rbind(dic$jasparTOtfclass,jasparTOtfclass)

# update dic 2: directly match genus name ---------------------------------------------------------
noAnno <- subset(jaspar, V2 %in% setdiff(jaspar$V2,dic$jasparTOtfclass$jaspar.id))
idx <- tolower(noAnno$V3) %in% tolower(tfclass$merge$genus.name) 
sum(idx) #15 more 
idx.2 <- tolower(tfclass$merge$genus.name)%in% tolower(noAnno$V3)
directMatch.2 <-  noAnno[idx,] %>% 
  mutate(V3=toupper(V3))%>%
  left_join(tfclass$merge[idx.2,]%>%mutate(genus.name=toupper(genus.name)),by = c("V3" = "genus.name"))%>%
  rename("V2"="jaspar.id","V3"="genus.name","tf.id"="ensembl.id")%>%
  separate_rows("ensembl.id") # 32 rows
  
dic$jasparTOtfclass <- rbind(dic$jasparTOtfclass,directMatch.2[,colnames(dic$jasparTOtfclass)])


# Update3: merge all  -----------------------------------------------------
dim(dic$jasparTOtfclass)
tmp.2 <- dic$jasparTOtfclass%>%
  group_by(jaspar.id)%>%
  summarise_all(funs(paste(unique(.),collapse = ";")))

tmp<-rename(jaspar,jaspar.id=V2,jaspar.name=V3) %>% 
  full_join(dic$jasparTOtfclass)
length(unique(jaspar$V2)) #579
length(unique(tmp$jaspar.id)) #579

tmp.1<- tmp %>% full_join(tfclass$merge%>%
               rename(ensembl.id=tf.id))%>%
  separate_rows(ensembl.id)
dim(tmp.1) #2099,11; 3389,11
length(unique(tmp.1$jaspar.id)) #580 (add NA)
length(unique(tfclass$merge$subfamily.id)) #396
length(unique(tmp.1$subfamily.id)) #397

dic$merged <- tmp.1

# update 4: remove duplicated rows ----------------------------------------
dic<- readRDS("./db/dic_jaspar_tfclass.rds")

# eg
dic$merged %>% filter(ensembl.id=="ENSG00000134954")
red <- lapply(unique(dic$merged$ensembl.id),function(x){
  tmp <- dic$merged %>% filter(ensembl.id==x)
  if(nrow(tmp[,-c(1:2)]%>% distinct())<nrow(tmp)){
    return(tmp[!is.na(tmp$jaspar.name),])
  }else{
    return(tmp)
  }
})
x <- "ENSG00000134954"
length(red)
red[[2]]
red <- do.call(rbind,red)
head(red)
red %>% filter(ensembl.id=="ENSG00000134954")

# handle NA i.e. rescures
# https://github.com/epigen-UCSD/atacMotif/blob/master/db/rescue_Jaspar.txt
sum(is.na(red$ensembl.id))
sum(is.na(unique(dic$merged$ensembl.id)))
dim(dic$merged[is.na(dic$merged$ensembl.id),]) #20
#write.csv(tmp,file="./db/rescue_Jaspar.csv",quote = F,row.names = F)
tmp <- read.csv("./db/rescue_Jaspar.csv",header = T,stringsAsFactors = F)
dic$merged <- rbind(red,tmp)


saveRDS(dic,file = "./db/dic_jaspar_tfclass.rds")


