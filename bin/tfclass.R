tfclass <- list()
tfclass$subfamily <- read.table("./data/subfamily_all.txt",sep = "\t",
                                col.names = c("id","name","seq"),fill = T,
                                stringsAsFactors = F)
tfclass$family <- read.table("./data/motifs/family.txt",sep = "\t",
                             col.names = c("id",'name'),
                             stringsAsFactors = F)

tfclass$genus <- read.table("./data/motifs/genus.txt",sep = "\t",
                             col.names = c("id",'name'),
                            stringsAsFactors = F)

# backup duplicated genes 
tfclass$genus.dup <-  tfclass$genus %>% filter(duplicated(name))
tfclass$genus <- tfclass$genus %>% filter(!duplicated(name)) %>% arrange(name)

# ensembl for tfs 
tfclass$tfs_species <- read.table("./data/motifs/tfs.txt",sep = "\t",
                            col.names = c("symbol",'ensembl'),
                            stringsAsFactors = F)
tfclass$tfs_species$symbol <- sub(" ","",tfclass$tfs_species$symbol)

# collapse ensemble id 
tfclass$tfs_species<-tfclass$tfs_species %>% group_by(symbol)%>%
  summarise(ensembl.id=paste(ensembl,collapse = ";"))

# merge   ----------------------------------------------------------------
# merge tf_species into genus 
# case 1: toupper 
tmp.genus <- toupper(as.character((tfclass$genus%>% arrange(name))$name)) #to upper
tmp.tfs <- toupper(as.character(tfclass$tfs_species$symbol)) #to upper
# case 2: rm - 
tmp.genus <- sub('-','',tmp.genus);tmp.tfs <- sub('-','',tmp.tfs) # remove -

# after 1st two case 
tmp.common <- intersect(tmp.tfs,tmp.genus) #1289
tmp.genus.s <- setdiff(tmp.genus,tmp.common) #164
tmp.tfs.s <- setdiff(tmp.tfs,tmp.common) #186


# original name  & merge 
genus.id.1 <- (tmp.genus %in% tmp.common)
genus.s <- tfclass$genus[!genus.id.1,] #original name,164 


tfs.id.1 <- (tmp.tfs %in% tmp.common)
tfs.s <- tfclass$tfs_species[!tfs.id.1,] #original name ,186


genus.ord.1 <- sapply(tmp.common,function(x) which(tmp.genus==x))
tfs.ord.1 <- sapply(tmp.common,function(x) which(tmp.tfs==x))
merge.1 <- cbind(tfclass$genus[genus.ord.1,],tfclass$tfs_species[tfs.ord.1,]) #1289

### merge 2 
tmp.grep <- unlist(sapply(tmp.tfs.s, function(x) if(length(grep(x,tmp.genus.s))==1) grep(x,tmp.genus.s)))
merge.2 <- data.frame(tfs=names(tmp.grep),genus=tmp.genus.s[tmp.grep])
tfs.id.2 <- sapply(names(tmp.grep), function(x) which(tmp.tfs.s == x))
merge.2 <- cbind(genus.s[tmp.grep,],tfs.s[tfs.id.2,]) #1289

### remaining tfs 
tfs.s.2 <- tfs.s[-tfs.id.2,] #34 
tmp.tfs.s.2 <- tmp.tfs.s[-as.numeric(tfs.id.2)]

# Merge 3
tmp.genus.s.1 <- sub(" ","",sub("\\(.*","",tmp.genus.s)) # before (
tmp.common.2 <- intersect(tmp.genus.s.1, tmp.tfs.s) #106 
tmp.tfs.s.2 <- setdiff(tmp.tfs.s,tmp.common.2) #80
tmp.genus.s.2 <- tmp.genus.s[!(tmp.genus.s.1 %in% tmp.common.2)] #58 

# grep -e "http://sybig.de/tfclass#Mus_musculus_ZNF41>" -A 6 tfclass.ttl | perl -ne '/(\d+.\d+.\d+.\d+.\d+)>/ && print "$1 \n"'
#for tf in ${tfs[@]}; do echo $tf; grep -e "http://sybig.de/tfclass#Mus_musculus_${tf}>" -A 6 tfclass.ttl | perl -ne '/(\d+.\d+.\d+.\d+.\d+)>/ && print "$1\n"'; done 
sink(file = "remainTFs.txt")
cat(as.character(tfs.s.2$symbol))
sink()
#  bash ./remainTFs.sh > remainTFs.res.tab 

tfs.s.2.res <- read.table("./data/remainTFs.res.tab",sep = "\t",stringsAsFactors = F,
                          col.names = c("symbol","id"))
tmp.dic <- tfclass$genus %>% column_to_rownames("id")
all.equal(tfs.s.2.res$symbol,as.character(tfs.s.2$symbol))
merge.3 <- data.frame(id=tfs.s.2.res[,2],
           name=tmp.dic[tfs.s.2.res$id,1],
           tfs.s.2,stringsAsFactors = F) 

# merge genus & tfs_species
tfclass$merge <- rbind(merge.1,merge.2,merge.3)
tfclass$tfs_species <- NULL 

# merge with subfamily

tfclass$subfamily <- tfclass$subfamily%>% mutate(seq = ifelse(grepl("http:",seq),NA,seq))
tmp.dic <- tfclass$subfamily %>% column_to_rownames("id")
tfclass$merge<- tfclass$merge %>%
  mutate(subfamily.id=sub("(\\.[0-9]+$)","",id))
tfclass$merge <- cbind(tfclass$merge,tmp.dic[tfclass$merge$subfamily.id,])
names(tfclass$merge) <- c("genus.id","genus.name","tf.symbol",'tf.id',"subfamily.id","subfamily.name","subfamily.seq")

# merge with family
tfclass$merge<- tfclass$merge %>%
  mutate(family.id=sub("(\\.[0-9]+$)","",subfamily.id))

tfclass$merge<- (data.frame(tfclass$merge,family.name=tmp.dic[tfclass$merge$family.id,]))

saveRDS(tfclass,file = "./data/tfclass.rds")


# add superfamily and class -----------------------------------------------
tfclass <- readRDS('./db/tfclass.rds')
tfclass$class <- fread('./db/class.txt',col.names = c('id','name','about'))
tfclass$superclass <- fread('./db/superclass.txt',col.names = c('id','name','about'))
saveRDS(tfclass,file = "./db/tfclass.rds")
