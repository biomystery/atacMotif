jaspar <- read.table("./db/jaspar2018_vertebrates.txt",skip = 1,stringsAsFactors = F)[,2:3] #579
tfclass <- readRDS("./db/tfclass.rds") 
dic<- readRDS("./db/dic_jaspar_tfclass.rds")


# rescue 20 jaspar motif --------------------------------------------------

# old 
a=nrow(dic$jasparTOtfclass %>% 
         dplyr::select(jaspar.id,genus.id) %>% 
         group_by(jaspar.id)%>%
         summarise(genus.id.2 = paste0(genus.id,collapse = ';')))
paste(a,nrow(jaspar),sep = "/")

head(dic$jasparTOtfclass)

# remove from jaspar
jaspar.rmlist <- c("MA0109.1","MA0619.1","MA0621.1")

