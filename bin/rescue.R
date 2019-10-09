jaspar <- read.table("../db/jaspar2018_vertebrates.txt",skip = 1,stringsAsFactors = F)[,2:3] #579
tfclass <- readRDS("../db/tfclass.rds")
dic<- readRDS("../db/dic_jaspar_tfclass.rds")


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

############################################################
## fix na in family.id in dic$merge
require(tidyverse)

(rowids <- dic$merged %>% rowid_to_column("id") %>% filter(!is.na(subfamily.id) &
                                                           is.na(family.id)) %>% pull(id))
## [1] 2846 2848 2849 2850 2851 2853 2854 2855 2857
dic$merged[rowids, ]
##########################################################################
##      jaspar.id  jaspar.name ensembl.id genus.id genus.name tf.symbol ##
## 2846  MA0019.1 Ddit3::Cebpa       <NA>     <NA>       <NA>      <NA> ##
## 2848  MA0149.1   EWSR1-FLI1       <NA>     <NA>       <NA>      <NA> ##
## 2849  MA0611.1          Dux       <NA>     <NA>       <NA>      <NA> ##
## 2850  MA0629.1       Rhox11       <NA>     <NA>       <NA>      <NA> ##
## 2851  MA0633.1       Twist2       <NA>     <NA>       <NA>      <NA> ##
## 2853  MA0643.1        Esrrg       <NA>     <NA>       <NA>      <NA> ##
## 2854  MA0659.1         MAFG       <NA>     <NA>       <NA>      <NA> ##
## 2855  MA0048.2        NHLH1       <NA>     <NA>       <NA>      <NA> ##
## 2857  MA0592.2        Esrra       <NA>     <NA>       <NA>      <NA> ##
##      subfamily.id subfamily.name subfamily.seq family.id family.name ##
## 2846      1.1.8.1           <NA>          <NA>      <NA>        <NA> ##
## 2848      3.5.2.1           <NA>          <NA>      <NA>        <NA> ##
## 2849      3.1.3.7           <NA>          <NA>      <NA>        <NA> ##
## 2850    3.1.3.23            <NA>          <NA>      <NA>        <NA> ##
## 2851      1.2.3.2           <NA>          <NA>      <NA>        <NA> ##
## 2853      2.1.1.2           <NA>          <NA>      <NA>        <NA> ##
## 2854      1.1.3.2           <NA>          <NA>      <NA>        <NA> ##
## 2855      1.2.3.1           <NA>          <NA>      <NA>        <NA> ##
## 2857      2.1.1.2           <NA>          <NA>      <NA>        <NA> ##
##########################################################################


(tmp<- dic$merged[rowids, ] %>% mutate(subfamily.id = trimws(subfamily.id)) %>%
    mutate(family.id = sub(".[0-9]+$", "", subfamily.id)) %>% left_join(tfclass$family %>%
                                                                        select(id, name), by = c(family.id = "id")) %>% select(-family.name) %>% rename(family.name = name))

########################################################################################
## + +   jaspar.id  jaspar.name ensembl.id genus.id genus.name tf.symbol subfamily.id ##
## 1  MA0019.1 Ddit3::Cebpa       <NA>     <NA>       <NA>      <NA>      1.1.8.1     ##
## 2  MA0149.1   EWSR1-FLI1       <NA>     <NA>       <NA>      <NA>      3.5.2.1     ##
## 3  MA0611.1          Dux       <NA>     <NA>       <NA>      <NA>      3.1.3.7     ##
## 4  MA0629.1       Rhox11       <NA>     <NA>       <NA>      <NA>     3.1.3.23     ##
## 5  MA0633.1       Twist2       <NA>     <NA>       <NA>      <NA>      1.2.3.2     ##
## 6  MA0643.1        Esrrg       <NA>     <NA>       <NA>      <NA>      2.1.1.2     ##
## 7  MA0659.1         MAFG       <NA>     <NA>       <NA>      <NA>      1.1.3.2     ##
## 8  MA0048.2        NHLH1       <NA>     <NA>       <NA>      <NA>      1.2.3.1     ##
## 9  MA0592.2        Esrra       <NA>     <NA>       <NA>      <NA>      2.1.1.2     ##
##   subfamily.name subfamily.seq family.id               family.name                 ##
## 1           <NA>          <NA>     1.1.8              CEBP-related                 ##
## 2           <NA>          <NA>     3.5.2               Ets-related                 ##
## 3           <NA>          <NA>     3.1.3         Paired-related HD                 ##
## 4           <NA>          <NA>                                                     ##
##     3.1.3         Paired-related HD                                                ##
## 5           <NA>          <NA>     1.2.3               Tal-related                 ##
## 6           <NA>          <NA>     2.1.1 Steroid hormone receptors                 ##
## 7           <NA>          <NA>     1.1.3               Maf-related                 ##
## 8           <NA>          <NA>     1.2.3               Tal-related                 ##
## 9           <NA>          <NA>     2.1.1 Steroid hormone receptors                 ##
########################################################################################

dic$merged[rowids, ] <- tmp
saveRDS(dic,'../db/dic_jaspar_tfclass.rds')
