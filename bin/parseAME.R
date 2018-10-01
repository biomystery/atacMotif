#! /bin/R

files <- list.files(pattern="*ame.tsv",recursive=T)

res <- lapply(files,read.table,header=T)
names(res) <- sub("_10.*","",files)

res.2 <- lapply(names(res),function(x){
    tmp <- res[[x]]
    data.frame(tmp[,-c(1,2)],sample=x,stringsAsFactors=F)
})

res <- do.call(rbind,res.2)

write.table(res,file="ame.all.res.txt",quote=F,sep="\t",col.names=NA)
