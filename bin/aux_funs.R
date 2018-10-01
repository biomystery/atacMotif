
require(DESeq2)
require(tidyverse)
require(rtracklayer)
require(RColorBrewer)
require(ggplot2)
require(BSgenome.Hsapiens.UCSC.hg38)
require(ChIPseeker)
require(ReactomePA)
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
require(cqn)
require(liftOver)
require(DT)
require(rGREAT)
require(TFBSTools)
require(JASPAR2018)
require(rvest)
require(magrittr)

col.map <- c("AEC"="darkgoldenrod","VEC"="dodgerblue","Mesoderm"="gray80")

plotAnnoPeak <- function(peak_list,splitToFile=F){
  #peak_list - bed file or GRange objects 
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  peakAnno <- annotatePeak(peak = peak_list,
                           tssRegion=c(-3000, 3000),level = "gene",
                           genomicAnnotationPriority=c("Promoter","Intron","Intergenic"),
                           TxDb=txdb)
  plotAnnoPie(peakAnno)
  
  # split peaks to different annotation 
  if(splitToFile){
    output.bed <- readLines(peak_list)
    con <- file("promoter.bed","w")
    writeLines(text = output.bed[peakAnno@detailGenomicAnnotation$Promoter],con = con)
    close(con)
  }
}


gcContent <- function(gr,bsgenome,...){
  reg.seqs <- getSeq(bsgenome, gr)
  base.frequency <- alphabetFrequency(reg.seqs)[,1:4]
  (base.frequency[,'G']+base.frequency[,'C'])/width(reg.seqs)
}

lift2hg19<-function(cur){
  
  path = list.files(path = "../", pattern = "hg38ToHg19.over.chain",recursive = T,full.names = T)
  ch = import.chain(path)
  seqlevelsStyle(cur) = "UCSC"  # necessary
  unlist(liftOver(cur, ch))
}

# convert tf_id to transcript_id 
getEnsemblFromJaspar<- function(jaspar.id="MA0114.3"){
  require(JASPAR2018)
  res <- getMatrixByID(JASPAR2018,jaspar.id)
  out<-sapply(sub("( +)","",res@tags$acc), function(x){
    res.get <- tryCatch(read_html(paste0("https://www.uniprot.org/uniprot/",x,".txt"))%>% 
      html_text(trim = T),error=function(cond) return(NA))
    pat.human <- ".*(ENSG([0-9]+)).*"
    pat.mouse <- ".*(ENSMUSG([0-9]+)).*"
    pat.mgi <- ".*(MGI:([0-9]+)).*"
    if(grepl(pat.mouse,res.get)) return(sub(pat.mouse,"\\1",res.get))
    if(grepl(pat.human,res.get)) return(sub(pat.human,"\\1",res.get))
    if(grepl(pat.mgi,res.get)) {
      require(biomaRt)
      mouse<- useMart("ensembl",dataset="mmusculus_gene_ensembl")
      getBM(attributes = "ensembl_gene_id",
            filters = "mgi_id",
            sub(pat.mgi,"\\1",res.get),
            mart = mouse)$ensembl_gene_id
      }
  })
  paste0(out,collapse = ";")
}

findSubfamily <- function(g.id,db.tfclass=tfclass$merge){
  ### 1. got enseml id 
  #g.id <-getEnsemblFromJaspar(jid)
  tf.id <- grep(g.id,db.tfclass$tf.id)
  subfamily.id<- paste0(db.tfclass[tf.id,"subfamily.id"],collapse = ";")
  #return(db.tfclass[db.tfclass$subfamily.id==subfamily.id,])
  #return(list(all=,
  #            subfamily.id= subfamily.id))
}

getExp <-  function(gg){
  require(RSQLite)
  sqlitePath <- "./data/zhang-rnaseq-db.sqlite"
  mydb <- dbConnect(SQLite(), sqlitePath)
  query <- sprintf("SELECT * FROM %s", "samples")
  metadata.samples <- dbGetQuery(mydb, query)
  
  data.frame(TPM= as.numeric(dbGetQuery(mydb,'SELECT * FROM genes_no_mt_tpm_rescale WHERE "gene_id" = :g',
                                        params=list(g=gg))[metadata.samples$id]),
             metadata.samples,Gene=gg)%>%dplyr::select(id,TPM)%>%spread(key = id,value = TPM)
}
