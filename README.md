jaspar stats
============

The `JASPAR motifs` were downloaded from [Jaspar](http://jaspar.genereg.net/downloads/) choosing `JASPAR2018_CORE_vertebrates_non-redundant` in `MEME` format. Then a customized function `getEnsemblFromJaspar` [code](./aux_funs.R) was used to get Ensembl ID from `Japsar id` via `uiportID` from Jaspar's record.

``` r
# number of TF motifs in Jaspar 
nrow(jaspar) 
```

    ## [1] 579

``` r
# How many of them are protein complex? 
length(grep(":",jaspar$V3)) 
```

    ## [1] 45

tfclass stats
=============

[TFclass database](http://tfclass.bioinf.med.uni-goettingen.de/about.jsf) were downloaded (2018-05-14), [parsed](./scripts/parseTfclassTTL.sh) and [extacted](./scripts/tfclass.R) the human and mouse genes' Ensembl ID.

``` r
# How many unique genus (i.e. TF genes) in parsed TFClass?
length(unique(tfclass$merge$genus.id)) 
```

    ## [1] 1440

``` r
# How many subfamilies in parsed TFClass?
length(unique(tfclass$merge$subfamily.id))
```

    ## [1] 396

``` r
# How many families in parsed TFClass?
length(unique(tfclass$merge$family.id))
```

    ## [1] 110

jaspar-tfclass dictionary stats
===============================

`JASPAR` and `TFClass` were majorly [mapped](./cmpDB.R) by matching `Ensembl ID` except `20` TFs by `TF Name` directly.

``` r
# how many Jasper motif have annotation? 
a=nrow(dic$jasparTOtfclass %>% 
  dplyr::select(jaspar.id,genus.id) %>% 
  group_by(jaspar.id)%>%
  summarise(genus.id.2 = paste0(genus.id,collapse = ';')))
paste(a,nrow(jaspar),sep = "/")
```

    ## [1] "559/579"

``` r
# which jaspar motif don't have annotation?
noAnno <- subset(jaspar, V2 %in% setdiff(jaspar$V2,dic$jasparTOtfclass$jaspar.id))
noAnno.2 <- dic$jasparTOensembl[setdiff(jaspar$V2,dic$jasparTOtfclass$jaspar.id)]
knitr::kable(data.frame(jaspar.id=noAnno$V2,jaspar.name=noAnno$V3,
           ensembl.id= noAnno.2[noAnno$V2],
           url= paste0("http://jaspar2018.genereg.net/matrix/",noAnno$V2))%>%
  rowid_to_column('idx'))
```

|  idx| jaspar.id | jaspar.name  | ensembl.id         | url                                             |
|----:|:----------|:-------------|:-------------------|:------------------------------------------------|
|    1| MA0019.1  | Ddit3::Cebpa | NULL;NULL          | <http://jaspar2018.genereg.net/matrix/MA0019.1> |
|    2| MA0089.1  | MAFG::NFE2L1 | NULL;NULL          | <http://jaspar2018.genereg.net/matrix/MA0089.1> |
|    3| MA0109.1  | HLTF         | NULL               | <http://jaspar2018.genereg.net/matrix/MA0109.1> |
|    4| MA0111.1  | Spz1         | ENSMUSG00000046957 | <http://jaspar2018.genereg.net/matrix/MA0111.1> |
|    5| MA0149.1  | EWSR1-FLI1   | NULL;NULL          | <http://jaspar2018.genereg.net/matrix/MA0149.1> |
|    6| MA0611.1  | Dux          | logical(0)         | <http://jaspar2018.genereg.net/matrix/MA0611.1> |
|    7| MA0619.1  | LIN54        | NULL               | <http://jaspar2018.genereg.net/matrix/MA0619.1> |
|    8| MA0621.1  | mix-a        | NULL               | <http://jaspar2018.genereg.net/matrix/MA0621.1> |
|    9| MA0629.1  | Rhox11       | ENSMUSG00000051038 | <http://jaspar2018.genereg.net/matrix/MA0629.1> |
|   10| MA0633.1  | Twist2       | ENSMUSG00000007805 | <http://jaspar2018.genereg.net/matrix/MA0633.1> |
|   11| MA0637.1  | CENPB        | ENSG00000125817    | <http://jaspar2018.genereg.net/matrix/MA0637.1> |
|   12| MA0643.1  | Esrrg        | ENSMUSG00000026610 | <http://jaspar2018.genereg.net/matrix/MA0643.1> |
|   13| MA0659.1  | MAFG         | ENSG00000197063    | <http://jaspar2018.genereg.net/matrix/MA0659.1> |
|   14| MA0048.2  | NHLH1        | ENSG00000171786    | <http://jaspar2018.genereg.net/matrix/MA0048.2> |
|   15| MA0690.1  | TBX21        | ENSG00000073861    | <http://jaspar2018.genereg.net/matrix/MA0690.1> |
|   16| MA0592.2  | Esrra        | ENSMUSG00000024955 | <http://jaspar2018.genereg.net/matrix/MA0592.2> |
|   17| MA0017.2  | NR2F1        | ENSG00000175745    | <http://jaspar2018.genereg.net/matrix/MA0017.2> |
|   18| MA0800.1  | EOMES        | ENSG00000163508    | <http://jaspar2018.genereg.net/matrix/MA0800.1> |
|   19| MA0802.1  | TBR1         | ENSG00000136535    | <http://jaspar2018.genereg.net/matrix/MA0802.1> |
|   20| MA1111.1  | NR2F2        | ENSG00000185551    | <http://jaspar2018.genereg.net/matrix/MA1111.1> |
