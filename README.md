jasper stats
============

``` r
# number of TF motifs in Jaspar 
nrow(jaspar) #579 TF motifs 
```

    ## [1] 579

``` r
# How many of them are protein complex? 
length(grep(":",jaspar$V3)) #45 complex 
```

    ## [1] 45

tfclass stats
=============

TFclass database were parsed and extacted the human and mouse genes.

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

``` r
# how many Jasper motif have annotation? 
a=nrow(dic$jasparTOtfclass %>% 
  dplyr::select(jaspar.id,genus.id) %>% 
  #distinct(jaspar.id,.keep_all = F)
  group_by(jaspar.id)%>%
  summarise(genus.id.2 = paste0(genus.id,collapse = ';')))
paste(a,nrow(jaspar),sep = "/")
```

    ## [1] "559/579"

``` r
# which jaspar motif don't have annotation?
noAnno <- subset(jaspar, V2 %in% setdiff(jaspar$V2,dic$jasparTOtfclass$jaspar.id))
noAnno.2 <- dic$jasparTOensembl[setdiff(jaspar$V2,dic$jasparTOtfclass$jaspar.id)]
nrow(noAnno) #20 out of 579 
```

    ## [1] 20

``` r
knitr::kable(data.frame(jaspar.id=noAnno$V2,jaspar.name=noAnno$V3,
           ensembl.id= noAnno.2[noAnno$V2])%>%
  rowid_to_column('idx'))
```

|  idx| jaspar.id | jaspar.name  | ensembl.id         |
|----:|:----------|:-------------|:-------------------|
|    1| MA0019.1  | Ddit3::Cebpa | NULL;NULL          |
|    2| MA0089.1  | MAFG::NFE2L1 | NULL;NULL          |
|    3| MA0109.1  | HLTF         | NULL               |
|    4| MA0111.1  | Spz1         | ENSMUSG00000046957 |
|    5| MA0149.1  | EWSR1-FLI1   | NULL;NULL          |
|    6| MA0611.1  | Dux          | logical(0)         |
|    7| MA0619.1  | LIN54        | NULL               |
|    8| MA0621.1  | mix-a        | NULL               |
|    9| MA0629.1  | Rhox11       | ENSMUSG00000051038 |
|   10| MA0633.1  | Twist2       | ENSMUSG00000007805 |
|   11| MA0637.1  | CENPB        | ENSG00000125817    |
|   12| MA0643.1  | Esrrg        | ENSMUSG00000026610 |
|   13| MA0659.1  | MAFG         | ENSG00000197063    |
|   14| MA0048.2  | NHLH1        | ENSG00000171786    |
|   15| MA0690.1  | TBX21        | ENSG00000073861    |
|   16| MA0592.2  | Esrra        | ENSMUSG00000024955 |
|   17| MA0017.2  | NR2F1        | ENSG00000175745    |
|   18| MA0800.1  | EOMES        | ENSG00000163508    |
|   19| MA0802.1  | TBR1         | ENSG00000136535    |
|   20| MA1111.1  | NR2F2        | ENSG00000185551    |

``` r
noAnno.2 # 20 Jaspers not found annotation 
```

    ##             MA0019.1             MA0089.1             MA0109.1 
    ##          "NULL;NULL"          "NULL;NULL"               "NULL" 
    ##             MA0111.1             MA0149.1             MA0611.1 
    ## "ENSMUSG00000046957"          "NULL;NULL"         "logical(0)" 
    ##             MA0619.1             MA0621.1             MA0629.1 
    ##               "NULL"               "NULL" "ENSMUSG00000051038" 
    ##             MA0633.1             MA0637.1             MA0643.1 
    ## "ENSMUSG00000007805"    "ENSG00000125817" "ENSMUSG00000026610" 
    ##             MA0659.1             MA0048.2             MA0690.1 
    ##    "ENSG00000197063"    "ENSG00000171786"    "ENSG00000073861" 
    ##             MA0592.2             MA0017.2             MA0800.1 
    ## "ENSMUSG00000024955"    "ENSG00000175745"    "ENSG00000163508" 
    ##             MA0802.1             MA1111.1 
    ##    "ENSG00000136535"    "ENSG00000185551"
