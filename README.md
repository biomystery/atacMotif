jaspar stats
============

The `JASPAR motifs` were downloaded from
[Jaspar](http://jaspar.genereg.net/downloads/) choosing
`JASPAR2018_CORE_vertebrates_non-redundant` in `MEME` format. Then a
customized function `getEnsemblFromJaspar` [code](./aux_funs.R) was used
to get Ensembl ID from `Japsar id` via `uiportID` from Jaspar's record.

    # number of TF motifs in Jaspar 
    nrow(jaspar) 

    ## [1] 579

    # How many of them are protein complex? 
    length(grep(":",jaspar$V3)) 

    ## [1] 45

tfclass stats
=============

[TFclass
database](http://tfclass.bioinf.med.uni-goettingen.de/about.jsf) were
downloaded (2018-05-14), [parsed](./scripts/parseTfclassTTL.sh) and
[extacted](./scripts/tfclass.R) the human and mouse genes' Ensembl ID.

    # How many unique genus (i.e. TF genes) in parsed TFClass?
    length(unique(tfclass$merge$genus.id)) 

    ## [1] 1440

    # How many subfamilies in parsed TFClass?
    length(unique(tfclass$merge$subfamily.id))

    ## [1] 396

    # How many families in parsed TFClass?
    length(unique(tfclass$merge$family.id))

    ## [1] 110

jaspar-tfclass dictionary stats
===============================

`JASPAR` and `TFClass` were majorly [mapped](./cmpDB.R) by matching
`Ensembl ID` except `20` TFs by `TF Name` directly.

    # how many Jasper motif have annotation? 
    a=nrow(dic$jasparTOtfclass %>% 
      dplyr::select(jaspar.id,genus.id) %>% 
      group_by(jaspar.id)%>%
      summarise(genus.id.2 = paste0(genus.id,collapse = ';')))
    paste(a,nrow(jaspar),sep = "/")

    ## [1] "559/579"

    # which jaspar motif don't have annotation?
    noAnno <- subset(jaspar, V2 %in% setdiff(jaspar$V2,dic$jasparTOtfclass$jaspar.id))
    noAnno.2 <- dic$jasparTOensembl[setdiff(jaspar$V2,dic$jasparTOtfclass$jaspar.id)]
    knitr::kable(data.frame(jaspar.id=noAnno$V2,jaspar.name=noAnno$V3,
               ensembl.id= noAnno.2[noAnno$V2],
               url= paste0("http://jaspar2018.genereg.net/matrix/",noAnno$V2))%>%
      rowid_to_column('idx'))

<table>
<thead>
<tr class="header">
<th align="right">idx</th>
<th align="left">jaspar.id</th>
<th align="left">jaspar.name</th>
<th align="left">ensembl.id</th>
<th align="left">url</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="right">1</td>
<td align="left">MA0019.1</td>
<td align="left">Ddit3::Cebpa</td>
<td align="left">NULL;NULL</td>
<td align="left"><a href="http://jaspar2018.genereg.net/matrix/MA0019.1" class="uri">http://jaspar2018.genereg.net/matrix/MA0019.1</a></td>
</tr>
<tr class="even">
<td align="right">2</td>
<td align="left">MA0089.1</td>
<td align="left">MAFG::NFE2L1</td>
<td align="left">NULL;NULL</td>
<td align="left"><a href="http://jaspar2018.genereg.net/matrix/MA0089.1" class="uri">http://jaspar2018.genereg.net/matrix/MA0089.1</a></td>
</tr>
<tr class="odd">
<td align="right">3</td>
<td align="left">MA0109.1</td>
<td align="left">HLTF</td>
<td align="left">NULL</td>
<td align="left"><a href="http://jaspar2018.genereg.net/matrix/MA0109.1" class="uri">http://jaspar2018.genereg.net/matrix/MA0109.1</a></td>
</tr>
<tr class="even">
<td align="right">4</td>
<td align="left">MA0111.1</td>
<td align="left">Spz1</td>
<td align="left">ENSMUSG00000046957</td>
<td align="left"><a href="http://jaspar2018.genereg.net/matrix/MA0111.1" class="uri">http://jaspar2018.genereg.net/matrix/MA0111.1</a></td>
</tr>
<tr class="odd">
<td align="right">5</td>
<td align="left">MA0149.1</td>
<td align="left">EWSR1-FLI1</td>
<td align="left">NULL;NULL</td>
<td align="left"><a href="http://jaspar2018.genereg.net/matrix/MA0149.1" class="uri">http://jaspar2018.genereg.net/matrix/MA0149.1</a></td>
</tr>
<tr class="even">
<td align="right">6</td>
<td align="left">MA0611.1</td>
<td align="left">Dux</td>
<td align="left">logical(0)</td>
<td align="left"><a href="http://jaspar2018.genereg.net/matrix/MA0611.1" class="uri">http://jaspar2018.genereg.net/matrix/MA0611.1</a></td>
</tr>
<tr class="odd">
<td align="right">7</td>
<td align="left">MA0619.1</td>
<td align="left">LIN54</td>
<td align="left">NULL</td>
<td align="left"><a href="http://jaspar2018.genereg.net/matrix/MA0619.1" class="uri">http://jaspar2018.genereg.net/matrix/MA0619.1</a></td>
</tr>
<tr class="even">
<td align="right">8</td>
<td align="left">MA0621.1</td>
<td align="left">mix-a</td>
<td align="left">NULL</td>
<td align="left"><a href="http://jaspar2018.genereg.net/matrix/MA0621.1" class="uri">http://jaspar2018.genereg.net/matrix/MA0621.1</a></td>
</tr>
<tr class="odd">
<td align="right">9</td>
<td align="left">MA0629.1</td>
<td align="left">Rhox11</td>
<td align="left">ENSMUSG00000051038</td>
<td align="left"><a href="http://jaspar2018.genereg.net/matrix/MA0629.1" class="uri">http://jaspar2018.genereg.net/matrix/MA0629.1</a></td>
</tr>
<tr class="even">
<td align="right">10</td>
<td align="left">MA0633.1</td>
<td align="left">Twist2</td>
<td align="left">ENSMUSG00000007805</td>
<td align="left"><a href="http://jaspar2018.genereg.net/matrix/MA0633.1" class="uri">http://jaspar2018.genereg.net/matrix/MA0633.1</a></td>
</tr>
<tr class="odd">
<td align="right">11</td>
<td align="left">MA0637.1</td>
<td align="left">CENPB</td>
<td align="left">ENSG00000125817</td>
<td align="left"><a href="http://jaspar2018.genereg.net/matrix/MA0637.1" class="uri">http://jaspar2018.genereg.net/matrix/MA0637.1</a></td>
</tr>
<tr class="even">
<td align="right">12</td>
<td align="left">MA0643.1</td>
<td align="left">Esrrg</td>
<td align="left">ENSMUSG00000026610</td>
<td align="left"><a href="http://jaspar2018.genereg.net/matrix/MA0643.1" class="uri">http://jaspar2018.genereg.net/matrix/MA0643.1</a></td>
</tr>
<tr class="odd">
<td align="right">13</td>
<td align="left">MA0659.1</td>
<td align="left">MAFG</td>
<td align="left">ENSG00000197063</td>
<td align="left"><a href="http://jaspar2018.genereg.net/matrix/MA0659.1" class="uri">http://jaspar2018.genereg.net/matrix/MA0659.1</a></td>
</tr>
<tr class="even">
<td align="right">14</td>
<td align="left">MA0048.2</td>
<td align="left">NHLH1</td>
<td align="left">ENSG00000171786</td>
<td align="left"><a href="http://jaspar2018.genereg.net/matrix/MA0048.2" class="uri">http://jaspar2018.genereg.net/matrix/MA0048.2</a></td>
</tr>
<tr class="odd">
<td align="right">15</td>
<td align="left">MA0690.1</td>
<td align="left">TBX21</td>
<td align="left">ENSG00000073861</td>
<td align="left"><a href="http://jaspar2018.genereg.net/matrix/MA0690.1" class="uri">http://jaspar2018.genereg.net/matrix/MA0690.1</a></td>
</tr>
<tr class="even">
<td align="right">16</td>
<td align="left">MA0592.2</td>
<td align="left">Esrra</td>
<td align="left">ENSMUSG00000024955</td>
<td align="left"><a href="http://jaspar2018.genereg.net/matrix/MA0592.2" class="uri">http://jaspar2018.genereg.net/matrix/MA0592.2</a></td>
</tr>
<tr class="odd">
<td align="right">17</td>
<td align="left">MA0017.2</td>
<td align="left">NR2F1</td>
<td align="left">ENSG00000175745</td>
<td align="left"><a href="http://jaspar2018.genereg.net/matrix/MA0017.2" class="uri">http://jaspar2018.genereg.net/matrix/MA0017.2</a></td>
</tr>
<tr class="even">
<td align="right">18</td>
<td align="left">MA0800.1</td>
<td align="left">EOMES</td>
<td align="left">ENSG00000163508</td>
<td align="left"><a href="http://jaspar2018.genereg.net/matrix/MA0800.1" class="uri">http://jaspar2018.genereg.net/matrix/MA0800.1</a></td>
</tr>
<tr class="odd">
<td align="right">19</td>
<td align="left">MA0802.1</td>
<td align="left">TBR1</td>
<td align="left">ENSG00000136535</td>
<td align="left"><a href="http://jaspar2018.genereg.net/matrix/MA0802.1" class="uri">http://jaspar2018.genereg.net/matrix/MA0802.1</a></td>
</tr>
<tr class="even">
<td align="right">20</td>
<td align="left">MA1111.1</td>
<td align="left">NR2F2</td>
<td align="left">ENSG00000185551</td>
<td align="left"><a href="http://jaspar2018.genereg.net/matrix/MA1111.1" class="uri">http://jaspar2018.genereg.net/matrix/MA1111.1</a></td>
</tr>
</tbody>
</table>
