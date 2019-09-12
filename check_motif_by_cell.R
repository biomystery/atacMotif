
# test NEUROD1 ------------------------------------------------------------
test_tf <- "NEUROD1";names(test_tf) <- "ENSG00000162992"
require(chromVAR)
require(chromVARmotifs)
jaspar_motifs <- getJasparMotifs() # default species is human

test_pwm <- human_pwms_v1[grep(names(test_tf),names(human_pwms_v1))]

# load islet matrix -------------------------------------------------------

# matrix 1: binary matrix peak by cell 
require(Matrix)
mat_peak_celltype <- read.table("./test/Islet_123.binary_annotation_matrix.txt",header = T)

apply(mat_peak_celltype[,-1],2,sum)
quantile(apply(mat_peak_celltype[,-1],1,sum))


