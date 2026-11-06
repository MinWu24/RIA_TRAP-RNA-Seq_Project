library(edgeR)

RIA_PA= read.table("RIA_OP_PA_matrix.csv", sep = ",", row.name = "Gene", header=TRUE)
str(RIA_PA)

RIA_PA = as.matrix(RIA_PA)
str(RIA_PA)

group =  factor(c(1,1,1,1,1,2,2,2,2,2))

DGEL.RIA_PA = DGEList(counts = RIA_PA, group = group)
str(DGEL.RIA_PA)

keep <- filterByExpr(DGEL.RIA_PA)
table(keep)
DGEL.RIA_PA <- DGEL.RIA_PA[keep,,keep.lib.size = FALSE]

cNorm.RIA_PA = calcNormFactors(DGEL.RIA_PA)
str(cNorm.RIA_PA)

ComDis.RIA_PA = estimateCommonDisp(cNorm.RIA_PA)
str(ComDis.RIA_PA)

TagDis.RIA_PA= estimateTagwiseDisp(ComDis.RIA_PA)
str(TagDis.RIA_PA)

et.RIA_PA = exactTest(TagDis.RIA_PA)
str(et.RIA_PA)

topTags(et.RIA_PA, n=100, sort.by="logFC")

RIA_PA.all = topTags(et.RIA_PA, n=90000, sort.by="logFC")
str(RIA_PA.all )

write.table(RIA_PA.all, file = "DE_RIA_OP_PA.csv", append = FALSE, quote = FALSE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE)

y.RIA_PA <- cpm(cNorm.RIA_PA, prior.count=2, log=TRUE)
head(y.RIA_PA)
