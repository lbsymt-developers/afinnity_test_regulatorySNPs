BiocManager::install("motifbreakR")
library(motifbreakR)

library(BSgenome)
available.SNPs()
snps <- c("rs1049673", "rs13032148", "rs1060743")
snps.mb <- motifbreakR::snps.from.rsid(rsid = snps,
                          dbSNP = SNPlocs.Hsapiens.dbSNP155.GRCh38::SNPlocs.Hsapiens.dbSNP155.GRCh38,
                          search.genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)

data(motifbreakR_motif)
motifbreakR_motif

results <- motifbreakR::motifbreakR(snpList = snps.mb, filterp = TRUE,
                       pwmList = motifbreakR_motif,
                       threshold = 1e-4,
                       method = "ic",
                       bkg = c(A=0.25,
                               C=0.25,
                               G=0.25,
                               T=0.25),
                       BPPARAM = BiocParallel::SerialParam())

rs13032148 <- results[names(results) %in% "rs13032148"]
rs13032148
rs13032148 <- calculatePvalue(rs13032148)
rs13032148

library(BSgenome.Hsapiens.UCSC.hg38)

plotMB(results = results, rsid = "rs13032148", effect = "strong")


rs1060743 <- results[names(results) %in% "rs1060743"]
rs1060743
rs1060743 <- calculatePvalue(rs1060743)
rs1060743

plotMB(results = results, rsid = "rs1060743", effect = "strong")


rs1049673 <- results[names(results) %in% "rs1049673"]
rs1049673
rs1049673 <- calculatePvalue(rs1049673)
rs1049673

plotMB(results = results, rsid = "rs1049673", effect = "strong")


