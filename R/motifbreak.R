BiocManager::install("motifbreakR")
library(motifbreakR)

library(BSgenome)
available.SNPs()
snps <- readr::read_csv("../Mapping_nonconding/data/SNPs_differentDataset_uniques.csv")
# require(GenomicRanges)
snps <- snps$RefSNP_id

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
                       BPPARAM = BiocParallel::MulticoreParam(7))

saveRDS(results, "results/motifs_all_AD_SNPs.rds")
results <- readRDS("results/motifs_all_AD_SNPs.rds")

find_critical_SNPs <- function(results){
  suppressMessages(library(dplyr))
  df <- results
  df <- df[df$effect == "strong", ]
  snps <- unique(names(df))

  future::plan("multisession")
  aa <- purrr::map(snps, purrr::safely(function(x){
    print(x)
    df_tmp <- df[names(df) %in% x,] %>%
      motifbreakR::calculatePvalue() %>%
      data.frame() %>%
      filter(Refpvalue <= 0.05)
  }), .progress = TRUE)



  df_snps <- data.frame()
  for(i in 1:length(snps)){
    print(i)
    df_tmp <- df[names(df) %in% snps[i],] %>%
      motifbreakR::calculatePvalue() %>%
      data.frame() %>%
      filter(Refpvalue <= 0.05)
    df_snps <- rbind(df_tmp, df_snps)
  }
  return(df_snps)
}


results_snps <- find_critical_SNPs(results)

r_strong <- results[results$effect=="strong"]
rs1001158 <- results[names(results) %in% "rs1001158:A"]
rs1001158
rs1001158 <- calculatePvalue(rs1001158)
rs1001158[rs1001158$Altpvalue<=0.05]


saveRDS(results_pvalue, "results/motifs_all_AD_SNPs_PVALUE.rds")
rs13032148 <- results[names(results) %in% "rs13032148"]

rs13032148 <- calculatePvalue(rs13032148)


library(BSgenome.Hsapiens.UCSC.hg38)

plotMB(results = results, rsid = "rs1001158:A", effect = "strong")


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


