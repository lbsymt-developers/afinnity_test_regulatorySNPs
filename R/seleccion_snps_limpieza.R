library(dplyr)

clean_snp_list <- function(list1){
  snp_tmp <- purrr::transpose(list1)
  snp_tmp <- snp_tmp[["result"]]

  df_dummy <- data.frame()
  for(i in 1:length(snp_tmp)){
    a <- snp_tmp[[i]]
    df_dummy <- rbind(df_dummy, a)
  }
  return(df_dummy)
}

snp1 <- readRDS("snps_1_1000.rds") %>%
  clean_snp_list()
snp2 <- readRDS("snps_1001_2000.rds") %>%
  clean_snp_list()
snp3 <- readRDS("snps_2001_3000.rds") %>%
  clean_snp_list()
snp4 <- readRDS("snps_3001_5000.rds") %>%
  clean_snp_list()
snp5 <- readRDS("snps_5001_6741.rds") %>%
  clean_snp_list()

snps <- rbind(snp1, snp2, snp3, snp4, snp5)
saveRDS(snps, "all_snps_motifs_strong.rds")

table(snps$geneSymbol)
a <- unique(snps$SNP_id)


