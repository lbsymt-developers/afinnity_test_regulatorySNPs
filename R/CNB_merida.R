load("data/resultados_mapeo_SNPs.rda")
snps <- readr::read_csv("data/LogFC_SNPs_celltype_grubman.csv")
snps <- snps[snps$celltype=="neuron",]

results <- data.frame(results)
ns <- c("n1", "n2", "n3", "n4", "n5", "n6")
newdf <- data.frame()
for(i in 1:length(ns)){
  pathway <- paste0("data/neuron_", ns[i],".csv")
  df <- readr::read_csv(pathway)
  df$subcluster <- ns[i]
  gene <- df[df$geneName %in% gen_sel, ]
  newdf <- rbind(newdf, gene)
}
sub_grubman <- dplyr::left_join(newdf, results,
                                 by = c("geneName" = "gene_name"))


morabito <- readr::read_csv("data/LogFC_SNPs_celltype_Morabito.csv")
morabito_ex <- morabito[morabito$celltype=="EX",]
morabito_inh <- morabito[morabito$celltype=="INH",]


# se seleccionaron los SNPs del gen SORL1

library(atSNP)
data(encode_library)
snp_info1 <- LoadSNPData(snpids = c("rs12272618", "rs2276412", "rs77819448"),
                         genome.lib ="BSgenome.Hsapiens.UCSC.hg38",
                         snp.lib = "SNPlocs.Hsapiens.dbSNP155.GRCh38", half.window.size = 30, default.par = TRUE, mutation = FALSE)
atsnp.scores <- ComputeMotifScore(encode_motif, snp_info1,
                                  ncores = 3)
atsnp.result <- ComputePValues(motif.lib = encode_motif,
                               snp.info = snp_info1,
                               motif.scores = atsnp.scores$motif.scores,
                               ncores = 3, testing.mc=TRUE)
library(qvalue)
qval_rank <- qvalue(atsnp.result$pval_rank, pi0=0.1)$qvalues
fdr <- p.adjust(atsnp.result$pval_snp)
atsnp.result <- cbind(atsnp.result, qval_rank, fdr)

match.subseq_result <- MatchSubsequence(snp.tbl = atsnp.scores$snp.tbl,
                                        motif.scores = atsnp.result, motif.lib = encode_motif,
                                        snpids = c("rs769446"),
                                        motifs = names(encode_motif), ncores = 1)
match.subseq_result[c("snpid", "motif", "IUPAC", "ref_match_seq", "snp_match_seq")]
match.seq <- dtMotifMatch(atsnp.scores$snp.tbl,
                          atsnp.scores$motif.scores,
                          snpids="rs769446", motifs="ELF1_disc3",
                          motif.lib = encode_motif)
match.seq <- match.seq[!duplicated(match.seq$snpid),]
plotMotifMatch(match.seq,  motif.lib = encode_motif)


match.seq_2 <- dtMotifMatch(atsnp.scores$snp.tbl,
                            atsnp.scores$motif.scores,
                            snpids="rs405509", motifs="POU2F3_2",
                            motif.lib = encode_motif)
match.seq <- match.seq[!duplicated(match.seq$snpid),]
plotMotifMatch(match.seq_2,  motif.lib = encode_motif)


