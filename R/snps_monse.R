library(atSNP)

data(encode_library)
length(encode_motif)
encode_motif[1]
encode_motif[[1]]
GetIUPACSequence(encode_motif[[1]])
length(encode_motifinfo)
head(encode_motifinfo)

encode_motifinfo[names(encode_motif[1])]
data(jaspar_library)
jaspar_motif[[1]]
jaspar_motifinfo[names(jaspar_motif[1])]

snp_info1 <- LoadSNPData(snpids = c("rs1049673", "rs13032148", "rs1060743"),
                         genome.lib ="BSgenome.Hsapiens.UCSC.hg38",
                         snp.lib = "SNPlocs.Hsapiens.dbSNP155.GRCh38", half.window.size = 30, default.par = TRUE, mutation = FALSE)

# install.packages("../../../Downloads/SNPlocs.Hsapiens.dbSNP155.GRCh38_0.99.21.tar",
#                  repos = NULL, type = "source")

str(snp_info1)
data(encode_library)
encode_motifinfo[names(encode_motif)]

atsnp.scores <- ComputeMotifScore(encode_motif, snp_info1,
                                  ncores = 3)
atsnp.scores$snp.tbl
atsnp.scores$motif.scores

atsnp.result <- ComputePValues(motif.lib = encode_motif,
                               snp.info = snp_info1,
                               motif.scores = atsnp.scores$motif.scores,
                               ncores = 3, testing.mc=TRUE)
atsnp.result
head(atsnp.result[order(atsnp.result$pval_rank), c("snpid", "motif", "pval_ref", "pval_snp", "pval_rank")])

library(qvalue)
qval_rank = qvalue(atsnp.result$pval_rank, pi0=0.1)$qvalues
fdr <- p.adjust(atsnp.result$pval_snp)


atsnp.result = cbind(atsnp.result, qval_rank, fdr)
selected <- head(atsnp.result[order(atsnp.result$pval_rank),
                  c("snpid", "motif", "pval_ref", "pval_snp", "pval_rank", "fdr")], 10000)
selected <- selected[selected$fdr<=0.05,]
selected <- selected[order(selected$fdr),]
readr::write_csv(selected, "data/motivos_union_Monse.csv")


match.subseq_result <- MatchSubsequence(snp.tbl = atsnp.scores$snp.tbl,
                                        motif.scores = atsnp.result, motif.lib = encode_motif,
                                        snpids = c("rs13032148"),
                                        motifs = names(encode_motif), ncores = 1)
match.subseq_result[c("snpid", "motif", "IUPAC", "ref_match_seq", "snp_match_seq")]
match.seq <- dtMotifMatch(atsnp.scores$snp.tbl,
                          atsnp.scores$motif.scores,
                          snpids="rs13032148", motifs="HEY1_disc2",
                          motif.lib = encode_motif)
match.seq <- match.seq[!duplicated(match.seq$snpid),]
plotMotifMatch(match.seq,  motif.lib = encode_motif)


match.seq_2 <- dtMotifMatch(atsnp.scores$snp.tbl,
                            atsnp.scores$motif.scores,
                            snpids="rs405509", motifs="POU2F3_2",
                            motif.lib = encode_motif)
match.seq <- match.seq[!duplicated(match.seq$snpid),]
plotMotifMatch(match.seq_2,  motif.lib = encode_motif)
