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

snp_info1 <- LoadSNPData(snpids = c("rs5050", "rs616488", "rs11249433", "rs182799", "rs12565013", "rs11208590"),
                         genome.lib ="BSgenome.Hsapiens.UCSC.hg38",
                         snp.lib = "SNPlocs.Hsapiens.dbSNP155.GRCh38", half.window.size = 30, default.par = TRUE, mutation = FALSE)

