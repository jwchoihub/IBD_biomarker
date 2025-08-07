# 02_differential_analysis.R
# Differential abundance (HC vs IBD; same for HC vs UC, HC vs CD, UC vs CD)

library(qiime2R)
library(ANCOMBC)
library(Maaslin2)

# 1) Import HI vs BD phyloseq
physeq.HCIBD <- qza_to_phyloseq(
  features = "table.qza",
  tree      = "rooted-tree-rep-seqs.qza",
  taxonomy  = "taxonomy_silva.qza",
  metadata  = "metadata.tsv"
)

# 2) LEfSe
lda.res <- ldamarker(
  physeq.HCIBD,
  group    = "Health_status",
  method   = "relative",
  pvalue   = 0.05,
  normalize= TRUE
)
lda.genus <- subset(na.omit(lda.res), rank=="Genus" & LDAscore>2)
write.csv(lda.genus, "lefse_genus_HCIBD.csv")

# 3) ANCOM-BC
ancom.res <- ancombc2(
  data       = physeq.HCIBD,
  fix_formula= "Health_status",
  p_adj_method="BH", prv_cut=0.2, struc_zero=TRUE
)
res.ancombc <- as.data.frame(ancom.res$res)
write.csv(res.ancombc, "ancombc_genus_HCIBD.csv")

# 4) MaAsLin2
fit_adj <- Maaslin2(
  input_data, input_meta, "maaslin_adj_HCIBD",
  transform      = "NONE",
  fixed_effects  = c("Health_status","Sex","Age"),
  reference      = "Health_status,HC",
  normalization  = "NONE",
  correction     = "BH",
  standardize    = FALSE
)
res_maaslin <- fit_adj$results
write.csv(res_maaslin, "maaslin_genus_HCIBD.csv")
