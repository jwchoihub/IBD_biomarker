# 01_diversity.R
# Calculate alpha & beta diversity (HC vs IBD; also works for HC vs UC, HC vs CD, UC vs CD)

# 1) Load packages
library(qiime2R)
library(phyloseq)
library(picante)
library(ggplot2)
library(dplyr)


# 2) Import phyloseq object from QIIME2 artifacts
physeq <- qza_to_phyloseq(
  features  = "table_subsampled.qza",
  tree      = "rooted-tree-rep-seqs_subsampled.qza",
  taxonomy  = "taxonomy_silva.qza",
  metadata  = "metadata.tsv"
)

# 3) Alpha diversity
evenness <- read_qza("core-metrics-results-subsampled/evenness_vector.qza")$data
adiv <- data.frame(
  Evenness     = evenness,
  Observed     = estimate_richness(physeq, measures="Observed"),
  Shannon      = estimate_richness(physeq, measures="Shannon"),
  Faith_PD     = pd(samp = as.data.frame(t(otu_table(physeq))), tree = phy_tree(physeq))[,1],
  Group        = sample_data(physeq)$Group
)

adiv_long <- adiv %>%
  tidyr::gather(metric, value, Evenness, Observed, Shannon, Faith_PD)

ggplot(adiv_long, aes(x=Group, y=value, fill=Group)) +
  theme_bw() +
  geom_boxplot(outlier.shape=16, outlier.size=1.2) +
  facet_wrap(~ metric, scales="free", ncol=4) +
  labs(x="", y="Alpha diversity") +
  stat_compare_means(method="wilcox.test", label="p.format") +
  theme(text = element_text(size=14))

# 4) Beta diversity (Unweighted UniFrac PCoA)
metadata <- read_q2metadata("metadata.tsv")
uwunifrac <- read_qza("core-metrics-results-subsampled/unweighted_unifrac_pcoa_results.qza")$data
pc_df <- uwunifrac$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata, by="SampleID")
expl <- round(100 * uwunifrac$ProportionExplained[c("PC1","PC2")])

ggplot(pc_df, aes(x=PC1, y=PC2, color=Group)) +
  theme_bw() +
  geom_point(size=3.5) +
  stat_ellipse(linetype=2) +
  labs(
    x = paste0("PC1: ", expl["PC1"], "%"),
    y = paste0("PC2: ", expl["PC2"], "%"),
    title = "Unweighted UniFrac"
  ) +
  theme(plot.title = element_text(size=19, face="bold"))
