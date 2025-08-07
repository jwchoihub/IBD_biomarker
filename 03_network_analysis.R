# 03_network_analysis.R
# Co-occurrence network (HC vs IBD; same for other group pairs)

library(qiime2R)
library(phyloseq)
library(igraph)
library(dplyr)
library(tidyr)

# 1) Load phyloseq
physeq <- qza_to_phyloseq(
  features = "table.qza",
  tree     = "rooted-tree-rep-seqs.qza",
  taxonomy = "taxonomy_silva.qza",
  metadata = "metadata.tsv"
)

# 2) Build relative-abundance matrix
otu_df <- as.data.frame(otu_table(physeq))
otu_rel <- sweep(otu_df, 2, colSums(otu_df), "/")
# filter by prevalence ≥1%
incid <- (otu_rel>0)*1
min_prev <- max(10, round(0.01*ncol(otu_rel)))
otu_filt <- otu_rel[rowSums(incid)>=min_prev, ]

# 3) Spearman correlation & network
cor_mat <- cor(t(otu_filt), method="spearman")
cor_mat[abs(cor_mat)<0.6] <- 0
net <- graph_from_adjacency_matrix(cor_mat, mode="lower", weighted=TRUE, diag=FALSE)
net <- delete_vertices(net, degree(net)==0)

# 4) Annotate nodes with taxonomy & abundance
tax_tbl <- as.data.frame(tax_table(physeq)) %>% rownames_to_column("ASV")
abund <- otu_rel %>%
  rownames_to_column("ASV") %>%
  pivot_longer(-ASV, names_to="SampleID", values_to="abund") %>%
  group_by(ASV) %>%
  summarise(rel_abund = sum(abund)) 

v_attr <- vertex_attr(net)$name %>%
  enframe(name=NULL, value="ASV") %>%
  left_join(tax_tbl, by="ASV") %>%
  left_join(abund, by="ASV")

net1 <- graph_from_data_frame(
  as_data_frame(net, what="edges"),
  directed=FALSE,
  vertices=as.data.frame(v_attr)
)
