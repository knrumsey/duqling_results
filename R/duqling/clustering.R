library(dbscan)
library(ggrepel)
load("data/results_paper.Rda")
load("data/results_paper_data.Rda")

duq1 <- process_sim_study(results)
duq2 <- process_sim_study(results_all)
duq <- join_sim_study(duq1, duq2)
df <- duq$df



# ---- Step 1. Define scenario ----
scenario_cols <- c("id", "n_train", "NSR", "design_type", "replication")
scenario_cols <- intersect(scenario_cols, names(df))
df$scenario_id <- apply(df[, scenario_cols, drop = FALSE], 1, paste, collapse = "_")

# ---- Step 2. Rank CRPS within each scenario ----
df_rank <- df %>%
  group_by(scenario_id) %>%
  mutate(rank_crps = rank(CRPS, ties.method = "min")) %>%
  ungroup()

# ---- Step 3. Emulator × scenario matrix ----
mat <- df_rank %>%
  dplyr::select(method, scenario_id, rank_crps) %>%
  tidyr::pivot_wider(names_from = scenario_id, values_from = rank_crps)

rownames_mat <- mat$method
mat <- as.matrix(mat[,-1])
rownames(mat) <- rownames_mat

# ---- Step 4. Distance matrix ----
cor_mat <- cor(t(mat), method = "spearman", use = "pairwise.complete.obs")
dist_mat <- as.dist(1 - cor_mat)

# ---- Step 5. 2D MDS projection ----
mds <- cmdscale(dist_mat, k = 2, eig = TRUE)
mds_df <- data.frame(
  method = rownames(mds$points),
  Dim1 = mds$points[,1],
  Dim2 = mds$points[,2]
)

# ---- Step 6. DBSCAN clustering ----
set.seed(123)
db <- dbscan(mds_df[,c("Dim1","Dim2")], eps = 0.22, minPts = 1)
mds_df$cluster <- factor(ifelse(db$cluster == 0, "Noise", paste0("C", db$cluster)))

# ---- Step 7. Plot with hidden legend + markers + non-overlapping labels ----
# Figure out how many clusters you got
n_clusters <- length(unique(mds_df$cluster))

# Predefine a vector of shape codes with enough variety
shape_codes <- c(16, 15, 17, 18, 3, 8, 7, 4, 3, 0)  # add more if needed

ggplot(mds_df, aes(x = Dim1, y = Dim2, color = cluster, shape = cluster)) +
  geom_point(size = 3, stroke = 1) +
  ggrepel::geom_text_repel(aes(label = method), size = 3, max.overlaps = Inf) +
  theme_minimal(base_size = 13) +
  guides(color = "none", shape = "none") +
  scale_shape_manual(values = shape_codes[seq_len(n_clusters)]) +
  labs(
    title = "Emulator clusters",
    x = "MDS Dimension 1", y = "MDS Dimension 2"
  )
ggsave(filename="figs/main/cluster_emu.eps", height=4.5, width=4.5)
