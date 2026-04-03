data("sim_study_testfuncs")
duq <- process_sim_study(sim_study_testfuncs, scale_CRPS = TRUE)

load("data/results_lagp.Rda")
duq_lagp <- process_sim_study(res_lagp)

load("data/results_lagp_smallnug.Rda")
duq_lagp_smallnug <- res_lagp
duq_lagp_smallnug$df$method <- paste0("small_", duq_lagp_smallnug$df$method)

duq <- join_sim_study(duq, duq_lagp)
duq <- join_sim_study(duq, duq_lagp_smallnug)

metrics <- c("CRPS_rel", "time_rel_log")
duq_sub <- filter_sim_study(duq, n_train = 1000, NSR = 0)
duq_sub <- ensure_metric(duq_sub, metrics, epsilon = c(0.001, 0), upper_bound = c(1000, Inf))

df <- duq_sub$df

# Compute mean metric per method
method_stats <- df %>%
  group_by(method) %>%
  summarise(across(all_of(metrics), mean, na.rm = TRUE), .groups = "drop")


# Identify laGP variants
criteria <- c("nn", "alcray", "alc", "mspe")
nbhd <- c(25, 50, 100, 200)
lagp_regex <- paste0("^(", paste(criteria, collapse = "|"), ")_(", paste(nbhd, collapse = "|"), ")$")

method_stats <- method_stats %>%
  mutate(
    nugget = if_else(str_detect(method, "^small_"), "small", "default"),
    method_clean = str_remove(method, "^small_"),
    is_lagp = str_detect(method_clean, lagp_regex),
    criterion = if_else(is_lagp, str_extract(method_clean, "^[^_]+"), NA_character_),
    nbhd_size = if_else(is_lagp, as.integer(str_extract(method_clean, "(?<=_)\\d+$")), NA_integer_),
    shape_group = if_else(is_lagp,
                          paste0(criterion, "_", nugget),
                          NA_character_)
  )

# Clean factor ordering
shape_levels <- c(
  "nn_default", "alcray_default", "alc_default", "mspe_default",
  "nn_small",   "alcray_small",   "alc_small",   "mspe_small"
)

method_stats$shape_group <- factor(method_stats$shape_group, levels = shape_levels)

# Pareto front (unchanged)
is_dominated <- function(i, df, xcol, ycol) {
  x_i <- df[[xcol]][i]
  y_i <- df[[ycol]][i]
  other <- df[-i, ]
  any(other[[xcol]] <= x_i & other[[ycol]] <= y_i &
        (other[[xcol]] < x_i | other[[ycol]] < y_i))
}

method_stats$pareto <- !sapply(seq_len(nrow(method_stats)), function(i) {
  is_dominated(i, method_stats, xcol = metrics[2], ycol = metrics[1])
})

pf <- method_stats %>% filter(pareto) %>% arrange(.data[[metrics[2]]])

# Plot
p <- ggplot(method_stats, aes_string(x = metrics[2], y = metrics[1])) +
  # Base layer
  geom_point(data = filter(method_stats, !is_lagp),
             color = "grey70", shape = 16, size = 2) +

  # laGP points
  geom_point(data = filter(method_stats, is_lagp),
             aes(shape = shape_group, color = factor(nbhd_size)),
             size = 3) +

  # Pareto front
  geom_path(data = pf,
            aes_string(x = metrics[2], y = metrics[1]),
            color = "black", linewidth = 0.6, linetype = "dashed") +

  scale_color_brewer(palette = "Set1", name = "Neighborhood Size") +

  scale_shape_manual(
    values = c(
      nn_default = 16,
      alcray_default = 17,
      alc_default = 15,
      mspe_default = 18,
      nn_small = 1,
      alcray_small = 2,
      alc_small = 0,
      mspe_small = 5
    ),
    labels = c(
      nn_default = "nn",
      alcray_default = "alcray",
      alc_default = "alc",
      mspe_default = "mspe",
      nn_small = "nn (dynamic)",
      alcray_small = "alcray (dynamic)",
      alc_small = "alc (dynamic)",
      mspe_small = "mspe (dynamic)"
    ),
    name = "Criterion"
  ) +

  labs(
    x = "log Relative Time\n← Better     Worse →",
    y = "Relative CRPS\n← Better     Worse →"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

p
ggsave(paste0("figs/main/pareto_lagp_combo.", extension), width=8, height=5)


