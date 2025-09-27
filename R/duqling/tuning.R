library(duqling)
library(tidyverse)

load("data/results_paper.Rda")
load("data/results_lagp.Rda")

duq <- process_sim_study(results)
duq_lagp <- process_sim_study(res_lagp)
duq <- join_sim_study(duq, duq_lagp)

#filter_sim_study(duq, n_train=1000, NSR=0) %>%
#  paretoplot_sim_study(metric=c("CRPS_rel", "time_rel_log"),
#                       epsilon=c(0.001, 0), upper_bound=c(100, Inf), show_legend = FALSE,
#                       title="CRPS & Speed (n=1000, NSR=0)", ylim=c(1, 70))


# Ensure metrics exist
metrics <- c("CRPS_rel", "time_rel_log")
duq_sub <- filter_sim_study(duq, n_train = 1000, NSR = 0)
duq_sub <- ensure_metric(duq_sub, metrics, epsilon = c(0.001, 0), upper_bound = c(100, Inf))

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
  mutate(is_lagp = str_detect(method, lagp_regex),
         criterion = if_else(is_lagp, str_extract(method, "^[^_]+"), NA_character_),
         nbhd_size = if_else(is_lagp, as.integer(str_extract(method, "(?<=_)\\d+$")), NA_integer_))

# Prepare colors and shapes
criterion_levels <- c("nn", "alcray", "alc", "mspe")
nbhd_levels <- sort(unique(nbhd))

# Compute Pareto front: Lower is better on both axes
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

method_stats$criterion <- factor(method_stats$criterion, levels=c("nn", "alcray", "alc", "mspe"))



# Extract Pareto front points
pf <- method_stats %>% filter(pareto) %>% arrange(.data[[metrics[2]]])

# Plot
p <- ggplot(method_stats, aes_string(x = metrics[2], y = metrics[1])) +
  # Base layer: grey dots
  geom_point(data = filter(method_stats, !is_lagp),
             color = "grey70", shape = 16, size = 2) +

  # Highlight laGP variants
  geom_point(data = filter(method_stats, is_lagp),
             aes(shape = criterion, color = factor(nbhd_size)),
             size = 3) +

  # Pareto front
  geom_path(data = pf, aes_string(x = metrics[2], y = metrics[1]),
            color = "black", linewidth = 0.6, linetype = "dashed") +

  scale_color_brewer(palette = "Set1", name = "Neighborhood Size") +
  scale_shape_manual(values = c(nn = 16, alcray = 17, alc = 15, mspe = 18),
                     name = "Criterion") +

  labs(
    x = "log Relative Time",
    y = "Relative CRPS"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

# Show the plot
p
ggsave("figs/main/pareto_lagp.eps", width=8, height=4)
