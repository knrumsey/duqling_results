load("data/results_paper.Rda") # -> results
duq <- process_sim_study(results, scale_CRPS = TRUE)

#================================================================
#    FIGURES
#================================================================
path <- "figs/main/Figure1.eps"
filter_sim_study(duq, n_train=1000, NSR=0) %>%
  rankplot_sim_study("CRPS", title="CRPS Rank (n=1000, NSR=0)")
ggsave(path, height=4, width=8)

path <- "figs/main/Figure2.eps"
filter_sim_study(duq, n_train=1000, NSR=0.1) %>%
  rankplot_sim_study("CRPS", title="CRPS Rank (n=1000, NSR=0.1)")
ggsave(path, height=4, width=8)

path <- "figs/main/Figure3.eps"
custom <- function(xx){
  log10(pmin(1, pmax(0.001, xx)))
}
filter_sim_study(duq, n_train=1000, NSR=0) %>%
  heatmap_sim_study(metric="CRPS_med", y_scale_fun = custom,
                    colorbar_labels = list(breaks=log10(c(0.001, 0.01,  0.1, 0.56, 1)),
                                           labels=c("<0.001","0.01", "0.1", "0.56", "")),
                    color_scale="mako",
                    title=FALSE)
ggsave(path, height=7.5, width=8)

path <- "figs/main/Figure4.eps"
custom <- function(xx){
  log10(pmin(1, pmax(0.015, xx)))
}
filter_sim_study(duq, n_train=1000, NSR=0.1) %>%
  heatmap_sim_study(metric="CRPS_med", y_scale_fun = custom,
                    colorbar_labels = list(breaks=log10(c(0.015,  0.1, 0.56, 1)),
                                           labels=c("<0.01", "0.1", "0.56", "")),
                    color_scale="mako",
                    title=FALSE)
ggsave(path, height=7.5, width=8)


path <- "figs/main/Figure5a.eps"
filter_sim_study(duq, n_train=1000, NSR=0) %>%
  paretoplot_sim_study(metric=c("CRPS_rel", "time_rel_log"),
                       epsilon=c(0.001, 0), upper_bound=c(100, Inf), show_legend = FALSE,
                       title="CRPS & Speed (n=1000, NSR=0)", ylim=c(1, 70))
ggsave(path, height=4.5, width=4.5)

path <- "figs/main/Figure5b.eps"
filter_sim_study(duq, n_train=1000, NSR=0.1) %>%
  paretoplot_sim_study(metric=c("CRPS_rel", "time_rel_log"),
                       epsilon=c(0.001, 0), upper_bound=c(100, Inf), show_legend = FALSE,
                       title="(n=1000, NSR=0.1)", ylim=c(1, 20))
ggsave(path, height=4.5, width=4.5)

path <- "figs/main/Figure6.eps"
filter_sim_study(duq, n_train=500, NSR=0) %>%
  rankplot_sim_study("CRPS", title="CRPS Rank (n=500, NSR=0)")
ggsave(path, height=4, width=8)

path <- "figs/main/Figure7.eps"
filter_sim_study(duq, n_train=5000, NSR=0) %>%
  rankplot_sim_study("CRPS", title="CRPS Rank (n=5000, NSR=0)")
ggsave(path, height=4, width=8)

path <- "figs/main/Figure6b.eps"
filter_sim_study(duq, n_train=500, NSR=0) %>%
  paretoplot_sim_study(metric=c("CRPS_rel", "time_rel_log"),
                       epsilon=c(0.001, 0), upper_bound=c(100, Inf), show_legend = FALSE,
                       title="CRPS & Speed (n=500, NSR=0)", ylim=c(1, 70))
ggsave(path, height=4.5, width=4.5)

path <- "figs/main/Figure7b.eps"
filter_sim_study(duq, n_train=5000, NSR=0) %>%
  paretoplot_sim_study(metric=c("CRPS_rel", "time_rel_log"),
                       epsilon=c(0.001, 0), upper_bound=c(100, Inf), show_legend = FALSE,
                       title="(n=5000, NSR=0)", ylim=c(1, 70))
ggsave(path, height=4.5, width=4.5)


path <- "figs/main/Figure8.eps"
load("data/results_paper_data.Rda")
duq_data <- process_sim_study(results_smalldata)
filter_sim_study(duq_data) %>%
  rankplot_sim_study("CRPS", title="CRPS Rank (Small Datasets)")
ggsave(path, height=4, width=8)

path <- "figs/main/Figure9.eps"
duq_data_big <- process_sim_study(results_bigdata)
filter_sim_study(duq_data_big) %>%
  rankplot_sim_study("CRPS_max", title="CRPS Rank (Big Datasets)")
ggsave(path, height=4, width=8)


path <- "figs/main/Figure10a.eps"
paretoplot_sim_study(duq_data, metric=c("CRPS_rel", "time_rel_log"),
                     epsilon=c(0.001, 0), upper_bound=c(100, Inf), show_legend = FALSE,
                     title="CRPS & Speed (Small Datasets)", ylim=c(1, 35))
ggsave(path, height=4.5, width=4.5)

path <- "figs/main/Figure10b.eps"
paretoplot_sim_study(duq_data_big, metric=c("CRPS_rel", "time_rel_log"),
                     epsilon=c(0.001, 0), upper_bound=c(100, Inf), show_legend = FALSE,
                     title="(Big Datasets)", ylim=c(1, 35))
ggsave(path, height=4.5, width=4.5)

path <- "figs/main/Figure11.eps"
custom <- function(xx){
  log10(pmin(1, pmax(0.001, xx)))
}
filter_sim_study(duq_data) %>%
  heatmap_sim_study(metric="CRPS_med", y_scale_fun = custom,
                    colorbar_labels = list(breaks=log10(c(0.001, 0.01,  0.1, 0.56, 1)),
                                           labels=c("<0.001", "0.01", "0.1", "0.56", "")),
                    color_scale="mako",
                    title="Median CRPS (Small Datasets)")
ggsave(path, height=5.5, width=8)

path <- "figs/main/Figure12.eps"
custom <- function(xx){
  log10(pmin(1, pmax(0.001, xx)))
}
filter_sim_study(duq_data_big) %>%
  heatmap_sim_study(metric="CRPS_med", y_scale_fun = custom,
                    colorbar_labels = list(breaks=log10(c(0.001, 0.01, 0.1, 0.56, 1)),
                                           labels=c("<0.001", "0.01", "0.1", "0.56", "")),
                    color_scale="mako",
                    title="Median CRPS (Big Datasets)")
ggsave(path, height=5, width=8)

