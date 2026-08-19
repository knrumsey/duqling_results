# From duqling
data("sim_study_testfuncs")
duq <- process_sim_study(sim_study_testfuncs, scale_CRPS = TRUE)

data("sim_study_realdata")
duq_data <- process_sim_study(sim_study_realdata, scale_CRPS=TRUE)

duq_data_small <- filter_sim_study(duq_data, expr = !big_data_flag )
duq_data_big   <- filter_sim_study(duq_data, expr = big_data_flag & fast_method_flag)

#================================================================
#    FIGURES
#================================================================
path <- paste0("figs/main/Figure1.", extension)
filter_sim_study(duq, n_train=1000, NSR=0) %>%
  rankplot_sim_study("CRPS", title="CRPS Rank (n=1000, NSR=0)")
ggsave(path, height=5, width=8)

path <- paste0("figs/main/Figure2.", extension)
filter_sim_study(duq, n_train=1000, NSR=0.1) %>%
  rankplot_sim_study("CRPS", title="CRPS Rank (n=1000, NSR=0.1)")
ggsave(path, height=5, width=8)

path <- paste0("figs/main/Figure3.", extension)
custom <- function(xx){
  log10(pmin(1, pmax(0.001, xx)))
}
filter_sim_study(duq, n_train=1000, NSR=0) %>%
  heatmap_sim_study(metric="CRPS_med", y_scale_fun = custom,
                    colorbar_labels = list(breaks=log10(c(0.001, 0.01,  0.1, 0.56, 1)),
                                           labels=c("<0.001","0.01", "0.1", "0.56", "")),
                    color_scale="mako",
                    title=FALSE)
ggsave(path, height=8.5, width=8)

path <- paste0("figs/main/Figure4.", extension)
custom <- function(xx){
  log10(pmin(1, pmax(0.015, xx)))
}
filter_sim_study(duq, n_train=1000, NSR=0.1) %>%
  heatmap_sim_study(metric="CRPS_med", y_scale_fun = custom,
                    colorbar_labels = list(breaks=log10(c(0.13,  0.3, 0.56, 1)),
                                           labels=c("0.1", "0.3", "0.56", "")),
                    color_scale="mako",
                    title=FALSE)
ggsave(path, height=8.5, width=8)

path <- paste0("figs/main/Figure5a.", extension)
filter_sim_study(duq, n_train=1000, NSR=0) %>%
  paretoplot_sim_study(metric=c("CRPS_rel", "time_rel_log"),
                       epsilon=c(0.001, 0), upper_bound=c(1000, Inf), show_legend = FALSE,
                       title="CRPS & Speed (n=1000, NSR=0)")
ggsave(path, height=4.5, width=4.5)

path <- paste0("figs/main/Figure5b.", extension)
filter_sim_study(duq, n_train=1000, NSR=0.1) %>%
  paretoplot_sim_study(metric=c("CRPS_rel", "time_rel_log"),
                       epsilon=c(0.001, 0), upper_bound=c(1000, Inf), show_legend = FALSE,
                       title="(n=1000, NSR=0.1)", ylim=c(1, 2.5))
ggsave(path, height=4.5, width=4.5)

path <- paste0("figs/main/Figure6.", extension)
filter_sim_study(duq, n_train=500, NSR=0) %>%
  rankplot_sim_study("CRPS", title="CRPS Rank (n=500, NSR=0)")
ggsave(path, height=5, width=8)

path <- paste0("figs/main/Figure7.", extension)
filter_sim_study(duq, n_train=5000, NSR=0) %>%
  rankplot_sim_study("CRPS", title="CRPS Rank (n=5000, NSR=0)")
ggsave(path, height=5, width=8)

path <- paste0("figs/main/Figure6b.", extension)
filter_sim_study(duq, n_train=500, NSR=0) %>%
  paretoplot_sim_study(metric=c("CRPS_rel", "time_rel_log"),
                       epsilon=c(0.001, 0), upper_bound=c(1000, Inf), show_legend = FALSE,
                       title="CRPS & Speed (n=500, NSR=0)")
ggsave(path, height=4.5, width=4.5)

path <- paste0("figs/main/Figure7b.", extension)
filter_sim_study(duq, n_train=5000, NSR=0) %>%
  paretoplot_sim_study(metric=c("CRPS_rel", "time_rel_log"),
                       epsilon=c(0.001, 0), upper_bound=c(1000, Inf), show_legend = FALSE,
                       title="(n=5000, NSR=0)", ylim=c(1, 300))
ggsave(path, height=4.5, width=4.5)

### REALDATA FIGURES

path <- paste0("figs/main/Figure8.", extension)
filter_sim_study(duq_data_small) %>%
  rankplot_sim_study("CRPS", title="CRPS Rank (Small Datasets)")
ggsave(path, height=5, width=8)

path <- paste0("figs/main/Figure9.", extension)
filter_sim_study(duq_data_big) %>%
  rankplot_sim_study("CRPS", title="CRPS Rank (Big Datasets)")
ggsave(path, height=5, width=8)

path <- paste0("figs/main/Figure8_new.", extension)
filter_sim_study(duq_data_small) %>%
  perfprofile_sim_study("CRPS", title="CRPS Profile (Small Datasets)")
ggsave(path, height=5, width=8)

path <- paste0("figs/main/Figure9_new.", extension)
filter_sim_study(duq_data_big) %>%
  perfprofile_sim_study("CRPS", title="CRPS Profile (Big Datasets)")
ggsave(path, height=5, width=8)


path <- paste0("figs/main/Figure10a.", extension)
paretoplot_sim_study(duq_data_small, metric=c("CRPS_rel", "time_rel_log"),
                     epsilon=c(0.001, 0), upper_bound=c(1000, Inf), show_legend = FALSE,
                     title="CRPS & Speed (Small Datasets)", ylim=c(1, 35))
ggsave(path, height=4.5, width=4.5)

path <- paste0("figs/main/Figure10b.", extension)
paretoplot_sim_study(duq_data_big, metric=c("CRPS_rel", "time_rel_log"),
                     epsilon=c(0.001, 0), upper_bound=c(1000, Inf), show_legend = FALSE,
                     title="(Big Datasets)", ylim=c(1,50))
ggsave(path, height=4.5, width=4.5)

path <- paste0("figs/main/Figure11.", extension)
custom <- function(xx){
  log10(pmin(1, pmax(0.001, xx)))
}
filter_sim_study(duq_data_small) %>%
  heatmap_sim_study(metric="CRPS_med", y_scale_fun = custom,
                    colorbar_labels = list(breaks=log10(c(0.001, 0.01,  0.1, 0.56, 1)),
                                           labels=c("<0.001", "0.01", "0.1", "0.56", "")),
                    color_scale="mako",
                    title="Median CRPS (Small Datasets)")
ggsave(path, height=5.5, width=8)

path <- paste0("figs/main/Figure12.", extension)
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

