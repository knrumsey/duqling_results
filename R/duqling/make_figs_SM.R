load("data/results_paper.Rda") # -> results
duq <- process_sim_study(results, scale_CRPS = TRUE)

## n = 1000, NSR=0
path <- "figs/SM/FigureSM1a.eps"
filter_sim_study(duq, n_train=1000, NSR=0) %>%
  rankplot_sim_study("RMSE", title="RMSE Rank (n=1000, NSR=0)")
ggsave(path, height=4, width=8)

path <- "figs/SM/FigureSM2a.eps"
filter_sim_study(duq, n_train=1000, NSR=0) %>%
  rankplot_sim_study("time", title="Time Rank (n=1000, NSR=0)")
ggsave(path, height=4, width=8)

path <- "figs/SM/FigureSM3a.eps"
custom <- function(xx){
  log10(pmin(1, pmax(0.001, xx)))
}
filter_sim_study(duq, n_train=1000, NSR=0) %>%
  heatmap_sim_study(metric="CRPS_Q3", y_scale_fun = custom,
                    colorbar_labels = list(breaks=log10(c(0.001, 0.01,  0.1, 0.56, 1)),
                                           labels=c("<0.001","0.01", "0.1", "0.56", "")),
                    color_scale="mako",
                    title="CRPS Q3 (n=1000, NSR=0)")
ggsave(path, height=7.5, width=8)

path <- "figs/SM/FigureSM4a.eps"
filter_sim_study(duq, n_train=1000, NSR=0) %>%
  heatmap_sim_study(metric="CRPS", y_scale_fun = custom,
                    colorbar_labels = list(breaks=log10(c(0.001, 0.01,  0.1, 0.56, 1)),
                                           labels=c("<0.001","0.01", "0.1", "0.56", "")),
                    color_scale="mako",
                    title="Mean CRPS (n=1000, NSR=0)")
ggsave(path, height=7.5, width=8)

path <- "figs/SM/FigureSM5a.eps"
filter_sim_study(duq, n_train=1000, NSR=0) %>%
  heatmap_sim_study(metric="FVU", y_scale_fun = custom,
                    colorbar_labels = list(breaks=log10(c(0.001, 0.01,  0.1, 0.5, 1)),
                                           labels=c("<0.001","0.01", "0.1", "0.5", "")),
                    color_scale="mako",
                    title="FVU (n=1000, NSR=0)")
ggsave(path, height=7.5, width=8)


## n = 1000, NSR=0.1
path <- "figs/SM/FigureSM1b.eps"
filter_sim_study(duq, n_train=1000, NSR=0.1) %>%
  rankplot_sim_study("RMSE", title="RMSE Rank (n=1000, NSR=0.1)")
ggsave(path, height=4, width=8)

path <- "figs/SM/FigureSM2b.eps"
filter_sim_study(duq, n_train=1000, NSR=0.1) %>%
  rankplot_sim_study("time", title="Time Rank (n=1000, NSR=0.1)")
ggsave(path, height=4, width=8)

path <- "figs/SM/FigureSM3b.eps"
custom <- function(xx){
  log10(pmin(1, pmax(0.001, xx)))
}
filter_sim_study(duq, n_train=1000, NSR=0.1) %>%
  heatmap_sim_study(metric="CRPS_Q3", y_scale_fun = custom,
                    colorbar_labels = list(breaks=log10(c(0.001, 0.01,  0.1, 0.56, 1)),
                                           labels=c("<0.001","0.01", "0.1", "0.56", "")),
                    color_scale="mako",
                    title="CRPS Q3 (n=1000, NSR=0.1)")
ggsave(path, height=7.5, width=8)

path <- "figs/SM/FigureSM4b.eps"
filter_sim_study(duq, n_train=1000, NSR=0.1) %>%
  heatmap_sim_study(metric="CRPS", y_scale_fun = custom,
                    colorbar_labels = list(breaks=log10(c(0.001, 0.01,  0.1, 0.56, 1)),
                                           labels=c("<0.001","0.01", "0.1", "0.56", "")),
                    color_scale="mako",
                    title="Mean CRPS (n=1000, NSR=0.1)")
ggsave(path, height=7.5, width=8)

path <- "figs/SM/FigureSM5b.eps"
filter_sim_study(duq, n_train=1000, NSR=0.1) %>%
  heatmap_sim_study(metric="FVU", y_scale_fun = custom,
                    colorbar_labels = list(breaks=log10(c(0.001, 0.01,  0.1, 0.5, 1)),
                                           labels=c("<0.001","0.01", "0.1", "0.5", "")),
                    color_scale="mako",
                    title="FVU (n=1000, NSR=0.1)")
ggsave(path, height=7.5, width=8)


## n = 500, NSR=0.1
path <- "figs/SM/FigureSM1c.eps"
filter_sim_study(duq, n_train=500, NSR=0.1) %>%
  rankplot_sim_study("RMSE", title="RMSE Rank (n=500, NSR=0.1)")
ggsave(path, height=4, width=8)

path <- "figs/SM/FigureSM2c.eps"
filter_sim_study(duq, n_train=500, NSR=0.1) %>%
  rankplot_sim_study("time", title="Time Rank (n=500, NSR=0.1)")
ggsave(path, height=4, width=8)

path <- "figs/SM/FigureSM3c.eps"
custom <- function(xx){
  log10(pmin(1, pmax(0.001, xx)))
}
filter_sim_study(duq, n_train=500, NSR=0.1) %>%
  heatmap_sim_study(metric="CRPS_Q3", y_scale_fun = custom,
                    colorbar_labels = list(breaks=log10(c(0.001, 0.01,  0.1, 0.56, 1)),
                                           labels=c("<0.001","0.01", "0.1", "0.56", "")),
                    color_scale="mako",
                    title="CRPS Q3 (n=500, NSR=0.1)")
ggsave(path, height=7.5, width=8)

path <- "figs/SM/FigureSM4c.eps"
filter_sim_study(duq, n_train=500, NSR=0.1) %>%
  heatmap_sim_study(metric="CRPS", y_scale_fun = custom,
                    colorbar_labels = list(breaks=log10(c(0.001, 0.01,  0.1, 0.56, 1)),
                                           labels=c("<0.001","0.01", "0.1", "0.56", "")),
                    color_scale="mako",
                    title="Mean CRPS (n=500, NSR=0.1)")
ggsave(path, height=7.5, width=8)

path <- "figs/SM/FigureSM5c.eps"
filter_sim_study(duq, n_train=500, NSR=0.1) %>%
  heatmap_sim_study(metric="FVU", y_scale_fun = custom,
                    colorbar_labels = list(breaks=log10(c(0.001, 0.01,  0.1, 0.5, 1)),
                                           labels=c("<0.001","0.01", "0.1", "0.5", "")),
                    color_scale="mako",
                    title="FVU (n=500, NSR=0.1)")
ggsave(path, height=7.5, width=8)

path <- "figs/SM/FigureSM6c.eps"
filter_sim_study(duq, n_train=500, NSR=0.1) %>%
  rankplot_sim_study("CRPS", title="RMSE Rank (n=500, NSR=0.1)")
ggsave(path, height=4, width=8)

path <- "figs/SM/FigureSM7c.eps"
filter_sim_study(duq, n_train=500, NSR=0.1) %>%
  heatmap_sim_study(metric="CRPS_med", y_scale_fun = custom,
                    colorbar_labels = list(breaks=log10(c(0.001, 0.01,  0.1, 0.5, 1)),
                                           labels=c("<0.001","0.01", "0.1", "0.5", "")),
                    color_scale="mako",
                    title="Median CRPS (n=500, NSR=0.1)")
ggsave(path, height=7.5, width=8)

## n = 5000, NSR=0
path <- "figs/SM/FigureSM1d.eps"
filter_sim_study(duq, n_train=5000, NSR=0) %>%
  rankplot_sim_study("RMSE", title="RMSE Rank (n=500, NSR=0)")
ggsave(path, height=4, width=8)

path <- "figs/SM/FigureSM2d.eps"
filter_sim_study(duq, n_train=5000, NSR=0) %>%
  rankplot_sim_study("time", title="Time Rank (n=500, NSR=0)")
ggsave(path, height=4, width=8)

path <- "figs/SM/FigureSM3d.eps"
custom <- function(xx){
  log10(pmin(1, pmax(0.001, xx)))
}
filter_sim_study(duq, n_train=5000, NSR=0) %>%
  heatmap_sim_study(metric="CRPS_Q3", y_scale_fun = custom,
                    colorbar_labels = list(breaks=log10(c(0.001, 0.01,  0.1, 0.56, 1)),
                                           labels=c("<0.001","0.01", "0.1", "0.56", "")),
                    color_scale="mako",
                    title="CRPS Q3 (n=500, NSR=0)")
ggsave(path, height=7.5, width=8)

path <- "figs/SM/FigureSM4d.eps"
filter_sim_study(duq, n_train=5000, NSR=0) %>%
  heatmap_sim_study(metric="CRPS", y_scale_fun = custom,
                    colorbar_labels = list(breaks=log10(c(0.001, 0.01,  0.1, 0.56, 1)),
                                           labels=c("<0.001","0.01", "0.1", "0.56", "")),
                    color_scale="mako",
                    title="Mean CRPS (n=500, NSR=0)")
ggsave(path, height=7.5, width=8)

path <- "figs/SM/FigureSM5d.eps"
filter_sim_study(duq, n_train=5000, NSR=0) %>%
  heatmap_sim_study(metric="FVU", y_scale_fun = custom,
                    colorbar_labels = list(breaks=log10(c(0.001, 0.01,  0.1, 0.5, 1)),
                                           labels=c("<0.001","0.01", "0.1", "0.5", "")),
                    color_scale="mako",
                    title="FVU (n=500, NSR=0)")
ggsave(path, height=7.5, width=8)
