#=================================
#  MAKE FIGURES
#=================================
res <- process_sim_study(res, scale_CRPS = TRUE)

heatmap_sim_study(filter_sim_study(res, NSR=0),
                  metric = "CRPS_rank",
                  title="CRPS Rank (NSR = 0)",
                  orientation = "horizontal", color_scale = "mako")
ggsave("figs/khaos/heat0.png", width=5, height=5, units="in", dpi=600)

heatmap_sim_study(filter_sim_study(res, NSR=0.5),
                  metric = "CRPS_rank",
                  title="    (NSR = 0.5)",
                  orientation = "horizontal", color_scale = "mako")
ggsave("figs/khaos/heat5.png", width=5, height=5, units="in", dpi=600)


#=================================
#  MAKE TABLES
#=================================
tab1 <- summarize_sim_study(filter_sim_study(res, NSR=0),
                            summarize=c("time", "CRPS"),
                            group_by = c("id"))

tab2 <- summarize_sim_study(filter_sim_study(res, NSR=0.5),
                            summarize=c("time", "CRPS"),
                            group_by = c("id"))


#=================================
#   MAKE SUPPLEMENTAL FIGURES
#=================================
filter_sim_study(res, NSR=0) %>%
  paretoplot_sim_study(metric=c("CRPS_rank", "time"),
                       show_legend=FALSE,
                       title="CRPS Rank vs Time (NSR=0)")
ggsave("figs/khaos/SM_pareto0.png", width=5, height=5, units="in", dpi=600)


filter_sim_study(res, NSR=0.5) %>%
  paretoplot_sim_study(metric=c("CRPS_rank", "time"),
                       show_legend=FALSE,
                       title="(NSR=0.5)")
ggsave("figs/khaos/SM_pareto5.png", width=5, height=5, units="in", dpi=600)

heatmap_sim_study(filter_sim_study(res, NSR=0),
                  metric = "RMSE_rank",
                  title="RMSE Rank (NSR = 0)",
                  orientation = "horizontal", color_scale = "mako")
ggsave("figs/khaos/SM_heat0_rmse.png", width=5, height=5, units="in", dpi=600)

heatmap_sim_study(filter_sim_study(res, NSR=0.5),
                  metric = "RMSE_rank",
                  title="    (NSR = 0.5)",
                  orientation = "horizontal", color_scale = "mako")
ggsave("figs/khaos/SM_heat5_rmse.png", width=5, height=5, units="in", dpi=600)


f_vec <- c("banana", "rabbits", "ishigami", "pollutant_uni", "friedman20")
for(f_curr in f_vec){
  filter_sim_study(res, id=f_curr, NSR=0) %>%
    boxplots_sim_study(y_scale_fun = log10,
                       title=paste0(f_curr, " (NSR=0)"), sort_by="alphabetical")
  ggsave(paste0("figs/khaos/SM_box0_", f_curr, ".png"),
         width=5, height=3, units="in", dpi=600)
  filter_sim_study(res, id=f_curr, NSR=0.5) %>%
    boxplots_sim_study(y_scale_fun = log10,
                       title=paste0(f_curr, " (NSR=0.5)"), sort_by="alphabetical")
  ggsave(paste0("figs/khaos/SM_box5_", f_curr, ".png"),
         width=5, height=3, units="in", dpi=600)
}






