load("data/results_paper.Rda")
duq <- process_sim_study(results, scale_CRPS = TRUE)

fname <- "ignition"
res_sub <- filter_sim_study(duq,
                            n_train=1000, NSR=0, id=fname)
boxplots_sim_study(res_sub,
                   ylim = c(min(res_sub$df$CRPS), 1),
                   y_scale_fun = log10,
                   title=paste(fname, "(n=1000, NSR=0)"))
ggsave(filename="figs/main/boxplots1.eps", height=3, width=5)

fname <- "banana"
res_sub <- filter_sim_study(duq,
                            n_train=1000, NSR=0, id=fname)
boxplots_sim_study(res_sub,
                   ylim = c(min(res_sub$df$CRPS), 1),
                   y_scale_fun = log10,
                   title=paste(fname, "(n=1000, NSR=0)"))
ggsave(filename="figs/main/boxplots2.eps", height=3, width=5)

fname <- "cube5"
res_sub <- filter_sim_study(duq,
                            n_train=1000, NSR=0.1, id=fname)
boxplots_sim_study(res_sub,
                   ylim = c(min(res_sub$df$CRPS), 1),
                   y_scale_fun = log10,
                   title=paste(fname, "(n=1000, NSR=0.1)"))
ggsave(filename="figs/main/boxplots3.eps", height=3, width=5)

fname <- "multivalley"
res_sub <- filter_sim_study(duq,
                            n_train=1000, NSR=0.1, id=fname)
boxplots_sim_study(res_sub,
                   ylim = c(min(res_sub$df$CRPS), 1),
                   y_scale_fun = log10,
                   title=paste(fname, "(n=1000, NSR=0.1)"))
ggsave(filename="figs/main/boxplots4.eps", height=3, width=5)


load("data/results_paper_data.Rda")
duq_data <- process_sim_study(results_all)
fname <- "Z_machine_max_vel1"
res_sub <- filter_sim_study(duq_data, id=fname)
boxplots_sim_study(res_sub,
                   ylim = c(min(res_sub$df$CRPS), 1),
                   y_scale_fun = log10,
                   title=paste(fname, ""))
ggsave(filename="figs/main/boxplots5.eps", height=3, width=5)


fname <- "nuclear_data"
res_sub <- filter_sim_study(duq_data, id=fname)
boxplots_sim_study(res_sub,
                   ylim = c(min(res_sub$df$CRPS), 1),
                   y_scale_fun = log10,
                   title=paste(fname, ""))
ggsave(filename="figs/main/boxplots6.eps", height=3, width=5)

## BPPRS vs Svecgp
sub_fnames <- c("rabbits", "park4", "friedman", "steel_column", "ebola", "borehole", "crater",
               "dms_harmonic", "dms_radial", "dms_additive", "dms_complicated", "dms_simple",
               "piston", "ishigami", "Gfunction", "borehole")
duq_sub <- filter_sim_study(duq, n_train=1000, method=c("bppr", "svecgp"), id=sub_fnames)

#Order funcs
lvls <- duq_sub$df %>%
  group_by(id) %>%
  summarise(mean_crps = mean(CRPS, na.rm = TRUE))
duq_sub$df$id <- factor(duq_sub$df$id, levels=lvls$id[order(lvls$mean_crps)])
duq_sub0 <- filter_sim_study(duq_sub, NSR=0)
duq_sub1 <- filter_sim_study(duq_sub, NSR=0.1)

# Example color palette for two methods
fill_colors <- c("bppr" = "#A6CEE3", "svecgp" = "#F8766D")
line_colors <- c("bppr" = "#1F78B4", "svecgp" = "#B22222")  # darker versions

ggplot(duq_sub0$df, aes(x = id, y = CRPS, fill = method, color = method)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_y_log10() +
  scale_fill_manual(values = fill_colors) +
  scale_color_manual(values = line_colors) +
  labs(x = "Test Function", y = "CRPS", fill = "Method", color = "Method") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  ggtitle("n=1000, NSR=0")
ggsave(filename="figs/main/boxplots_comp0.eps", height=3, width=8)

ggplot(duq_sub1$df, aes(x = id, y = CRPS, fill = method, color = method)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_y_log10() +
  scale_fill_manual(values = fill_colors) +
  scale_color_manual(values = line_colors) +
  labs(x = "Test Function", y = "CRPS", fill = "Method", color = "Method") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  ggtitle("n=1000, NSR=0.1")
ggsave(filename="figs/main/boxplots_comp1.eps", height=3, width=8)




