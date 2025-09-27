#===========================================================
#        REAL DATA APPLICATIONS (USING DUQLING)
#===========================================================
data2 <- duqling::get_emulation_data("pbx9501_gold")
XX <- apply(data2$X, 2, BASS:::scale_range)
yy <- (data2$y - mean(data2$y))/sd(data2$y)

fit_pbx <- fit2(XX, yy)
sob2 <- sobol_khaos(fit_pbx)
TT <- cbind(sob2$T, sob2$leftover)
colnames(TT) <- c(c("r0", "r1", "r2", "a", "b", "w"), "epsilon")

original_order <- colnames(TT)
TT_long <- as.data.frame(TT) %>%
  mutate(draw = 1:n()) %>%
  pivot_longer(-draw, names_to = "variable", values_to = "sobol") %>%
  mutate(variable = factor(variable, levels = original_order))

pretty_labels <- c(
  r0 = "italic(r)[0]",
  r1 = "italic(r)[1]",
  r2 = "italic(r)[2]",
  a = "italic(a)",
  b = "italic(b)",
  w = "italic(w)",
  epsilon = "epsilon"
)

ggplot(TT_long, aes(x = variable, y = sobol)) +
  geom_violin(fill = "#3b9f9f", color = "#0a0c32", scale = "width", width = 0.8) +
  scale_x_discrete(labels = function(x) parse(text = pretty_labels[x])) +
  labs(x = NULL, y = "Total Sobol Index") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, family = "serif", size = 8),
    axis.title.y = element_text(family = "serif", size = 14),
    panel.grid.major.x = element_blank()
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  coord_cartesian(ylim = c(0, 1)) +
  ggtitle("Cylinder Experiments")
ggsave("figs/khaos/sobol_cylex.png", width=5, height=5, units="in", dpi=600)


### WINE DATA
data <- read.csv("data/winequality-white.csv", sep=";")
XX <- apply(data[,-12], 2, BASS:::scale_range)
yy <- data[,12] - 2
fit_wine <- ordinal_khaos(XX, yy)
sob2 <- sobol_khaos(fit_wine)
TT <- cbind(sob2$T, sob2$leftover)
colnames(TT) <- c(colnames(XX), "epsilon")

original_order <- colnames(TT)
TT_long <- as.data.frame(TT) %>%
  mutate(draw = 1:n()) %>%
  pivot_longer(-draw, names_to = "variable", values_to = "sobol") %>%
  mutate(variable = factor(variable, levels = original_order))

pretty_labels <- c(
  fixed.acidity         = "italic(F)[acid]",
  volatile.acidity      = "italic(V)[acid]",
  citric.acid           = "italic(C)[acid]",
  residual.sugar        = "italic(R)[sugar]",
  chlorides             = "italic(Cl)",
  free.sulfur.dioxide   = "italic(S)[free]",
  total.sulfur.dioxide  = "italic(S)[total]",
  density               = "italic(rho)",
  pH                    = "italic(pH)",
  sulphates             = "italic(SO)[4]^'-'",
  alcohol               = "italic(A)",
  epsilon               = "epsilon"
)


ggplot(TT_long, aes(x = variable, y = sobol)) +
  geom_violin(fill = "#3b9f9f", color = "#0a0c32", scale = "width", width = 0.8) +
  scale_x_discrete(labels = function(x) parse(text = pretty_labels[x])) +
  labs(x = NULL, y = "Total Sobol Index") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, family = "serif", size = 8),
    axis.title.y = element_text(family = "serif", size = 14),
    panel.grid.major.x = element_blank()
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  coord_cartesian(ylim = c(0, 1)) +
  ggtitle("Wine Quality (Ordinal)")
ggsave("figs/khaos/sobol_wine.png", width=5, height=5, units="in", dpi=600)
