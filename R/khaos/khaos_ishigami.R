set.seed(12345)
p <- 4
sob_true <- c(0.56, 0.44, 0.24, 0)
n_vec <- c(50, 100, 200, 300, 400, 500)
#n_vec <- c(50, 200, 500)
reps <- 5
estimates <- array(NA, dim = c(reps, length(n_vec), 1000, p))
for(i in seq_along(n_vec)){
  for(j in 1:reps){
    X <- lhs::randomLHS(n_vec[i], p)
    y <- apply(X, 1, duqling::ishigami)
    fit <- adaptive_khaos(X, y, legacy=FALSE)
    sob <- sobol_khaos(fit, plot_it=FALSE)

    boxplot(sob$T, ylim=c(0, 1), main=paste(n_vec[i], ", j=", j))
    points(1:4, sob_true, pch=16, cex=2, col='orange')

    estimates[j,i,,] <- as.matrix(sob$T)
  }
}

# Make plot
png(filename="figs/khaos/convergence.png", height=5, width=7, units="in", res=300)
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
cols <- RColorBrewer::brewer.pal(5, "Set1")
n_vals <- n_vec        # Sample sizes
true_vals <- sob_true  # True Sobol index values
for (v in 1:p) {
  mat <- apply(estimates[,,,v], c(1, 2), mean)
  avg <- colMeans(mat)
  plot(range(n_vals), c(0, 1), type = "n",
       xlab = "Sample size (n)", ylab = "Total Sobol Index",
       main = paste0("x", v))

  # Add replication lines
  for (r in 1:reps) {
    for(rr in 1:1000){
      lines(n_vals, estimates[r,,rr,v], col=adjustcolor(cols[r], alpha.f=0.05), lwd=1)
    }
    adj<- 0.8
    lines(n_vals, mat[r, ], col = adjustcolor(cols[r], red.f=adj, green.f=adj, blue.f=adj), lwd = 2)
    points(n_vals, mat[r,], col= adjustcolor(cols[r], red.f=adj, green.f=adj, blue.f=adj), pch=16)
  }

  # Add average lines
  lines(n_vals, avg, col = "black", lwd = 2.5, lty=2)
  points(max(n_vals), true_vals[v], pch = 8, col = "orange", cex = 1.5, lwd=2)
  abline(h = true_vals[v], col = "orange", lty = 3, lwd=3)

  # Add legend
  if (v == p) {
    legend("topright",
           legend = c("MCMC draws",
                      "Posterior mean (per replication)",
                      "Overall mean",
                      "True value"),
           col = c(adjustcolor("gray50", alpha.f = 0.3),
                   adjustcolor("gray10"),
                   "black",
                   "orange"),
           lty = c(1, 1, 2, 3),
           lwd = c(1, 2, 1.5, 1.5),
           pch = c(NA, 16, NA, 8),
           pt.cex = c(1, 1.2, 1, 1),
           bty = "n", cex = 0.8)
  }
}
dev.off()

