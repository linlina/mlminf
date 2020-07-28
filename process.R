library(xtable)

# Make Figures 1 and 2 using below function.
p_list <- c(300, 600)
q_list <- c(1,2) 
s_list <- c(5, 10)
beta_min_list <- c(0.5, 1)
model_list <- 1:2 

test_figures <- function(model, filename_suffix, pdf_name, group_test = FALSE) {
  tab <- c()
  pow_track <- c()
  type1_track <- c()
  for (p_index in 1:length(p_list)){
    for (q_index in 1:length(q_list)){
      temp <- c()
      for (beta_min_index in 1:length(beta_min_list)){
        for (s_index in 1:length(s_list)){
          p <- p_list[p_index]
          q <- q_list[q_index]
          beta_min <- beta_min_list[beta_min_index]
          s <- s_list[s_index]
          
          filename <- paste(paste("Output/Ridge_lmm", "p", p, "q", q, "model", model, "s", s, "beta.min", beta.min, filename_suffix, sep = "_"), ".RData", sep = "")
          load(filename)
          
          if (!group_test) {
            pow <- mean(colMeans(test.tracker)[1:s])
            type1 <- mean(colMeans(test.tracker)[-c(1:s)])
          } else {
            pow <- colMeans(test.tracker)[1]
            type1 <- colMeans(test.tracker)[2]
          }
          temp <- cbind(temp, paste("(", round(type1, 3), ",", " ", round(pow, 3), ")", sep = ""))
          pow_track <- c(pow_track, pow)
          type1_track <- c(type1_track, type1)
        }
      }
      tab <- rbind(tab, temp)
    }
  }
  
  pdf(file = pdf_name, height = 5, width = 5)
  plot(type1_track, pow_track, xlim = c(1e-8, 1e-1), ylim = c(0, 1), ylab = "Avg. power", log = "x", xlab = "Avg. type I error", bg = "#FF023440", col = "#FF0234", pch = 21, cex = 1)
  abline(v = 0.05, lty = 2, col = "darkgrey")
  dev.off()
}


# Make Figure 3.
pdf(file = "conf_int_coverage.pdf", height = 5, width = 6)
plot(tracker, ylim = c(0.78, 1), cex = 1, bg = "#FF023440", col = "#FF0234", pch = 21, type = "n", ylab = "Confidence interval coverage", xlab = "Index")
abline(0.95, 0, col = "#5b5b5b")
abline(0.95, 0, col = "#5b5b5b50", lwd = 5)
points(lmm_tracker, cex = 1, bg = "#FF023440", col = "#FF0234", pch = 21)
points(ridge_tracker, cex = 1, bg = "#1a3b9130", col = "#1a3b9150", pch = 21)
points(lasso_tracker, cex = 1, bg = "#e8b10030", col = "#e8b10050", pch = 21)
dev.off()
