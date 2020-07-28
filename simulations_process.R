#processing simulations

load(file = "ridge_henderson_simulation_v2_lm_unknown_fixed_sc.RData")

ridge_henderson_tracker <- tracker

load(file = "ridge_indep_simulation_v2_fixed.RData")

ridge_indep_tracker <- tracker

load(file = "lasso_indep_simulation_v2_fixed.RData")

lasso_indep_tracker <- tracker


pdf(file = "conf_int_coverage.pdf", height = 5, width = 6)
plot(ridge_henderson_tracker, ylim = c(0.78, 1), cex = 1, bg = "#FF023440", col = "#FF0234", pch = 21, type = "n", ylab = "Confidence interval coverage", xlab = "Index")
abline(0.95, 0, col = "#5b5b5b")
abline(0.95, 0, col = "#5b5b5b50", lwd = 5)
points(ridge_henderson_tracker, cex = 1, bg = "#FF023440", col = "#FF0234", pch = 21)

points(ridge_indep_tracker, cex = 1, bg = "#1a3b9130", col = "#1a3b9150", pch = 21)

points(lasso_indep_tracker, cex = 1, bg = "#e8b10030", col = "#e8b10050", pch = 21)
dev.off()

pdf(file = "conf_int_coverage_support.pdf", height = 5, width = 5)
plot(ridge_henderson_tracker[1:5], ylim = c(0.78, 1), cex = 0.8, bg = "#FF023440", col = "#FF0234", pch = 21, type = "n", ylab = "CI coverage", xlab = "Coordinate")
abline(0.95, 0, col = "#5b5b5b")
abline(0.95, 0, col = "#5b5b5b50", lwd = 5)
points(ridge_henderson_tracker[1:5], cex = 0.8, bg = "#FF023440", col = "#FF0234", pch = 21)

points(ridge_indep_tracker[1:5], cex = 0.8, bg = "#1a3b9130", col = "#1a3b91", pch = 21)

points(lasso_indep_tracker[1:5], cex = 0.8, bg = "#e8b10030", col = "#e8b100", pch = 21)
dev.off()

pdf(file = "conf_int_printer_friendly.pdf", height = 5, width = 6)
load("Output/Ridge_lmm_conf_int_fixed_no_delta.RData")
load("Output/Ridge_lmm_conf_int_fixed_ridgedelta.RData")
plot(lmm.tracker, ylim = c(0.78, 1), cex = 1, bg = "#FF023440", col = "#FF0234", pch = 21, type = "n", ylab = "Confidence interval coverage", xlab = "Index")
abline(0.95, 0, col = "#5b5b5b")
abline(0.95, 0, col = "#5b5b5b50", lwd = 5)
points(lmm.tracker, cex = 1, bg = "#FF023450", col = "#FF0234", pch = 16)
points(ridge.tracker, cex = 1, bg = "#1a3b9190", col = "#1a3b9190", pch = 0)
points(lasso.tracker, cex = 1, bg = "#e8b10090", col = "#e8b10090", pch = 4)
legend("bottomright", c("lmm", "ridge", "lasso"), col = c("#FF0234", "#1a3b91", "#e8b100"), pch = c(16, 0, 4))
dev.off()