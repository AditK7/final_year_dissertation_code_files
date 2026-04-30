# ============================================
# Binding Affinity Distribution - Both Screens
# Overlaid histogram comparing beta-lactamase vs PBP3
# Updated with two-tier threshold scheme
# ============================================
library(ggplot2)

# ---- Load data ----
beta <- read.csv("C:\\Users\\Lenovo\\OneDrive\\Documents\\Final year project\\htvs_ranked_results.csv")
pbp3 <- read.csv("C:\\Users\\Lenovo\\OneDrive\\Documents\\Final year project\\htvs_ranked_results_5df9.csv")

beta$target <- "S70G Beta-Lactamase (3Q07)"
pbp3$target <- "PBP3 (5DF9)"
colnames(beta)[colnames(beta) == "vina_score_kcal_mol"] <- "score"
colnames(pbp3)[colnames(pbp3) == "vina_score_kcal_mol"] <- "score"

combined <- rbind(
  beta[, c("score", "target")],
  pbp3[, c("score", "target")]
)

# ---- Reference scores and half-strength thresholds ----
pip_beta <- -9.617
pip_pbp3 <- -10.76
RT_LN2 <- 0.41
half_beta <- pip_beta + RT_LN2   # -9.207
half_pbp3 <- pip_pbp3 + RT_LN2   # -10.350

# ---- Plot ----
p <- ggplot(combined, aes(x = score, fill = target)) +
  geom_histogram(binwidth = 0.5, alpha = 1, position = "stack", colour = "white", linewidth = 0.2) +
  
  # Reference (strong-hit) thresholds — dashed, thicker
  geom_vline(xintercept = pip_beta, linetype = "solid", colour = "#e63946", linewidth = 0.7) +
  geom_vline(xintercept = pip_pbp3, linetype = "solid", colour = "#2a9d8f", linewidth = 0.7) +
  
  # Half-strength (reasonable-hit) thresholds — dotted, thinner
  geom_vline(xintercept = half_beta, linetype = "dashed", colour = "#e63946", linewidth = 0.6) +
  geom_vline(xintercept = half_pbp3, linetype = "dashed", colour = "#2a9d8f", linewidth = 0.6) +
  
  # 3Q07 reference label (top row)
  annotate("text", x = pip_beta - 0.55, y = Inf, vjust = 2, hjust = 1,
           label = paste0("3Q07 redock\n(", pip_beta, ")"),
           size = 2.6, colour = "#e63946", fontface = "italic") +
  
  # 3Q07 half-strength label (second row down)
  annotate("text", x = half_beta + 0.6, y = Inf, vjust = 4, hjust = 0,
           label = paste0("3Q07 ½-strength\n(", half_beta, ")"),
           size = 2.6, colour = "#e63946", fontface = "italic") +
  
  # 5DF9 reference label (third row down)
  annotate("text", x = pip_pbp3 - 0.55, y = Inf, vjust = 6, hjust = 1,
           label = paste0("5DF9 redock\n(", pip_pbp3, ")"),
           size = 2.6, colour = "#2a9d8f", fontface = "italic") +
  
  # 5DF9 half-strength label (fourth row down)
  annotate("text", x = half_pbp3 + 0.6, y = Inf, vjust = 9, hjust = 0,
           label = paste0("5DF9 ½-strength\n(", half_pbp3, ")"),
           size = 2.6, colour = "#2a9d8f", fontface = "italic") +
  
  scale_x_reverse() +
  scale_fill_manual(values = c(
    "S70G Beta-Lactamase (3Q07)" = "#e63946",
    "PBP3 (5DF9)" = "#2a9d8f"
  )) +
  labs(
    title = "Affinity Distribution: Beta-Lactamase vs PBP3 Screens",
    subtitle = paste0("Beta-lactamase: n=", nrow(beta), " | PBP3: n=", nrow(pbp3),
                      " | Solid = redocking reference, Dashed = half-strength threshold"),
    x = expression(paste("Vina ", Delta, "G (kcal/mol)")),
    y = "Number of Compounds",
    fill = "Target"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, colour = "grey40"),
    legend.position = "bottom"
  )

p

ggsave("C:\\Users\\Lenovo\\OneDrive\\Documents\\Final year project\\affinity_distribution.png", p, width = 10, height = 6, dpi = 300)
cat("Plot saved to: affinity_distribution.png\n")

# ---- Summary stats ----
cat("\n--- Beta-Lactamase ---\n")
cat(sprintf("  n = %d\n", nrow(beta)))
cat(sprintf("  Mean: %.3f kcal/mol\n", mean(beta$score)))
cat(sprintf("  Median: %.3f kcal/mol\n", median(beta$score)))
cat(sprintf("  SD: %.3f\n", sd(beta$score)))

cat("\n--- PBP3 ---\n")
cat(sprintf("  n = %d\n", nrow(pbp3)))
cat(sprintf("  Mean: %.3f kcal/mol\n", mean(pbp3$score)))
cat(sprintf("  Median: %.3f kcal/mol\n", median(pbp3$score)))
cat(sprintf("  SD: %.3f\n", sd(pbp3$score)))