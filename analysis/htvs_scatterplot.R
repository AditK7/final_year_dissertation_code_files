# ============================================
# Dual-Target Scatter Plot v4
# ALL compounds, two-tier threshold scheme
# Size encodes hit category for visual hierarchy
# ============================================
library(ggplot2)
library(ggrepel)

# ---- Load data ----
beta <- read.csv("C:\\Users\\Lenovo\\OneDrive\\Documents\\Final year project\\htvs_ranked_results.csv")
pbp3 <- read.csv("C:\\Users\\Lenovo\\OneDrive\\Documents\\Final year project\\htvs_ranked_results_5df9.csv")

colnames(beta)[colnames(beta) == "vina_score_kcal_mol"] <- "beta_score"
colnames(pbp3)[colnames(pbp3) == "vina_score_kcal_mol"] <- "pbp3_score"

# Merge — automatically keeps only compounds with scores for both targets
plot_data <- merge(
  beta[, c("ligand_number", "beta_score", "chembl_id", "molecular_weight")],
  pbp3[, c("ligand_number", "pbp3_score")],
  by = "ligand_number"
)
cat(sprintf("Plotting %d compounds with scores for both targets\n", nrow(plot_data)))

# ---- Reference scores and half-strength thresholds ----
pip_beta <- -9.617
pip_pbp3 <- -10.76
RT_LN2 <- 0.41
half_beta <- pip_beta + RT_LN2
half_pbp3 <- pip_pbp3 + RT_LN2

# ---- Classify compounds (two-tier scheme) ----
plot_data$category <- with(plot_data, ifelse(
  beta_score < pip_beta & pbp3_score < pip_pbp3, "Strong dual hit",
  ifelse(
    beta_score < half_beta & pbp3_score < half_pbp3, "Reasonable dual hit",
    ifelse(beta_score < pip_beta, "Beta-lactamase hit only",
           ifelse(pbp3_score < pip_pbp3, "PBP3 hit only",
                  "Below both thresholds"))
  )
))

plot_data$category <- factor(
  plot_data$category,
  levels = c("Strong dual hit", "Reasonable dual hit",
             "Beta-lactamase hit only", "PBP3 hit only",
             "Below both thresholds")
)

# ---- Size grouping (three tiers) ----
plot_data$size_group <- with(plot_data, ifelse(
  category %in% c("Strong dual hit", "Reasonable dual hit"), "dual",
  ifelse(category %in% c("Beta-lactamase hit only", "PBP3 hit only"), "single",
         "below")
))
plot_data$size_group <- factor(plot_data$size_group,
                               levels = c("dual", "single", "below"))

# ---- Label all four dual hits ----
label_these <- c("CHEMBL8432", "CHEMBL10890", "CHEMBL268785", "CHEMBL11008")
plot_data$label <- ifelse(plot_data$chembl_id %in% label_these,
                          as.character(plot_data$chembl_id), "")

# ---- Plot ----
# Order so background points are drawn first, hits drawn on top
plot_data <- plot_data[order(plot_data$size_group, decreasing = TRUE), ]

p <- ggplot(plot_data, aes(x = beta_score, y = pbp3_score,
                           colour = category, size = size_group, alpha = size_group)) +
  geom_point() +
  geom_text_repel(
    aes(label = label),
    size = 2.8,
    max.overlaps = 30,
    segment.size = 0.3,
    show.legend = FALSE,
    box.padding = 0.4,
    point.padding = 0.3,
    inherit.aes = TRUE
  ) +
  
  geom_vline(xintercept = pip_beta, linetype = "solid", colour = "red", alpha = 0.5) +
  geom_hline(yintercept = pip_pbp3, linetype = "solid", colour = "red", alpha = 0.5) +
  geom_vline(xintercept = half_beta, linetype = "dashed", colour = "red", alpha = 0.4) +
  geom_hline(yintercept = half_pbp3, linetype = "dashed", colour = "red", alpha = 0.4) +
  
  scale_x_reverse() +
  scale_y_reverse() +
  scale_colour_manual(values = c(
    "Strong dual hit" = "#e63946",
    "Reasonable dual hit" = "#f4a261",
    "Beta-lactamase hit only" = "#457b9d",
    "PBP3 hit only" = "#2a9d8f",
    "Below both thresholds" = "#8d99ae"
  )) +
  scale_size_manual(values = c("dual" = 4, "single" = 2.5, "below" = 1),
                    guide = "none") +
  scale_alpha_manual(values = c("dual" = 1, "single" = 0.8, "below" = 0.4),
                     guide = "none") +
  labs(
    title = "Dual-Target Virtual Screening: S70G Beta-Lactamase vs PBP3",
    subtitle = paste0(nrow(plot_data), " compounds docked against both targets. ",
                      "Solid lines = redocking reference; dashed lines = half-strength threshold."),
    x = expression(paste("S70G Beta-Lactamase Binding ", Delta, "G (kcal/mol)")),
    y = expression(paste("PBP3 Binding ", Delta, "G (kcal/mol)")),
    colour = "Classification"
  ) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, colour = "grey40"),
    legend.position = "bottom"
  )

p
ggsave("C:\\Users\\Lenovo\\OneDrive\\Documents\\Final year project\\dual_target_scatter.png",
       p, width = 10, height = 8, dpi = 300)