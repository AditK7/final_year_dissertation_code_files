# ============================================
# Full Heavy-Atom RMSD Calculation
# 5DF9 Crystal vs Docked Ligand (59J)
# Matching by atom name
# ============================================

# Crystal structure coordinates (38 heavy atoms)
crystal <- data.frame(
  name = c("C","O","C","N","O","S","C","N","O","C","N","O","C","N","O","C","N","O","C","O","C","O","C","O","C","C","C","C","C","C","C","C","C","C","C","C","C","C"),
  x = c(-42.525,-39.757,-43.713,-42.625,-43.653,-43.364,-43.943,-43.016,-44.694,-43.857,-43.161,-41.0,-42.852,-42.990,-41.0,-42.937,-44.246,-41.144,-43.507,-45.114,-42.324,-46.163,-43.716,-43.974,-41.389,-41.823,-42.359,-42.144,-42.856,-42.322,-40.719,-44.399,-45.002,-44.831,-41.650,-39.994,-40.439,-45.705),
  y = c(-8.769,-12.934,-12.054,-12.202,-8.458,-13.215,-10.626,-10.110,-9.117,-9.700,-7.715,-12.636,-14.667,-5.661,-14.807,-14.714,-3.410,-6.743,-9.329,-6.423,-13.547,-4.663,-15.593,-10.152,-13.672,-9.862,-6.713,-4.658,-3.284,-10.294,-10.449,-5.597,-4.547,-2.335,-11.290,-11.529,-11.932,-1.448),
  z = c(-39.336,-42.293,-35.717,-34.714,-35.190,-37.015,-36.181,-37.189,-38.466,-35.012,-40.159,-32.390,-34.775,-41.628,-32.591,-36.241,-42.912,-40.714,-38.324,-41.226,-34.162,-42.974,-33.983,-33.855,-32.981,-40.153,-40.822,-42.300,-42.438,-41.512,-39.612,-41.776,-42.584,-43.730,-42.201,-40.364,-41.606,-42.870)
)

# Docked ligand coordinates (38 heavy atoms)
docked <- data.frame(
  name = c("N","N","N","N","N","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","O","O","O","O","O","O","O","O","O","S"),
  x = c(-42.927,-43.162,-43.037,-44.256,-42.435,-43.467,-42.534,-42.412,-42.184,-44.448,-42.837,-45.045,-44.835,-44.781,-41.978,-42.592,-40.897,-42.049,-40.312,-40.863,-43.814,-43.570,-42.115,-42.837,-42.679,-43.509,-41.124,-43.677,-44.656,-41.246,-45.169,-46.231,-40.311,-40.779,-40.643,-43.533,-43.693,-43.292),
  y = c(-10.550,-8.564,-6.271,-3.931,-12.503,-10.175,-9.670,-7.372,-5.290,-6.138,-3.882,-5.041,-2.809,-1.544,-10.797,-11.152,-11.487,-12.183,-12.604,-12.938,-11.003,-12.404,-13.817,-15.096,-14.966,-15.837,-13.882,-10.013,-10.259,-7.279,-6.945,-5.093,-13.976,-14.994,-12.818,-8.776,-10.406,-13.641),
  z = c(-37.238,-40.389,-41.453,-42.608,-34.676,-38.545,-39.630,-40.710,-42.148,-41.517,-42.212,-42.264,-43.364,-42.536,-40.510,-41.843,-40.050,-42.592,-40.867,-42.086,-36.166,-35.631,-34.062,-36.037,-34.582,-33.697,-32.923,-35.055,-38.773,-40.362,-40.947,-42.588,-42.834,-32.438,-32.462,-35.312,-33.871,-36.875)
)

cat("Crystal atoms:", nrow(crystal), "\n")
cat("Docked atoms: ", nrow(docked), "\n\n")

# Simple greedy nearest-neighbour matching within each element type
match_by_element <- function(crystal_df, docked_df) {
  pairs <- data.frame(
    crystal_idx = integer(),
    docked_idx = integer(),
    element = character(),
    dist = numeric(),
    stringsAsFactors = FALSE
  )
  
  elements <- unique(crystal_df$name)
  
  for (elem in elements) {
    c_idx <- which(crystal_df$name == elem)
    d_idx <- which(docked_df$name == elem)
    
    cat(sprintf("Element %s: %d crystal, %d docked\n", elem, length(c_idx), length(d_idx)))
    
    if (length(c_idx) != length(d_idx)) {
      cat(sprintf("  WARNING: count mismatch for %s!\n", elem))
      # Use the smaller count
      n <- min(length(c_idx), length(d_idx))
    } else {
      n <- length(c_idx)
    }
    
    # For small groups, try all permutations if <= 8
    # Otherwise use greedy nearest neighbour
    if (n <= 8) {
      # Try all permutations to find optimal matching
      best_total <- Inf
      best_perm <- NULL
      
      if (n == 1) {
        perms <- matrix(1, nrow = 1, ncol = 1)
      } else {
        # Generate permutations
        perm_list <- combinat::permn(n)
        perms <- do.call(rbind, perm_list)
      }
      
      for (p in 1:nrow(perms)) {
        perm <- perms[p, ]
        total <- 0
        for (k in 1:n) {
          ci <- c_idx[k]
          di <- d_idx[perm[k]]
          total <- total + (crystal_df$x[ci] - docked_df$x[di])^2 +
            (crystal_df$y[ci] - docked_df$y[di])^2 +
            (crystal_df$z[ci] - docked_df$z[di])^2
        }
        if (total < best_total) {
          best_total <- total
          best_perm <- perm
        }
      }
      
      for (k in 1:n) {
        ci <- c_idx[k]
        di <- d_idx[best_perm[k]]
        d <- sqrt((crystal_df$x[ci] - docked_df$x[di])^2 +
                    (crystal_df$y[ci] - docked_df$y[di])^2 +
                    (crystal_df$z[ci] - docked_df$z[di])^2)
        pairs <- rbind(pairs, data.frame(
          crystal_idx = ci, docked_idx = di,
          element = elem, dist = d, stringsAsFactors = FALSE
        ))
      }
    } else {
      # Greedy nearest neighbour for larger groups (e.g., 20 carbons)
      available <- d_idx
      for (k in 1:length(c_idx)) {
        ci <- c_idx[k]
        best_d <- Inf
        best_di <- NA
        for (di in available) {
          d <- sqrt((crystal_df$x[ci] - docked_df$x[di])^2 +
                      (crystal_df$y[ci] - docked_df$y[di])^2 +
                      (crystal_df$z[ci] - docked_df$z[di])^2)
          if (d < best_d) {
            best_d <- d
            best_di <- di
          }
        }
        pairs <- rbind(pairs, data.frame(
          crystal_idx = ci, docked_idx = best_di,
          element = elem, dist = best_d, stringsAsFactors = FALSE
        ))
        available <- available[available != best_di]
      }
    }
  }
  
  return(pairs)
}

# Check if combinat is available, install if not
if (!requireNamespace("combinat", quietly = TRUE)) {
  cat("Installing 'combinat' package...\n")
  install.packages("combinat", repos = "https://cloud.r-project.org")
}
library(combinat)

cat("\nMatching atoms by element type...\n\n")
pairs <- match_by_element(crystal, docked)

# Calculate RMSD
sq_dists <- pairs$dist^2
rmsd <- sqrt(mean(sq_dists))

# Output
cat("\n")
cat("=" , rep("=", 55), "\n", sep = "")
cat("ATOM PAIR DISTANCES\n")
cat("=" , rep("=", 55), "\n", sep = "")
cat(sprintf("%-4s  %-12s  %-12s  %s\n", "Elem", "Crystal idx", "Docked idx", "Distance (A)"))
cat(rep("-", 50), "\n", sep = "")
  
for (i in 1:nrow(pairs)) {
  cat(sprintf("%-4s  %-12d  %-12d  %.4f\n",
              pairs$element[i], pairs$crystal_idx[i],
              pairs$docked_idx[i], pairs$dist[i]))
}

cat("\n")
cat("=" , rep("=", 55), "\n", sep = "")
cat(sprintf("FULL HEAVY-ATOM RMSD: %.4f Angstroms\n", rmsd))
cat(sprintf("Number of atom pairs: %d\n", nrow(pairs)))
cat(sprintf("Mean per-atom distance: %.4f A\n", mean(pairs$dist)))
cat(sprintf("Max per-atom distance:  %.4f A\n", max(pairs$dist)))
cat(sprintf("Min per-atom distance:  %.4f A\n", min(pairs$dist)))
cat("=" , rep("=", 55), "\n", sep = "")
