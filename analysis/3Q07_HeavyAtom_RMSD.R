# ============================================
# Full Heavy-Atom RMSD Calculation
# 3Q07 Crystal vs Docked Piperacillin (WPP)
# Matching by element type + nearest neighbour
# ============================================

library(combinat)

crystal <- data.frame(
  name = c("N","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","O","N","N","C","N","N","O","O","O","O","O","O","S"),
  x = c(-6.262,-5.638,-9.098,-3.461,-4.392,-9.812,-3.799,-4.663,-11.463,-7.572,-15.023,-8.078,-11.898,-1.961,-13.252,-14.756,-12.318,-8.081,-7.122,-7.145,-8.139,-9.126,-4.199,-3.602,-9.408,-11.046,-5.264,-13.432,-3.914,-8.353,-1.435,-1.268,-9.071,-14.173,-11.621,-5.911),
  y = c(-16.369,-16.708,-17.813,-18.450,-19.632,-14.372,-20.948,-19.665,-12.756,-16.460,-11.247,-15.903,-14.410,-18.670,-13.845,-12.063,-11.826,-16.872,-16.766,-17.641,-18.607,-18.689,-16.517,-15.496,-15.322,-13.850,-18.213,-12.648,-17.791,-16.877,-18.459,-18.977,-13.998,-14.477,-15.370,-19.311),
  z = c(12.994,14.275,10.231,12.705,12.339,10.818,12.829,10.840,10.051,12.775,12.484,11.479,11.823,12.552,11.888,11.236,10.893,10.346,9.343,8.259,8.163,9.143,14.055,13.762,11.640,10.934,14.401,11.356,13.938,13.598,11.431,13.551,9.923,12.402,12.524,13.199)
)

docked <- data.frame(
  name = c("N","N","N","N","N","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","O","O","O","O","O","O","O","S"),
  x = c(-6.724,-4.264,-9.749,-11.257,-13.846,-7.997,-6.148,-4.741,-5.565,-3.634,-4.352,-3.607,-4.510,-2.111,-8.389,-8.239,-7.350,-8.913,-7.176,-8.745,-7.885,-10.003,-11.659,-12.093,-13.094,-13.438,-15.095,-14.773,-8.799,-4.283,-1.518,-1.477,-9.097,-11.779,-14.137,-5.958),
  y = c(-16.630,-17.671,-15.390,-13.806,-13.025,-16.291,-16.852,-16.452,-18.282,-18.274,-19.587,-20.792,-19.683,-18.270,-15.919,-17.015,-16.881,-18.218,-17.928,-19.262,-19.115,-14.287,-12.681,-14.334,-12.950,-13.751,-12.285,-10.839,-16.187,-15.358,-18.396,-18.220,-13.728,-15.218,-13.925,-19.483),
  z = c(12.880,14.015,11.247,10.494,10.513,12.686,14.208,14.080,14.398,12.832,12.435,12.997,10.923,12.785,11.288,10.289,9.229,10.477,8.320,9.570,8.488,10.567,11.344,9.571,11.762,9.486,10.469,10.152,13.585,13.801,11.687,13.868,9.971,8.790,8.489,13.184)
)

cat("Crystal atoms:", nrow(crystal), "\n")
cat("Docked atoms: ", nrow(docked), "\n\n")

# Match atoms by element type with optimal/greedy assignment
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
    
    n <- min(length(c_idx), length(d_idx))
    
    if (n <= 8) {
      best_total <- Inf
      best_perm <- NULL
      
      perm_list <- combinat::permn(n)
      
      for (p in seq_along(perm_list)) {
        perm <- perm_list[[p]]
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

cat("\nMatching atoms by element type...\n\n")
pairs <- match_by_element(crystal, docked)

sq_dists <- pairs$dist^2
rmsd <- sqrt(mean(sq_dists))

cat("\n")
cat("=========================================================\n")
cat("ATOM PAIR DISTANCES\n")
cat("=========================================================\n")
cat(sprintf("%-4s  %-12s  %-12s  %s\n", "Elem", "Crystal idx", "Docked idx", "Distance (A)"))
cat(rep("-", 50), "\n", sep = "")

for (i in 1:nrow(pairs)) {
  cat(sprintf("%-4s  %-12d  %-12d  %.4f\n",
              pairs$element[i], pairs$crystal_idx[i],
              pairs$docked_idx[i], pairs$dist[i]))
}

cat("\n")
cat("=========================================================\n")
cat(sprintf("FULL HEAVY-ATOM RMSD: %.4f Angstroms\n", rmsd))
cat(sprintf("Number of atom pairs: %d\n", nrow(pairs)))
cat(sprintf("Mean per-atom distance: %.4f A\n", mean(pairs$dist)))
cat(sprintf("Max per-atom distance:  %.4f A\n", max(pairs$dist)))
cat(sprintf("Min per-atom distance:  %.4f A\n", min(pairs$dist)))
cat("=========================================================\n")