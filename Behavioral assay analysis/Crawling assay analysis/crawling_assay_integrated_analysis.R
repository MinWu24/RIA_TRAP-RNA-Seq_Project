df = read.csv("Combined.csv", header=TRUE)
str(df)

# Core data analysis
library(dplyr)
library(tidyr)     

# 1. Calculate response_score (weighted Z-sums)
df$response_score <- (scale(df$max_turn) * 0.10) +
  (scale(df$turn_95) * 0.1) +
  (scale(df$sharp_turn_proportion) * 0.1) +
  (scale(df$max_turn_accel) * 0.10) +
  (scale(df$tortuosity) * 0.1) +
  (scale(df$traveling_distance) * 0.5)

# 2. Ensure response_score is positive (shift if needed)
df$response_score <- df$response_score - min(df$response_score, na.rm=TRUE) + 0.1

# 3. Calculate UALS = Response
#df$UALS <- df$response_score 
df$UALS <- df$response_score * (1 - df$navigation_index)

# 4. Scale to 0-1 for visualization
df$UALS_scaled <- (df$UALS - min(df$UALS, na.rm=TRUE)) / 
  (max(df$UALS, na.rm=TRUE) - min(df$UALS, na.rm=TRUE))

# ============================================
# STEP 2: CALCULATE NAIVE MEANS (METHOD 1 - SAFE)
# ============================================

# Calculate mean UALS for naive animals within each strain
naive_means <- df %>%
  filter(group == "naive") %>%
  group_by(strain) %>%
  summarise(
    naive_mean_UALS = mean(UALS_scaled, na.rm = TRUE),
    naive_sd_UALS = sd(UALS_scaled, na.rm = TRUE),
    naive_n = n(),
    .groups = 'drop'
  )

# View the naive means
cat("\n=== Naive Means by Strain ===\n")
print(naive_means)

# ============================================
# STEP 3: JOIN NAIVE MEANS TO MAIN DATAFRAME
# ============================================

# Method A: Using left_join
df <- df %>%
  left_join(naive_means, by = "strain")

# Check that the join worked
cat("\n=== Column names after join ===\n")
print(names(df))

# Check if naive_mean_UALS exists
if("naive_mean_UALS" %in% names(df)) {
  cat("\n✓ naive_mean_UALS exists in the dataframe!\n")
} else {
  cat("\n✗ naive_mean_UALS NOT found. Trying alternative method...\n")
  
  # Method B: Manual merge (if left_join fails)
  df <- merge(df, naive_means, by = "strain", all.x = TRUE)
}

# ============================================
# STEP 4: CALCULATE LEARNING INDEX
# ============================================

# Now calculate Learning Index
df <- df %>%
  mutate(
    learning_index = UALS_scaled - naive_mean_UALS
  )

# Verify the calculation
cat("\n=== Learning Index Summary ===\n")
print(summary(df$learning_index))

# ============================================
# STEP 5: VERIFY NAIVE ANIMALS HAVE LI ~ 0
# ============================================

naive_check <- df %>%
  filter(group == "naive") %>%
  group_by(strain) %>%
  summarise(
    mean_LI = mean(learning_index, na.rm = TRUE),
    sd_LI = sd(learning_index, na.rm = TRUE),
    n = n(),
    .groups = 'drop'
  )

cat("\n=== Learning Index for Naive Animals (should be ~0) ===\n")
print(naive_check)
    
