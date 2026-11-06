
library(signal)

# ============================================================
# 1. GAP DETECTION FUNCTIONS
# ============================================================

detect_gaps_aggressive <- function(x, y, max_gap = 200) {
  n <- length(x)
  
  if (n < 2) {
    return(list(gap_indices = integer(0), has_gaps = FALSE, distances = numeric(0)))
  }
  
  # Calculate distances
  distances <- numeric(n-1)
  for (i in 1:(n-1)) {
    distances[i] <- sqrt((x[i+1] - x[i])^2 + (y[i+1] - y[i])^2)
  }
  
  # Method 1: Absolute large gaps (> max_gap)
  absolute_gaps <- which(distances > max_gap)
  
  # Method 2: Adaptive threshold (10x median)
  median_dist <- median(distances, na.rm = TRUE)
  adaptive_threshold <- median_dist * 10 + 10
  adaptive_gaps <- which(distances > adaptive_threshold)
  
  # Method 3: Look for large jumps in consecutive points
  far_points <- c()
  for (i in 2:(n-1)) {
    dist_to_prev <- sqrt((x[i] - x[i-1])^2 + (y[i] - y[i-1])^2)
    dist_to_next <- sqrt((x[i+1] - x[i])^2 + (y[i+1] - y[i])^2)
    if (dist_to_prev > max_gap && dist_to_next > max_gap) {
      far_points <- c(far_points, i)
    }
  }
  
  # Combine all gap indices
  all_gaps <- unique(sort(c(absolute_gaps, adaptive_gaps)))
  
  # If we have far points, convert to gap indices
  if (length(far_points) > 0) {
    for (fp in far_points) {
      if (fp > 1 && fp <= n) {
        all_gaps <- c(all_gaps, fp - 1)
      }
    }
    all_gaps <- unique(sort(all_gaps))
  }
  
  return(list(
    gap_indices = all_gaps,
    has_gaps = length(all_gaps) > 0,
    distances = distances,
    n_gaps = length(all_gaps)
  ))
}

split_path_at_gaps <- function(x, y, gap_indices) {
  segments <- list()
  start_idx <- 1
  
  if (length(gap_indices) == 0) {
    return(list(list(x = x, y = y)))
  }
  
  for (gap in gap_indices) {
    if (gap > start_idx) {
      segments[[length(segments) + 1]] <- list(
        x = x[start_idx:gap],
        y = y[start_idx:gap]
      )
    }
    start_idx <- gap + 1
  }
  
  # Add last segment
  if (start_idx <= length(x)) {
    segments[[length(segments) + 1]] <- list(
      x = x[start_idx:length(x)],
      y = y[start_idx:length(y)]
    )
  }
  
  return(segments)
}

get_best_segment_aggressive <- function(segments, min_points = 5) {
  if (length(segments) == 0) return(NULL)
  
  # Find all segments with enough points
  valid_segments <- list()
  for (i in 1:length(segments)) {
    if (length(segments[[i]]$x) >= min_points) {
      valid_segments[[length(valid_segments) + 1]] <- segments[[i]]
    }
  }
  
  if (length(valid_segments) == 0) {
    # No valid segments - try to combine all
    combined_x <- c()
    combined_y <- c()
    for (i in 1:length(segments)) {
      combined_x <- c(combined_x, segments[[i]]$x)
      combined_y <- c(combined_y, segments[[i]]$y)
    }
    if (length(combined_x) >= min_points) {
      return(list(x = combined_x, y = combined_y, combined = TRUE))
    }
    return(NULL)
  }
  
  # Return the largest valid segment
  max_len <- 0
  best_idx <- 1
  for (i in 1:length(valid_segments)) {
    if (length(valid_segments[[i]]$x) > max_len) {
      max_len <- length(valid_segments[[i]]$x)
      best_idx <- i
    }
  }
  
  return(valid_segments[[best_idx]])
}

# ============================================================
# 2. EXTREME FALLBACK FUNCTION
# ============================================================

calculate_sharpness_extreme_fallback <- function(x, y, verbose = FALSE) {
  n <- length(x)
  
  if (n < 3) {
    if (verbose) cat("  EXTREME FALLBACK: Too few points\n")
    return(NA)
  }
  
  # Use VERY heavy smoothing
  window_size <- min(51, floor(n/2))
  if (window_size %% 2 == 0) window_size <- window_size + 1
  if (window_size < 5) window_size <- 5
  
  if (verbose) cat("  EXTREME FALLBACK: Using window size", window_size, "\n")
  
  # Simple moving average with padding
  kernel <- rep(1/window_size, window_size)
  
  x_pad <- c(rep(x[1], window_size), x, rep(x[n], window_size))
  y_pad <- c(rep(y[1], window_size), y, rep(y[n], window_size))
  
  x_smooth <- as.numeric(stats::filter(x_pad, kernel, sides = 2))
  y_smooth <- as.numeric(stats::filter(y_pad, kernel, sides = 2))
  
  x_smooth <- x_smooth[!is.na(x_smooth)]
  y_smooth <- y_smooth[!is.na(y_smooth)]
  
  # Trim to original length
  start_idx <- ceiling(window_size/2) + window_size
  end_idx <- start_idx + n - 1
  
  if (end_idx <= length(x_smooth)) {
    x_smooth <- x_smooth[start_idx:end_idx]
    y_smooth <- y_smooth[start_idx:end_idx]
  } else {
    x_smooth <- x_smooth[1:min(n, length(x_smooth))]
    y_smooth <- y_smooth[1:min(n, length(y_smooth))]
  }
  
  if (length(x_smooth) < 3) {
    if (verbose) cat("  EXTREME FALLBACK: Smoothing failed\n")
    return(NA)
  }
  
  # Calculate approximate heading
  n_smooth <- length(x_smooth)
  heading <- numeric(n_smooth)
  
  for (i in 2:(n_smooth-1)) {
    dx <- x_smooth[min(i+1, n_smooth)] - x_smooth[max(i-1, 1)]
    dy <- y_smooth[min(i+1, n_smooth)] - y_smooth[max(i-1, 1)]
    heading[i] <- atan2(dy, dx)
  }
  heading[1] <- heading[2]
  heading[n_smooth] <- heading[n_smooth-1]
  
  if (any(is.na(heading)) || any(is.infinite(heading))) {
    # Very simple fallback: use overall direction
    overall_direction <- atan2(y[n] - y[1], x[n] - x[1])
    heading <- rep(overall_direction, n_smooth)
  }
  
  # Calculate turn angles
  turn_angles <- diff(heading)
  turn_angles <- turn_angles[!is.na(turn_angles)]
  
  if (length(turn_angles) < 2) {
    return(list(
      sharpness_score = 0.5,
      max_turn = 0.5,
      turn_95 = 0.5,
      sharp_turn_proportion = 0.5,
      tortuosity = 2.0,
      max_turn_accel = 0,
      n_sharp_turns = 0,
      smoothed_x = x_smooth,
      smoothed_y = y_smooth,
      heading = heading,
      turn_angles = turn_angles,
      fallback = TRUE,
      extreme_fallback = TRUE
    ))
  }
  
  # Simple metrics
  max_turn <- max(abs(turn_angles), na.rm = TRUE)
  turn_95 <- quantile(abs(turn_angles), 0.95, na.rm = TRUE)
  n_sharp <- sum(abs(turn_angles) > 0.3, na.rm = TRUE)
  sharp_prop <- n_sharp / length(turn_angles)
  
  # Tortuosity
  path_length <- sum(sqrt(diff(x_smooth)^2 + diff(y_smooth)^2))
  straight_dist <- sqrt((x_smooth[n_smooth] - x_smooth[1])^2 + 
                          (y_smooth[n_smooth] - y_smooth[1])^2)
  if (straight_dist < 1) straight_dist <- 1
  tortuosity <- path_length / straight_dist
  
  # Component scores
  max_turn_score <- exp(-max_turn * 1.5)
  turn_95_score <- exp(-turn_95 * 1.5)
  sharp_prop_score <- exp(-sharp_prop * 10)
  tortuosity_score <- exp(-(tortuosity - 1) * 1.5)
  
  final_score <- 0.35 * max_turn_score + 0.25 * turn_95_score + 
    0.20 * sharp_prop_score + 0.20 * tortuosity_score
  
  final_score <- max(0, min(1, final_score))
  
  if (verbose) cat("  EXTREME FALLBACK: Score =", round(final_score, 4), "\n")
  
  return(list(
    sharpness_score = final_score,
    max_turn = max_turn,
    turn_95 = turn_95,
    sharp_turn_proportion = sharp_prop,
    tortuosity = tortuosity,
    max_turn_accel = 0,
    n_sharp_turns = n_sharp,
    smoothed_x = x_smooth,
    smoothed_y = y_smooth,
    heading = heading,
    turn_angles = turn_angles,
    fallback = TRUE,
    extreme_fallback = TRUE
  ))
}

# ============================================================
# 3. MAIN SHARPNESS CALCULATION FUNCTION
# ============================================================

calculate_turn_sharpness_robust <- function(x, y, 
                                            smooth_window = 51,
                                            turn_threshold = 0.3,
                                            verbose = FALSE) {
  
  n <- length(x)
  
  # --- Check 1: Minimum points ---
  if (n < 3) {
    if (verbose) cat("  ERROR: Less than 3 points (n=", n, ")\n")
    return(NA)
  }
  
  # --- Check if all points are identical ---
  if (all(x == x[1]) && all(y == y[1])) {
    if (verbose) cat("  All points are identical\n")
    return(list(
      sharpness_score = 1.0,
      max_turn = 0,
      turn_95 = 0,
      n_sharp_turns = 0,
      sharp_turn_proportion = 0,
      max_turn_accel = 0,
      tortuosity = 1.0,
      smoothed_x = x,
      smoothed_y = y,
      heading = rep(0, n),
      turn_angles = rep(0, n-1),
      has_gaps = FALSE,
      n_gaps = 0,
      fallback = TRUE,
      extreme_fallback = TRUE
    ))
  }
  
  # --- Check if points are collinear (straight line) ---
  if (n > 2) {
    cross_products <- numeric(n-2)
    for (i in 1:(n-2)) {
      v1 <- c(x[i+1] - x[i], y[i+1] - y[i])
      v2 <- c(x[i+2] - x[i+1], y[i+2] - y[i+1])
      cross_products[i] <- v1[1]*v2[2] - v1[2]*v2[1]
    }
    
    if (max(abs(cross_products)) < 0.001) {
      if (verbose) cat("  Points are collinear (straight line)\n")
      return(list(
        sharpness_score = 1.0,
        max_turn = 0,
        turn_95 = 0,
        n_sharp_turns = 0,
        sharp_turn_proportion = 0,
        max_turn_accel = 0,
        tortuosity = 1.0,
        smoothed_x = x,
        smoothed_y = y,
        heading = rep(atan2(y[n] - y[1], x[n] - x[1]), n),
        turn_angles = rep(0, n-1),
        has_gaps = FALSE,
        n_gaps = 0,
        fallback = TRUE,
        extreme_fallback = TRUE
      ))
    }
  }
  
  # --- Detect gaps ---
  gap_info <- detect_gaps_aggressive(x, y, max_gap = 1000)
  
  if (gap_info$has_gaps) {
    if (verbose) cat("  Found", gap_info$n_gaps, "gaps\n")
    
    segments <- split_path_at_gaps(x, y, gap_info$gap_indices)
    if (verbose) cat("  Split into", length(segments), "segments\n")
    
    best_segment <- get_best_segment_aggressive(segments, min_points = 5)
    
    if (!is.null(best_segment) && length(best_segment$x) >= 5) {
      x_clean <- best_segment$x
      y_clean <- best_segment$y
      if (verbose) cat("  Using segment with", length(x_clean), "points\n")
    } else {
      if (verbose) cat("  No good segment, using extreme fallback\n")
      return(calculate_sharpness_extreme_fallback(x, y, verbose = verbose))
    }
  } else {
    x_clean <- x
    y_clean <- y
    if (verbose) cat("  No gaps detected\n")
  }
  
  # --- If still too few points ---
  if (length(x_clean) < 5) {
    if (verbose) cat("  Too few points, using extreme fallback\n")
    return(calculate_sharpness_extreme_fallback(x, y, verbose = verbose))
  }
  
  # --- Smoothing ---
  n_clean <- length(x_clean)
  
  if (smooth_window %% 2 == 0) smooth_window <- smooth_window + 1
  if (smooth_window > n_clean) smooth_window <- n_clean - (n_clean %% 2 == 0)
  if (smooth_window < 5) smooth_window <- min(5, n_clean)
  
  if (verbose) cat("  Smooth window:", smooth_window, "\n")
  
  x_smooth <- NULL
  y_smooth <- NULL
  
  # Method 1: Savitzky-Golay
  tryCatch({
    x_smooth <- sgolayfilt(x_clean, p = 3, n = smooth_window)
    y_smooth <- sgolayfilt(y_clean, p = 3, n = smooth_window)
    if (verbose) cat("  Savitzky-Golay succeeded\n")
  }, error = function(e) {
    if (verbose) cat("  Savitzky-Golay failed\n")
    x_smooth <<- NULL
    y_smooth <<- NULL
  })
  
  # Method 2: Moving average
  if (is.null(x_smooth) || is.null(y_smooth) || 
      length(x_smooth) < 5 || any(is.na(x_smooth))) {
    if (verbose) cat("  Trying moving average\n")
    kernel <- rep(1/smooth_window, smooth_window)
    
    x_pad <- c(rep(x_clean[1], smooth_window), x_clean, rep(x_clean[n_clean], smooth_window))
    y_pad <- c(rep(y_clean[1], smooth_window), y_clean, rep(y_clean[n_clean], smooth_window))
    
    x_temp <- as.numeric(stats::filter(x_pad, kernel, sides = 2))
    y_temp <- as.numeric(stats::filter(y_pad, kernel, sides = 2))
    
    x_temp <- x_temp[!is.na(x_temp)]
    y_temp <- y_temp[!is.na(y_temp)]
    
    start_idx <- ceiling(smooth_window/2) + smooth_window
    end_idx <- start_idx + n_clean - 1
    
    if (end_idx <= length(x_temp)) {
      x_smooth <- x_temp[start_idx:end_idx]
      y_smooth <- y_temp[start_idx:end_idx]
    } else {
      x_smooth <- x_temp[1:min(n_clean, length(x_temp))]
      y_smooth <- y_temp[1:min(n_clean, length(y_temp))]
    }
  }
  
  # --- If smoothing failed ---
  if (is.null(x_smooth) || is.null(y_smooth) || 
      length(x_smooth) < 5 || length(y_smooth) < 5 ||
      any(is.na(x_smooth)) || any(is.na(y_smooth))) {
    if (verbose) cat("  All smoothing failed, using extreme fallback\n")
    return(calculate_sharpness_extreme_fallback(x, y, verbose = verbose))
  }
  
  # --- Calculate metrics ---
  n_smooth <- length(x_smooth)
  
  # Heading
  heading <- numeric(n_smooth)
  for (i in 2:(n_smooth-1)) {
    dx <- x_smooth[min(i+1, n_smooth)] - x_smooth[max(i-1, 1)]
    dy <- y_smooth[min(i+1, n_smooth)] - y_smooth[max(i-1, 1)]
    heading[i] <- atan2(dy, dx)
  }
  heading[1] <- heading[2]
  heading[n_smooth] <- heading[n_smooth-1]
  
  if (any(is.na(heading)) || any(is.infinite(heading))) {
    if (verbose) cat("  Invalid headings, using extreme fallback\n")
    return(calculate_sharpness_extreme_fallback(x, y, verbose = verbose))
  }
  
  # Turn angles
  turn_angles <- diff(heading)
  
  # Unwrap
  turn_angles_unwrapped <- turn_angles
  for (i in 2:length(turn_angles_unwrapped)) {
    if (!is.na(turn_angles_unwrapped[i]) && !is.na(turn_angles_unwrapped[i-1])) {
      diff_angle <- turn_angles_unwrapped[i] - turn_angles_unwrapped[i-1]
      if (diff_angle > pi) {
        turn_angles_unwrapped[i] <- turn_angles_unwrapped[i] - 2*pi
      } else if (diff_angle < -pi) {
        turn_angles_unwrapped[i] <- turn_angles_unwrapped[i] + 2*pi
      }
    }
  }
  
  turn_angles_unwrapped <- turn_angles_unwrapped[!is.na(turn_angles_unwrapped)]
  
  if (length(turn_angles_unwrapped) < 3) {
    if (verbose) cat("  Too few turn angles, using extreme fallback\n")
    return(calculate_sharpness_extreme_fallback(x, y, verbose = verbose))
  }
  
  # Metrics
  max_turn <- max(abs(turn_angles_unwrapped), na.rm = TRUE)
  turn_95 <- quantile(abs(turn_angles_unwrapped), 0.95, na.rm = TRUE)
  n_sharp_turns <- sum(abs(turn_angles_unwrapped) > turn_threshold, na.rm = TRUE)
  sharp_turn_proportion <- n_sharp_turns / length(turn_angles_unwrapped)
  
  turn_acceleration <- diff(turn_angles_unwrapped)
  max_turn_accel <- max(abs(turn_acceleration), na.rm = TRUE)
  
  path_length <- sum(sqrt(diff(x_smooth)^2 + diff(y_smooth)^2))
  straight_distance <- sqrt((x_smooth[n_smooth] - x_smooth[1])^2 + 
                              (y_smooth[n_smooth] - y_smooth[1])^2)
  
  if (straight_distance < 1e-10) {
    straight_distance <- 1e-10
  }
  
  tortuosity <- path_length / straight_distance
  
  # Component scores
  max_turn_score <- exp(-max_turn * 1.5)
  turn_95_score <- exp(-turn_95 * 1.5)
  sharp_prop_score <- exp(-sharp_turn_proportion * 10)
  accel_score <- exp(-max_turn_accel * 0.5)
  tortuosity_score <- exp(-(tortuosity - 1) * 1.5)
  
  final_score <- (
    0.35 * max_turn_score +      
      0.25 * turn_95_score +        
      0.20 * sharp_prop_score +     
      0.10 * accel_score +          
      0.10 * tortuosity_score       
  )
  
  final_score <- max(0, min(1, final_score))
  
  return(list(
    sharpness_score = final_score,
    max_turn = max_turn,
    turn_95 = turn_95,
    n_sharp_turns = n_sharp_turns,
    sharp_turn_proportion = sharp_turn_proportion,
    max_turn_accel = max_turn_accel,
    tortuosity = tortuosity,
    smoothed_x = x_smooth,
    smoothed_y = y_smooth,
    heading = heading,
    turn_angles = turn_angles_unwrapped,
    has_gaps = gap_info$has_gaps,
    n_gaps = gap_info$n_gaps,
    fallback = FALSE,
    extreme_fallback = FALSE
  ))
}

# ============================================================
process_worm_file_fixed <- function(position_file, verbose = FALSE) {
  # Validate input file exists
  if (!file.exists(position_file)) {
    warning(paste("Position file does not exist:", position_file))
    return(NULL)
  }
  
  # Extract base name
  base_name <- sub("_[pP]osition\\.csv$", "", position_file)
  
  # Find end file
  end_file_lower <- paste0(base_name, "_end.csv")
  end_file_upper <- paste0(base_name, "_End.csv")
  
  if (file.exists(end_file_lower)) {
    end_file <- end_file_lower
  } else if (file.exists(end_file_upper)) {
    end_file <- end_file_upper
  } else {
    end_file_caps <- paste0(base_name, "_END.csv")
    if (file.exists(end_file_caps)) {
      end_file <- end_file_caps
    } else {
      warning(paste("No corresponding end file found for:", position_file))
      return(NULL)
    }
  }
  
  tryCatch({
    # Read files
    worm <- read.csv(position_file, sep = ",", skip = 4, header = TRUE, stringsAsFactors = FALSE)
    end_data <- read.csv(end_file, sep = ",", skip = 4, header = TRUE, stringsAsFactors = FALSE)
    
    # Clean data - REMOVE NA VALUES
    colnames(worm) <- c("Frame", "Time", "x", "y")
    
    # Remove rows with NA in x or y
    worm_clean <- worm[!is.na(worm$x) & !is.na(worm$y), ]
    
    # If too many NAs, warn but continue
    if (nrow(worm_clean) < nrow(worm)) {
      warning(paste("Removed", nrow(worm) - nrow(worm_clean), "rows with NA values"))
    }
    
    worm_data <- data.frame(x = worm_clean$x, y = worm_clean$y)
    
    # Get center from end file
    if (ncol(end_data) >= 4) {
      colnames(end_data) <- c("Frame", "Time", "x", "y")
      # Also clean end data
      end_clean <- end_data[!is.na(end_data$x) & !is.na(end_data$y), ]
      if (nrow(end_clean) > 0) {
        center <- data.frame(x = end_clean$x[1], y = end_clean$y[1])
      } else {
        warning(paste("End file has NA values, using last valid worm point"))
        if (nrow(worm_data) > 0) {
          center <- tail(worm_data, n = 1)
        } else {
          return(NULL)
        }
      }
    } else {
      warning(paste("End file format unexpected for:", end_file))
      if (nrow(worm_data) > 0) {
        center <- tail(worm_data, n = 1)
      } else {
        return(NULL)
      }
    }
    
    # Check if we have enough points
    if (nrow(worm_data) < 3) {
      warning(paste("Too few valid points after cleaning:", nrow(worm_data)))
      return(NULL)
    }
    
    # ============================================================
    # CALCULATE TURN SHARPNESS
    # ============================================================
    
    sharpness_results <- list()
    for (window in c(31, 41, 51, 61)) {
      if (verbose) cat("  Trying window:", window, "\n")
      
      result <- tryCatch({
        calc_result <- calculate_turn_sharpness_robust(
          worm_data$x, worm_data$y, 
          smooth_window = window, 
          verbose = verbose
        )
        if (is.list(calc_result) && !is.null(calc_result$sharpness_score)) {
          calc_result
        } else {
          NULL
        }
      }, error = function(e) {
        if (verbose) cat("  ERROR in window", window, ":", e$message, "\n")
        NULL
      })
      
      if (!is.null(result) && is.list(result)) {
        if (!is.na(result$sharpness_score) && is.finite(result$sharpness_score)) {
          sharpness_results[[as.character(window)]] <- result
          if (verbose) cat("  Window", window, "successful. Score:", 
                           round(result$sharpness_score, 4), "\n")
        }
      }
    }
    
    # Use median across successful windows
    if (length(sharpness_results) > 0) {
      scores <- sapply(sharpness_results, function(r) {
        if (is.list(r) && !is.null(r$sharpness_score)) {
          r$sharpness_score
        } else {
          NA
        }
      })
      scores <- scores[!is.na(scores)]
      
      if (length(scores) > 0) {
        turn_sharpness <- median(scores, na.rm = TRUE)
      } else {
        turn_sharpness <- NA
      }
      
      detailed_result <- sharpness_results[[1]]
      if (is.list(detailed_result)) {
        max_turn <- ifelse(!is.null(detailed_result$max_turn), detailed_result$max_turn, NA)
        turn_95 <- ifelse(!is.null(detailed_result$turn_95), detailed_result$turn_95, NA)
        sharp_turn_proportion <- ifelse(!is.null(detailed_result$sharp_turn_proportion), 
                                        detailed_result$sharp_turn_proportion, NA)
        max_turn_accel <- ifelse(!is.null(detailed_result$max_turn_accel),  
                                 detailed_result$max_turn_accel, NA)  
        tortuosity <- ifelse(!is.null(detailed_result$tortuosity), detailed_result$tortuosity, NA)
        has_gaps <- ifelse(!is.null(detailed_result$has_gaps), detailed_result$has_gaps, FALSE)
        n_gaps <- ifelse(!is.null(detailed_result$n_gaps), detailed_result$n_gaps, 0)
        fallback <- ifelse(!is.null(detailed_result$fallback), detailed_result$fallback, FALSE)
        extreme_fallback <- ifelse(!is.null(detailed_result$extreme_fallback), 
                                   detailed_result$extreme_fallback, FALSE)
      } else {
        max_turn <- NA
        turn_95 <- NA
        sharp_turn_proportion <- NA
        max_turn_accel <- NA  
        tortuosity <- NA
        has_gaps <- FALSE
        n_gaps <- 0
        fallback <- FALSE
        extreme_fallback <- FALSE
      }
    } else {
      turn_sharpness <- NA
      max_turn <- NA
      turn_95 <- NA
      sharp_turn_proportion <- NA
      max_turn_accel <- NA 
      tortuosity <- NA
      detailed_result <- NULL
      has_gaps <- FALSE
      n_gaps <- 0
      fallback <- TRUE
      extreme_fallback <- TRUE
    }
    
    # ============================================================
    # NAVIGATION INDEX CALCULATION
    # ============================================================
    
    worm_data_with_center <- worm_data
    worm_data_with_center <- rbind(worm_data_with_center, center)
    
    distance <- numeric()
    n <- nrow(worm_data_with_center)
    for (i in 1:(n-1)) {
      distance[i] <- sqrt((worm_data_with_center[i+1, 1] - worm_data_with_center[i, 1])^2 + 
                            (worm_data_with_center[i+1, 2] - worm_data_with_center[i, 2])^2)
    }
    
    sec <- seq(1, n, by = 14)
    worm_2sec <- worm_data_with_center[sec, ]
    m <- nrow(worm_2sec)
    
    dis <- numeric()
    distance2 <- numeric()
    a <- numeric()
    v_radial <- numeric()
    v_real <- numeric()
    ni <- numeric()
    
    for (j in 1:(m-1)) {
      dis[j] <- sqrt((center[1, 1] - worm_2sec[j, 1])^2 + (center[1, 2] - worm_2sec[j, 2])^2)
      distance2[j] <- sqrt((worm_2sec[j+1, 1] - worm_2sec[j, 1])^2 + (worm_2sec[j+1, 2] - worm_2sec[j, 2])^2)
      dis[j+1] <- sqrt((center[1, 1] - worm_2sec[j+1, 1])^2 + (center[1, 2] - worm_2sec[j+1, 2])^2)
      
      cos_val <- (dis[j]^2 + distance2[j]^2 - dis[j+1]^2) / (2 * dis[j] * distance2[j] + 0.0001)
      cos_val <- pmax(pmin(cos_val, 1), -1)
      
      a[j] <- acos(cos_val)
      v_radial[j] <- distance2[j] * cos(a[j])
      v_real[j] <- sum(distance[(14*j-13):(14*j)], na.rm = TRUE)
      ni[j] <- v_radial[j] / (v_real[j] + 0.0001)
    }
    
    ni <- ni[is.finite(ni)]
    
    if (length(ni) > 0) {
      navigation_index <- mean(ni, na.rm = TRUE)
    } else {
      navigation_index <- NA
    }
    
    traveling_distance <- sum(distance, na.rm = TRUE) / 1000
    
    # ============================================================
    # EXTRACT INFORMATION FROM FILENAME
    # ============================================================
    
    date_match <- regmatches(position_file, regexpr("20\\d{6}", position_file))
    date <- ifelse(length(date_match) > 0, date_match[1], NA)
    
    if (grepl("naive", position_file, ignore.case = TRUE)) {
      group <- "naive"
    } else if (grepl("trained", position_file, ignore.case = TRUE)) {
      group <- "trained"
    } else {
      group <- "unknown"
    }
    
    filename_lower <- tolower(position_file)
    
    if (grepl("n2", filename_lower)) {
      strain <- "N2"
    } else if (grepl("tm2332", filename_lower)) {
      if (grepl("tat-2b", filename_lower) && !grepl("tat-2e", filename_lower)) {
        strain <- "Pglr-3s::tat-2b"
      } else if (grepl("tat-2e", filename_lower)) {
        strain <- "Pglr-3s::tat-2e"
      } else {
        strain <- "tm2332"
      }
    } else {
      strain <- "unknown"
    }
    
    # ============================================================
    # PLOT WITH GAP HANDLING
    # ============================================================
    
    n_points <- nrow(worm_data)
    colors <- rainbow(n_points, start = 0, end = 0.75)
    
    # Detect gaps for plotting
    gap_info <- detect_gaps_aggressive(worm_data$x, worm_data$y, max_gap = 1000)
    
    title_text <- paste0(basename(position_file), 
                         "\nNI = ", round(navigation_index, 3), 
                         ", Sharpness = ", round(turn_sharpness, 3),
                         ifelse(has_gaps, paste0(" [", n_gaps, " gaps]"), ""),
                         ifelse(fallback || extreme_fallback, " [fallback]", ""))
    
    # Plot ranges
    all_x <- c(worm_data$x, center$x)
    all_y <- c(worm_data$y, center$y)
    
    x_range <- range(all_x, na.rm = TRUE)
    y_range <- range(all_y, na.rm = TRUE)
    x_pad <- diff(x_range) * 0.1
    y_pad <- diff(y_range) * 0.1
    
    # Create plot
    plot(worm_data$x, worm_data$y, pch = 19, cex = 0.5, 
         main = title_text, 
         xlab = "X Position", ylab = "Y Position",
         xlim = c(x_range[1] - x_pad, x_range[2] + x_pad),
         ylim = c(y_range[1] - y_pad, y_range[2] + y_pad),
         col = colors,
         type = "p")
    
    # Draw segments, skipping gaps
    for (i in 1:(n_points-1)) {
      is_gap <- FALSE
      if (gap_info$has_gaps) {
        if (i %in% gap_info$gap_indices) {
          is_gap <- TRUE
        }
      }
      
      if (!is_gap) {
        segments(worm_data$x[i], worm_data$y[i], 
                 worm_data$x[i+1], worm_data$y[i+1],
                 col = colors[i], lwd = 0.5)
      }
    }
    
    # Mark gaps
    if (gap_info$has_gaps) {
      for (gap_idx in gap_info$gap_indices) {
        if (gap_idx < n_points) {
          points(worm_data$x[gap_idx], worm_data$y[gap_idx], 
                 col = "red", pch = 4, cex = 1.5)
        }
      }
    }
    
    # Start point
    points(worm_data$x[1], worm_data$y[1], 
           col = colors[1], cex = 1.5, pch = 18, lwd = 1)
    
    # End point
    if (n_points > 1) {
      points(worm_data$x[n_points], worm_data$y[n_points], 
             col = colors[n_points], cex = 1.2, pch = 15, lwd = 1)
    }
    
    # Destination
    points(center$x, center$y, 
           col = "black", cex = 3, pch = 17, lwd = 2)
    
    # Add smoothed path
    if (!is.null(detailed_result) && is.list(detailed_result) && 
        !is.null(detailed_result$smoothed_x) && !is.null(detailed_result$smoothed_y)) {
      lines(detailed_result$smoothed_x, 
            detailed_result$smoothed_y, 
            col = "darkred", lwd = 2, lty = 2)
    }
    
    # Legend
    legend("topright", 
           legend = c("Start", "Path", "End", "Destination", "Gap", "Smoothed"),
           col = c(colors[1], "gray", if(n_points>1) colors[n_points] else colors[1], 
                   "black", "red", "darkred"),
           pch = c(18, 19, 15, 17, 4, NA),
           lty = c(NA, NA, NA, NA, NA, 2),
           lwd = c(NA, NA, NA, NA, NA, 2),
           pt.cex = c(1, 0.5, 1, 1, 1, NA),
           bty = "n")
    
    # ============================================================
    # RETURN RESULTS - NOW INCLUDING max_turn_accel
    # ============================================================
    
    return(data.frame(
      position_file = basename(position_file),
      end_file = basename(end_file),
      date = date,
      group = group,
      strain = strain,
      navigation_index = navigation_index,
      turn_sharpness = turn_sharpness,
      max_turn = max_turn,
      turn_95 = turn_95,
      sharp_turn_proportion = sharp_turn_proportion,
      max_turn_accel = max_turn_accel,  # ← Now defined
      tortuosity = tortuosity,
      traveling_distance = traveling_distance,
      n_points = nrow(worm_data),
      n_points_original = nrow(worm),
      n_na_removed = nrow(worm) - nrow(worm_data),
      n_points_with_center = n,
      valid_ni_points = length(ni),
      has_gaps = has_gaps,
      n_gaps = n_gaps,
      fallback_used = fallback || extreme_fallback,
      center_x = center$x,
      center_y = center$y,
      stringsAsFactors = FALSE
    ))
    
  }, error = function(e) {
    warning(paste("Error processing file:", position_file, "-", e$message))
    return(NULL)
  })
}
# ============================================================
# 5. RUN THE ANALYSIS
# ============================================================

data_directory <- "./positions"

if (!dir.exists(data_directory)) {
  stop(paste("Directory does not exist:", data_directory))
}

cat("Looking for files in directory:", data_directory, "\n")

position_files <- list.files(path = data_directory, 
                             pattern = "*[pP]osition\\.csv$", 
                             full.names = TRUE)
position_files <- unique(position_files)

if (length(position_files) == 0) {
  stop(paste("No position files found in:", data_directory))
}

cat("Found", length(position_files), "position files to process\n")

# Create PDF
n_plots <- length(position_files)
n_pages <- ceiling(n_plots / 10)
cat("Creating PDF with", n_pages, "page(s) for", n_plots, "plots\n")

pdf("worm_tracking_plots_fixed.pdf", 
    width = 15, height = 9, 
    onefile = TRUE)

par(mfrow = c(2, 5))
par(mar = c(3, 3, 2, 1))
par(oma = c(1, 1, 1, 1))

results_list <- list()
error_count <- 0
plot_count <- 0

for (file in position_files) {
  tryCatch({
    result <- process_worm_file_fixed(file, verbose = FALSE)
    if (!is.null(result)) {
      results_list[[file]] <- result
      plot_count <- plot_count + 1
      cat(sprintf("Processed and plotted (%d/%d): %s\n", 
                  plot_count, n_plots, basename(file)))
    }
  }, error = function(e) {
    error_count <- error_count + 1
    warning(paste("Error processing file:", file, "-", e$message))
  })
}

dev.off()

cat("\n========================================\n")
cat("PDF saved to: worm_tracking_plots_fixed.pdf\n")
cat("========================================\n")
cat("Total plots generated:", plot_count, "\n")
cat("Files with errors:", error_count, "\n")

# Combine results
if (length(results_list) > 0) {
  final_results <- do.call(rbind, results_list)
  rownames(final_results) <- NULL
  
  cat("\n========================================\n")
  cat("PROCESSING COMPLETE\n")
  cat("========================================\n")
  cat("Files processed successfully:", length(results_list), "\n")
  
  # Check NA removal statistics
  if ("n_na_removed" %in% colnames(final_results)) {
    total_na_removed <- sum(final_results$n_na_removed, na.rm = TRUE)
    cat("Total NA rows removed across all files:", total_na_removed, "\n")
  }
  
  # Check fallbacks
  if ("fallback_used" %in% colnames(final_results)) {
    n_fallback <- sum(final_results$fallback_used, na.rm = TRUE)
    cat("Files using fallback:", n_fallback, "\n")
  }
  
  write.csv(final_results, "worm_analysis_fixed.csv", row.names = FALSE)
  cat("Results saved to: worm_analysis_fixed.csv\n")
  
  cat("\nSummary of results:\n")
  print(summary(final_results))
}

