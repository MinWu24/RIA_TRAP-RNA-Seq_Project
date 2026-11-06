library(signal)

# ============================================================
# 1. GAP DETECTION FUNCTIONS
# ============================================================

detect_gaps_aggressive <- function(x, y, max_gap = 200) {
  n <- length(x)
  
  if (n < 2) {
    return(list(gap_indices = integer(0), has_gaps = FALSE, distances = numeric(0)))
  }
  
  distances <- numeric(n-1)
  for (i in 1:(n-1)) {
    distances[i] <- sqrt((x[i+1] - x[i])^2 + (y[i+1] - y[i])^2)
  }
  
  absolute_gaps <- which(distances > max_gap)
  
  median_dist <- median(distances, na.rm = TRUE)
  adaptive_threshold <- median_dist * 10 + 10
  adaptive_gaps <- which(distances > adaptive_threshold)
  
  far_points <- c()
  for (i in 2:(n-1)) {
    dist_to_prev <- sqrt((x[i] - x[i-1])^2 + (y[i] - y[i-1])^2)
    dist_to_next <- sqrt((x[i+1] - x[i])^2 + (y[i+1] - y[i])^2)
    if (dist_to_prev > max_gap && dist_to_next > max_gap) {
      far_points <- c(far_points, i)
    }
  }
  
  all_gaps <- unique(sort(c(absolute_gaps, adaptive_gaps)))
  
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

# ============================================================
# 2. MAIN PROCESSING FUNCTION
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
    
    # Clean data
    colnames(worm) <- c("Frame", "Time", "x", "y")
    worm_clean <- worm[!is.na(worm$x) & !is.na(worm$y), ]
    
    if (nrow(worm_clean) < nrow(worm)) {
      warning(paste("Removed", nrow(worm) - nrow(worm_clean), "rows with NA values"))
    }
    
    worm_data <- data.frame(x = worm_clean$x, y = worm_clean$y)
    
    # Get center from end file
    if (ncol(end_data) >= 4) {
      colnames(end_data) <- c("Frame", "Time", "x", "y")
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
    # PLOT
    # ============================================================
    
    n_points <- nrow(worm_data)
    colors <- rainbow(n_points, start = 0, end = 0.75)
    
    title_text <- paste0(basename(position_file), 
                         "\nNI = ", round(navigation_index, 3))
    
    all_x <- c(worm_data$x, center$x)
    all_y <- c(worm_data$y, center$y)
    
    x_range <- range(all_x, na.rm = TRUE)
    y_range <- range(all_y, na.rm = TRUE)
    x_pad <- diff(x_range) * 0.1
    y_pad <- diff(y_range) * 0.1
    
    plot(worm_data$x, worm_data$y, pch = 19, cex = 0.5, 
         main = title_text, 
         xlab = "X Position", ylab = "Y Position",
         xlim = c(x_range[1] - x_pad, x_range[2] + x_pad),
         ylim = c(y_range[1] - y_pad, y_range[2] + y_pad),
         col = colors,
         type = "p")
    
    # Draw path segments
    for (i in 1:(n_points-1)) {
      segments(worm_data$x[i], worm_data$y[i], 
               worm_data$x[i+1], worm_data$y[i+1],
               col = colors[i], lwd = 0.5)
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
    
    # Legend
    legend("topright", 
           legend = c("Start", "Path", "End", "Destination"),
           col = c(colors[1], "gray", if(n_points>1) colors[n_points] else colors[1], "black"),
           pch = c(18, 19, 15, 17),
           pt.cex = c(1, 0.5, 1, 1),
           bty = "n")
    
    # ============================================================
    # RETURN RESULTS
    # ============================================================
    
    return(data.frame(
      position_file = basename(position_file),
      end_file = basename(end_file),
      date = date,
      group = group,
      strain = strain,
      navigation_index = navigation_index,
      traveling_distance = traveling_distance,
      n_points = nrow(worm_data),
      n_points_original = nrow(worm),
      n_na_removed = nrow(worm) - nrow(worm_data),
      n_points_with_center = n,
      valid_ni_points = length(ni),
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
# 3. RUN THE ANALYSIS
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
  
  if ("n_na_removed" %in% colnames(final_results)) {
    total_na_removed <- sum(final_results$n_na_removed, na.rm = TRUE)
    cat("Total NA rows removed across all files:", total_na_removed, "\n")
  }
  
  write.csv(final_results, "worm_analysis_fixed.csv", row.names = FALSE)
  cat("Results saved to: worm_analysis_fixed.csv\n")
  
  cat("\nSummary of results:\n")
  print(summary(final_results))
}
