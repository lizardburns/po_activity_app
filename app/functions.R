######## Phenoloxidase activity app functions =================================

calculate_vmax <- function(well, df){
  out <- list()
  absorbance <- df[[well]]
  times <- df[[2]]
  
  # at what point does the absorbance start increasing consistently?
  # looks for 2 consecutive positive gradients (3 points)
  gradient_check <- tibble(diff = diff(absorbance)) %>% 
    mutate(direction = if_else(diff > 0, 1, 0)) %>% 
    mutate(lag1 = lag(direction),
           direction_diff = if_else(direction == lag1, "same", "switched")) %>% 
    mutate(seq_len = if_else(direction == 1 & lag1 == 1, 2, 0))
  
  # what's the first point we won't throw out?
  # if the line is jumping up and down repeatedly set it to NA in order to throw out this well
  switch_tally <- gradient_check %>% 
    count(direction_diff) %>% 
    filter(!is.na(direction_diff))
  
  switch_count <- switch_tally %>% filter(direction_diff == "switched") %>% pull(n)
  # steady_count <- switch_tally %>% filter(direction_diff == "same") %>% pull(n)
  
  if(length(switch_count) > 0){
    if(switch_count / sum(switch_tally$n, na.rm = TRUE) > .75){
      cut_point <- NA
    } else {cut_point <- which((gradient_check %>% slice(1:ceiling(nrow(.) * 2 / 3)))$seq_len == 2)[1] - 1}
  } else {
    # use slice as don't want to consider late ascents
    cut_point <- which((gradient_check %>% slice(1:ceiling(nrow(.) * 2 / 3)))$seq_len == 2)[1] - 1
  }
  
  # get the cleaned up data
  if(is.na(cut_point)){
    out$vmax <- NA_real_
    out$p <- ggplot() + 
      geom_point(data = df, 
                 aes(x = df[[2]], y = df[[well]]), 
                 alpha = .5) +
      labs(title = paste("Well", well), x = "Time /s", y = "Absorbance") +
      theme_bw()
  } else {
    times <- times[cut_point:length(times)]
    absorbance <- absorbance[cut_point:length(absorbance)]
    
    # fit a smooth line to the remaining data
    span <- .6
    # don't fit line to end 3 time points to prevent vmax corresponding to short sharp climb towards end of run
    predict_times <- data.frame(times = min(times):floor(times[length(times)-3]))
    # fit loess
    fitted_vals <- suppressWarnings(predict(loess(absorbance ~ times, span = span), 
                                            newdata = predict_times))
    slopes <- diff(fitted_vals)#/diff(predict_times$times)
    
    out$vmax <- max(slopes, na.rm = TRUE)
    
    vmax_time <- predict_times$times[which.max(slopes %>% unname())] + .5
    
    out$p <- ggplot() + 
      geom_point(data = df, 
                 aes(x = df[[2]], y = df[[well]]), alpha = .5) +
      geom_line(aes(x = predict_times$times, y = fitted_vals), 
                inherit.aes = FALSE, na.rm = TRUE) +
      geom_vline(aes(xintercept = vmax_time),
                 colour = "#5ab4ac") +
      labs(title = paste("Well", well), x = "Time /s", y = "Absorbance") +
      theme_bw()
  }
  
  return(out)
}

# run the vmax calculation and plotting for all reactions on plate
calculate_dabs_plots <- function(df){
  reactions <- df %>% 
    select_if(str_detect(names(.), "[A-Z]{1}[0-9]+")) %>% 
    names()
  
  reactions %>% 
    map(~calculate_vmax(.x, df))
}

calculate_results_summary <- function(lambda475, lambda600, lambda490, protein_table, blank475, blank490){
  n_samples <- length(lambda475) / 2
  
  results <- tibble(well = rep(paste0(rep(LETTERS[1:8], each = 12), rep(1:12, times = 8)), times = 3),
                    wavelength = rep(c(475, 600, 490), each = 96),
                    vmax = c(lambda475, lambda600, lambda490)) %>% 
    mutate(dabs_min = vmax * 60) %>% 
    pivot_wider(names_from = wavelength, values_from = c(vmax, dabs_min)) %>% 
    mutate(dabs_min_475_600 = dabs_min_475 - dabs_min_600,
           sample_index = rep(1:n_samples, each = 2), 
           replicate_id = rep(c("rep1", "rep2"), times = n_samples)) %>% 
    pivot_wider(names_from = replicate_id, 
                values_from = c(well, vmax_475, vmax_600, vmax_490, dabs_min_475, dabs_min_600, dabs_min_490, dabs_min_475_600)) %>% 
    mutate(replicates_flag_475 = dplyr::if_else(dabs_min_475_rep1 / dabs_min_475_rep2 > 2 | dabs_min_475_rep1 / dabs_min_475_rep2 < .5, "warning", NA_character_),
           replicates_flag_600 = dplyr::if_else(dabs_min_600_rep1 / dabs_min_600_rep2 > 2 | dabs_min_600_rep1 / dabs_min_600_rep2 < .5, "warning", NA_character_),
           replicates_flag_490 = dplyr::if_else(dabs_min_490_rep1 / dabs_min_490_rep2 > 2 | dabs_min_490_rep1 / dabs_min_490_rep2 < .5, "warning", NA_character_)) %>% 
    left_join(protein_table, ., by = c("well" = "well_rep1")) %>% 
    mutate(mean_dabs_min_475_600 = (dabs_min_475_600_rep1 + dabs_min_475_600_rep2) / 2,
           mean_dabs_min_490 = (dabs_min_490_rep1 + dabs_min_490_rep2) / 2)
  
  if(str_detect(str_to_lower(paste(results$sampleID, collapse = " ")), "blank")){
    tmp_blank475 <- filter(results, str_to_lower(sampleID) == "blank") %>% pull(mean_dabs_min_475_600)
    tmp_blank490 <- filter(results, str_to_lower(sampleID) == "blank") %>% pull(mean_dabs_min_490)
    
    if(is.na(tmp_blank475)){
      blank_val_475_600 <- 0
    } else {
      blank_val_475_600 <- tmp_blank475
    }
    
    if(is.na(tmp_blank490)){
      blank_val_490 <- 0
    } else {
      blank_val_490 <- tmp_blank490
    }
    
  } else {
    blank_val_475_600 <- blank475
    blank_val_490 <- blank490
  }
  
  results %>% mutate(`mean_dabs_min_475_600 - blank` = mean_dabs_min_475_600 - blank_val_475_600,
                     `mean_dabs_min_490 - blank` = mean_dabs_min_490 - blank_val_490,
                     mean_dabs_min_475_600_mgprotein = `mean_dabs_min_475_600 - blank` / mg_protein,
                     mean_dabs_min_490_mgprotein = `mean_dabs_min_490 - blank` / mg_protein) %>% 
    select(-sample_index)
}

calculate_490_summary <- function(lambda490, protein_table, blank490){
  n_samples <- length(lambda490) / 2
  
  results <- tibble(well = paste0(rep(LETTERS[1:8], each = 12), rep(1:12, times = 8)),
                    # wavelength = rep(490, 96),
                    vmax = lambda490) %>% 
    mutate(dabs_min = vmax * 60,
           sample_index = rep(1:n_samples, each = 2), 
           replicate_id = rep(c("rep1", "rep2"), times = n_samples)) %>% 
    pivot_wider(names_from = replicate_id, 
                values_from = c(well, vmax, dabs_min)) %>% 
    mutate(replicates_flag_490 = dplyr::if_else(dabs_min_rep1 / dabs_min_rep2 > 2 | dabs_min_rep1 / dabs_min_rep2 < .5, "warning", NA_character_)) %>% 
    left_join(protein_table, ., by = c("well" = "well_rep1")) %>% 
    mutate(mean_dabs_min_490 = (dabs_min_rep1 + dabs_min_rep2) / 2)
  
  if(str_detect(str_to_lower(paste(results$sampleID, collapse = " ")), "blank")){
    blank_val_490 <- filter(results, str_to_lower(sampleID) == "blank") %>% pull(mean_dabs_min_490)
  } else {
    blank_val_490 <- blank490
  }
  
  results %>% mutate(`mean_dabs_min_490 - blank` = mean_dabs_min_490 - blank_val_490,
                     mean_dabs_min_490_mgprotein = `mean_dabs_min_490 - blank` / mg_protein) %>% 
    select(-sample_index)
}

# transpose_tibble <- function(df){
#   df %>% 
#     tibble::rownames_to_column() %>% 
#     tidyr::pivot_longer(-rowname) %>% 
#     tidyr::pivot_wider(names_from = rowname, values_from = value) %>% 
#     dplyr::select(-name)
# }
