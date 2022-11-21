

mean_norm <- function(x, na.rm = TRUE) {
  return((x)/(mean(x)))
}

mean_center <- function(x, na.rm = TRUE) {
  return((x)-(mean(x)))
}


# fair_zscore <- function(x,vars,group_var){
#    
#   x <- x %>% 
#     mutate(across(vars,asinh))
#   
#   mean_sd <- x %>%
#     group_by_(group_var) %>% 
#     summarise(across(vars,sd)) %>% ungroup() %>% 
#     summarise(across(vars,mean,.names = '{col}_mean_sd')) 
#     
#   
#   x %>% 
#     bind_cols(mean_sd) %>% 
#     group_by_(group_var) %>% 
#     mutate(across(.cols = vars, .fns = ~ (.x - mean(.x))/mean_sd(.x,)))
# 
#   
# }

calc_mean_ir <- function(data, iridium_channels, q){
  
  
  perc <- data %>% 
    select(iridium_channels) %>% 
    mutate_all(asinh) %>% 
    pivot_longer(names_to = "iridium_channels", values_to = "values", everything()) %>% 
    group_by(iridium_channels) %>% 
    summarise(percentile = quantile(values, probs = q)) %>% 
    ungroup() %>% 
    mutate(sf = percentile[iridium_channels == iridium_channels[1]]/percentile) %>% #first is the reference
    select(iridium_channels, sf)
  
  
  mean_ir <- data %>% 
    select(id,iridium_channels) %>% 
    mutate_at(vars(iridium_channels), asinh) %>%
    pivot_longer(names_to = "iridium_channels", values_to = "value", -id) %>% 
    left_join(perc, by = "iridium_channels") %>% 
    mutate(value = value*sf) %>% 
    group_by(id) %>% 
    summarize(mean_ir = mean(value)) %>%  
    ungroup() %>% 
    mutate(mean_ir = sinh(mean_ir))
  
  return(mean_ir)
}

calc_mean_used_bc <- function(data, bc_channels, bc_col, n_used_bc, q){
  
  
  
perc <- data %>% 
    select(bc_channels) %>% 
    mutate_all(asinh) %>% 
    pivot_longer(names_to = "bc_channels", values_to = "values", everything()) %>% 
    group_by(bc_channels) %>% 
    summarise(percentile = quantile(values, probs = q)) %>% 
    ungroup() %>% 
    mutate(sf = percentile[bc_channels == bc_channels[1]]/percentile) %>% #first is the reference
    select(bc_channels, sf)
    
  
usedBC <- data %>% 
    select(bc_col, bc_channels) %>%
    mutate_at(vars(bc_channels), asinh) %>%
    pivot_longer(names_to  = "bc_channels", values_to = "value", -c(bc_col)) %>% 
    left_join(perc, by = "bc_channels") %>% # percentile normalization
    mutate(value = value*sf) %>% 
    select(-sf) %>% 
    group_by_(bc_col,"bc_channels") %>%
    summarise_at(vars(value), mean) %>% 
    ungroup() %>% 
    group_by_(bc_col) %>% 
    mutate(ranks = rank(-value)) %>% 
    mutate(use = ifelse(ranks %in% c(1:n_used_bc),1,0)) %>% 
    ungroup() %>% 
    select(bc_col,bc_channels,use) 
  
  mean_used_bc <- data %>% 
    select(id,bc_col, bc_channels) %>% 
    mutate_at(vars(bc_channels), asinh) %>% 
    pivot_longer(names_to = "bc_channels", values_to  ="value", -c(bc_col,id)) %>% 
    left_join(usedBC, by = c(bc_col,"bc_channels")) %>% 
    left_join(perc, by = "bc_channels") %>% # percentile normalization
    mutate(value = value*sf) %>% 
    select(-sf) %>% 
    filter(use == 1) %>% # we only care about used barcodes
    group_by(id) %>% 
    summarise(mean_used_bc = mean(value)) %>% 
    ungroup() %>% 
    mutate(mean_used_bc = sinh(mean_used_bc))
  
  
  return(mean_used_bc)
  
}
