library(fastDummies)

# RUC function

regress_surrogates <- function(data, markers, surrogates, trans_surr, interaction) {
  # line column should be called "line"
  
  # data := data frame with cells in rows, features in columns. CyTOF measurements must be in linear scale.
  # The data set must contain information from ONE CyTOF run.
  
  # markers := vector with marker names (same names must be in the data frame "data") for which we want to apply the regression (y).
  
  # surrogates := vector with surrogates names (same names must be in the data frame "data") that serve as a proxy of the cell size or unwanted covariance (x).
  
  # trans_surr := "per_line" if you want to regress out unwanted covariance within but not between lines -> more conservative
  #    "across_line" if you want to regress-out unwanted covariance within and between the lines -> bold
  
  
  #slope = "per_line", or "across_lines"
  
  #positivecoef <- as.logical(positivecoef)
  interaction <- as.logical(interaction)
  
  if (trans_surr == "per_line") {
    
    data_reg <- data %>%
      mutate_at(vars(markers,surrogates),asinh) %>%
      group_by(line) %>%
      mutate_at(vars(surrogates),mean_center) 
    
  } else if (trans_surr == "across_line") {
    
    data_reg <- data %>%
      mutate_at(vars(markers,surrogates),asinh) %>%
      mutate_at(vars(surrogates),mean_center) %>%
      ungroup()
    
  } else {
    print("Please specify argument 'trans_surr'")
  }
  
  
  # add dummy variables for the lines
  data_reg <- data_reg %>% dummy_cols(select_columns = "line", remove_first_dummy = TRUE) 
  dummy_line_var <- data_reg %>% select(contains("line_")) %>% colnames()
  n_dummy_var <- length(dummy_line_var)
  dummy_values <- data_reg %>% group_by(line) %>% summarise_at(vars(dummy_line_var), max)
  
  
  #model function
  
  if(interaction == TRUE){
    
    intercept_dummy <- paste(dummy_line_var, collapse = " + ")
    
    slope_dummy <- paste(paste0(expand.grid(surrogates,dummy_line_var)[,1],":",expand.grid(surrogates,dummy_line_var)[,2]), collapse =  " + ")
    
    modelfunction <- paste("y",
                           paste0(intercept_dummy, " + ", paste(surrogates, collapse = " + "), " + ", 
                                  slope_dummy), sep = " ~ ")
    
    #lowerbounds <- c(rep(-Inf, length(dummy_line_var)+1), #intercepts
    # rep(0,length(surrogates) + dim(expand.grid(surrogates,dummy_line_var))[1]))
    
    n_coeff <- 1 + length(dummy_line_var) + length(surrogates)+ dim(expand.grid(surrogates,dummy_line_var))[1]
    
    model_coefficients <- matrix(NA,nrow = length(markers),
                                 ncol = n_coeff)
    
    
  }else{
    
    modelfunction <- paste("y",paste(c(surrogates,colnames(data_reg)[str_detect(colnames(data_reg),"line_")]), collapse = " + "),sep = " ~ ")
    
    #lowerbounds <- c(-Inf,rep(0,length(surrogates)), rep(-Inf, n_dummy_var))
    
    n_coeff <- 1 + length(surrogates) + length(dummy_line_var) 
    
    model_coefficients <- matrix(NA,nrow = length(markers),
                                 ncol = n_coeff)
    
  }
  
  # if(positivecoef == FALSE){
  # lowerbounds[lowerbounds==0] <- -Inf
  # }
  
  
  # model coef
  
  
  model_residuals <- matrix(NA, nrow = dim(data_reg)[1], ncol = 1+length(markers))
  adjr2 <- matrix(NA, nrow = length(markers), ncol = 1)
  
  
  
  for (i in 1:length(markers)){
    
    print(paste0("Fitting ", markers[i], " for interaction = ", as.character(interaction), " & trans_surr = ", trans_surr))
    
    #model <- colf::colf_nls(as.formula(gsub("y",markers[i],modelfunction)), 
    # data = data_reg, 
    # lower = lowerbounds, start = rep(0, length(lowerbounds)))
    
    model <- lm(formula = as.formula(gsub("y",markers[i],modelfunction)),
                data = data_reg)
    
    model_coefficients[i,] <- t(as.matrix(coef(model)))
    
    model_residuals[,i+1] <- as.matrix(resid(model))
    
    SSres <- sum(residuals(model)^2)
    SStot <- sum((data_reg[,markers[i]] - mean(as.matrix(data_reg[,markers[i]])))^2)
    model_rsquared <- 1-(SSres/SStot)
    n <- dim(data_reg)[1]
    p <- length(coef(model))-1
    adjr2[i,] <- 1-(1-model_rsquared)*(n-1)/(n-p-1)
    
  }
  
  adjr2 <- as.data.frame(adjr2)
  colnames(adjr2) <- "adj_r_squared"
  adjr2 %>% mutate(marker = markers) -> adjr2
  
  model_residuals[,1] <- data_reg$id
  model_coefficients <- as.data.frame(model_coefficients)
  colnames(model_coefficients) <- colnames(t(as.matrix(coef(model))))
  
  model_coefficients %>% mutate(marker = markers) -> model_coefficients
  
  model_residuals <- as.data.frame(model_residuals)
  colnames(model_residuals) <- c("id",as.vector(markers))
  
  #### standardized regression coefficients
  
  if (interaction == TRUE){
    
    all_coeffs <- model_coefficients %>%
      pivot_longer(names_to = "coef_key", values_to = "value", -marker) %>%
      mutate(coef = str_extract_all(string = coef_key,pattern = paste(surrogates,collapse = "|"), simplify = TRUE)[,1]) %>%
      mutate(coef = ifelse(coef %in% surrogates, coef, "intercept")) %>%
      mutate(coef_type = ifelse(coef %in% surrogates, "surrogate", "intercept")) %>%
      mutate(line = str_extract(string = coef_key, pattern = paste(str_remove(dummy_line_var,"line_"), collapse = "|"))) %>%
      mutate(line = ifelse(is.na(line),as.character(dummy_values$line[1]),line)) %>% select(-coef_key)
    
    intercept_coef <- all_coeffs %>% 
      filter(coef_type == "intercept") %>% 
      select(-coef, -coef_type) %>% 
      group_by(marker) %>% 
      mutate(added_value = ifelse(line ==as.character(dummy_values$line[1]), 0, value[line == as.character(dummy_values$line[1])])) %>% 
      ungroup() %>%  mutate(eff_value = value + added_value)
    
    slope_coef <-  all_coeffs %>% 
      filter(coef_type == "surrogate") %>% 
      select(-coef_type) %>% group_by(marker,coef) %>% 
      mutate(added_value = ifelse(line == as.character(dummy_values$line[1]), 0, value[line == as.character(dummy_values$line[1])])) %>% 
      ungroup() %>%  mutate(eff_value = value + added_value)
    
    sd_markers <- data_reg %>% group_by(line) %>% 
      summarise_at(vars(markers,surrogates),sd) %>% ungroup() %>% pivot_longer(names_to = "marker",values_to = "std",-line)
    
    slope_coef <- sd_markers %>% filter(marker %in% markers) %>% mutate(std_y = std) %>% select(-std) %>% left_join(
      sd_markers %>% filter(marker %in% surrogates) %>% mutate(std_x = std, coef = marker) %>% select(-std,-marker), by = "line") %>% left_join(slope_coef, by = c("line", "marker", "coef")) %>% mutate(stand_value = eff_value*std_x/std_y)
    
    
  }else{
    
    all_coeffs <- model_coefficients %>%
      pivot_longer(names_to = "coef_key", values_to = "value", -marker) %>%
      mutate(coef = str_extract_all(string = coef_key,pattern = paste(surrogates,collapse = "|"), simplify = TRUE)[,1]) %>%
      mutate(coef = ifelse(coef %in% surrogates, coef, "intercept")) %>%
      mutate(coef_type = ifelse(coef %in% surrogates, "surrogate", "intercept")) %>%
      mutate(line = str_extract(string = coef_key, pattern = paste(str_remove(dummy_line_var,"line_"), collapse = "|"))) %>%
      mutate(line = ifelse(coef_type == "intercept" & is.na(line),
                           as.character(dummy_values$line[1]),line)) %>% select(-coef_key)
    
    intercept_coef <- all_coeffs %>% 
      filter(coef_type == "intercept") %>% 
      select(-coef, -coef_type) %>% 
      group_by(marker) %>% 
      mutate(added_value = ifelse(line == as.character(dummy_values$line[1]), 0, value[line == as.character(dummy_values$line[1])])) %>% 
      ungroup() %>%  
      mutate(eff_value = value + added_value)
    
    slope_coef <-  all_coeffs %>% filter(coef_type == "surrogate") %>% select(-coef_type)  %>%
      mutate(eff_value = value) %>% select(marker, value, coef, eff_value)
    
    sd_markers <- data_reg %>% summarise_at(vars(markers,surrogates),sd) %>% pivot_longer(names_to = "marker",values_to = "std", everything())
    
    slope_coef <- slope_coef %>% 
      left_join(sd_markers %>% 
                  filter(marker %in% markers) %>% 
                  mutate(std_y = std) %>% 
                  select(-std), 
                by = "marker") %>%
      left_join(sd_markers %>% 
                  filter(marker %in% surrogates) %>% 
                  mutate(std_x = std) %>% 
                  select(-std) %>% 
                  set_colnames(c("coef","std_x")), 
                by = "coef") %>% 
      mutate(stand_value = eff_value*std_x/std_y)
    
  }
  
  
  ######## New values
  new_values <- model_residuals
  
  for (m in (as.vector(markers))){
    new_values[m] <- as.numeric(as.vector(unlist(model_residuals[m]))) + #residuals
      model_coefficients[model_coefficients$marker == m,1]   # intercept for all
    param_line_offset <- model_coefficients %>% select(dummy_line_var) %>% colnames()
    
    for (i in param_line_offset){
      new_values[m] <- new_values[m] + 
        model_coefficients[model_coefficients$marker == m,i]*data_reg[,i]
    }
  }
  
  data_reg <- data_reg %>% select(-markers) %>% left_join(x = ., y = new_values, by = "id") %>% mutate_at(vars(markers),sinh)
  
  
  data_reg <- data_reg %>% select(-surrogates) %>% left_join(data %>% select(id,surrogates)) 
  
  out_reg <- list(data_reg,markers, surrogates, trans_surr, interaction,modelfunction,model_coefficients,slope_coef, intercept_coef, model_residuals,adjr2)
  names(out_reg) <- c("data_reg","markers", "surrogates", "trans_surr","interaction", "modelfunction","model_coefficients","slope_coef", "intercept_coef",  "model_residuals","adjr2")
  
  return(out_reg)
}
