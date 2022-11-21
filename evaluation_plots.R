
filter <- dplyr::filter

# compare before/after
#pearson correlation
plot_pearson_corr <- function(data,surrogates,markers,out_reg,dir, clust, width, height){
  
  p <- list()
  
  for (l in unique(data$line)){
    
    
    a <- data %>% 
      filter(line == l) %>% 
      select(surrogates,markers) %>% 
      mutate_all(asinh) %>% 
      cor(method = "pearson") %>% 
      ggcorrplot(type = "upper",outline.col = "white", show.legend = TRUE, title = "Before", hc.order = clust) 
    
    
    b <- out_reg$data_reg %>% 
      filter(line == l) %>% 
      select(surrogates,markers) %>% 
      mutate_all(asinh) %>% 
      cor(method = "pearson") %>% 
      ggcorrplot(type = "upper",outline.col = "white", show.legend = TRUE, title = "After", hc.order = clust) 
    
    
    p[[l]] <- ggarrange(plotlist = list(a,b), ncol = 2)
    
    if(clust == TRUE){
      c = "_clust"
    }else{
      c = ""
    }
    
    ggsave(filename = paste0(dir,"pearsoncorr",c,"/",l,".png"), plot = p[[l]],
           width = width, height = height)
    
    
  }
  
  return(p)
  
}

# scatter plots with regression line

plot_scatter_fit <- function(data,markers,surrogates,out_reg,frac, dir, width, height){
  
  
  if (out_reg$trans_surr == "per_line") {
    
    dt_before <- data %>%
      mutate_at(vars(markers,surrogates),asinh) %>%
      group_by(line) %>%
      mutate_at(vars(surrogates),mean_center) #### improve and use fair_scale!!!!
    
    
  } else if (out_reg$trans_surr == "across_line") {
    
    dt_before <- data %>%
      mutate_at(vars(markers,surrogates),asinh) %>%
      mutate_at(vars(surrogates),mean_center) %>%
      ungroup()
    
  }
  
  dt_after <- out_reg$data_reg %>% 
    mutate_at(vars(markers,surrogates),asinh) %>%
    mutate_at(vars(surrogates),mean_center) %>%
    ungroup()
  
  p <- list()
  ids_plot <- dt_before %>% group_by(line) %>% 
    sample_frac(frac) %>% pull(id)
  
  for (m in markers){
    
    
    if (out_reg$interaction == TRUE){
      line_dt <- data.frame(line = out_reg$slope_coef %>% filter(marker == m) %>% pull(line),
                            surrogates = out_reg$slope_coef %>% filter(marker == m) %>% pull(coef),
                            slope = out_reg$slope_coef %>% filter(marker == m) %>% pull(eff_value)) %>% 
        left_join(out_reg$intercept_coef %>% filter(marker == m) %>% mutate(intercept = eff_value) %>% select(line, intercept), by = "line")
      
    }else{
      
      line_dt <- data.frame(line = out_reg$intercept_coef %>% filter(marker == m) %>% pull(line),
                            intercept = out_reg$intercept_coef %>% filter(marker == m) %>% pull(eff_value)) %>% 
        left_join(
          data.frame(surrogates = rep(out_reg$slope_coef %>% filter(marker == m) %>% pull(coef),length(unique(out_reg$intercept_coef$line))),
                     slope = rep(out_reg$slope_coef %>% filter(marker == m) %>% pull(eff_value),length(unique(out_reg$intercept_coef$line)))) %>% mutate(line = rep(unique(out_reg$intercept_coef$line),times = rep(length(unique(out_reg$slope_coef$coef)),length(unique(out_reg$intercept_coef$line))))), by = "line")
      
      
    }
    
    
    
    a <- dt_before %>% 
      filter(id %in% ids_plot) %>% 
      ungroup() %>% 
      select(m,surrogates, line) %>% 
      pivot_longer(names_to = "surrogates",values_to = "surr_value",surrogates) %>% 
      ggplot(aes_string(x = "surr_value", y = m, fill = "line", color = "line")) + 
      geom_point(alpha = 0.5, size = 0.1)  + 
      facet_wrap(~surrogates, scale = "free_x", ncol = 1) +
      theme_minimal() + 
      geom_vline(xintercept = 0) +
      #geom_hline(yintercept = 0) +
      geom_abline(mapping = aes(intercept = intercept, slope = slope, col = line), 
                  data = line_dt) + ggtitle("Before")
    
    b <- dt_after %>% 
      group_by(line) %>% 
      filter(id %in% ids_plot) %>% 
      select(m,surrogates, line) %>% 
      pivot_longer(names_to = "surrogates",values_to = "surr_value",surrogates) %>% 
      ggplot(aes_string(x = "surr_value", y = m, fill = "line", color = "line")) + 
      geom_point(alpha = 0.5, size = 0.1)  + 
      facet_wrap(~surrogates, scale = "free_x", ncol = 1)+
      theme_minimal() + 
      geom_vline(xintercept = 0) +
      guides(colour = guide_legend(override.aes = list(size=3, alpha = 1))) + 
      #geom_hline(yintercept = 0) +
      ggtitle("After")
    
    p[[m]] <-  ggarrange(plotlist = list(a,b), ncol = 2)
    
    ggsave(filename = paste0(dir,"scatter_fit/",m,".png"), plot =  p[[m]],
           width = width, height = height)
    
    
  }
  
  return(p)
  
}


# scatter plots with regression line

plot_residuals <- function(data,markers,surrogates,out_reg,frac,dir, width, height){
  
  
  dt <- data %>%
    mutate_at(vars(markers),asinh) %>% 
    select(id,markers,line) %>% 
    pivot_longer(names_to = "marker", values_to = "original_values", markers) %>% 
    left_join(out_reg$model_residuals %>% pivot_longer(names_to = "marker", values_to = "residuals", markers), by = c("id","marker"))
  
  
  
  ids_plot <- data %>% group_by(line) %>% 
    sample_frac(frac) %>% pull(id)
  
  p <-  dt %>% filter(id %in% ids_plot ) %>% 
    ggplot(aes(x = original_values, y = residuals, color = line)) + 
    geom_point(alpha = 0.5, size = 0.1)  + 
    facet_wrap(~marker) +
    geom_hline(yintercept = 0)
  
  ggsave(filename = paste0(dir,"original_residuals.png"), plot = p, width = width, height = height)
  
  
  return(p)
  
}


# model coefficients


plot_model_coeffs <- function(out_reg,dir,width,height){
  
  p <- out_reg$intercept_coef %>% 
    ggplot(aes(x = line, y = eff_value, fill = line)) + 
    geom_col() + 
    facet_wrap(~marker, scale = "free") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  
  
  ggsave(filename = paste0(dir,"coeffs_intercepts.png"), plot = p, width = width, height = height)
  
  if(out_reg$interaction==FALSE){
    
    p <- out_reg$slope_coef %>% 
      ggplot(aes(x = coef, y = stand_value)) + 
      geom_col() + 
      facet_wrap(~marker) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    ggsave(filename = paste0(dir,"coeffs_slope.png"), plot = p,  width = width, height = height)
    
  }else{
    
    p <- out_reg$slope_coef %>% 
      ggplot(aes(x = coef, y = stand_value, fill = line)) + 
      geom_col(position = "dodge") + 
      facet_wrap(~marker) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    ggsave(filename = paste0(dir,"coeffs_slope.png"), plot = p,  width = width, height = height)
    
  }
  
  
  return(p)
  
}


# adjusted R2: goodness of fit

plot_adjR2 <- function(list_out_reg,list_names,dir, width, height){
  
  dt <- list_out_reg[[1]]$adjr2 %>% mutate(data = list_names[[1]]) 
  
  if(length(list_out_reg) > 1){
    for (i in 2:length(list_out_reg)){
      dt <- dt %>% rbind(list_out_reg[[i]]$adjr2 %>% mutate(data = list_names[[i]]))
    }
  }

  
  p <- dt %>% 
    ggplot(aes(x = data, y = adj_r_squared, fill = data, label = round(adj_r_squared,2))) + 
    geom_col() + 
    facet_wrap(~marker) + 
    geom_label() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  ggsave(filename = paste0(dir,"adjR2.png"), plot = p, width = width, height = height)
  return(p)
  
}



### Compare different models

#distributions

plot_distributions <- function(data,markers,list_out_reg,list_names,dir, width, height){
  
  dt <- data %>% select(markers,line) %>% mutate(data = "original") 
  
  for (i in 1:length(list_out_reg)){
    dt <- dt %>% rbind( list_out_reg[[i]]$data_reg %>% select(markers,line) %>% mutate(data = list_names[[i]])) 
  }
  
  p <- list() 
  
  for (m in markers){
    
    p[[m]] <- dt %>% 
      mutate_at(vars(m),asinh) %>% 
      mutate(data = factor(data, levels = c("original",list_names))) %>% 
      ggplot(aes_string(x = m, y = "line", color = "data")) + 
      geom_density_ridges(alpha = 0)
    
    
    ggsave(filename = paste0(dir,"distributions/",m,".png"), plot = p[[m]], width = width, height = height)
    
  }
  
  return(p)
  
  
}

# mean values

plot_meanvalues <- function(data,markers,list_out_reg,list_names,dir, width, height){
  
  
  dt <- data %>% select(markers,line) %>% mutate(data = "original") 
  
  for (i in 1:length(list_out_reg)){
    dt <- dt %>% rbind(list_out_reg[[i]]$data_reg %>% select(markers,line) %>% mutate(data = list_names[[i]])) 
  }
  
  
  p <- dt %>% 
    mutate_at(vars(markers),asinh) %>% 
    group_by(line,data) %>% 
    summarize_at(vars(markers),mean) %>% 
    pivot_longer(names_to = "marker", values_to = "mean_value", markers) %>% 
    mutate(data = factor(data, levels = c("original",list_names))) %>% 
    mutate(marker = factor(marker, levels = markers)) %>% 
    ggplot(aes(x = line, y = mean_value, fill  = data)) +
    geom_col(position = "dodge") + 
    facet_wrap(~marker, scale = "free") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  ggsave(paste0(dir,"meanvalues.png"), plot = p, width = width, height = height)
  
  return(p)
  
}


