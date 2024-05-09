library(lubridate)

### Organizando o código

# Prepare the data by time and delay
create_data <- function(data, D = 15, return_data_real = F){
  # input: data.frame with 3 columns (onset date, digitation date, disease)
  # output: a data.frame with the structure by time and delay
  
  if(return_data_real){
    # Ver se será últil, e quando será útil
    # (MODIFICAR DEPOIS PARA SER USADO COM OUTROS DADOS)
    data_real <- dados |>
      filter(dt_sin_pri > ymd("2017-01-01")) |>
      select(dt_digita, dt_sin_pri, arbo) |>
      tidyr::drop_na(dt_digita) |>
      mutate(delay = as.numeric(floor(difftime(dt_digita, dt_sin_pri, units = "weeks")))) |>
      group_by(dt_sin_pri, delay, arbo)|>
      summarise(y = n()) |>
      ungroup() |>
      mutate(epiweek = epiweek(dt_sin_pri), year = year(dt_sin_pri)) %>% 
      group_by(year, epiweek, arbo) |>
      summarise(y = sum(y), date = min(dt_sin_pri)) 
    
    return(data_real)
  }
  
  data <- data |> 
    janitor::clean_names() |> 
    tidyr::drop_na(dt_digita) |>
    dplyr::mutate(delay = as.numeric(floor(difftime(dt_digita, dt_sin_pri, units = "weeks")))) |>
    dplyr::group_by(dt_sin_pri, delay, arbo)|>
    dplyr::summarise(y = n()) |>
    dplyr::filter(delay <= D) 
  
  data <- data |>
    dplyr::ungroup() |>
    dplyr::mutate(epiweek = epiweek(dt_sin_pri), year = year(dt_sin_pri)) |>
    dplyr::group_by(epiweek, year, delay, arbo) |>
    dplyr::summarise(y = sum(y)) |> 
    tidyr::pivot_wider(names_from = delay, values_from = y) |>
    dplyr::mutate_if(~is.numeric(.), ~coalesce(.,0)) |>
    tidyr::pivot_longer(cols = !(epiweek:arbo), names_to = "delay", values_to = "y") |>
    dplyr::mutate(delay = as.numeric(delay)) |>
    dplyr::arrange(year, epiweek, delay) |> 
    tidyr::pivot_wider(names_from = arbo, values_from = y)
  
  return(data)
}

# Filter data and prepare it to enter the model
model_data <- function(data, max_year = 2019, max_epiweek = 20){
  # input: output from the create_data function
  # output: data with the required structure for the model
  data %>% 
    dplyr::ungroup() %>% 
    dplyr::filter(year < max_year | epiweek <= max_epiweek) %>% 
    dplyr::mutate(chik_obs = ifelse(epiweek + delay > max_epiweek & year == max_year, NA, Chik),
                  den_obs = ifelse(epiweek + delay > max_epiweek & year == max_year, NA, Den)) -> data_model
  
  data_model %>% 
    dplyr::distinct(epiweek, year) %>% 
    dplyr::mutate(time = 1:n()) %>% 
    dplyr::inner_join(data_model, by = c("epiweek", "year")) -> data_model
  
  y_new <- data_model %>% 
    dplyr::mutate(report_time = time + delay) %>% 
    dplyr::select(time, delay, report_time, y1_obs = chik_obs, 
                  y2_obs = den_obs, y1 = Chik, y2 = Den)
  
  return(y_new)
}



run_inla_model <- function(data, evol_time = "rw2", evol_delay = "rw1",
                           family = "nbinomial", model = NULL, 
                           data_inla = NULL){
  # Future plans: adapt this function to include regressor and higher dimensions
  if(is.null(data_inla)){
    
  
  y_df <- data
  
  n <- nrow(y_df)
  
  y_bvar <- matrix(NA, nrow = 2 * n, ncol = 2)
  y_bvar[1:n, 1] <- y_df$y1_obs
  y_bvar[n + 1:n, 2] <- y_df$y2_obs
  
  I <- matrix(NA, nrow = 2 * n, ncol = 2)
  I[1:n, 1] <- 1
  I[n + 1:n, 2] <- 1
  
  
  fixed_covariate <- data.frame(alphai = as.factor(c(rep(1, n), rep(2, n))),
                                id_time1 = c(rep(1,n), rep(NA, n)),
                                id_delay1 = c(rep(1,n), rep(NA, n)),
                                id_time2 = c(rep(NA, n), rep(1, n)),
                                id_delay2 = c(rep(NA, n), rep(1, n)),
                                idrep = rep(1:2, each = n))
  
  
  random_covariate <- list(time1  = c(y_df$time, rep(NA, n)),
                           delay1 = c(y_df$delay, rep(NA, n)),
                           time2  = c(rep(NA, n), y_df$time),
                           delay2 = c(rep(NA, n), y_df$delay),
                           time = rep(y_df$time, 2),
                           delay = rep(y_df$delay, 2))
  
  
  Y <<- y_bvar
  joint_data <<- cbind(Y, fixed_covariate, random_covariate)}
  
  
  joint_model <- Y ~ -1 + alphai +
    + f(time,  model = evol_time, replicate = idrep,
        hyper = list(prec = list(param = c(0.001, 0.001)))) +
    + f(delay, model = evol_delay, replicate = idrep) +
    + f(time1, model = evol_time, constr = F) +
    + f(time2, copy  = "time1", fixed = TRUE) + 
    + f(delay1, model = evol_delay, constr = F) +
    + f(delay2, copy  = "delay1", fixed = TRUE) 
  
  if(!is.null(model)){joint_model <- model}
  if(!is.null(data_inla)){joint_data <- data_inla}
  
  if(family == "poisson"){
    control_family_list <- list()
  }else if(family == "nbinomial"){
    control_family_list <- list(
      list(hyper = list("theta" = list(prior = "loggamma", param = c(0.001, 0.001)))),
      list(hyper = list("theta" = list(prior = "loggamma", param = c(0.001, 0.0001))))
    )
  }
  
  tryCatch({
    output0 <- INLA::inla(joint_model, family = rep(family, 2),
                          data = joint_data,
                          control.predictor = list(link = 1, compute = T),
                          control.compute = list(config = T, dic=T, waic = T,
                                                 cpo = T),
                          control.family = control_family_list
    )  
  }, error = function(cond){
    message(cond)
    return(NA)
  })
  # output0 <- INLA::inla(joint_model, family = rep(family, 2),
  #                       data = joint_data,
  #                       control.predictor = list(link = 1, compute = T),
  #                       control.compute = list(config = T, dic=T, waic = T,
  #                                              cpo = T),
  #                       control.family = control_family_list
  # )
  
  return(list(inla_output = output0, 
              joint_data = joint_data))
}


gg.age <- function(x, dados.gg, idx, n = 1984){
  idx_1 <- idx[idx <= n]
  idx_2 <- idx[idx > n]
  data.aux <- dados.gg
  Tmin <- min(dados.gg$time[idx])
  data.aux$`1`[idx_1] <- x[,1]
  data.aux$`2`[idx_2] <- x[,2]
  data.aggregated <- data.aux |>
    ## Selecionando apenas os dias faltantes a partir
    ## do domingo da respectiva ultima epiweek
    ## com dados faltantes
    dplyr::filter(time >= Tmin  ) |>
    dplyr::group_by(time) |>
    dplyr::summarise(
      Y1 = sum(`1`, na.rm = T), 
      Y2 = sum(`2`, na.rm = T), .groups = "keep" 
    )
  data.aggregated
}

## AINDA RESTA ADAPTAR PARA POISSON
sample_predictive <- function(inla_output, joint_data, family = "nbinomial", 
                              return_posterior = F){
  
  # Sample from the posterior 
  samples <- INLA::inla.posterior.sample(n = 1000, inla_output)
  
  # Which cells haven't been observed yet
  index_missing <- which(is.na(joint_data$`1`) & is.na(joint_data$`2`))
  
  # Sample from the "missing" Y_{t,d} given the parameters sampled from the
  # posterior
  vector_samples_bvar <- lapply(X = samples,
                                FUN = function(x, idx = index_missing, 
                                               n = nrow(joint_data)/2){
                                  idx_1 <- idx[idx <= n]
                                  idx_2 <- idx[idx > n]
                                  unif.log = 1
                                  
                                  if(family == "nbinomial"){
                                    y1 <- stats::rnbinom(n = idx_1,
                                                         mu = exp(x$latent[idx_1]),
                                                         size = x$hyperpar[1]
                                    ) * unif.log
                                    
                                    y2 <- stats::rnbinom(n = idx_2,
                                                         mu = exp(x$latent[idx_2]),
                                                         size = x$hyperpar[2]
                                    ) * unif.log  
                                  }else if(family == "poisson"){
                                    y1 <- stats::rpois(n = idx_1,
                                                       lambda = exp(x$latent[idx_1])
                                    ) * unif.log
                                    
                                    y2 <- stats::rpois(n = idx_2,
                                                       lambda = exp(x$latent[idx_2])
                                    ) * unif.log
                                  }
                                  
                                  
                                  return(cbind(y1,y2))
                                } )
  
  
  tibble.samples.0 <- lapply(X = vector_samples_bvar,
                             FUN = gg.age,
                             dados.gg = joint_data,
                             idx = index_missing,
                             n = nrow(joint_data)/2)
  
  srag.pred.0 <- dplyr::bind_rows(tibble.samples.0, .id = "sample")
  
  if(return_posterior){return(list(posterior_samples = samples,
                                   predictive = srag.pred.0))}
  return(srag.pred.0)
}


bivar_nowcasting_summary <- function(predictive_samples){
  total.summy <- predictive_samples |>
    dplyr::group_by(time, sample) |>
    dplyr::summarise(Y1 = sum(Y1, na.rm = T),
                     Y2 = sum(Y2, na.rm = T)) |>
    dplyr::group_by(time) |>
    dplyr::summarise(Median_y1 = stats::median(Y1, na.rm = T),
                     LI_y1  = stats::quantile(Y1, probs = 0.025, na.rm = T),
                     LS_y1  = stats::quantile(Y1, probs = 0.975, na.rm = T),
                     LIb_y1 = stats::quantile(Y1, probs = 0.25, na.rm = T),
                     LSb_y1 = stats::quantile(Y1, probs = 0.75, na.rm = T),
                     Median_y2 = stats::median(Y2, na.rm = T),
                     LI_y2  = stats::quantile(Y2, probs = 0.025, na.rm = T),
                     LS_y2  = stats::quantile(Y2, probs = 0.975, na.rm = T),
                     LIb_y2 = stats::quantile(Y2, probs = 0.25, na.rm = T),
                     LSb_y2 = stats::quantile(Y2, probs = 0.75, na.rm = T),
                     .groups = "drop")
  
  summary_bivar <- total.summy
  return(summary_bivar)
}


weighted_interval_score <- function(pred_samples, data_real, alpha = 5:95/100){
  
  pred_time = sort(unique(pred_samples$time))
  y_real <- data_real %>% filter(time %in% pred_time) %>% 
    select(time, y1, y2)
  
  
  
  pred_samples$Y1
  for(t in 1:length(unique(pred_samples$time))){
    y_t <- pred_samples |> filter(time == t) #|> pull(matches("y"))
    quants_t <- sapply(alpha, \(alpha){
      quantile(y_t, c(alpha/2, 1- alpha/2))
      }) |> 
      t()
    wis_t <- interval_score(true_values = y_real,lower = quants_t[,1],
                             upper = quants_t[,2],
                             weigh = T,
                             interval_range = 100*(1-alpha)) |>
      mean() 
      

      
  }
  sapply(alpha, \(x){quantile(pred_samples, prob = c(x, 1-x))})
  IS_95 = interval_score(true_values = obs_chik, 
                         lower = LI_y1, upper = LS_y1, 
                         interval_range = 95)
}


scores_posterior <- function(pred_samples, data_real){
  require(scoringRules)
  require(magrittr)
  
  pred_time = sort(unique(pred_samples$time))
  y_real <- data_real %>% filter(time %in% pred_time)
  
  energy_score <- sapply(1:length(pred_time), 
                             function(i){
                               es_sample(y_real[i,2:3] %>%
                                           as.numeric(), 
                                         pred_samples %>%
                                           filter(time == pred_time[i]) %>% 
                                           .[,c(3,4)] %>%
                                           as.matrix() %>% t())})
  
  
  
  variogram_score <- sapply(1:length(pred_time), 
                            function(i){
                              vs_sample(y_real[i,2:3] %>% as.numeric(),
                                        pred_samples %>%
                                          filter(time == pred_time[i]) %>%
                                          .[,c(3,4)] %>%
                                          as.matrix() %>% t())})
  
  
  
  return(cbind.data.frame(energy_score, variogram_score))
  
}

selection_criteria <- function(summary_pred, data_real, inla_output){
  ## Criteria: waic, dic, length of the last interval, mean interval score, cpo?
  ## mse, mae, mape
  require(scoringutils)
  require(metrica)
  
  obs_chik <- data_real %>% 
    filter(time %in% summary_pred$time) %>% 
    pull(y1)
  
  obs_dengue <- data_real %>% 
    filter(time %in% summary_pred$time) %>% 
    pull(y2)
  
  
  metrics_chik <- metrica::metrics_summary(obs = obs_chik, 
                                           pred = summary_pred$Median_y1,
                                           metrics_list = c("MAPE", "MAE", "MSE"),
                                           type = "regression") |> 
    tidyr::pivot_wider(names_from = Metric, values_from = Score)
  
  metrics_dengue <- metrica::metrics_summary(obs = obs_dengue, 
                                             pred = summary_pred$Median_y2, 
                                             metrics_list = c("MAPE", "MAE", "MSE"),
                                             type = "regression") |> 
    tidyr::pivot_wider(names_from = Metric, values_from = Score)
  
  interval_metrics_chik <- summary_pred |> 
    arrange(time) |> 
    mutate(IS_95 = interval_score(true_values = obs_chik, 
                                  lower = LI_y1, upper = LS_y1, 
                                  interval_range = 95),
           IS_50 = interval_score(true_values = obs_chik, 
                                  lower = LIb_y1, upper = LSb_y1,
                                  interval_range = 50),
           range_90 = LS_y1 - LI_y1,
           range_50 = LSb_y1 - LIb_y1
    ) |> 
    select(IS_95:range_50) |> 
    summarise_all(list(last = ~last(.),
                       mean = ~mean(.))) |> 
    select(!(starts_with("range") & ends_with("mean")))
  
  interval_metrics_dengue <- summary_pred |> 
    arrange(time) |> 
    mutate(IS_95 = interval_score(true_values = obs_dengue, 
                                  lower = LI_y2, upper = LS_y2, 
                                  interval_range = 95),
           IS_50 = interval_score(true_values = obs_dengue, 
                                  lower = LIb_y2, upper = LSb_y2,
                                  interval_range = 50),
           range_90 = LS_y2 - LI_y2,
           range_50 = LSb_y2 - LIb_y2
    ) |> 
    select(IS_95:range_50) |> 
    summarise_all(list(last = ~last(.),
                       mean = ~mean(.))) |> 
    select(!(starts_with("range") & ends_with("mean")))
  
  if(is.null(inla_output)){
    information_metrics <- c(0,0,0)
  }else{
    information_metrics <- c(
      "waic" = inla_output$waic$waic,
      "dic" = inla_output$dic$dic,
      "CPO" = -sum(inla_output$cpo$cpo, na.rm = T)
    )
  }

  
  
  results_vector <- c(
    information_metrics, 
    metrics_chik[1,],
    interval_metrics_chik[1,],
    metrics_dengue[1,],
    interval_metrics_dengue[1,]
  ) |> 
    as.numeric()
  
  names(results_vector) <- c("waic", "dic",  "CPO",  "MAE_y1",  "MAPE_y1", "MSE_y1", 
                             "IS_95_last_y1", "IS_50_last_y1", "range_90_last_y1", "range_50_last_y1",
                             "IS_95_mean_y1", "IS_50_mean_y1","MAE_y2","MAPE_y2","MSE_y2", 
                             "IS_95_last_y2","IS_50_last_y2", "range_90_last_y2", 
                             "range_50_last_y2", "IS_95_mean_y2",
                             "IS_50_mean_y2")
  
  return(results_vector)
  
}


selection_criteria2 <- function(summary_pred, data_real, inla_output){
  ## Criteria: waic, dic, length of the last interval, mean interval score, cpo?
  ## mse, mae, mape
  require(scoringutils)
  require(metrica)
  
  obs_chik <- data_real %>% 
    filter(time %in% summary_pred$time) %>% 
    pull(y1)
  
  obs_dengue <- data_real %>% 
    filter(time %in% summary_pred$time) %>% 
    pull(y2)
  
  
  metrics_chik <- metrica::metrics_summary(obs = obs_chik, 
                                           pred = summary_pred$Median_y1,
                                           metrics_list = c("MAPE", "MAE", "MSE"),
                                           type = "regression") |> 
    tidyr::pivot_wider(names_from = Metric, values_from = Score)
  
  metrics_dengue <- metrica::metrics_summary(obs = obs_dengue, 
                                             pred = summary_pred$Median_y2, 
                                             metrics_list = c("MAPE", "MAE", "MSE"),
                                             type = "regression") |> 
    tidyr::pivot_wider(names_from = Metric, values_from = Score)
  
  interval_metrics_chik <- summary_pred |> 
    arrange(time) |> 
    mutate(IS_95 = interval_score(true_values = obs_chik, 
                                  lower = LI_y1, upper = LS_y1, 
                                  interval_range = 95),
           IS_50 = interval_score(true_values = obs_chik, 
                                  lower = LIb_y1, upper = LSb_y1,
                                  interval_range = 50),
           range_90 = LS_y1 - LI_y1,
           range_50 = LSb_y1 - LIb_y1
    ) |> 
    select(IS_95:range_50) |> 
    summarise_all(list(last = ~last(.),
                       mean = ~mean(.))) |> 
    select(!(starts_with("range") & ends_with("mean")))
  
  interval_metrics_dengue <- summary_pred |> 
    arrange(time) |> 
    mutate(IS_95 = interval_score(true_values = obs_dengue, 
                                  lower = LI_y2, upper = LS_y2, 
                                  interval_range = 95),
           IS_50 = interval_score(true_values = obs_dengue, 
                                  lower = LIb_y2, upper = LSb_y2,
                                  interval_range = 50),
           range_90 = LS_y2 - LI_y2,
           range_50 = LSb_y2 - LIb_y2
    ) |> 
    select(IS_95:range_50) |> 
    summarise_all(list(last = ~last(.),
                       mean = ~mean(.))) |> 
    select(!(starts_with("range") & ends_with("mean")))
  
  
  interval_metrics_joint <- summary_pred |> 
    arrange(time) |> 
    mutate(IS_95_rel = interval_score(true_values = obs_dengue, 
                                  lower = LI_y2, upper = LS_y2, 
                                  interval_range = 95)/obs_dengue + 
             interval_score(true_values = obs_chik, 
                            lower = LI_y1, upper = LS_y1, 
                            interval_range = 95)/obs_chik ,
           IS_50_rel = interval_score(true_values = obs_dengue, 
                                      lower = LI_y2, upper = LS_y2, 
                                      interval_range = 50)/obs_dengue + 
             interval_score(true_values = obs_chik, 
                            lower = LI_y1, upper = LS_y1, 
                            interval_range = 50)/obs_chik,
           range_90_rel = (LS_y2 - LI_y2)/obs_dengue + (LS_y1 - LI_y1)/obs_chik ,
           range_50_rel = (LSb_y2 - LIb_y2)/obs_dengue + (LSb_y1 - LIb_y1)/obs_chik,
           
    ) |> 
    select(IS_95_rel:range_50_rel) |> 
    summarise_all(list(last = ~last(.),
                       mean = ~mean(.))) #|> 
    #select(!(starts_with("range") & ends_with("mean")))
  
  
  if(is.null(inla_output)){
    information_metrics <- c(0,0,0)
  }else{
    information_metrics <- c(
      "waic" = inla_output$waic$waic,
      "dic" = inla_output$dic$dic,
      "CPO" = -sum(inla_output$cpo$cpo, na.rm = T)
    )
  }
  
  
  results_vector <- c(
    information_metrics, 
    metrics_chik[1,],
    interval_metrics_chik[1,],
    metrics_dengue[1,],
    interval_metrics_dengue[1,],
    interval_metrics_chik[1,]
  ) |> 
    as.numeric()
  
  # names(results_vector) <- c("waic", "dic",  "CPO",  "MAE_y1",  "MAPE_y1", "MSE_y1", 
  #                            "IS_95_last_y1", "IS_50_last_y1", "range_90_last_y1", "range_50_last_y1",
  #                            "IS_95_mean_y1", "IS_50_mean_y1","MAE_y2","MAPE_y2","MSE_y2", 
  #                            "IS_95_last_y2","IS_50_last_y2", "range_90_last_y2", 
  #                            "range_50_last_y2", "IS_95_mean_y2",
  #                            "IS_50_mean_y2")
  
  return(list(results_vector,
         interval_metrics_joint))
  
}



information_metrics <- function(inla_output){
  information_metrics <- c(
    "waic" = inla_output$waic$waic,
    "dic" = inla_output$dic$dic,
    "CPO" = -sum(inla_output$cpo$cpo, na.rm = T)
  )
  return(information_metrics)
}


selection_criteria_summary <- function(summary_pred, data_real){

  require(scoringutils)
  require(metrica)
  
  y1_real <- data_real %>% 
    arrange(time) %>%
    filter(time %in% summary_pred$time) %>% 
    pull(y1)
  
  y2_real <- data_real %>% 
    arrange(time) %>%
    filter(time %in% summary_pred$time) %>% 
    pull(y2)
  
  summary_pred |> 
    arrange(time) |> 
    mutate(y1_real = y1_real,
           y2_real = y2_real,
           perc_error_y1     = abs((Median_y1 - y1_real)/y1_real),
           squared_error_y1  = (Median_y1 - y1_real)^2,
           absolute_error_y1 = abs(Median_y1 - y1_real),
           perc_error_y2     = abs((Median_y2 - y2_real)/y2_real),
           squared_error_y2  = (Median_y2 - y2_real)^2,
           absolute_error_y2 = abs(Median_y2 - y2_real),
           
           width_95_y1 = LS_y1 - LI_y1,
           width_50_y1 = LSb_y1 - LIb_y1,
           is_95_y1 = interval_score(true_values = y1_real, 
                                     lower = LI_y1, upper = LS_y1, 
                                     interval_range = 95, weigh = F),
           is_50_y1 = interval_score(true_values = y1_real, 
                                      lower = LIb_y1, upper = LSb_y1, 
                                      interval_range = 50, weigh = F),
           coverage_95_y1 = (y1_real <= LS_y1 & y1_real >= LI_y1),
           coverage_50_y1 = (y1_real <= LSb_y1 & y1_real >= LIb_y1),
           width_95_y2 = LS_y2 - LI_y2,
           width_50_y2 = LSb_y2 - LIb_y2,
           is_95_y2 = interval_score(true_values = y2_real, 
                                     lower = LI_y2, upper = LS_y2, 
                                     interval_range = 95, weigh = F),
           is_50_y2 = interval_score(true_values = y2_real, 
                                     lower = LIb_y2, upper = LSb_y2, 
                                     interval_range = 50, weigh = F),
           coverage_95_y2 = (y2_real <= LS_y2 & y2_real >= LI_y2),
           coverage_50_y2 = (y2_real <= LSb_y2 & y2_real >= LIb_y2),
           coverage_95_both = coverage_95_y1 & coverage_95_y2,
           coverage_50_both = coverage_50_y1 & coverage_50_y2,
           coverage_95_atleastone = coverage_95_y1 | coverage_95_y2,
           coverage_50_atleastone = coverage_50_y1 | coverage_50_y2,
           relative_width_95 = width_95_y1/y1_real + width_95_y2/y2_real,
           relative_width_50 = width_50_y1/y1_real + width_50_y2/y2_real,
           relative_is_95 = is_95_y1/y1_real + is_95_y2/y2_real,
           relative_is_50 = is_50_y1/y1_real + is_50_y2/y2_real
           ) -> results
  
  return(results)
  
}





