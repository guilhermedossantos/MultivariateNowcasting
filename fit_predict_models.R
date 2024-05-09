### Funções para ajuste e predição com os modelos bivariado e univariado

## gg.age para as contagens com base na preditiva ao longo do tempo
gg.age.bivar <- function(x, dados.gg, idx, n = 1984){
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
} #N = 1984 VAI DAR ERRO

gg.age <- function(x, dados.gg, idx){
  data.aux <- dados.gg
  Tmin <- min(dados.gg$time[idx])
  data.aux$Y[idx] <- x
  data.aggregated <- data.aux |>
    ## Selecionando apenas os dias faltantes a partir
    ## do domingo da respectiva ultima epiweek
    ## com dados faltantes
    dplyr::filter(time >= Tmin  ) |>
    dplyr::group_by(time) |>
    dplyr::summarise(
      Y = sum(Y), .groups = "keep"
    )
  data.aggregated
}


fit_bivar <- function(data_model, evol_time = "rw2", 
                      evol_delay = "rw1"){
  
  data_model %>% 
    distinct(epiweek, year) %>% 
    mutate(time = 1:n()) %>% 
    inner_join(data_model, by = c("epiweek", "year")) -> data_model
  
  y_new <- data_model %>% 
    mutate(report_time = time + delay) %>% 
    select(time, delay, report_time, y1_obs = chik_obs, y2_obs = den_obs, y1 = Chik, y2 = Den)
  
  
  y_df <- y_new
  
  #### MODELO BIVARIADO ####
  n <- nrow(y_df)
  
  y_bvar <- matrix(NA, nrow = 2 * n, ncol = 2)
  y_bvar[1:n, 1] <- y_new$y1_obs
  y_bvar[n + 1:n, 2] <- y_new$y2_obs
  
  I <- matrix(NA, nrow = 2 * n, ncol = 2)
  I[1:n, 1] <- 1
  I[n + 1:n, 2] <- 1
  
  
  fixed_covariate = data.frame(alphai = as.factor(c(rep(1, n), rep(2, n))),
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
  
  
  Y = y_bvar
  joint_data <- cbind(Y, fixed_covariate, random_covariate)
  
  
  joint_model <- Y ~ -1 + alphai +
    + f(time,  model = evol_time, replicate = idrep,
        hyper = list(prec = list(param = c(0.001, 0.001)))) +
    + f(delay, model = evol_delay, replicate = idrep) +
    + f(time1, model = evol_time) +
    + f(time2, copy  = "time1", fixed = TRUE) + 
    + f(delay1, model = evol_delay) +
    + f(delay2, copy  = "delay1", fixed = TRUE) 
  
  
  output0 <- INLA::inla(joint_model, family = rep("nbinomial", 2),
                        data = joint_data,
                        control.predictor = list(link = 1, compute = T),
                        control.compute = list( config = T, waic=F, dic=F),
                        control.family = list(
                          list(hyper = list("theta" = list(prior = "loggamma", param = c(0.001, 0.001)))),
                          list(hyper = list("theta" = list(prior = "loggamma", param = c(0.001, 0.0001))))
                        )
                        
  )
  
  
  #### preditivas - bivariado ####
  samples <- INLA::inla.posterior.sample(n = 1000, output0)
  
  
  
  index_missing <- which(is.na(joint_data$`1`) & is.na(joint_data$`2`))
  
  
  vector_samples_bvar <- lapply(X = samples,
                                FUN = function(x, idx = index_missing, n = nrow(y_df)){
                                  idx_1 <- idx[idx <= n]
                                  idx_2 <- idx[idx > n]
                                  unif.log = 1
                                  
                                  y1 <- stats::rnbinom(n = idx_1,
                                                       mu = exp(x$latent[idx_1]),
                                                       size = x$hyperpar[1]
                                  ) * unif.log
                                  
                                  y2 <- stats::rnbinom(n = idx_2,
                                                       mu = exp(x$latent[idx_2]),
                                                       size = x$hyperpar[2]
                                  ) * unif.log
                                  
                                  return(cbind(y1,y2))
                                } )
  
  
  
  ## Step 4: Applying the age aggregation on each posterior
  tibble.samples.0 <- lapply(X = vector_samples_bvar,
                             FUN = gg.age.bivar,
                             dados = joint_data,
                             idx = index_missing,
                             n = nrow(y_df))
  
  
  srag.pred.0 <- dplyr::bind_rows(tibble.samples.0, .id = "sample")
  
  total.summy <- srag.pred.0 |>
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

fit_uni <- function(data_model, disease = "dengue", evol_time = "rw2", 
                    evol_delay = "rw1"){
  data_model %>% 
    distinct(epiweek, year) %>% 
    mutate(time = 1:n()) %>% 
    inner_join(data_model, by = c("epiweek", "year")) -> data_model
  
  y_new <- data_model %>% 
    mutate(report_time = time + delay) %>% 
    select(time, delay, report_time, y1_obs = chik_obs, y2_obs = den_obs, y1 = Chik, y2 = Den)
  
  if(disease %in% c("chik","chikungunya")){
    dados_uni <- y_new %>% 
      select(time, delay, report_time, Y = y1_obs, count = y1)
  }else if(disease %in% c("den", "dengue")){
    dados_uni <- y_new %>% 
      select(time, delay, report_time, Y = y2_obs, count = y2)
  }
  
  model <- Y ~ 1 +
    f(time,
      model = evol_time,
      hyper = list("prec" = list(prior = "loggamma",
                                 param = c(0.001, 0.001))
      )) +
    f(delay, model = evol_delay,
      hyper = list("prec" = list(prior = "loggamma",
                                 param = c(0.001, 0.001)))
    )
  
  output_uni <- INLA::inla(model, family = "nbinomial",
                            data = dados_uni,
                            control.predictor = list(link = 1, compute = T),
                            control.compute = list(config = T, waic=F, dic=F),
                            control.family = list(
                              hyper = list("theta" = list(prior = "loggamma",
                                                          param = c(0.001,
                                                                    0.001))
                              )
                            )
  )
  
  
  index.missing <- which(is.na(dados_uni$Y))
  
  
  ## Step 1: Sampling from the approximate posterior distribution using INLA
  srag.samples0.list <- INLA::inla.posterior.sample(n = 1000, output_uni)
  
  ## Give a parameter to trajectories, TO-DO
  
  ## Step 2: Sampling the missing triangle from the likelihood using INLA estimates
  vector.samples0 <- lapply(X = srag.samples0.list,
                            FUN = function(x, idx = index.missing){
                              # if(zeroinflated){
                              #   unif.log <- as.numeric(runif(idx,0,1) < x$hyperpar[2])
                              # }else{
                              unif.log = 1
                              # }
                              stats::rnbinom(n = idx,
                                             mu = exp(x$latent[idx]),
                                             size = x$hyperpar[1]
                              ) * unif.log
                            } )
  
  ## Step 3: Calculate N_{a,t} for each triangle sample {N_{t,a} : t=Tactual-Dmax+1,...Tactual}
  
  
  ## Step 4: Applying the age aggregation on each posterior
  tibble.samples.0 <- lapply( X = vector.samples0,
                              FUN = gg.age,
                              dados = dados_uni,
                              idx = index.missing)
  
  srag.pred.0 <- dplyr::bind_rows(tibble.samples.0, .id = "sample")
  
  
  total.summy <- srag.pred.0 |>
    dplyr::group_by(time, sample) |>
    dplyr::summarise(Y = sum(Y, na.rm = T)) |>
    dplyr::group_by(time) |>
    dplyr::summarise(Median = stats::median(Y, na.rm = T),
                     LI = stats::quantile(Y, probs = 0.025, na.rm = T),
                     LS = stats::quantile(Y, probs = 0.975, na.rm = T),
                     LIb = stats::quantile(Y, probs = 0.25, na.rm = T),
                     LSb = stats::quantile(Y, probs = 0.75, na.rm = T),
                     .groups = "drop")
  

  output <- total.summy
  
  return(output)
}

calculate_score <- function(output, data_real, disease = "dengue", bivar = F){
  require(scoringutils)
  
  if(disease %in% c("dengue", "den")){
    obs_data <- data_real %>% 
      filter(time > nrow(data_real) - nrow(output)) %>% 
      pull(y2)
    
    if(bivar){
      output <- output |> 
        rename(LI = LI_y2,
               LS = LS_y2,
               LIb = LIb_y2,
               LSb = LSb_y2)
    }
      
  }else if(disease %in% c("chik","chikungunya")){
    obs_data <- data_real %>% 
      filter(time > nrow(data_real) - nrow(output)) %>% 
      pull(y1)
    if(bivar){
      output <- output |> 
        rename(LI = LI_y1,
               LS = LS_y1,
               LIb = LIb_y1,
               LSb = LSb_y1)
      }
  }
  
  output |> 
    mutate(score_95 = interval_score(true_values = obs_data, 
                                     lower = LI, upper = LS, interval_range = 95),
           range_95 = LS- LI,
           score_50 = interval_score(true_values = obs_data, 
                                     lower = LIb, upper = LSb, interval_range = 50),
           range_50 = LSb - LIb) 
  
}




