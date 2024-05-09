library(INLA)
library(ggplot2)
library(dplyr)
library(nowcaster)
library(lubridate)
library(patchwork)
library(scoringutils)
library(purrr)

## Funções
source("functions_multivar_nowcasting.R")

## Dados
dados <- readRDS("Dados\\den_chik_2017_2019_RJ.RData")
dados <- dados %>% janitor::clean_names()

info_dengue_data <- read.csv("info_dengue_rj.csv") |> 
  janitor::clean_names()

dengue_rj_google_trends <- read_csv("dengue_rj_google_trends.csv", skip = 2)
dengue_rj_google_trends <- janitor::clean_names(dengue_rj_google_trends)

M3 <- Y ~ -1 + alphai +
  + f(time,  model = "rw2", replicate = idrep,
      hyper = list(prec = list(param = c(0.001, 0.001)))) +
  + f(delay, model = "rw1", replicate = idrep,
      hyper = list(prec = list(param = c(0.001, 0.001)))) +
  + f(time1, model = "rw2",
      hyper = list(prec = list(param = c(0.001, 0.001)))) +
  + f(time2, copy  = "time1", fixed = TRUE,
      hyper = list(prec = list(param = c(0.001, 0.001)))) + 
  + f(delay1, model = "rw1",
      hyper = list(prec = list(param = c(0.001, 0.001)))) +
  + f(delay2, copy  = "delay1", fixed = TRUE,
      hyper = list(prec = list(param = c(0.001, 0.001)))) 

M0 <- Y ~ -1 + alphai +
  + f(time,  model = "rw2", replicate = idrep,
      hyper = list(prec = list(param = c(0.001, 0.001)))) +
  + f(delay, model = "rw1", replicate = idrep,
      hyper = list(prec = list(param = c(0.001, 0.001)))) 


Inicio_geral <- proc.time()
suppressMessages(
for(max_delay in c(15)){
  dados |>
    filter(dt_sin_pri > ymd("2017-01-01")) |>
    select(dt_digita, dt_sin_pri, arbo) |> 
    create_data(D = max_delay) |> 
    model_data(max_year = 2019, max_epiweek = 52) -> data_model
  
  ## Dados reais para comparação final
  dados |>
    filter(dt_sin_pri > ymd("2017-01-01")) |>
    select(dt_digita, dt_sin_pri, arbo) |> 
    create_data(D = max_delay, return_data_real = T) -> data_real
  
  
  for(window_size in c(70)){
    
    maxtime <- 116:140 #seq(60, max(data_model$time), by = 4)
    
    #### Janela window_size ####
    data_model_list <- purrr::map(maxtime, ~{
      maxdate <- ymd("2017-01-01") + weeks(.x-1)
      dados |>
        filter(dt_sin_pri > ymd("2017-01-01")) |>
        select(dt_digita, dt_sin_pri, arbo) |> 
        create_data(D = max_delay) |> 
        model_data(max_year = year(maxdate), 
                   max_epiweek = max(epiweek(maxdate), week(maxdate))) |> 
        filter(time <= .x & time > .x - window_size) #|>
        #mutate(time = time - min(time) + 1)
    })
    
    data_real_list <- purrr::map(maxtime, ~{
      # maxdate <- ymd("2017-01-01") + weeks(.x)
      data_real |> 
        ungroup() |> 
        mutate(time = rep(1:(nrow(data_real)/2), each = 2)) |> 
        filter(time <= .x & time > .x - window_size) %>% 
        tidyr::pivot_wider(names_from = "arbo", values_from = "y") %>% 
        rename(y1 = Chik, y2 = Den)
    })
    
    cat(paste("Ajustando modelos com delay", max_delay, "e janela", window_size, 
              "\nInício desse passo:", Sys.time()))
    inic <- proc.time()
    bivariate_fit <- map(data_model_list, 
                         ~{run_inla_model(data = .x, model = M0)})
    (tempo_bivariate_fit_30 = proc.time() - inic)
    
    
    file_suffix <- paste0("delay_", max_delay,"_window_", window_size, "_25runs_independente.RData")
    save(bivariate_fit, file = paste0("bivariate_fit_", file_suffix))
    
    cat(paste("Amostrando de modelos com delay", max_delay, "e janela", window_size, 
              "\nInício desse passo:", Sys.time()))
    inic <- proc.time()
    bivariate_samples <- map(bivariate_fit, 
                             ~{sample_predictive(.x$inla_output, 
                                                 joint_data = .x$joint_data)})
    (tempo_bivariate_samples_30 <- proc.time() - inic)
    
    
    # save(univariate_fit, file = "univariate_sliding_window_30_delay_10_final.RData")
    
    save(tempo_bivariate_fit_30, 
         tempo_bivariate_samples_30, file = paste0("tempos_bivariate", file_suffix))
    
    save(bivariate_samples, file = paste0("bivariate_samples_", file_suffix))
  }

})
(tempo_geral <- proc.time() - Inicio_geral)
# save(tempo_geral, "tempo_so_biv_geral.RData")



## Amanhã rodar janela univariado pro que der melhor desses entre os bivariados
## Já tenho né, posso só rodar os 4 e concatenar 
# 
# times_fit <- list() -> times_samples
# window_sizes <- c(30, 40, 50, 60, 70, 116)
# for(i in 1:6){
#   file_suffix <-  paste0("delay_", 20,"_window_", window_sizes[i], "_25runs.RData")
#   load(paste0("tempos_bivariate", file_suffix))
#   times_fit[[i]] <- tempo_bivariate_fit_30
#   times_samples[[i]] <- tempo_bivariate_samples_30
# }
# 
# sapply(times_fit, \(x) x[3]/60)
# sapply(times_samples, \(x) x[3]/60)
# 
# cbind(sapply(times_fit, \(x) x[3]/60),
#       sapply(times_samples, \(x) x[3]/60))





M4 <- Y ~ -1 + alphai +
  + f(time,  model = "rw2", replicate = idrep,
      hyper = list(prec = list(param = c(0.001, 0.001)))) +
  + f(delay, model = "rw1", replicate = idrep,
      hyper = list(prec = list(param = c(0.001, 0.001)))) +
  + f(time1, model = "rw2",
      hyper = list(prec = list(param = c(0.001, 0.001)))) +
  + f(time2, copy  = "time1", fixed = TRUE) +
  + f(delay1, model = "rw1",
      hyper = list(prec = list(param = c(0.001, 0.001)))) +
  + f(delay2, copy  = "delay1", fixed = TRUE) +
  + tweets_ef_dengue +
  + tweets_ef_chik +
  + f(id_tweets, tweets, model = "rw1", replicate = idrep,
      hyper = list(prec = list(param = c(0.001, 0.001))))

M5 <- Y ~ -1 + alphai +
  + f(time,  model = "rw2", replicate = idrep,
      hyper = list(prec = list(param = c(0.001, 0.001)))) +
  + f(delay, model = "rw1", replicate = idrep,
      hyper = list(prec = list(param = c(0.001, 0.001)))) +
  + f(time1, model = "rw2",
      hyper = list(prec = list(param = c(0.001, 0.001)))) +
  + f(time2, copy  = "time1", fixed = TRUE) +
  + f(delay1, model = "rw1",
      hyper = list(prec = list(param = c(0.001, 0.001)))) +
  + f(delay2, copy  = "delay1", fixed = TRUE) +
  + tweets_ef_dengue +
  + tweets_ef_chik 

M0 <- Y ~ -1 + alphai +
  + f(time,  model = "rw2", replicate = idrep,
      hyper = list(prec = list(param = c(0.001, 0.001)))) +
  + f(delay, model = "rw1", replicate = idrep) 






Inicio_geral <- proc.time()
suppressMessages(
  for(max_delay in c(15)){
    dados |>
      filter(dt_sin_pri > ymd("2017-01-01")) |>
      select(dt_digita, dt_sin_pri, arbo) |> 
      create_data(D = max_delay) |> 
      model_data(max_year = 2019, max_epiweek = 52) -> data_model
    
    ## Dados reais para comparação final
    dados |>
      filter(dt_sin_pri > ymd("2017-01-01")) |>
      select(dt_digita, dt_sin_pri, arbo) |> 
      create_data(D = max_delay, return_data_real = T) -> data_real
    
    
    for(window_size in c(139)){
      
      maxtime <- 140 #116:140
      
      #### Janela window_size ####
      data_model_list <- purrr::map(maxtime, ~{
        maxdate <- ymd("2017-01-01") + weeks(.x-1)
        dados |>
          filter(dt_sin_pri > ymd("2017-01-01")) |>
          select(dt_digita, dt_sin_pri, arbo) |> 
          create_data(D = max_delay) |> 
          model_data(max_year = year(maxdate), 
                     max_epiweek = max(epiweek(maxdate), week(maxdate))) |> 
          filter(time <= .x & time > .x - window_size) #|>
          #mutate(time = time - min(time) + 1)
      })
      
      data_real_list <- purrr::map(maxtime, ~{
        # maxdate <- ymd("2017-01-01") + weeks(.x)
        data_real |> 
          ungroup() |> 
          mutate(time = rep(1:(nrow(data_real)/2), each = 2)) |> 
          filter(time <= .x & time > .x - window_size) %>% 
          tidyr::pivot_wider(names_from = "arbo", values_from = "y") %>% 
          rename(y1 = Chik, y2 = Den)
      })
      
      tweets_list <- purrr::map(maxtime, ~{
        maxdate <- ymd("2017-01-01") + weeks(.x-1)
        info_dengue_data |> 
          arrange(ymd(data_ini_se)) |> 
          filter(data_ini_se <= maxdate) |> 
          pull(tweet) %>%
          .[(.x - window_size + 1):.x]
      })
      
      joint_data_list <- purrr::map2(data_model_list, tweets_list ,~{
        tweets <- log(.y + 1)
        
        
        y_df <- .x
        n <- nrow(y_df)
        
        y_bvar <- matrix(NA, nrow = 2 * n, ncol = 2)
        y_bvar[1:n, 1] <- y_df$y1_obs
        y_bvar[n + 1:n, 2] <- y_df$y2_obs
        
        I <- matrix(NA, nrow = 2 * n, ncol = 2)
        I[1:n, 1] <- 1
        I[n + 1:n, 2] <- 1
        
        
        fixed_covariate = data.frame(alphai = as.factor(c(rep(1, n), rep(2, n))),
                                     id_time1 = c(rep(1,n), rep(NA, n)),
                                     id_delay1 = c(rep(1,n), rep(NA, n)),
                                     id_time2 = c(rep(NA, n), rep(1, n)),
                                     id_delay2 = c(rep(NA, n), rep(1, n)),
                                     idrep = rep(1:2, each = n),
                                     id_tweets = c(y_df$time, y_df$time))
        
        
        random_covariate <- list(time1  = c(y_df$time, rep(NA, n)),
                                 delay1 = c(y_df$delay, rep(NA, n)),
                                 time2  = c(rep(NA, n), y_df$time),
                                 delay2 = c(rep(NA, n), y_df$delay),
                                 time = rep(y_df$time, 2),
                                 delay = rep(y_df$delay, 2),
                                 tweets = c(rep(tweets, each = 16), 
                                            rep(tweets, each = 16)),
                                 tweets_ef_dengue = c(rep(NA, n), 
                                                      rep(tweets, each = 16)),
                                 tweets_ef_chik = c(rep(tweets, each = 16), 
                                                    rep(NA, n)))
        
        
        Y = y_bvar
        joint_data <- cbind(Y, fixed_covariate, random_covariate)
      })
      
      
      cat(paste("Ajustando modelos com delay", max_delay, "e janela", window_size, 
                "\nInício desse passo:", Sys.time()))
      inic <- proc.time()
      bivariate_fit <- map2(data_model_list, joint_data_list, 
                           ~{run_inla_model(data = .x, model = M4,
                                            data_inla = .y)})
      (tempo_bivariate_fit_30 = proc.time() - inic)
      
      
      # file_suffix <- paste0("delay_", max_delay,"_window_", window_size, "_25runs_regressora_fixo.RData")
      # save(bivariate_fit, file = paste0("bivariate_fit_", file_suffix))
      # 
      # cat(paste("Amostrando de modelos com delay", max_delay, "e janela", window_size, 
      #           "\nInício desse passo:", Sys.time()))
      # inic <- proc.time()
      bivariate_samples <- map(bivariate_fit, 
                               ~{sample_predictive(.x$inla_output, 
                                                   joint_data = .x$joint_data)})
      # (tempo_bivariate_samples_30 <- proc.time() - inic)
      # 
      # 
      # # save(univariate_fit, file = "univariate_sliding_window_30_delay_10_final.RData")
      # 
      # save(tempo_bivariate_fit_30, 
      #      tempo_bivariate_samples_30, file = paste0("tempos_bivariate", file_suffix))
      # 
      # save(bivariate_samples, file = paste0("bivariate_samples_", file_suffix))
    }
    
  })
(tempo_geral <- proc.time() - Inicio_geral)






