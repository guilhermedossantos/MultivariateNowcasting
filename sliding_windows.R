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

dados |>
  filter(dt_sin_pri > ymd("2017-01-01")) |>
  select(dt_digita, dt_sin_pri, arbo) |> 
  create_data(D = 10) |> 
  model_data(max_year = 2019, max_epiweek = 52) -> data_model

## Dados reais para comparação final
dados |>
  filter(dt_sin_pri > ymd("2017-01-01")) |>
  select(dt_digita, dt_sin_pri, arbo) |> 
  create_data(D = 10, return_data_real = T) -> data_real

maxtime <- 114:134 #seq(60, max(data_model$time), by = 4)

#### Janela 30 ####
data_model_list <- purrr::map(maxtime, ~{
  maxdate <- ymd("2017-01-01") + weeks(.x)
  dados |>
    filter(dt_sin_pri > ymd("2017-01-01")) |>
    select(dt_digita, dt_sin_pri, arbo) |> 
    create_data(D = 10) |> 
    model_data(max_year = year(maxdate), 
               max_epiweek = max(epiweek(maxdate), week(maxdate))) |> 
    filter(time <= .x & time > .x - 30)
})

data_real_list <- purrr::map(maxtime, ~{
  # maxdate <- ymd("2017-01-01") + weeks(.x)
  data_real |> 
    ungroup() |> 
    mutate(time = rep(1:(nrow(data_real)/2), each = 2)) |> 
    filter(time <= .x & time > .x - 30) %>% 
    tidyr::pivot_wider(names_from = "arbo", values_from = "y") %>% 
    rename(y1 = Chik, y2 = Den)
})


M3 <- Y ~ -1 + alphai +
  + f(time,  model = "rw2", replicate = idrep,
      hyper = list(prec = list(param = c(0.001, 0.001)))) +
  + f(delay, model = "rw1", replicate = idrep) +
  + f(time1, model = "rw2") +
  + f(time2, copy  = "time1", fixed = TRUE) + 
  + f(delay1, model = "rw1") +
  + f(delay2, copy  = "delay1", fixed = TRUE) 

M0 <- Y ~ -1 + alphai +
  + f(time,  model = "rw2", replicate = idrep,
      hyper = list(prec = list(param = c(0.001, 0.001)))) +
  + f(delay, model = "rw1", replicate = idrep) 


inic <- proc.time()
bivariate_fit <- map(data_model_list, 
                     ~{run_inla_model(data = .x, model = M3)},
                     .progress = TRUE)
(tempo_bivariate_fit_30 = proc.time() - inic)


inic_uni <- proc.time()
univariate_fit <- map(data_model_list, 
                      ~{run_inla_model(.x, model = M0)}, 
                      .progress = TRUE)
(tempo_univariate_fit_30 = proc.time() - inic_uni)


save(bivariate_fit, file = "bivariate_sliding_window_30_delay_10_final.RData")
save(univariate_fit, file = "univariate_sliding_window_30_delay_10_final.RData")



# load("bivariate_sliding_window.RData")
# load("univariate_sliding_window.RData")

inic <- proc.time()
bivariate_samples <- map(bivariate_fit, 
                         ~{sample_predictive(.x$inla_output, 
                                             joint_data = .x$joint_data)},
                         .progress = TRUE)
(tempo_bivariate_samples_30 <- proc.time() - inic)

inic2 <- proc.time()
univariate_samples <- map(univariate_fit, 
                          ~{sample_predictive(.x$inla_output, 
                                              joint_data = .x$joint_data)},
                          .progress = TRUE)
tempo_univariate_samples_30 = proc.time() - inic2


save(tempo_bivariate_fit_30, 
     tempo_univariate_fit_30, 
     tempo_bivariate_samples_30,
     tempo_univariate_samples_30, file = "tempos_sliding_delay_10_window__30.RData")

save(bivariate_samples, file = "bivariate_samples_sliding_window_30_delay_10_final.RData")
save(univariate_samples, file = "univariate_samples_sliding_window_30_delay_10_final.RData")



#### Janela 40 ####
data_model_list <- purrr::map(maxtime, ~{
  maxdate <- ymd("2017-01-01") + weeks(.x)
  dados |>
    filter(dt_sin_pri > ymd("2017-01-01")) |>
    select(dt_digita, dt_sin_pri, arbo) |> 
    create_data(D = 10) |> 
    model_data(max_year = year(maxdate), 
               max_epiweek = max(epiweek(maxdate), week(maxdate))) |> 
    filter(time <= .x & time > .x - 40)
})

data_real_list <- purrr::map(maxtime, ~{
  # maxdate <- ymd("2017-01-01") + weeks(.x)
  data_real |> 
    ungroup() |> 
    mutate(time = rep(1:(nrow(data_real)/2), each = 2)) |> 
    filter(time <= .x & time > .x - 40) %>% 
    tidyr::pivot_wider(names_from = "arbo", values_from = "y") %>% 
    rename(y1 = Chik, y2 = Den)
})


inic <- proc.time()
bivariate_fit <- map(data_model_list, 
                     ~{run_inla_model(data = .x, model = M3)},
                     .progress = TRUE)
(tempo_bivariate_fit_40 = proc.time() - inic)


inic_uni <- proc.time()
univariate_fit <- map(data_model_list, 
                      ~{run_inla_model(.x, model = M0)}, 
                      .progress = TRUE)
(tempo_univariate_fit_40 = proc.time() - inic_uni)


save(bivariate_fit, file = "bivariate_sliding_window_40_delay_10_final.RData")
save(univariate_fit, file = "univariate_sliding_window_40_delay_10_final.RData")



# load("bivariate_sliding_window.RData")
# load("univariate_sliding_window.RData")

inic <- proc.time()
bivariate_samples <- map(bivariate_fit, 
                     ~{sample_predictive(.x$inla_output, 
                                         joint_data = .x$joint_data)},
                     .progress = TRUE)
(tempo_bivariate_samples_40 <- proc.time() - inic)

inic2 <- proc.time()
univariate_samples <- map(univariate_fit, 
                          ~{sample_predictive(.x$inla_output, 
                                              joint_data = .x$joint_data)},
                          .progress = TRUE)
tempo_univariate_samples_40 = proc.time() - inic2


save(tempo_bivariate_fit_40, 
     tempo_univariate_fit_40, 
     tempo_bivariate_samples_40,
     tempo_univariate_samples_40, file = "tempos_sliding_delay_10_window__40.RData")

save(bivariate_samples, file = "bivariate_samples_sliding_window_40_delay_10_final.RData")
save(univariate_samples, file = "univariate_samples_sliding_window_40_delay_10_final.RData")




#### Janela 50 ####

data_model_list <- purrr::map(maxtime, ~{
  maxdate <- ymd("2017-01-01") + weeks(.x)
  dados |>
    filter(dt_sin_pri > ymd("2017-01-01")) |>
    select(dt_digita, dt_sin_pri, arbo) |> 
    create_data(D = 10) |> 
    model_data(max_year = year(maxdate), 
               max_epiweek = max(epiweek(maxdate), week(maxdate))) |> 
    filter(time <= .x & time > .x - 50)
})

data_real_list <- purrr::map(maxtime, ~{
  # maxdate <- ymd("2017-01-01") + weeks(.x)
  data_real |> 
    ungroup() |> 
    mutate(time = rep(1:(nrow(data_real)/2), each = 2)) |> 
    filter(time <= .x & time > .x - 50) %>% 
    tidyr::pivot_wider(names_from = "arbo", values_from = "y") %>% 
    rename(y1 = Chik, y2 = Den)
})



inic <- proc.time()
bivariate_fit <- map(data_model_list, 
                     ~{run_inla_model(data = .x, model = M3)},
                     .progress = TRUE)
(tempo_bivariate_fit_50 = proc.time() - inic)


inic_uni <- proc.time()
univariate_fit <- map(data_model_list, 
                      ~{run_inla_model(.x, model = M0)}, 
                      .progress = TRUE)
(tempo_univariate_fit_50 = proc.time() - inic_uni)


save(bivariate_fit, file = "bivariate_sliding_window_50_delay_10_final.RData")
save(univariate_fit, file = "univariate_sliding_window_50_delay_10_final.RData")



# load("bivariate_sliding_window.RData")
# load("univariate_sliding_window.RData")

inic <- proc.time()
bivariate_samples <- map(bivariate_fit, 
                         ~{sample_predictive(.x$inla_output, 
                                             joint_data = .x$joint_data)},
                         .progress = TRUE)
(tempo_bivariate_samples_50 <- proc.time() - inic)

inic2 <- proc.time()
univariate_samples <- map(univariate_fit, 
                          ~{sample_predictive(.x$inla_output, 
                                              joint_data = .x$joint_data)},
                          .progress = TRUE)
tempo_univariate_samples_50 = proc.time() - inic2


save(tempo_bivariate_fit_50, 
     tempo_univariate_fit_50, 
     tempo_bivariate_samples_50,
     tempo_univariate_samples_50, file = "tempos_sliding_delay_10_window__50.RData")

save(bivariate_samples, file = "bivariate_samples_sliding_window_50_delay_10_final.RData")
save(univariate_samples, file = "univariate_samples_sliding_window_50_delay_10_final.RData")




#### janela 60 ####
data_model_list <- purrr::map(maxtime, ~{
  maxdate <- ymd("2017-01-01") + weeks(.x)
  dados |>
    filter(dt_sin_pri > ymd("2017-01-01")) |>
    select(dt_digita, dt_sin_pri, arbo) |> 
    create_data(D = 10) |> 
    model_data(max_year = year(maxdate), 
               max_epiweek = max(epiweek(maxdate), week(maxdate))) |> 
    filter(time <= .x & time > .x - 60)
})

data_real_list <- purrr::map(maxtime, ~{
  # maxdate <- ymd("2017-01-01") + weeks(.x)
  data_real |> 
    ungroup() |> 
    mutate(time = rep(1:(nrow(data_real)/2), each = 2)) |> 
    filter(time <= .x & time > .x - 60) %>% 
    tidyr::pivot_wider(names_from = "arbo", values_from = "y") %>% 
    rename(y1 = Chik, y2 = Den)
})



inic <- proc.time()
bivariate_fit <- map(data_model_list, 
                     ~{run_inla_model(data = .x, model = M3)},
                     .progress = TRUE)
(tempo_bivariate_fit_60 = proc.time() - inic)


inic_uni <- proc.time()
univariate_fit <- map(data_model_list, 
                      ~{run_inla_model(.x, model = M0)}, 
                      .progress = TRUE)
(tempo_univariate_fit_60 = proc.time() - inic_uni)


save(bivariate_fit, file = "bivariate_sliding_window_60_delay_10_final.RData")
save(univariate_fit, file = "univariate_sliding_window_60_delay_10_final.RData")


inic <- proc.time()
bivariate_samples <- map(bivariate_fit, 
                         ~{sample_predictive(.x$inla_output, 
                                             joint_data = .x$joint_data)},
                         .progress = TRUE)
(tempo_bivariate_samples_60 <- proc.time() - inic)

inic2 <- proc.time()
univariate_samples <- map(univariate_fit, 
                          ~{sample_predictive(.x$inla_output, 
                                              joint_data = .x$joint_data)},
                          .progress = TRUE)
tempo_univariate_samples_60 = proc.time() - inic2


save(tempo_bivariate_fit_60, 
     tempo_univariate_fit_60, 
     tempo_bivariate_samples_60,
     tempo_univariate_samples_60, file = "tempos_sliding_delay_10_window__60.RData")

save(bivariate_samples, file = "bivariate_samples_sliding_window_60_delay_10_final.RData")
save(univariate_samples, file = "univariate_samples_sliding_window_60_delay_10_final.RData")



#### Full data ####
data_model_list <- purrr::map(maxtime, ~{
  maxdate <- ymd("2017-01-01") + weeks(.x)
  dados |>
    filter(dt_sin_pri > ymd("2017-01-01")) |>
    select(dt_digita, dt_sin_pri, arbo) |> 
    create_data(D = 10) |> 
    model_data(max_year = year(maxdate), 
               max_epiweek = max(epiweek(maxdate), week(maxdate))) |> 
    filter(time <= .x)
})

data_real_list <- purrr::map(maxtime, ~{
  # maxdate <- ymd("2017-01-01") + weeks(.x)
  data_real |> 
    ungroup() |> 
    mutate(time = rep(1:(nrow(data_real)/2), each = 2)) |> 
    filter(time <= .x ) %>% 
    tidyr::pivot_wider(names_from = "arbo", values_from = "y") %>% 
    rename(y1 = Chik, y2 = Den)
})



inic <- proc.time()
bivariate_fit <- map(data_model_list, 
                     ~{run_inla_model(data = .x, model = M3)},
                     .progress = TRUE)
(tempo_bivariate_fit_full_data = proc.time() - inic)


inic_uni <- proc.time()
univariate_fit <- map(data_model_list, 
                      ~{run_inla_model(.x, model = M0)}, 
                      .progress = TRUE)
(tempo_univariate_fit_full_data = proc.time() - inic_uni)


save(bivariate_fit, file = "bivariate_sliding_window_full_data_delay_10_final.RData")
save(univariate_fit, file = "univariate_sliding_window_full_data_delay_10_final.RData")


inic <- proc.time()
bivariate_samples <- map(bivariate_fit, 
                         ~{sample_predictive(.x$inla_output, 
                                             joint_data = .x$joint_data)},
                         .progress = TRUE)
(tempo_bivariate_samples_full_data <- proc.time() - inic)

inic2 <- proc.time()
univariate_samples <- map(univariate_fit, 
                          ~{sample_predictive(.x$inla_output, 
                                              joint_data = .x$joint_data)},
                          .progress = TRUE)
tempo_univariate_samples_full_data = proc.time() - inic2


save(tempo_bivariate_fit_full_data, 
     tempo_univariate_fit_full_data, 
     tempo_bivariate_samples_full_data,
     tempo_univariate_samples_full_data, file = "tempos_sliding_delay_10_window__full_data.RData")

save(bivariate_samples, file = "bivariate_samples_sliding_window_full_data_delay_10_final.RData")
save(univariate_samples, file = "univariate_samples_sliding_window_full_data_delay_10_final.RData")


### FULL DATA DELAY 15 ####
data_model_list <- purrr::map(maxtime, ~{
  maxdate <- ymd("2017-01-01") + weeks(.x)
  dados |>
    filter(dt_sin_pri > ymd("2017-01-01")) |>
    select(dt_digita, dt_sin_pri, arbo) |> 
    create_data(D = 15) |> 
    model_data(max_year = year(maxdate), 
               max_epiweek = max(epiweek(maxdate), week(maxdate))) |> 
    filter(time <= .x)
})

data_real_list <- purrr::map(maxtime, ~{
  # maxdate <- ymd("2017-01-01") + weeks(.x)
  data_real |> 
    ungroup() |> 
    mutate(time = rep(1:(nrow(data_real)/2), each = 2)) |> 
    filter(time <= .x ) %>% 
    tidyr::pivot_wider(names_from = "arbo", values_from = "y") %>% 
    rename(y1 = Chik, y2 = Den)
})



inic <- proc.time()
bivariate_fit <- map(data_model_list, 
                     ~{run_inla_model(data = .x, model = M3)},
                     .progress = TRUE)
(tempo_bivariate_fit_full_data = proc.time() - inic)


inic_uni <- proc.time()
univariate_fit <- map(data_model_list, 
                      ~{run_inla_model(.x, model = M0)}, 
                      .progress = TRUE)
(tempo_univariate_fit_full_data = proc.time() - inic_uni)


save(bivariate_fit, file = "bivariate_sliding_window_full_data_delay_15_final.RData")
save(univariate_fit, file = "univariate_sliding_window_full_data_delay_15_final.RData")


inic <- proc.time()
bivariate_samples <- map(bivariate_fit, 
                         ~{sample_predictive(.x$inla_output, 
                                             joint_data = .x$joint_data)},
                         .progress = TRUE)
(tempo_bivariate_samples_full_data <- proc.time() - inic)

inic2 <- proc.time()
univariate_samples <- map(univariate_fit, 
                          ~{sample_predictive(.x$inla_output, 
                                              joint_data = .x$joint_data)},
                          .progress = TRUE)
tempo_univariate_samples_full_data = proc.time() - inic2


save(tempo_bivariate_fit_full_data, 
     tempo_univariate_fit_full_data, 
     tempo_bivariate_samples_full_data,
     tempo_univariate_samples_full_data, file = "tempos_sliding_delay_15_window__full_data.RData")

save(bivariate_samples, file = "bivariate_samples_sliding_window_full_data_delay_15_final.RData")
save(univariate_samples, file = "univariate_samples_sliding_window_full_data_delay_15_final.RData")


#### full data delay 20 ####

data_model_list <- purrr::map(maxtime, ~{
  maxdate <- ymd("2017-01-01") + weeks(.x)
  dados |>
    filter(dt_sin_pri > ymd("2017-01-01")) |>
    select(dt_digita, dt_sin_pri, arbo) |> 
    create_data(D = 20) |> 
    model_data(max_year = year(maxdate), 
               max_epiweek = max(epiweek(maxdate), week(maxdate))) |> 
    filter(time <= .x)
})

data_real_list <- purrr::map(maxtime, ~{
  # maxdate <- ymd("2017-01-01") + weeks(.x)
  data_real |> 
    ungroup() |> 
    mutate(time = rep(1:(nrow(data_real)/2), each = 2)) |> 
    filter(time <= .x ) %>% 
    tidyr::pivot_wider(names_from = "arbo", values_from = "y") %>% 
    rename(y1 = Chik, y2 = Den)
})



inic <- proc.time()
bivariate_fit <- map(data_model_list, 
                     ~{run_inla_model(data = .x, model = M3)},
                     .progress = TRUE)
(tempo_bivariate_fit_full_data = proc.time() - inic)


inic_uni <- proc.time()
univariate_fit <- map(data_model_list, 
                      ~{run_inla_model(.x, model = M0)}, 
                      .progress = TRUE)
(tempo_univariate_fit_full_data = proc.time() - inic_uni)


save(bivariate_fit, file = "bivariate_sliding_window_full_data_delay_20_final.RData")
save(univariate_fit, file = "univariate_sliding_window_full_data_delay_20_final.RData")


inic <- proc.time()
bivariate_samples <- map(bivariate_fit, 
                         ~{sample_predictive(.x$inla_output, 
                                             joint_data = .x$joint_data)},
                         .progress = TRUE)
(tempo_bivariate_samples_full_data <- proc.time() - inic)

inic2 <- proc.time()
univariate_samples <- map(univariate_fit, 
                          ~{sample_predictive(.x$inla_output, 
                                              joint_data = .x$joint_data)},
                          .progress = TRUE)
tempo_univariate_samples_full_data = proc.time() - inic2


save(tempo_bivariate_fit_full_data, 
     tempo_univariate_fit_full_data, 
     tempo_bivariate_samples_full_data,
     tempo_univariate_samples_full_data, file = "tempos_sliding_delay_20_window__full_data.RData")

save(bivariate_samples, file = "bivariate_samples_sliding_window_full_data_delay_20_final.RData")
save(univariate_samples, file = "univariate_samples_sliding_window_full_data_delay_20_final.RData")








##### graficos #####




# load("bivariate_samples_sliding_window.RData")
# load("univariate_samples_sliding_window.RData")

# inic = proc.time()
# bivariate_samples[[12]] <- sample_predictive(bivariate_fit$inla_output, 
#                                              joint_data = bivariate_fit$joint_data)
# univariate_samples[[12]] <- sample_predictive(uniivariate_fit$inla_output, 
#                                                joint_data = univariate_fit$joint_data)
# proc.time() - inic


# bivariate_summary <- map(bivariate_samples, bivar_nowcasting_summary)
# univariate_summary <- map(univariate_samples, bivar_nowcasting_summary)
# 
# 
# biv_scores_samples <- map2(bivariate_samples, data_real_list, scores_posterior) 
# uni_scores_samples <- map2(univariate_samples, data_real_list, scores_posterior) 
# 
# biv_scores_samples |> 
#   map(~summarise_all(.x, mean)) |> 
#   reduce(rbind) -> biv_mean_scores
# 
# uni_scores_samples |> 
#   map(~summarise_all(.x, mean)) |> 
#   reduce(rbind) -> uni_mean_scores
# 
# ts.plot(biv_mean_scores$energy_score)
# lines(uni_mean_scores$energy_score, col = 2)
# legend("topleft", col = 1:2, lty = 1, legend = c("Bivariate", "Univariate"))
# 
# ts.plot(biv_mean_scores$variogram_score)
# lines(uni_mean_scores$variogram_score, col = 2)
# legend("topleft", col = 1:2, lty = 1, legend = c("Bivariate", "Univariate"))
# 
# 
# 
# criteria_biv <- mapply(\(x,y,z){selection_criteria(x,y,z$inla_output)},
#                        bivariate_summary,
#                        data_real_list,
#                        bivariate_fit) |> t()
# 
# criteria_uni <- mapply(\(x,y,z){selection_criteria(x,y,z$inla_output)},
#                        univariate_summary,
#                        data_real_list,
#                        univariate_fit) |> t() #map2(selection_criteria())
# 
# 
# 
# criteria_biv |> 
#   as.data.frame() -> criteria_biv
# 
# criteria_uni |> 
#   as.data.frame() -> criteria_uni
# 
# 
# ts.plot(criteria_biv$IS_95_mean_y1)
# lines(criteria_uni$IS_95_mean_y1, col = 2)
# 
# ts.plot(criteria_biv$range_90_last_y1, ylim = c(0, 6750))
# lines(criteria_uni$range_90_last_y1, col = 2)
# 
# y1_last <- sapply(data_real_list,\(x){last(x$y1)})
# y2_last <- sapply(data_real_list,\(x){last(x$y2)})
# 
# ts.plot(criteria_biv$range_90_last_y1/y1_last + criteria_biv$range_90_last_y2/y2_last)
# lines(criteria_uni$range_90_last_y1/y1_last + criteria_uni$range_90_last_y2/y2_last, col = 2)
# 
# ts.plot(criteria_biv$range_90_last_y1, ylim = c(0, 6750))
# lines(criteria_uni$range_90_last_y1, col = 2)
# 
# 
# ts.plot(criteria_biv$MAPE_y1)
# lines(criteria_uni$MAPE_y1, col = 2)
# 
# ts.plot(criteria_biv$MSE_y1)
# lines(criteria_uni$MSE_y1, col = 2)
# 
# 
# # bivariate_fit <- future_map(data_model_list, ~{
# #   run_inla_model(data_model = .x, model = M3)
# # }, seed = TRUE)
# # 
# # univariate_fit <- future_map(data_model_list, ~{
# #   run_inla_model(.x, model = M0)
# # }, .progress = T)
# 
# 
# bind_rows(bivariate_summary[[12]], univariate_summary[[12]]) |> 
#   mutate(model = rep(c("Bivariate", "Univariate"), each = 20))
#   mutate(limit = substr(limit, 1, nchar(limit) - 1)) |> 
#   filter(limit %in% c("LI", "LS"),
#          name == "y2") |> 
#   pivot_wider(names_from = limit, values_from = value) |> 
#   mutate(range_90 = LS - LI) |> 
#   ggplot() +
#   aes(x = time, y = range_90, colour = model) +
#   geom_line() +
#   labs(title = "Dengue") +
#   theme_bw() 
# 
# 
# 
# 
# 



