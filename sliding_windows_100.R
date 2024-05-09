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
  create_data() |> 
  model_data(max_year = 2019, max_epiweek = 52) -> data_model

## Dados reais para comparação final
dados |>
  filter(dt_sin_pri > ymd("2017-01-01")) |>
  select(dt_digita, dt_sin_pri, arbo) |> 
  create_data(return_data_real = T) -> data_real

maxtime <- seq(100, max(data_model$time), by = 4)

data_model_list <- purrr::map(maxtime, ~{
  maxdate <- ymd("2017-01-01") + weeks(.x)
  dados |>
    filter(dt_sin_pri > ymd("2017-01-01")) |>
    select(dt_digita, dt_sin_pri, arbo) |> 
    create_data() |> 
    model_data(max_year = year(maxdate), 
               max_epiweek = max(epiweek(maxdate), week(maxdate))) |> 
    filter(time <= .x & time > .x - 100)
})

data_real_list <- purrr::map(maxtime, ~{
  # maxdate <- ymd("2017-01-01") + weeks(.x)
  data_real |> 
    ungroup() |> 
    mutate(time = rep(1:(nrow(data_real)/2), each = 2)) |> 
    filter(time <= .x & time > .x - 100) %>% 
    tidyr::pivot_wider(names_from = "arbo", values_from = "y") %>% 
    rename(y1 = Chik, y2 = Den)
})

# #### Joint data para inclusão de tweets ####
# info_dengue_data <- read.csv("info_dengue_rj.csv") |> 
#   janitor::clean_names()
# 
# info_dengue_data |> 
#   arrange(ymd(data_ini_se)) |> 
#   filter(se <= "201920") -> info_dengue_data
# 
# 
# tweets <- log(info_dengue_data$tweet + 1)
# 
# y_df <- data_model
# n <- nrow(y_df)
# 
# y_bvar <- matrix(NA, nrow = 2 * n, ncol = 2)
# y_bvar[1:n, 1] <- y_df$y1_obs
# y_bvar[n + 1:n, 2] <- y_df$y2_obs
# 
# I <- matrix(NA, nrow = 2 * n, ncol = 2)
# I[1:n, 1] <- 1
# I[n + 1:n, 2] <- 1
# 
# 
# fixed_covariate = data.frame(alphai = as.factor(c(rep(1, n), rep(2, n))),
#                              id_time1 = c(rep(1,n), rep(NA, n)),
#                              id_delay1 = c(rep(1,n), rep(NA, n)),
#                              id_time2 = c(rep(NA, n), rep(1, n)),
#                              id_delay2 = c(rep(NA, n), rep(1, n)),
#                              idrep = rep(1:2, each = n),
#                              id_tweets = c(y_df$time, y_df$time))
# 
# 
# random_covariate <- list(time1  = c(y_df$time, rep(NA, n)),
#                          delay1 = c(y_df$delay, rep(NA, n)),
#                          time2  = c(rep(NA, n), y_df$time),
#                          delay2 = c(rep(NA, n), y_df$delay),
#                          time = rep(y_df$time, 2),
#                          delay = rep(y_df$delay, 2),
#                          tweets = c(rep(tweets, each = 16), 
#                                     rep(tweets, each = 16)),
#                          tweets_ef_dengue = c(rep(NA, n), 
#                                               rep(tweets, each = 16)),
#                          tweets_ef_chik = c(rep(tweets, each = 16), 
#                                             rep(NA, n)))
# 
# 
# Y = y_bvar
# joint_data <- cbind(Y, fixed_covariate, random_covariate)

# library(furrr)
# plan(multisession, workers = 3)

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


univariate_fit <- map(data_model_list, 
                      ~{run_inla_model(.x, model = M0)}, 
                      .progress = TRUE)
(tempo = proc.time() - inic)


save(bivariate_fit, file = "bivariate_sliding_window100.RData")
save(univariate_fit, file = "univariate_sliding_window100.RData")

# load("bivariate_sliding_window100.RData")
# load("univariate_sliding_window100.RData")

inic <- proc.time()
bivariate_samples <- map(bivariate_fit, 
                         ~{sample_predictive(.x$inla_output, 
                                             joint_data = .x$joint_data)},
                         .progress = TRUE)


univariate_samples <- map(univariate_fit, 
                          ~{sample_predictive(.x$inla_output, 
                                              joint_data = .x$joint_data)},
                          .progress = TRUE)
tempo_100 = proc.time() - inic



save(bivariate_samples, file = "bivariate_samples_sliding_window100.RData")
save(univariate_samples, file = "univariate_samples_sliding_window100.RData")

## 63 minutos para amostrar

# load("bivariate_samples_sliding_window100.RData")
# load("univariate_samples_sliding_window100.RData")

# inic = proc.time()
# bivariate_samples[[12]] <- sample_predictive(bivariate_fit$inla_output, 
#                                              joint_data = bivariate_fit$joint_data)
# univariate_samples[[12]] <- sample_predictive(uniivariate_fit$inla_output, 
#                                                joint_data = univariate_fit$joint_data)
# proc.time() - inic


bivariate_summary <- map(bivariate_samples, bivar_nowcasting_summary)
univariate_summary <- map(univariate_samples, bivar_nowcasting_summary)


biv_scores_samples <- map2(bivariate_samples, data_real_list[1:15], scores_posterior) 
uni_scores_samples <- map2(univariate_samples, data_real_list[1:15], scores_posterior) 

biv_scores_samples |> 
  map(~summarise_all(.x, mean)) |> 
  reduce(rbind) -> biv_mean_scores

uni_scores_samples |> 
  map(~summarise_all(.x, mean)) |> 
  reduce(rbind) -> uni_mean_scores

ts.plot(biv_mean_scores$energy_score)
lines(uni_mean_scores$energy_score, col = 2)
legend("topleft", col = 1:2, lty = 1, legend = c("Bivariate", "Univariate"))

ts.plot(biv_mean_scores$variogram_score)
lines(uni_mean_scores$variogram_score, col = 2)
legend("topleft", col = 1:2, lty = 1, legend = c("Bivariate", "Univariate"))



criteria_biv <- mapply(\(x,y,z){selection_criteria(x,y,z$inla_output)},
                       bivariate_summary,
                       data_real_list,
                       bivariate_fit) |> t()

criteria_uni <- mapply(\(x,y,z){selection_criteria(x,y,z$inla_output)},
                       univariate_summary,
                       data_real_list,
                       univariate_fit) |> t() #map2(selection_criteria())



criteria_biv |> 
  as.data.frame() -> criteria_biv

criteria_uni |> 
  as.data.frame() -> criteria_uni


ts.plot(criteria_biv$IS_95_mean)
lines(criteria_uni$IS_95_mean, col = 2)

ts.plot(criteria_biv$range_90_last, ylim = c(0, 6750))
lines(criteria_uni$range_90_last, col = 2)

ts.plot(criteria_biv$MAPE)
lines(criteria_uni$MAPE, col = 2)

ts.plot(criteria_biv$MSE)
lines(criteria_uni$MSE, col = 2)


# bivariate_fit <- future_map(data_model_list, ~{
#   run_inla_model(data_model = .x, model = M3)
# }, seed = TRUE)
# 
# univariate_fit <- future_map(data_model_list, ~{
#   run_inla_model(.x, model = M0)
# }, .progress = T)











