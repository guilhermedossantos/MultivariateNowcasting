### Análise de sensibilidade
library(INLA)
library(ggplot2)
library(dplyr)
library(nowcaster)
library(lubridate)
library(patchwork)
library(scoringutils)
library(purrr)
library(tidyr)

## Funções
source("functions_multivar_nowcasting.R")


# Dados -------------------------------------------------------------------

## Dados
dados <- readRDS("Dados\\den_chik_2017_2019_RJ.RData")
dados <- dados %>% janitor::clean_names()

info_dengue_data <- read.csv("info_dengue_rj.csv") |> 
  janitor::clean_names()

dengue_rj_google_trends <- readr::read_csv("dengue_rj_google_trends.csv", skip = 2)
dengue_rj_google_trends <- dengue_rj_google_trends |> 
  janitor::clean_names()

# dados |>
#   filter(dt_sin_pri > ymd("2017-01-01")) |>
#   select(dt_digita, dt_sin_pri, arbo) |> 
#   create_data() |> 
#   model_data() -> data_model

## Dados reais para comparação final
dados |>
  filter(dt_sin_pri > ymd("2017-01-01")) |>
  select(dt_digita, dt_sin_pri, arbo) |> 
  create_data(return_data_real = T) -> data_real

dados |>
  filter(dt_sin_pri > ymd("2017-01-01")) |>
  select(dt_digita, dt_sin_pri, arbo) |> 
  create_data(D = 15) |> 
  model_data(max_year = 2019, max_epiweek = 36) -> data_model

## Dados reais para comparação final
dados |>
  filter(dt_sin_pri > ymd("2017-01-01")) |>
  select(dt_digita, dt_sin_pri, arbo) |> 
  create_data(D = max_delay, return_data_real = T) -> data_real




info_dengue_data |> 
  arrange(ymd(data_ini_se)) |> 
  filter(data_ini_se <= ymd("2017-01-01") + weeks(max(data_model$time)) - 1) |> 
  pull(tweet) -> tweets




# Declaração dos modelos --------------------------------------------------





# joint_model <- Y ~ -1 + alphai +
#   + f(time,  model = "rw1", replicate = idrep,
#       hyper = list(prec = list(param = c(0.001, 0.001)))) +
#   + f(delay, model = "rw1", replicate = idrep) +
#   + f(time1, model = "rw1") +
#   + f(time2, copy  = "time1", fixed = TRUE) + 
#   + f(delay1, model = "rw1") +
#   + f(delay2, copy  = "delay1", fixed = TRUE)



# Joint data --------------------------------------------------------------

tweets1 <- tweets  
tweets <- log(tweets + 1)


y_df <- data_model
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


# 
# fit_biv_uni <- run_inla_model(data = data_model, model = M5)
# samples_biv_uni <- sample_predictive(inla_output = fit_biv_uni$inla_output,
#                                      joint_data = fit_biv_uni$joint_data)
# 
# summary_biv_uni <- bivar_nowcasting_summary(samples_biv_uni)





# Prioris -----------------------------------------------------------------

HN.prior = "expression:
  tau0 = 0.001;
  sigma = exp(-theta/2);
  log_dens = log(2) - 0.5 * log(2 * pi) + 0.5 * log(tau0);
  log_dens = log_dens - 0.5 * tau0 * sigma^2;
  log_dens = log_dens - log(2) - theta / 2;
  return(log_dens);  
"

HC.prior  = "expression:
  sigma = exp(-theta/2);
  gamma = 25;
  log_dens = log(2) - log(pi) - log(gamma);
  log_dens = log_dens - log(1 + (sigma / gamma)^2);
  log_dens = log_dens - log(2) - theta / 2;
  return(log_dens);
"

HT.prior = "expression:
  sigma = exp(-theta/2);
  nu = 3;
  log_dens = 0 - 0.5 * log(nu * pi) - (-0.1207822);
  log_dens = log_dens - 0.5 * (nu + 1) * log(1 + sigma * sigma);
  log_dens = log_dens - log(2) - theta / 2;
  return(log_dens);
"

UN.prior = "expression:
  log_dens = 0 - log(2) - theta / 2;
  return(log_dens);
"



prior.list = list(
  default = list(prec = list(prior = "loggamma", param = c(1, 0.00005))),
  used = list(prec = list(prior = "loggamma", param = c(0.001, 0.001))),
  # gamma2 = list(prec = list(prior = "loggamma", param = c(0.0001, 0.0001))),
  half.normal = list(prec = list(prior = HN.prior)),
  half.cauchy = list(prec = list(prior = HC.prior)),
  h.t = list(prec = list(prior = HT.prior)),
  uniform = list(prec = list(prior = UN.prior))#,
  # pc.prec = list(prec = list(prior = "pc.prec", param = c(5, 0.01)))
  
) 



# Execução para cada priori -----------------------------------------------


models <- lapply(prior.list, function(tau.prior) {
  
  joint_model <- Y ~ -1 + alphai +
    + f(time,  model = "rw2", replicate = idrep,
        hyper = tau.prior) +
    + f(delay, model = "rw1", replicate = idrep,
        hyper = tau.prior) +
    + f(time1, model = "rw2",
        hyper = tau.prior) +
    + f(time2, copy  = "time1", fixed = TRUE) +
    + f(delay1, model = "rw1",
        hyper = tau.prior) +
    + f(delay2, copy  = "delay1", fixed = TRUE) +
    + tweets_ef_dengue +
    + tweets_ef_chik 
  
  control_family_list <- list(
    list(hyper = list("theta" = tau.prior)),
    list(hyper = list("theta" = tau.prior))
  )
  
  output0 <- INLA::inla(joint_model, family = rep("nbinomial", 2),
                        data = joint_data,
                        control.predictor = list(link = 1, compute = T),
                        control.compute = list(config = T, dic=T, waic = T,
                                               cpo = T),
                        control.family = control_family_list
  ) 
  
  
  return(output0)
})



# Gráficos ----------------------------------------------------------------
prior_names <- names(models)

## LABELS ERRADOS
#IDEIA: MOSTRAR ESTIMATIVA SSEGUNDO CADA UM, E FAZER O NOWCASTING OCM PELO MENOS
# 3 PRA VER SE TEM DIFERENÇA 

## phi1
lapply(1:length(prior.list), \(i) models[[i]]$marginals.hyperpar$`size for the nbinomial observations (1/overdispersion)` |> 
         as.data.frame() |> 
         mutate(prior = prior_names[i])) |> 
  reduce(rbind) |> 
  ggplot() +
  geom_line(aes(x = x, y = y, colour = prior), size = 1) +
  labs(x = expression(phi[1]),
       colour = "Prior",
       y = expression(paste(p, "(", phi[1], " | ", bold(n), ")"))) +
  scale_color_hue(h = c(0,360), labels =  c("IG(1, 5*10^-5)", 
                                            "Half-t",
                                            "Half-Cauchy",
                                            "Half-Normal",
                                            "Uniform",
                                            "IG(0.001, 0.001)")) +
  theme_bw(base_size = 13) -> posterior_phi1

## phi2
lapply(1:length(prior.list), \(i) models[[i]]$marginals.hyperpar$`size for the nbinomial observations (1/overdispersion)[2]` |> 
         as.data.frame() |> 
         mutate(prior = prior_names[i])) |> 
  reduce(rbind) |> 
  ggplot() +
  geom_line(aes(x = x, y = y, colour = prior), size = 1) +
  labs(x = expression(phi[2]),
       colour = "Prior",
       y = expression(paste(p, "(", phi[2], " | ", bold(n), ")"))) +
  scale_color_hue(h = c(0,360), labels =  c("IG(1, 5*10^-5)", 
                                           "Half-t",
                                           "Half-Cauchy",
                                           "Half-Normal",
                                           "Uniform",
                                           "IG(0.001, 0.001)")) +
  theme_bw(base_size = 13) -> posterior_phi2


## beta_t
lapply(1:length(prior.list), \(i) models[[i]]$marginals.hyperpar$`Precision for time` |> 
         as.data.frame() |> 
         mutate(prior = prior_names[i])) |> 
  reduce(rbind) |> 
  ggplot() +
  geom_line(aes(x = x, y = y, colour = prior), size = 1) +
  labs(x = expression(sigma[beta]^"-2"),
       colour = "Prior",
       y = expression(paste(p, "(", sigma[beta]^"-2", " | ", bold(n), ")"))) +
  scale_color_hue(h = c(0,360), labels =  c("IG(1, 5*10^-5)", 
                                            "Half-t",
                                            "Half-Cauchy",
                                            "Half-Normal",
                                            "Uniform",
                                            "IG(0.001, 0.001)")) +
  theme_bw(base_size = 13) -> posterior_beta

## delta_t
lapply(1:length(prior.list), \(i) models[[i]]$marginals.hyperpar$`Precision for time1` |> 
         as.data.frame() |> 
         mutate(prior = prior_names[i])) |> 
  reduce(rbind) |> 
  ggplot() +
  geom_line(aes(x = x, y = y, colour = prior), size = 1) +
  labs(x = expression(sigma[delta]^"-2"),
       colour = "Prior",
       y = expression(paste(p, "(",sigma[delta]^"-2", " | ", bold(n), ")"))) +
  # scale_x_continuous(breaks = c(0, 250, 500, 750, 1000)) +
  scale_color_hue(h = c(0,360), labels =  c("IG(1, 5*10^-5)", 
                                            "Half-t",
                                            "Half-Cauchy",
                                            "Half-Normal",
                                            "Uniform",
                                            "IG(0.001, 0.001)")) +
  theme_bw(base_size = 13) +
  theme(plot.margin = unit(c(.2, .5, .2, .5), "cm")) -> posterior_delta


## gamma
lapply(1:length(prior.list), \(i) models[[i]]$marginals.hyperpar$`Precision for delay` |> 
         as.data.frame() |> 
         mutate(prior = prior_names[i])) |> 
  reduce(rbind) |> 
  ggplot() +
  geom_line(aes(x = x, y = y, colour = prior), size = 1) +
  labs(x = expression(sigma[gamma]^"-2"),
       colour = "Prior",
       y = expression(paste(p, "(",sigma[gamma]^"-2", " | ", bold(n), ")"))) +
  scale_color_hue(h = c(0,360), labels =  c("IG(1, 5*10^-5)", 
                                            "Half-t",
                                            "Half-Cauchy",
                                            "Half-Normal",
                                            "Uniform",
                                            "IG(0.001, 0.001)")) +
  theme_bw(base_size = 13) -> posterior_gamma

## psi_d
lapply(1:length(prior.list), \(i) models[[i]]$marginals.hyperpar$`Precision for delay1` |> 
         as.data.frame() |> 
         mutate(prior = prior_names[i])) |> 
  reduce(rbind) |> 
  ggplot() +
  geom_line(aes(x = x, y = y, colour = prior), size = 1) +
  labs(x = expression(sigma[psi]^"-2"),
       colour = "Prior",
       y = expression(paste(p, "(",sigma[psi]^"-2", " | ", bold(n), ")"))) +
  scale_color_hue(h = c(0,360), labels =  c("IG(1, 5*10^-5)", 
                                            "Half-t",
                                            "Half-Cauchy",
                                            "Half-Normal",
                                            "Uniform",
                                            "IG(0.001, 0.001)")) +
  theme_bw(base_size = 13) -> posterior_psi


library(ggpubr)

graphs_list <- mget(ls(pattern = "posterior_")) 

ggarrange(plotlist = graphs_list[c(4,5,1,2,3,6)], ncol = 2, nrow = 3, common.legend = T,
          legend = "bottom", align = "hv")
ggsave("prior_sensitivity_graphs_legenda_sd.pdf", height = 9, width = 7.5)



# Estimativas -------------------------------------------------------------



lapply(1:length(prior.list), \(i) models[[i]]$summary.hyperpar |> 
       select(`0.5quant`, `0.025quant`,  `0.975quant`) |> 
         as.data.frame() |> 
         mutate(prior = prior_names[i])) -> estimates_priors

save(estimates_priors, file = "estimates_prior.RData")

# Nowcasting --------------------------------------------------------------
pred_samples <- lapply(models, \(x) sample_predictive(x, joint_data = joint_data))
summaries <- lapply(pred_samples, bivar_nowcasting_summary)


Sys.setlocale(locale = "English")

lapply(1:length(summaries), \(i){
  summaries[[i]] |> 
    mutate(Model = prior_names[i])
}) |> 
  purrr::reduce(rbind) |> 
  mutate(date = ymd("2017-01-01") + weeks(time - 1)) |>  
  ggplot() +
  aes(y = Median_y1, x = date, colour = Model, ymax = LS_y1, ymin = LI_y1,
      fill = Model) +
  geom_line() +
  geom_ribbon(alpha = 0.1) +
  labs(y = "Predictions for Chikungunya",
       x = "",
       colour = "Prior",
       fill = "Prior") +
  # scale_x_date(date_breaks = "1 month", date_labels ) + 
  scale_color_hue(h = c(0,360), labels =  c("Gamma(1, 5*10^-5)", 
                                            "Half-t",
                                            "Half-Cauchy",
                                            "Half-Normal",
                                            "Uniform",
                                            "Gamma(0.001, 0.001)")) +
  scale_fill_hue(h = c(0,360), labels =  c("Gamma(1, 5*10^-5)", 
                                           "Half-t",
                                           "Half-Cauchy",
                                           "Half-Normal",
                                           "Uniform",
                                           "Gamma(0.001, 0.001)")) +
  theme_bw(base_size = 13) -> pred_chik_priors

lapply(1:length(summaries), \(i){
  summaries[[i]] |> 
    mutate(Model = prior_names[i])
}) |> 
  purrr::reduce(rbind) |> 
  mutate(date = ymd("2017-01-01") + weeks(time - 1)) |> 
  ggplot() +
  aes(y = Median_y2, x = date, colour = Model, ymax = LS_y2, ymin = LI_y2,
      fill = Model) +
  geom_line() +
  geom_ribbon(alpha = 0.1) +
  labs(y = "Predictions for Dengue",
       x = "",
       colour = "Prior",
       fill = "Prior") +
  scale_color_hue(h = c(0,360), labels =  c("Gamma(1, 5*10^-5)", 
                                            "Half-t",
                                            "Half-Cauchy",
                                            "Half-Normal",
                                            "Uniform",
                                            "Gamma(0.001, 0.001)")) +
  scale_fill_hue(h = c(0,360), labels =  c("Gamma(1, 5*10^-5)", 
                                           "Half-t",
                                           "Half-Cauchy",
                                           "Half-Normal",
                                           "Uniform",
                                           "Gamma(0.001, 0.001)")) +
  theme_bw(base_size = 13) -> pred_dengue_priors




lapply(1:length(summaries), \(i){
  summaries[[i]] |> 
    mutate(Model = prior_names[i],
           range_95_y1 = LS_y1 - LI_y1,
           range_95_y2 = LS_y2 - LI_y2)
}) |> 
  purrr::reduce(rbind) |> 
  mutate(date = ymd("2017-01-01") + weeks(time - 1)) |> 
  ggplot() +
  aes(y = range_95_y1, x = date, colour = Model, fill = Model) +
  geom_line() +
  labs(y = "Width of the 95% interval (Chikungunya)",
       x = "",
       colour = "Prior") +
  scale_color_hue(h = c(0,360), labels =  c("Gamma(1, 5*10^-5)", 
                                            "Half-t",
                                            "Half-Cauchy",
                                            "Half-Normal",
                                            "Uniform",
                                            "Gamma(0.001, 0.001)")) +
  theme_bw(base_size = 13) -> ranges_chik_priors


lapply(1:length(summaries), \(i){
  summaries[[i]] |> 
    mutate(Model = prior_names[i],
           range_95_y1 = LS_y1 - LI_y1,
           range_95_y2 = LS_y2 - LI_y2)
}) |> 
  purrr::reduce(rbind) |> 
  mutate(date = ymd("2017-01-01") + weeks(time - 1)) |> 
  ggplot() +
  aes(y = range_95_y2, x = date, colour = Model, fill = Model) +
  geom_line() +
  labs(y = "Width of the 95% interval (Dengue)",
       x = "",
       colour = "Prior") +
  scale_color_hue(h = c(0,360), labels =  c("Gamma(1, 5*10^-5)", 
                                            "Half-t",
                                            "Half-Cauchy",
                                            "Half-Normal",
                                            "Uniform",
                                            "Gamma(0.001, 0.001)")) +
  theme_bw(base_size = 13) -> ranges_dengue_prios


ggarrange(pred_chik_priors, pred_dengue_priors,
          ranges_chik_priors, ranges_dengue_prios,
          ncol = 2, nrow = 2, common.legend = TRUE, 
          legend = "bottom")
ggsave("sensitivity_predictions_primeiro.pdf", height = 7, width = 8.5)



# lapply(1:length(summaries), \(i){
#   summaries[[i]] |> 
#     mutate(Model = prior_names[i])
# }) |> 
#   purrr::reduce(rbind) |> 
#   ggplot() +
#   aes(y = Median_y1, x = time, colour = Model, ymax = LS_y1, ymin = LI_y1,
#       fill = Model) +
#   geom_line() #+
  #geom_ribbon(alpha = 0.2)


# summaries |> 
#   lapply(\(summy) {
#     summy |> 
#       names()
#       ggplot() +
#       geom_line(aes(y = Median_y1, x = time)) +
#       geom_ribbon(aes(x = time, ymin = LI_y1, ymax = LS_y1), 
#                   fill = "gray50", alpha = 0.5)
#   })



# Nowcasting em um momento diferente --------------------------------------
prior_names <- names(prior.list)

dados |>
  filter(dt_sin_pri > ymd("2017-01-01")) |>
  select(dt_digita, dt_sin_pri, arbo) |> 
  create_data(D = 15) |> 
  model_data(max_year = 2019, max_epiweek = 22) -> data_model

## Dados reais para comparação final
dados |>
  filter(dt_sin_pri > ymd("2017-01-01")) |>
  select(dt_digita, dt_sin_pri, arbo) |> 
  create_data(D = max_delay, return_data_real = T) -> data_real


inic <- proc.time()
Tmax <- 126
tweets1 <- tweets[1:Tmax]  
tweets <- log(tweets[1:Tmax] + 1)


y_df <- data_model |> filter(time <= Tmax)
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


models2 <- lapply(prior.list, function(tau.prior) {
  
  joint_model <- Y ~ -1 + alphai +
    + f(time,  model = "rw2", replicate = idrep,
        hyper = tau.prior) +
    + f(delay, model = "rw1", replicate = idrep,
        hyper = tau.prior) +
    + f(time1, model = "rw2",
        hyper = tau.prior) +
    + f(time2, copy  = "time1", fixed = TRUE) +
    + f(delay1, model = "rw1",
        hyper = tau.prior) +
    + f(delay2, copy  = "delay1", fixed = TRUE) +
    + tweets_ef_dengue +
    + tweets_ef_chik 
  
  control_family_list <- list(
    list(hyper = list("theta" = tau.prior)),
    list(hyper = list("theta" = tau.prior))
  )
  
  output0 <- INLA::inla(joint_model, family = rep("nbinomial", 2),
                        data = joint_data,
                        control.predictor = list(link = 1, compute = T),
                        control.compute = list(config = T, dic=T, waic = T,
                                               cpo = T),
                        control.family = control_family_list
  ) 
  
  
  return(output0)
})


pred_samples2 <- lapply(models, \(x) sample_predictive(x, joint_data = joint_data))
summaries2 <- lapply(pred_samples, bivar_nowcasting_summary)

Sys.setlocale(locale = "English")

lapply(1:length(summaries2), \(i){
  summaries2[[i]] |> 
    mutate(Model = prior_names[i])
}) |> 
  purrr::reduce(rbind) |> 
  mutate(date = ymd("2017-01-01") + weeks(time - 1)) |>  
  ggplot() +
  aes(y = Median_y1, x = date, colour = Model, ymax = LS_y1, ymin = LI_y1,
      fill = Model) +
  geom_line() +
  geom_ribbon(alpha = 0.1) +
  labs(y = "Predictions for Chikungunya",
       x = "",
       colour = "Prior",
       fill = "Prior") +
  # scale_x_date(date_breaks = "1 month", date_labels ) + 
  scale_color_hue(h = c(0,360), labels =  c("IG(1, 5*10^-5)", 
                                            "Half-t",
                                            "Half-Cauchy",
                                            "Half-Normal",
                                            "Uniform",
                                            "IG(0.001, 0.001)")) +
  scale_fill_hue(h = c(0,360), labels =  c("IG(1, 5*10^-5)", 
                                           "Half-t",
                                           "Half-Cauchy",
                                           "Half-Normal",
                                           "Uniform",
                                           "IG(0.001, 0.001)")) +
  theme_bw(base_size = 13) -> pred_chik_priors

lapply(1:length(summaries2), \(i){
  summaries2[[i]] |> 
    mutate(Model = prior_names[i])
}) |> 
  purrr::reduce(rbind) |> 
  mutate(date = ymd("2017-01-01") + weeks(time - 1)) |> 
  ggplot() +
  aes(y = Median_y2, x = date, colour = Model, ymax = LS_y2, ymin = LI_y2,
      fill = Model) +
  geom_line() +
  geom_ribbon(alpha = 0.1) +
  labs(y = "Predictions for Dengue",
       x = "",
       colour = "Prior",
       fill = "Prior") +
  scale_color_hue(h = c(0,360), labels =  c("IG(1, 5*10^-5)", 
                                            "Half-t",
                                            "Half-Cauchy",
                                            "Half-Normal",
                                            "Uniform",
                                            "IG(0.001, 0.001)")) +
  scale_fill_hue(h = c(0,360), labels =  c("IG(1, 5*10^-5)", 
                                           "Half-t",
                                           "Half-Cauchy",
                                           "Half-Normal",
                                           "Uniform",
                                           "IG(0.001, 0.001)")) +
  theme_bw(base_size = 13) -> pred_dengue_priors




lapply(1:length(summaries2), \(i){
  summaries2[[i]] |> 
    mutate(Model = prior_names[i],
           range_95_y1 = LS_y1 - LI_y1,
           range_95_y2 = LS_y2 - LI_y2)
}) |> 
  purrr::reduce(rbind) |> 
  mutate(date = ymd("2017-01-01") + weeks(time - 1)) |> 
  ggplot() +
  aes(y = range_95_y1, x = date, colour = Model, fill = Model) +
  geom_line() +
  labs(y = "Width of the 95% interval (Chikungunya)",
       x = "",
       colour = "Prior") +
  scale_color_hue(h = c(0,360), labels =  c("IG(1, 5*10^-5)", 
                                            "Half-t",
                                            "Half-Cauchy",
                                            "Half-Normal",
                                            "Uniform",
                                            "IG(0.001, 0.001)")) +
  theme_bw(base_size = 13) -> ranges_chik_priors


lapply(1:length(summaries2), \(i){
  summaries2[[i]] |> 
    mutate(Model = prior_names[i],
           range_95_y1 = LS_y1 - LI_y1,
           range_95_y2 = LS_y2 - LI_y2)
}) |> 
  purrr::reduce(rbind) |> 
  mutate(date = ymd("2017-01-01") + weeks(time - 1)) |> 
  ggplot() +
  aes(y = range_95_y2, x = date, colour = Model, fill = Model) +
  geom_line() +
  labs(y = "Width of the 95% interval (Dengue)",
       x = "",
       colour = "Prior") +
  scale_color_hue(h = c(0,360), labels =  c("IG(1, 5*10^-5)", 
                                            "Half-t",
                                            "Half-Cauchy",
                                            "Half-Normal",
                                            "Uniform",
                                            "IG(0.001, 0.001)")) +
  theme_bw(base_size = 13) -> ranges_dengue_prios

# (tempo = proc.time() - inic)

ggarrange(pred_chik_priors, pred_dengue_priors,
          ranges_chik_priors, ranges_dengue_prios,
          ncol = 2, nrow = 2, common.legend = TRUE, 
          legend = "bottom")
ggsave("sensitivity_predictions_legenda_sd.pdf", height = 10, width = 12)



data_real |> 
  ungroup() |> 
  mutate(time = rep(1:(nrow(data_real)/2), each = 2)) |>
  tidyr::pivot_wider(names_from = "arbo", values_from = "y") %>% 
  rename(y1 = Chik, y2 = Den) -> data_real

lapply(summaries, selection_criteria_summary, data_real) |> 
  reduce(rbind) |> 
  mutate(prior = rep(prior_names, 
                     each = 15)) -> selection_summy


pred_samples |> 
  map(~scores_posterior(.x, data_real)) |> 
  reduce(rbind) |> 
  cbind(selection_summy)|>  
  group_by(prior) |> 
  summarise_all(.funs = list(mean = mean,
                             last = last)) |> 
  mutate(perc_error_y1_mean = perc_error_y1_mean*100,
         perc_error_y2_mean = perc_error_y2_mean*100) |> 
  select("MAPE (N_1)" = perc_error_y1_mean,
         "MAPE (N_2)" = perc_error_y2_mean,
         "Coverage (N_1)" = coverage_95_y1_mean,
         "Coverage (N_2)" = coverage_95_y2_mean,
         "Coverage (both)" = coverage_95_both_mean,
         "Last coverage (both)" = coverage_95_both_last,
         "mean relWidth" = relative_width_95_mean,
         "mean relIS" = relative_is_95_mean,
         "mean last relWidth" = relative_width_95_last,
         "mean last relIS" = relative_is_95_last,
         "mean Energy Score" = energy_score_mean,
         "mean Variogram Score" = variogram_score_mean,
         "last Energy Score" = energy_score_last,
         "last Variogram Score" = variogram_score_last
  ) |> 
  t() |> View()
kableExtra::kbl(format = "latex", booktabs = T, digits = 3) |> 
  add_header_above(c(" " = 1,
                     "Bivariate (NB)" = 1,
                     "Bivariate (poisson)" = 1,
                     "Univariate" = 1)) |>
  add_header_above(c(" " = 1, "Model" = 3))













