library(INLA)
library(ggplot2)
library(dplyr)
library(nowcaster)
library(lubridate)
library(patchwork)
library(scoringutils)

## Funções
source("functions_multivar_nowcasting.R")

## Dados
dados <- readRDS("Dados\\den_chik_2017_2019_RJ.RData")
dados <- dados %>% janitor::clean_names()

dados |>
  filter(dt_sin_pri > ymd("2017-01-01")) |>
  select(dt_digita, dt_sin_pri, arbo) |> 
  create_data() |> 
  model_data() -> data_model

## Dados reais para comparação final
dados |>
  filter(dt_sin_pri > ymd("2017-01-01")) |>
  select(dt_digita, dt_sin_pri, arbo) |> 
  create_data(return_data_real = T) -> data_real

data_real %>% 
  ungroup() %>% 
  filter(year < 2019 | epiweek <= 20) %>% 
  mutate(time = rep(unique(data_model$time), each = 2)) %>% 
  tidyr::pivot_wider(names_from = "arbo", values_from = "y") %>% 
  rename(y1 = Chik, y2 = Den) -> data_real

#### Joint data para inclusão de tweets ####
info_dengue_data <- read.csv("info_dengue_rj.csv") |> 
  janitor::clean_names()

info_dengue_data |> 
  arrange(ymd(data_ini_se)) |> 
    filter(se <= "201920") -> info_dengue_data


tweets <- log(info_dengue_data$tweet + 1)

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

 
## Modelos investigados ####

M1 <- Y ~ -1 + alphai +
  + f(time,  model = "rw2", replicate = idrep,
      hyper = list(prec = list(param = c(0.001, 0.001)))) +
  + f(delay, model = "rw1", replicate = idrep) +
  + f(time1, model = "rw2") +
  + f(time2, copy  = "time1", fixed = TRUE) +
  + f(delay1, model = "rw1") +
  + f(delay2, copy  = "delay1", fixed = TRUE) +
  + tweets_ef_dengue +
  + tweets_ef_chik +
  + f(id_tweets, tweets, model = "rw1", replicate = idrep)

M2 <- Y ~ -1 + alphai +
  + f(time,  model = "rw2", replicate = idrep,
      hyper = list(prec = list(param = c(0.001, 0.001)))) +
  + f(delay, model = "rw1", replicate = idrep) +
  + f(time1, model = "rw2") +
  + f(time2, copy  = "time1", fixed = TRUE) +
  + f(delay1, model = "rw1") +
  + f(delay2, copy  = "delay1", fixed = TRUE) +
  + tweets_ef_dengue +
  + tweets_ef_chik


M3 <- Y ~ -1 + alphai +
  + f(time,  model = "rw2", replicate = idrep,
      hyper = list(prec = list(param = c(0.001, 0.001)))) +
  + f(delay, model = "rw1", replicate = idrep) +
  + f(time1, model = "rw2") +
  + f(time2, copy  = "time1", fixed = TRUE) + 
  + f(delay1, model = "rw1") +
  + f(delay2, copy  = "delay1", fixed = TRUE) 

M4 <- Y ~ -1 + alphai +
  + f(time,  model = "rw2", replicate = idrep,
      hyper = list(prec = list(param = c(0.001, 0.001)))) +
  #+ f(delay, model = "rw1", replicate = idrep) +
  + f(time1, model = "rw2") +
  + f(time2, copy  = "time1", fixed = TRUE) + 
  + f(delay1, model = "rw1") +
  + f(delay2, copy  = "delay1", fixed = TRUE)

M5 <- Y ~ -1 + alphai +
  + f(time,  model = "rw2", replicate = idrep,
      hyper = list(prec = list(param = c(0.001, 0.001)))) +
  + f(delay, model = "rw1", replicate = idrep) 

M6 <- Y ~ -1 + alphai +
  + f(time1, model = "rw2") +
  + f(time2, copy  = "time1", fixed = TRUE) + 
  + f(delay1, model = "rw1") +
  + f(delay2, copy  = "delay1", fixed = TRUE)



#### Estimação ####
inic <- proc.time()
models <- mget(ls(pattern = "M"))

model_runs <- lapply(models, \(x) run_inla_model(data_model, model = x,
                                                 data_inla = joint_data))
predictive_samples <- lapply(model_runs, \(x) sample_predictive(x$inla_output, 
                                                                joint_data = joint_data))

tempo <- proc.time() - inic


load("results\\model_runs.RData")
load("results\\predictive_samples.RData")

total.summy <- lapply(predictive_samples, bivar_nowcasting_summary)

### Métricas ####
mapply(function(x,y)selection_criteria(x, data_real = data_real, y$inla_output), 
       total.summy,
       model_runs) |> 
  t() |> 
  round(2) -> metrics

colnames(metrics)[-c(3,4,6,7,8,10,12,13,15,16,17,19,21)] <- c("WAIC", "DIC", "MAPE", 
                                                              "Range of last 95% CI",
                                                              "MIS",
                                                              "MAPE",
                                                              "Range of last 95% CI","MIS")


kable(metrics[,-c(3,4,6,7,8,10,12,13,15,16,17,19,21)], "latex", booktabs = T) |> 
  kable_styling(latex_options = c("scale_down")) |> 
  add_header_above(c(" " = 3 , "Chikungunya correction" = 3, 
                     "Dengue correction" = 3)) # %>% 
  # column_spec(2, bold = )

data_no_delay <- data_model |> 
   group_by(time) |> 
   summarise(y1_obs = sum(y1_obs, na.rm = T),
             y2_obs = sum(y2_obs, na.rm = T)) |> 
  left_join(data_real |> select(time, date))

#### Gráficos ####
total.summy[3:4] |>
  purrr::map(~.x |> left_join(data_real, by = "time") |>
               ggplot(aes(x = date, y = Median_y1, col = 'Nowcasting')) +
               geom_line(data = data_real %>% filter(time > 100),
                         aes(date, y = y1, col = 'Observed'), size = 0.9) +
               geom_line(data = data_no_delay %>% filter(time > 100),
                         aes(date, y = y1_obs, col = 'Observed without delay'), size = 0.9) +
               geom_ribbon(aes(ymin = LI_y1, ymax = LS_y1, col = NA), alpha = 0.2, show.legend = F) +
               geom_line(size = 0.9) +
               labs(title = "Chikungunya") +
               theme_bw() +
               #theme(legend.position = "bottom", axis.text.x = element_text(angle = 90)) +
               scale_color_manual(values = c('grey50', 'black', 'darkred'), name = '')+
               scale_x_date(date_breaks = '2 weeks', date_labels = '%V/%y', name = 'Date in Weeks') +
               scale_y_continuous(breaks = seq(0, 12000, by = 2000), limits = c(0, 13000)) +
               labs(x = '', y = 'Nº Cases')+
               theme(legend.position = "bottom",
                     axis.title.x = element_text(size = 12),
                     axis.title.y = element_text(size = 12),
                     axis.text.x = element_text(size = 10,  angle = 90),
                     axis.text.y = element_text(size = 10),
                     legend.title = element_text(size = 12),
                     legend.text = element_text(size = 11)  
               )) -> grafs_chik 


total.summy[3:4] |>
  purrr::map(~.x |> left_join(data_real, by = "time") |>
               ggplot(aes(x = date, y = Median_y2, col = 'Nowcasting')) +
               geom_line(data = data_real %>% filter(time > 100),
                         aes(date, y = y2, col = 'Observed'), size = 0.9) +
               geom_line(data = data_no_delay %>% filter(time > 100),
                         aes(date, y = y2_obs, col = 'Observed without delay'), , size = 0.9) + 
               geom_ribbon(aes(ymin = LI_y2, ymax = LS_y2, col = NA), alpha = 0.2, show.legend = F) +
               geom_line(size = 0.9) +
               labs(title = "Dengue") +
               theme_bw() +
               #theme(legend.position = "bottom", axis.text.x = element_text(angle = 90)) +
               scale_color_manual(values = c('grey50', 'black', 'darkred'), name = '')+
               scale_x_date(date_breaks = '2 weeks', date_labels = '%V/%y', name = 'Date in Weeks')+
               scale_y_continuous(breaks = seq(0,6000, 1500), limits = c(0,7000)) +
               labs(x = '', y = 'Nº Cases') +
               theme(legend.position = "bottom",
                     axis.title.x = element_text(size = 12),
                     axis.title.y = element_text(size = 12),
                     axis.text.x = element_text(size = 10, angle = 90),
                     axis.text.y = element_text(size = 10),
                     legend.title = element_text(size = 12),
                     legend.text = element_text(size = 11)  
               )) -> grafs_dengue

library(latex2exp)
grafs_chik[[1]] <- grafs_chik[[1]] +
  ggtitle("Nowcasting Chikungunya cases with M3") +
  annotate(geom = "text", 
           label = TeX("Linear predictor = $\\alpha_i + \\beta_{t,i} + \\gamma_{d,i} + \\delta_t + \\psi_d$",
                       output = "character"),
           x = ymd("2019-02-01"), y= 12000,
           parse = T, size = 4.5)

grafs_chik[[2]] <- grafs_chik[[2]] +
  ggtitle("Nowcasting Chikungunya cases with M4") +
  annotate(geom = "text", 
           label = TeX("Linear predictor = $\\alpha_i + \\beta_{t,i}  + \\delta_t + \\psi_d $",
                       output = "character"),
           x = ymd("2019-02-01"), y= 12000,
           parse = T, size = 4.5)


grafs_dengue[[1]] <- grafs_dengue[[1]] +
  ggtitle("Nowcasting Dengue cases with M3") +
  annotate(geom = "text", 
           label = TeX("Linear predictor = $\\alpha_i + \\beta_{t,i} + \\gamma_{d,i} + \\delta_t + \\psi_d$",
                       output = "character"),
           x = ymd("2019-02-01"), y= 5000,
           parse = T, size = 4.5)

grafs_dengue[[2]] <- grafs_dengue[[2]] +
  ggtitle("Nowcasting Dengue cases with M4") +
  annotate(geom = "text", 
           label = TeX("Linear predictor = $\\alpha_i + \\beta_{t,i}  + \\delta_t + \\psi_d $",
                       output = "character"),
           x = ymd("2019-02-01"), y= 5000,
           parse = T, size = 4.5)







gridExtra::grid.arrange(grobs = c(grafs_chik, grafs_dengue))



gridExtra::grid.arrange(grobs = grafs_dengue)


output <- model_runs$M3$inla_output


# random_effects <- output$summary.random |> 
#   purrr::map2(names(output$summary.random), ~ .x |> mutate(effect = .y)) |> 
#   bind_rows()
# 
# 
# 
# 
# time_effects |> 
#   filter(stringr::str_detect(effect, "time")) |> 
#   select("ID", "0.5quant", "effect") |> 
#   tidyr::pivot_wider(values_from = `0.5quant`, names_from = effect) |> 
#   mutate()
#   ggplot() +
#   aes(x = time, 
#       y = effect,
#       col = disease) +
#   geom_line()

  
vectime = 1:124
effects_time <- data.frame(time = 1:124,
                           date = unique(data_real$date)[1:124],
                           time_chik = output$summary.random$time$`0.5quant`[vectime] + output$summary.random$time1$`0.5quant`,
                           time_dengue = output$summary.random$time$`0.5quant`[vectime + 124] + output$summary.random$time2$`0.5quant`,
                           time_uni = output$summary.random$time1$`0.5quant`)

effects_time |> 
  ggplot() +
  geom_line(aes(x = time, y = time_chik, col = "Chikungunya"), size = 1) +
  geom_line(aes(x = time, y = time_dengue, col = "Dengue"), size = 1) +
  # geom_line(aes(x = time, y = time_uni, col = "common", linetype = "Common"), 
  #           size = 1, linetype = 2) +
  labs(x = "Time", y = "Effect", colour = "Disease") + 
  theme_light() +
  theme(legend.position = "bottom",
         axis.title.x = element_text(size = 12),
         axis.title.y = element_text(size = 12),
         axis.text.x = element_text(size = 11),
         axis.text.y = element_text(size = 11),
         legend.title = element_text(size = 12),
         legend.text = element_text(size = 11)  
  ) -> gtime




effects_delay <- data.frame(delay = 1:16,
                            delay_chik = output$summary.random$delay$`0.5quant`[1:16] + output$summary.random$delay1$`0.5quant`,
                            delay_dengue = output$summary.random$delay$`0.5quant`[1:16 + 16] + output$summary.random$delay2$`0.5quant`)

effects_delay |> 
  ggplot() +
  geom_line(aes(x = delay, y = delay_chik, col = "Chikungunya"), size = 1) +
  geom_line(aes(x = delay, y = delay_dengue, col = "Dengue"), size = 1) +
  scale_x_continuous(breaks = seq(1, 16, by = 2), labels = seq(0, 15, by = 2)) + 
  labs(x = "Delay", y = "Effect", colour = "Disease") + 
  theme_light() +
  theme(legend.position = "bottom",
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11)  
  ) -> gdelay


library(ggpubr)
ggarrange(gtime, gdelay, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")



info_dengue_data <- read.csv("info_dengue_rj.csv") |> 
  janitor::clean_names()

ggplot(data = info_dengue_data |> filter(se <= "201952")) +
  geom_line(aes(x = ymd(data_ini_se), y = (tweet), colour = "Tweets"), size = 1) +
  labs(x = "Date in weeks", y = "Number of tweets") +
  scale_color_manual(values = "#65B0EC") +
  scale_x_date(breaks = "6 months") + 
  theme_light() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11))

ts.plot(aaam1$inla_output$summary.random$id_tweets$`0.5quant`)


