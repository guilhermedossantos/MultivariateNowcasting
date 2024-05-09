library(INLA)
library(ggplot2)
library(dplyr)
library(nowcaster)
library(lubridate)
library(patchwork)
# library(scoringutils)

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


## Modelos 
M1 <- Y ~ -1 + alphai +
  + f(time,  model = "rw2", replicate = idrep,
      hyper = list(prec = list(param = c(0.001, 0.001)))) +
  + f(delay, model = "rw1", replicate = idrep) +
  + f(time1, model = "rw2") +
  + f(time2, copy  = "time1", fixed = TRUE) + 
  + f(delay1, model = "rw1") +
  + f(delay2, copy  = "delay1", fixed = TRUE) 

M2 <- Y ~ -1 + alphai +
  + f(time,  model = "rw2", replicate = idrep,
      hyper = list(prec = list(param = c(0.001, 0.001)))) +
  + f(delay, model = "rw1", replicate = idrep) 

M3 <- Y ~ -1 + alphai +
  + f(time,  model = "rw2", replicate = idrep,
      hyper = list(prec = list(param = c(0.001, 0.001)))) +
  #+ f(delay, model = "rw1", replicate = idrep) +
  + f(time1, model = "rw2") +
  + f(time2, copy  = "time1", fixed = TRUE) + 
  + f(delay1, model = "rw1") +
  + f(delay2, copy  = "delay1", fixed = TRUE)

M4 <- Y ~ -1 + alphai +
  + f(time1, model = "rw2") +
  + f(time2, copy  = "time1", fixed = TRUE) + 
  + f(delay1, model = "rw1") +
  + f(delay2, copy  = "delay1", fixed = TRUE)

# M5 <- Y ~ -1 + alphai +
#   + f(time,  model = "rw2", replicate = idrep,
#       hyper = list(prec = list(param = c(0.001, 0.001)))) +
#   + f(delay, model = "rw1", replicate = idrep) +
#   + f(time1, model = "rw2") +
#   + f(time2, copy  = "time1", fixed = TRUE) + 
#   + f(delay1, model = "rw1") +
#   + f(delay2, copy  = "delay1", fixed = TRUE) #+
#   + tweets

inic <- proc.time()
models <- mget(ls(pattern = "M"))

model_runs <- lapply(models, \(x) run_inla_model(data_model, model = x))
predictive_samples <- lapply(model_runs, \(x) sample_predictive(x$inla_output, joint_data = joint_data))

tempo <- proc.time() - inic


total.summy <- lapply(predictive_samples, bivar_nowcasting_summary)

# dados |>
#   filter(dt_sin_pri > ymd("2017-01-01")) |>
#   select(dt_digita, dt_sin_pri, arbo) |>
#   create_data(return_data_real = T) -> data_real
# 
# data_real %>%
#   ungroup() %>%
#   filter(year < 2019 | epiweek <= 20) %>%
#   mutate(time = rep(unique(data_model$time), each = 2)) %>%
#   tidyr::pivot_wider(names_from = "arbo", values_from = "y") %>%
#   rename(y1 = Chik, y2 = Den) -> data_real


total.summy |>
  purrr::map(~.x |> left_join(data_real, by = "time") |>
        ggplot(aes(x = date, y = Median_y1, col = 'Nowcasting')) +
        geom_line(data = data_real %>% filter(time > 100),
                  aes(date, y = y1, col = 'Observed')) +
        geom_ribbon(aes(ymin = LI_y1, ymax = LS_y1, col = NA), alpha = 0.2, show.legend = F) +
        geom_line() +
        labs(title = "Chikungunya") +
        theme_bw() +
        theme(legend.position = "bottom", axis.text.x = element_text(angle = 90)) +
        scale_color_manual(values = c('grey50', 'black'), name = '')+
        scale_x_date(date_breaks = '2 weeks', date_labels = '%V/%y', name = 'Date in Weeks')+
        labs(x = '', y = 'Nº Cases')) -> grafs_chik #|>
  #patchwork::`+`

total.summy |>
  purrr::map(~.x |> left_join(data_real, by = "time") |>
               ggplot(aes(x = date, y = Median_y2, col = 'Nowcasting')) +
               geom_line(data = data_real %>% filter(time > 100),
                         aes(date, y = y2, col = 'Observed')) +
               geom_ribbon(aes(ymin = LI_y2, ymax = LS_y2, col = NA), alpha = 0.2, show.legend = F) +
               geom_line() +
               labs(title = "Dengue") +
               theme_bw() +
               theme(legend.position = "bottom", axis.text.x = element_text(angle = 90)) +
               scale_color_manual(values = c('grey50', 'black'), name = '')+
               scale_x_date(date_breaks = '2 weeks', date_labels = '%V/%y', name = 'Date in Weeks')+
               labs(x = '', y = 'Nº Cases')) -> grafs_dengue

gridExtra::grid.arrange(grobs = grafs_chik)
gridExtra::grid.arrange(grobs = grafs_dengue)

load("results\\model_runs.RData")
load("results\\predictive_samples.RData")


# summary(model_runs$M1$inla_output$waic)


info_dengue_data <- read.csv("info_dengue_rj.csv") |> 
  janitor::clean_names()

info_dengue_data |> 
  filter(se <= "201920") -> info_dengue_data


tweets <- log(info_dengue_data$tweet + 1)

#y_df <- data_model
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

M5 <- Y ~ -1 + alphai +
  + f(time,  model = "rw2", replicate = idrep,
      hyper = list(prec = list(param = c(0.001, 0.001)))) +
  + f(delay, model = "rw1", replicate = idrep) +
  + f(time1, model = "rw2") +
  + f(time2, copy  = "time1", fixed = TRUE) +
  + f(delay1, model = "rw1") +
  + f(delay2, copy  = "delay1", fixed = TRUE) +
  + tweets_ef_dengue +
  + tweets_ef_chik

M6 <- Y ~ -1 + alphai +
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

model_runs$M5 <- run_inla_model(data_model, model = M5, data_inla = joint_data)
model_runs$M6 <- run_inla_model(data_model, model = M6, data_inla = joint_data)


predictive_samples$M5 <- sample_predictive(model_runs$M5$inla_output, joint_data = joint_data)
predictive_samples$M6 <- sample_predictive(model_runs$M6$inla_output, joint_data = joint_data)


inic <- proc.time()
aaatempo <- run_inla_model(data_model, model = M5, data_inla = joint_data)
aaatempopred <- sample_predictive(aaatempo$inla_output, joint_data = joint_data)
tempo <- proc.time() - inic

### building the table 

model_runs_outputs <- lapply(model_runs, \(x) x$inla_output)
model_runs_outputs <- model_runs_outputs[c(6,5,1,3,2,4)]
total.summy_metrics <- total.summy[c(6,5,1,3,2,4)]
mapply(function(x,y)selection_criteria(x, data_real = data_real, y), 
       total.summy_metrics,
       model_runs_outputs) |> 
  t() |> 
  round(2) -> metrics

metrics[-c()] |> 
  xtable::xtable()

library(kableExtra)

colnames(metrics)[-c(3,4,6,7,8,10,12,13,15,16,17,19,21)] <- c("WAIC", "DIC", "MAPE", 
"Range of last 95% CI",
"MIS",
"MAPE",
"Range of last 95% CI","MIS")


 # |> 
 #  as.data.frame() |> 
 #  lapply(\(x){ifelse(x == min(x), cell_spec(round(x,2), bold = TRUE), round(x,2))}) |> 

kable(metrics[,-c(3,4,6,7,8,10,12,13,15,16,17,19,21)], "latex", booktabs = T) |> 
  kable_styling(latex_options = c("scale_down")) |> 
  add_header_above(c(" " = 2 , "Chikungunya correction" = 3, 
                     "Dengue correction" = 3)) %>% 
  column_spec(2, bold = )


View(predictive_samples$M3)


model_runs$M3$inla_output$summary.hyperpar |> 
  select(`0.5quant`) |> 
  round(2) |> 
  xtable::xtable()
  






