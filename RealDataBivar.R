
library(INLA)
library(ggplot2)
library(dplyr)
library(nowcaster)
library(lubridate)
#library(truncnorm)


# load("Dados\\den_chik_2017_2019_RJ.RData")
dados <- readRDS("Dados\\den_chik_2017_2019_RJ.RData")
dados <- dados %>% janitor::clean_names()

dados %>% 
  filter(dt_sin_pri > ymd("2017-01-01")) %>% 
  group_by(year(dt_sin_pri), epiweek(dt_sin_pri), arbo) %>% 
  summarise(y = n(), date = min(dt_sin_pri)) %>%  
  ggplot() +
  aes(x = date, y = y, colour = arbo) +
  geom_line() 


data <- dados |>
  filter(dt_sin_pri > ymd("2017-01-01")) |>
  select(dt_digita, dt_sin_pri, arbo) |>
  tidyr::drop_na(dt_digita) |>
  mutate(delay = as.numeric(floor(difftime(dt_digita, dt_sin_pri, units = "weeks")))) |>
  group_by(dt_sin_pri, delay, arbo)|>
  summarise(y = n()) |>
  filter(delay <= 15) 

data |> 
  ungroup() |>
  filter(arbo == "Den") |>
  mutate(epiweek = epiweek(dt_sin_pri), year = year(dt_sin_pri)) %>% 
  group_by(year, epiweek, delay) |>
  summarise(y = sum(y), date = min(dt_sin_pri)) |>
  group_by(year, epiweek) |>
  mutate(date = min(date)) |>
  group_by(date) |>
  summarise(y_real = cumsum(y), delay = delay) |>
  #filter(delay %in% c(0,7,15)) |>
ggplot() +
  aes(x = date, y = y_real, col = factor(delay)) +
  geom_line() +
  scale_x_date(breaks = "4 weeks") +
  theme(axis.text.x = element_text(angle = 45))


data |>
  ungroup() |>
  mutate(epiweek = epiweek(dt_sin_pri), year = year(dt_sin_pri)) |>
  group_by(epiweek, year, delay, arbo) |>
  summarise(y = sum(y)) |> 
  tidyr::pivot_wider(names_from = delay, values_from = y) |>
  mutate_if(~is.numeric(.),~coalesce(.,0)) |>
  tidyr::pivot_longer(cols = 4:19, names_to = "delay", values_to = "y") |>
  mutate(delay = as.numeric(delay)) |>
  arrange(year, epiweek, delay) |> 
  tidyr::pivot_wider(names_from = arbo, values_from = y) -> data #|>
  #mutate_at(.vars = vars(Chik:Den), ~coalesce(., 0)) -> data


data %>% 
  ungroup() %>% 
  filter(year < 2019 | epiweek <= 20) %>% 
  mutate(chik_obs = ifelse(epiweek + delay > 20 & year == 2019, NA, Chik),
         den_obs = ifelse(epiweek + delay > 20 & year == 2019, NA, Den)) -> data_model
  # mutate_at(.vars = vars(Chik:Den), 
  #         ~ifelse(epiweek + delay > 20 & year == 2019, NA, .x)) -> data_model
  
data_model %>% 
  distinct(epiweek, year) %>% 
  mutate(time = 1:n()) %>% 
  inner_join(data_model, by = c("epiweek", "year")) -> data_model


colnames(y_new)

y_new <- data_model %>% 
  mutate(report_time = time + delay) %>% 
  select(time, delay, report_time, y1_obs = chik_obs, y2_obs = den_obs, y1 = Chik, y2 = Den)
  

y_df <- y_new

#### Preparando dados para o modelo
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
  + f(time,  model = "rw2", replicate = idrep,
      hyper = list(prec = list(param = c(0.001, 0.001)))) +
  + f(delay, model = "rw1", replicate = idrep) +
  + f(time1, model = "rw2") +
  + f(time2, copy  = "time1", fixed = FALSE) + 
  + f(delay1, model = "rw1") +
  + f(delay2, copy  = "delay1", fixed = FALSE) 


output0 <- INLA::inla(joint_model, family = rep("nbinomial", 2),
                      data = joint_data,
                      control.predictor = list(link = 1, compute = T),
                      control.compute = list( config = T, waic=F, dic=F),
                      control.family = list(
                        list(hyper = list("theta" = list(prior = "loggamma", param = c(0.001, 0.001)))),
                        list(hyper = list("theta" = list(prior = "loggamma", param = c(0.001, 0.0001))))
                      )
                      
)


vectime <- unique(y_new$time)
ts.plot(output0$summary.random$time$`0.5quant`[vectime])
ts.plot(output0$summary.random$time$`0.5quant`[vectime + length(vectime)])
ts.plot(output0$summary.random$time1$`0.5quant`)
ts.plot(output0$summary.random$time2$`0.5quant`)

ts.plot(output0$summary.random$delay$`0.5quant`[1:16])
ts.plot(output0$summary.random$delay$`0.5quant`[17:32])
ts.plot(output0$summary.random$delay1$`0.5quant`)
ts.plot(output0$summary.random$delay2$`0.5quant`)

ts.plot(output0$summary.random$delay1$`0.5quant` + output0$summary.random$delay$`0.5quant`[1:16])
# lines(gamma_d[,1] + psi, col = 2)
ts.plot(output0$summary.random$delay2$`0.5quant` + output0$summary.random$delay$`0.5quant`[17:32])
# lines(gamma_d[,2] + psi, col = 2)


ts.plot(output0$summary.random$time1$`0.5quant` + output0$summary.random$time$`0.5quant`[1:124], ylim = c(-6, 9))
lines(beta[,1] + delta, col = 2)

ts.plot(output0$summary.random$time2$`0.5quant` + output0$summary.random$time$`0.5quant`[125:248])
lines(beta[,2] + delta, col = 2)

ts.plot(c(LP1, LP2))
lines(output0$summary.linear.predictor$`0.5quant`, col = 2)



samples <- INLA::inla.posterior.sample(n = 1000, output0)



index_missing <- which(is.na(joint_data$`1`) & is.na(joint_data$`2`))


vector_samples_bvar <- lapply(X = samples,
                              FUN = function(x, idx = index_missing, n = 1984){
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

## Step 4: Applying the age aggregation on each posterior
tibble.samples.0 <- lapply(X = vector_samples_bvar,
                           FUN = gg.age,
                           dados = joint_data,
                           idx = index_missing)


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



dados_by_week <- y_new |>
  group_by(time) |> 
  summarise(y1 = sum(y1, na.rm = T),
            y2 = sum(y2, na.rm = T))


  



# dados %>% 
#   filter(dt_sin_pri > ymd("2017-01-01")) %>% 
#   group_by(year(dt_sin_pri), epiweek(dt_sin_pri), arbo) %>% 
#   summarise(y = n(), date = min(dt_sin_pri)) %>% 
#   tidyr::pivot_wider(names_from = "arbo", values_from = "y") -> aaa
# 
# dados_by_week$y1 <- aaa$Den[1:124]
# dados_by_week$y2 <- aaa$Chik[1:124]




dados_by_week |> 
  ggplot()+
  geom_line(data = dados_by_week, aes(time, y = y1, col = 'Observed'))+
  theme_bw()+
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90)) +
  scale_color_manual(values = c('grey50', 'black'), name = '')+
  labs(x = '', y = 'Nº Cases')


total.summy |> 
  ggplot(aes(x = time, y = Median_y1, col = 'Nowcasting')) +
  geom_line(data = dados_by_week %>% filter(time > 100), aes(time, y = y1, col = 'Observed')) +
  geom_ribbon(aes(ymin = LI_y1, ymax = LS_y1, col = NA), alpha = 0.2, show.legend = F) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90)) +
  # scale_color_manual(values = c('grey50', 'black'), name = '')+
  # scale_x_date(date_breaks = '2 weeks', date_labels = '%V/%y', name = 'Date in Weeks')+
  labs(x = '', y = 'Nº Cases')


total.summy |> 
  ggplot(aes(x = time, y = Median_y2, col = 'Nowcasting')) +
  geom_line(data = dados_by_week %>% filter(time > 100), aes(time, y = y2, col = 'Observed')) +
  geom_ribbon(aes(ymin = LI_y2, ymax = LS_y2, col = NA), alpha = 0.2, show.legend = F) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90)) +
  # scale_color_manual(values = c('grey50', 'black'), name = '')+
  #scale_x_date(date_breaks = '2 weeks', date_labels = '%V/%y', name = 'Date in Weeks')+
  labs(x = '', y = 'Nº Cases')


#### Abordagem univariada







