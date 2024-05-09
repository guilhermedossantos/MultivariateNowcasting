library(ggplot2)
library(dplyr)


## Funções
source("functions_multivar_nowcasting.R")



## Dados
dados <- readRDS("Dados\\den_chik_2017_2019_RJ.RData")
dados <- dados %>% janitor::clean_names()


### Gráfico delay dengue
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


data <- dados |>
  filter(dt_sin_pri > ymd("2017-01-01")) |>
  select(dt_digita, dt_sin_pri, arbo) |>
  tidyr::drop_na(dt_digita) |>
  mutate(delay = as.numeric(floor(difftime(dt_digita, dt_sin_pri, units = "weeks")))) |>
  group_by(dt_sin_pri, delay, arbo)|>
  summarise(y = n()) #|>
#filter(delay <= 15) 

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
  filter(delay %in% c(0, 01, 05, 10, 15, max(delay))) |>
  mutate(delay = ifelse(delay == max(delay), "Total", delay)) |> 
  ggplot() +
  aes(x = date, y = y_real, col = factor(delay)) +
  geom_line() +
  scale_x_date(breaks = "8 weeks") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 


### Proporção

data <- dados  |>
  filter(dt_sin_pri > ymd("2017-01-01")) |> 
  select(dt_digita, dt_sin_pri, arbo) |> 
  tidyr::drop_na(dt_digita) |>
  dplyr::mutate(delay = as.numeric(floor(difftime(dt_digita, dt_sin_pri, units = "weeks")))) |>
  dplyr::group_by(dt_sin_pri, delay, arbo)|>
  dplyr::summarise(y = n())


data |> 
  group_by(dt_sin_pri, arbo) |> 
  summarise(y = sum(y)) -> y_time

data |> 
  left_join(y_time, 
            by = c("dt_sin_pri", "arbo"),
            suffix = c("obs", "real")) |> 
  mutate(prop_obs = yobs/yreal) -> data_prop_obs


data_prop_obs |> 
  arrange(dt_sin_pri, delay) |> 
  group_by(dt_sin_pri, arbo) |> 
  mutate(prop_obs_cum = cumsum(prop_obs)) -> data_eda

data_eda |> 
  filter(arbo == "Den") |> 
  pull(yobs) |> 
  sum() -> reported_dengue

data_eda |> 
  filter(arbo == "Chik") |> 
  pull(yobs) |> 
  sum() -> reported_chik

data_eda |> 
  filter(arbo == "Den", delay <= 15) |> 
  pull(yobs) |> 
  sum() -> delay_less_15_dengue


delay_less_15_dengue/reported_dengue

data_eda |> 
  filter(arbo == "Chik", delay <= 15) |> 
  pull(yobs) |> 
  sum() -> delay_less_15_chik

delay_less_15_chik/reported_chik

# uns 95% reportados em menos de 15 semanas 
# 89 e 86% em 10 semanas



sapply(0:25, function(D){
  data_eda |> 
    filter(arbo == "Den", delay <= D) |> 
    pull(yobs) |> 
    sum() -> delay_less_15_dengue
  
  
  delay_less_15_dengue/reported_dengue
  
  data_eda |> 
    filter(arbo == "Chik", delay <= D) |> 
    pull(yobs) |> 
    sum() -> delay_less_15_chik
  
  delay_less_15_chik/reported_chik
  
  return(c("D" = D, 
           "Dengue" = delay_less_15_dengue/reported_dengue, 
           "Chikungunya" = delay_less_15_chik/reported_chik))
}) |> t() -> delays_prop
delays_prop
xtable::xtable(delays_prop)

delays_prop |> 
  as_tibble() |> 
  mutate(Dengue = c(Dengue[1],diff(Dengue)),
         Chikungunya = c(Chikungunya[1],diff(Chikungunya))) |> 
  tidyr::pivot_longer(cols = -D,names_to = "Disease", values_to = "Proportion") |>
  ggplot() +
  geom_col(aes(x = D, y = Proportion, fill = Disease), position = "dodge") +
  geom_hline(yintercept = 0.01, linetype = "dashed") +
  labs(x = "Number of weeks of delay", y = "Percentage of reported cases") + 
  scale_x_continuous(breaks = seq(0,20,5)) + 
  scale_y_continuous(labels = scales::percent) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

ggsave("perc_reported_delay.pdf")


### DELAY MÉDIO E MEDIANO
data |> 
  group_by(arbo) |> 
  summarise(mean_delay = weighted.mean(delay, y),
            median_delay = median(rep(delay, times = y)))

# data_eda |> 
#   filter(arbo == "Chik") |> 
#   ggplot() +
#   aes(x = dt_sin_pri, y = yobs, fill = factor(delay)) +
#   geom_col(position = "fill") #delay mudando no tempo?

# 
# data |> 
#   filter(delay <= 50) |> 
ggplot(data = NULL) +
  geom_density(aes(x = rep(data$delay, times = data$y), 
                   fill = rep(data$arbo, times = data$y)), 
               alpha = 0.6, bw = 2) +
  scale_x_continuous(limits = c(0,25)) +
  labs(x = "Delay",
       y = "Density",
       fill = "Disease") +
  theme_bw(base_size = 13)




