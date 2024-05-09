## Análise sliding windows 
library(stringr)
library(INLA)
library(dplyr)
library(ggplot2)

## Funções
source("functions_multivar_nowcasting.R")


dados <- readRDS("Dados\\den_chik_2017_2019_RJ.RData")
dados <- dados %>% janitor::clean_names()

dados |>
  filter(dt_sin_pri > ymd("2017-01-01")) |>
  select(dt_digita, dt_sin_pri, arbo) |> 
  create_data(D = 15) |> 
  model_data(max_year = 2019, max_epiweek = 52) -> data_model

## Dados reais para comparação final
dados |>
  filter(dt_sin_pri > ymd("2017-01-01")) |>
  select(dt_digita, dt_sin_pri, arbo) |> 
  create_data(D = 15, return_data_real = T) -> data_real

maxtime <- 116:140 #seq(60, max(data_model$time), by = 4)

data_real_list <- purrr::map(maxtime, ~{
  # maxdate <- ymd("2017-01-01") + weeks(.x)
  data_real |> 
    ungroup() |> 
    mutate(time = rep(1:(nrow(data_real)/2), each = 2)) |> 
    filter(time <= .x & time > .x - 70) %>% 
    tidyr::pivot_wider(names_from = "arbo", values_from = "y") %>% 
    rename(y1 = Chik, y2 = Den)
})

data_model_list <- purrr::map(maxtime, ~{
  maxdate <- ymd("2017-01-01") + weeks(.x)
  dados |>
    filter(dt_sin_pri > ymd("2017-01-01")) |>
    select(dt_digita, dt_sin_pri, arbo) |> 
    create_data(D = 15) |> 
    model_data(max_year = year(maxdate), 
               max_epiweek = max(epiweek(maxdate), week(maxdate))) |> 
    filter(time <= .x & time > .x - 70)
})



## Delay 15
tempo_files <- list.files(pattern = "tempos")
tempo_delay15 <- tempo_files[(str_detect(tempo_files, "delay_15.*25runs"))]
# tempo_delay15 <- tempo_delay15[c(2:6,1)]

tempos_bivariate_fit <- list() -> tempos_bivariate_samples

for(i in 1:6){
  load(tempo_delay15[i])  
  tempos_bivariate_fit[[i]] <- tempo_bivariate_fit_30
  tempos_bivariate_samples[[i]] <- tempo_bivariate_samples_30
}

tempos_bivariate_fit <- tempos_bivariate_fit |> 
  sapply(\(x) round(x[3]/60,2)) 

tempos_bivariate_samples <- tempos_bivariate_samples |> 
  sapply(\(x) round(x[3]/60,2)) 

tempos_biv_15 <- cbind(tempos_bivariate_fit, tempos_bivariate_samples)

## Delay 10 
tempo_files <- list.files(pattern = "tempos")
tempo_delay10 <- tempo_files[(str_detect(tempo_files, "delay_10.*25runs"))]
tempo_delay10 <- tempo_delay10[c(2:6,1)]

tempos_bivariate_fit <- list() -> tempos_bivariate_samples

for(i in 1:6){
  load(tempo_delay10[i])  
  tempos_bivariate_fit[[i]] <- tempo_bivariate_fit_30
  tempos_bivariate_samples[[i]] <- tempo_bivariate_samples_30
}

tempos_bivariate_fit <- tempos_bivariate_fit |> 
  sapply(\(x) round(x[3]/60,2)) 

tempos_bivariate_samples <- tempos_bivariate_samples |> 
  sapply(\(x) round(x[3]/60,2)) 

(tempos_biv_10 <- cbind(tempos_bivariate_fit, tempos_bivariate_samples))

## Delay 20 
tempo_files <- list.files(pattern = "tempos")
tempo_delay20 <- tempo_files[(str_detect(tempo_files, "delay_20.*25runs"))]
tempo_delay20 <- tempo_delay20[c(2:6,1)]

tempos_bivariate_fit <- list() -> tempos_bivariate_samples

for(i in 1:6){
  load(tempo_delay20[i])  
  tempos_bivariate_fit[[i]] <- tempo_bivariate_fit_30
  tempos_bivariate_samples[[i]] <- tempo_bivariate_samples_30
}

tempos_bivariate_fit <- tempos_bivariate_fit |> 
  sapply(\(x) round(x[3]/60,2)) 

tempos_bivariate_samples <- tempos_bivariate_samples |> 
  sapply(\(x) round(x[3]/60,2)) 

(tempos_biv_20 <- cbind(tempos_bivariate_fit, tempos_bivariate_samples))


### Delay15
# fit_files <- list.files(pattern = "(bivariate_sliding.*_final.RData)")
# 
# fit_delay15 <- fit_files[(!str_detect(fit_files,"delay")| str_detect(fit_files,"delay_15"))]
# 
# fits_delay15 <- list()
# 
# for(i in 1:length(fit_delay15)){
#   load(fit_delay15[i])
#   fits_delay15[[i]] <- lapply(bivariate_fit, \(x) x$inla_output)
# }



# sample_files <- list.files(pattern = "(bivariate_samples.*_final.RData)")
# sample_delay15 <- sample_files[(!str_detect(sample_files,"delay")| str_detect(sample_files,"delay_15"))]

sample_files <- list.files(pattern = "(bivariate_samples.*_25runs.RData)")
sample_delay15 <- sample_files[str_detect(sample_files,"delay_15")]
sample_delay15 <- c(sample_delay15[2:6], sample_delay15[1])




samples_list <- list()
summaries <- list()


for(i in 1:length(sample_delay15)){
  load(sample_delay15[i])
  samples_list[[i]] <- bivariate_samples
  summaries[[i]] <- suppressMessages(lapply(bivariate_samples, 
                           bivar_nowcasting_summary))
}

summaries_test <- list()
summaries_vs_es <- list()
coverage_mean <- list() -> coverage_mean2 -> coverage_mean3
mpi <- list() -> mis
for(i in 1:6){
  summaries_test[[i]] <- purrr::map2(summaries[[i]], 
                                     data_real_list,
                                     selection_criteria_summary)
  
  summaries_vs_es[[i]] <- purrr::map2(samples_list[[i]], 
                                      data_real_list,
                                      scores_posterior)
  
  summaries_test[[i]] <- purrr::map2(summaries_test[[i]], 
                                     summaries_vs_es[[i]],
                                     cbind)
  
  
  
  coverage_mean[[i]] <- lapply(summaries_test[[i]], \(x) mean(x$coverage_95_y1 & x$coverage_95_y2)) |> 
    unlist() |> 
    mean()
  
  coverage_mean2[[i]] <- lapply(summaries_test[[i]], \(x) x$coverage_95_y1 | x$coverage_95_y2) |> 
    unlist() |> 
    mean()
  
  coverage_mean3[[i]] <- lapply(summaries_test[[i]], \(x) x$coverage_95_y1) |> 
    unlist() |> 
    mean()
  
  mpi[[i]] <- lapply(summaries_test[[i]], \(x) (x$relative_width_95)) |>
    unlist() |>
  mean()
  
  mis[[i]] <- lapply(summaries_test[[i]], \(x) (x$relative_is_95)) |>
    unlist() |>
    mean()
  
  
}

window_sizes = c(3:7*10, "Full data")
summaries_metrics_df <- list()
for(i in 1:6){
  summaries_metrics_df[[i]] <- lapply(summaries_test[[i]], function(summy_df){
    summy_df %>%
      summarise_at(vars(perc_error_y1:variogram_score), .funs = list(mean = mean,
                                                                    last = last))
  }) %>% 
    
    do.call(rbind, .) %>%
    mutate(window_size = window_sizes[i],
           window_index = 1:nrow(.)) %>%
    select(window_size, everything())
  
  
}

summaries_metrics_df <- do.call(rbind, summaries_metrics_df)

ggplot(summaries_metrics_df) +
  geom_line(aes(x= window_index, y = perc_error_y1_mean, col = window_size))
ggplot(summaries_metrics_df) +
  geom_line(aes(x= window_index, y = perc_error_y2_mean, col = window_size))

ggplot(summaries_metrics_df) +
  geom_line(aes(x= window_index, y = perc_error_y1_mean + perc_error_y2_mean, col = window_size))


ggplot(summaries_metrics_df) +
  geom_line(aes(x= window_index, y = relative_width_95_mean, col = window_size))

ggplot(summaries_metrics_df) +
  geom_line(aes(x= window_index, y = relative_is_95_mean, col = window_size))

ggplot(summaries_metrics_df) +
  geom_line(aes(x= window_index, y = relative_is_95_last, col = window_size))

ggplot(summaries_metrics_df) +
  geom_line(aes(x= window_index, y = relative_width_95_last, col = window_size))

ggplot(summaries_metrics_df) +
  geom_line(aes(x= window_index, y = coverage_95_y1_mean , col = window_size))

ggplot(summaries_metrics_df) +
  geom_line(aes(x= window_index, y = coverage_95_y2_mean , col = window_size))


ggplot(summaries_metrics_df) +
  geom_line(aes(x= window_index, y = coverage_95_both_mean, col = window_size))

ggplot(summaries_metrics_df) +
  geom_line(aes(x= window_index, y = coverage_95_atleastone_mean, col = window_size))

ggplot(summaries_metrics_df) +
  geom_line(aes(x= window_index, y = variogram_score_mean, col = window_size))

ggplot(summaries_metrics_df) +
  geom_line(aes(x= window_index, y = energy_score_mean, col = window_size))

ggplot(summaries_metrics_df) +
  geom_line(aes(x= window_index, y = variogram_score_last, col = window_size))





summaries_metrics_df |> 
  select(window_size, MAPE_y1 = perc_error_y1_mean, MAPE_y2 = perc_error_y2_mean,
         coverage_95_y1_mean, coverage_95_y2_mean, coverage_95_both_mean,
         coverage_95_both_last,
         width_95_y1_mean, width_95_y2_mean,
         width_95_y1_last, width_95_y2_last,
         relative_width_95_mean, relative_is_95_mean,
         relative_width_95_last, relative_is_95_last,
         energy_score_mean, energy_score_last,
         variogram_score_mean, variogram_score_last
         ) |> 
  group_by(window_size) |> 
  summarise_all(mean) -> metrics_by_window_size


library(kableExtra)
metrics_by_window_size |> 
  mutate(fit_time = tempos_biv_15[,1],
         sample_time = tempos_biv_15[,2],
         MAPE_y1 = 100*MAPE_y1,
         MAPE_y2 = 100*MAPE_y2) |>
  #kableExtra::add_header_above("Window size") |> 
  select("Fit time" = fit_time,
         "Sampling time" = sample_time,
         "MAPE (Chikungunya)" = MAPE_y1,
         "MAPE (Dengue)" = MAPE_y2, 
         "Coverage (Chikungunya)" = coverage_95_y1_mean,
         "Coverage (Dengue)" = coverage_95_y2_mean,
         "Coverage (both)" = coverage_95_both_mean,
         "Last coverage (both)" = coverage_95_both_last,
         "mean relWidth" = relative_width_95_mean,
         "mean relIS" = relative_is_95_mean,
         "mean last relWidth" = relative_width_95_last,
         "mean last relIS" = relative_is_95_last,
         "mean Energy Score" = energy_score_mean,
         "mean Variogram Score" = variogram_score_mean
         ) |> 
  t() |> 
  kableExtra::kbl(format = "latex", booktabs = T, digits = 3) |> 
  add_header_above(c(" " = 1,
                     "30" = 1,
                     "40" = 1,
                     "50" = 1,
                     "60" = 1,
                     "70" = 1,
                     "Full data" = 1)) |>
  add_header_above(c(" " = 1, "Window Size" = 6))

  knitr::kable()  |> 
  
  xtable::xtable()
  


# selection_criteria_summary(summaries[[5]][[1]], 
#                    data_real = data_real_list[[1]]) |> 
#   pull(relative_width_95) |> 
#   ts.plot()



summary_list_60 <- summaries[[5]]

Sys.setlocale("LC_ALL", "English")
g_list_y1 <- list() -> g_list_y2
for (i in 1:length(summary_list_60)){
  data_no_delay <- data_model_list[[i]] |> 
    mutate(date = ymd("2017-01-01") + weeks(time - 1),
           y1_obs = as.numeric(y1_obs)) |> 
    group_by(date, time) |>
    summarise_at(vars(y1_obs:y2), sum, na.rm = T)
  
  data_real |> 
    filter(arbo == "Chik") |> 
    ggplot() +
    geom_line(aes(x = date, y = y, 
                  colour = "Real observed",
                  linetype = "Real observed"), size = 0.6) +
    geom_line(mapping = aes(x = date, y = y1_obs, 
                            colour = "Observed at the moment",
                            linetype = "Observed at the moment"), data = data_no_delay,
              size = 0.6) +
    geom_line(mapping = aes(x = date, y = Median_y1, 
                            colour = "Estimated", 
                            linetype = "Estimated"),
              size = 0.6, data = summary_list_60[[i]] |> inner_join(data_no_delay, by = "time")) +
    geom_ribbon(aes(ymin = LI_y1, ymax = LS_y1, x = date), fill = "black", alpha = 0.3,
                data = summary_list_60[[i]] |> inner_join(data_no_delay, by = "time")) +
    annotate("text", x = ymd("2019-01-25"), y = Inf, hjust = 0.5, vjust = 1.2,
             label = bquote(W[.(i)]), size = 5) +
    scale_x_date(limits = c(ymd("2019-01-01"), ymd("2019-10-01"))) +
    # labs(y = ifelse(i %% 5 == 1,"Chikungunya cases",""), x= ifelse(i > 20, "Date","")) +
    labs(y = "", x= "") +
    scale_linetype_manual(values = c("dashed", "dotted", "solid")) +
    scale_color_manual(values = c("red", "black", "black")) +
    theme_bw() +
    guides(color = guide_legend(override.aes = list(linetype = c("dashed", "dotted", "solid")), 
                                title = "Legend"),
           linetype = FALSE)  -> g_list_y1[[i]]
    
    data_real |> 
      filter(arbo == "Den") |> 
      ggplot() +
      geom_line(aes(x = date, y = y, 
                    colour = "Real observed",
                    linetype = "Real observed"), size = 0.9) +
      geom_line(mapping = aes(x = date, y = y2_obs, 
                              colour = "Observed at the moment",
                              linetype = "Observed at the moment"), data = data_no_delay,
                size = 0.9) +
      geom_line(mapping = aes(x = date, y = Median_y2, 
                              colour = "Estimated", 
                              linetype = "Estimated"),
                size = 0.9, data = summary_list_60[[i]] |> inner_join(data_no_delay, by = "time")) +
      geom_ribbon(aes(ymin = LI_y2, ymax = LS_y2, x = date), fill = "black", alpha = 0.3,
                  data = summary_list_60[[i]] |> inner_join(data_no_delay, by = "time")) +
      scale_x_date(limits = c(ymd("2019-01-01"), ymd("2019-10-01"))) +
      annotate("text", x = ymd("2019-01-25"), y = Inf, hjust = 0.5, vjust = 1.2, label = bquote(W[.(i)]), size = 5) +
      #labs(y = ifelse(i %% 5 == 1,"Dengue cases",""), x= ifelse(i > 20, "Date","")) +
      labs(y = "", x= "") +
      scale_linetype_manual(values = c("dashed", "dotted", "solid")) +
      scale_color_manual(values = c("red", "black", "black")) +
      theme_bw() +
      guides(color = guide_legend(override.aes = list(linetype = c("dashed", "dotted", "solid")), 
                                  title = "Legend"),
             linetype = FALSE) -> g_list_y2[[i]]
  
}
# 
# leg1 <- cowplot::get_legend(g_list_y1[[1]])
# 
# g_list_y1 <- lapply(g_list_y1, \(x) x + theme(legend.position = "none"))

library(ggpubr)
library(grid)
ggpubr::ggarrange(plotlist = g_list_y1,
                  common.legend = T, legend= "top") -> figure

annotate_figure(figure, 
                left = textGrob("Chikungunya cases", rot = 90, vjust = 1, gp = gpar(cex = 1.1)),
                bottom = textGrob("Date", gp = gpar(cex = 1.1)))
# ggsave("sliding_windows_corrections_chik_v2.pdf", height = 10, width = 12)


ggpubr::ggarrange(plotlist = g_list_y2,
                  common.legend = T, legend= "top") -> figure

annotate_figure(figure, 
                left = textGrob("Dengue cases", rot = 90, vjust = 1, gp = gpar(cex = 1.1)),
                bottom = textGrob("Date", gp = gpar(cex = 1.1)))
ggsave("sliding_windows_corrections_dengue_v2.pdf", height = 10, width = 12)


gridExtra::grid.arrange(grobs = g_list_y1, ncol = 5, bottom = leg1)
gridExtra::grid.arrange(grobs = g_list_y2, ncol = 5)

data_no_delay <- data_model_list[[1]] |> 
  mutate(date = ymd("2017-01-01") + weeks(time - 1),
         y1_obs = as.numeric(y1_obs)) |> 
  group_by(date, time) |>
  summarise_at(vars(y1_obs:y2), sum, na.rm = T)

data_real |> 
  filter(arbo == "Chik") |> 
ggplot() +
  geom_line(aes(x = date, y = y)) +
  geom_line(mapping = aes(x = date, y = y1_obs), data = data_no_delay, linetype = "dashed") +
  geom_line(mapping = aes(x = date, y = Median_y1), colour = "red", linetype = "dashed",
            data = summary_list_60[[1]] |> inner_join(data_no_delay, by = "time")) +
  geom_ribbon(aes(ymin = LI_y1, ymax = LS_y1, x = date), fill = "black", alpha = 0.1,
              data = summary_list_60[[1]] |> inner_join(data_no_delay, by = "time")) +
  scale_x_continuous(limits = c(ymd("2019-01-01"), ymd("2019-12-31"))) + 
  theme_bw()





#### Analise entre modelos####

sample_files <- list.files(pattern = "(samples.*_25runs_)|(samples.*delay_15_window_70.*_25runs)")
sample_delay15 <- sample_files[str_detect(sample_files,"delay_15")]
sample_delay15 <- c(sample_delay15[2:4], sample_delay15[1])


samples_list <- list()
summaries <- list()


for(i in 1:length(sample_delay15)){
  load(sample_delay15[i])
  if(str_detect(sample_delay15[i], "25runs_")){
    
    samples_list[[i]] <- purrr::map2(bivariate_samples,
                                     1:25,
                                     ~{
                                       .x |> 
                                         mutate(time = time + 70 - 25 + .y)
                                     })
  }else{
    samples_list[[i]] <- bivariate_samples
  }
  summaries[[i]] <- suppressMessages(lapply(samples_list[[i]], 
                                            bivar_nowcasting_summary))
}

summaries_test <- list()
summaries_vs_es <- list()
coverage_mean <- list() -> coverage_mean2 -> coverage_mean3
mpi <- list() -> mis
#summaries <- summaries1
for(i in seq_along(summaries)){
  
  summaries_test[[i]] <- purrr::map2(summaries[[i]], 
                                     data_real_list,
                                     selection_criteria_summary)
  
  summaries_vs_es[[i]] <- purrr::map2(samples_list[[i]], 
                                      data_real_list,
                                      scores_posterior)
  
  summaries_test[[i]] <- purrr::map2(summaries_test[[i]], 
                                     summaries_vs_es[[i]],
                                     cbind)
  
  
  
  coverage_mean[[i]] <- lapply(summaries_test[[i]], \(x) mean(x$coverage_95_y1 & x$coverage_95_y2)) |> 
    unlist() |> 
    mean()
  
  coverage_mean2[[i]] <- lapply(summaries_test[[i]], \(x) x$coverage_95_y1 | x$coverage_95_y2) |> 
    unlist() |> 
    mean()
  
  coverage_mean3[[i]] <- lapply(summaries_test[[i]], \(x) x$coverage_95_y1) |> 
    unlist() |> 
    mean()
  
  mpi[[i]] <- lapply(summaries_test[[i]], \(x) (x$relative_width_95)) |>
    unlist() |>
    mean()
  
  mis[[i]] <- lapply(summaries_test[[i]], \(x) (x$relative_is_95)) |>
    unlist() |>
    mean()
  
  
}

models = c("M0", "M1", "M2", "M3")
summaries_metrics_df <- list()
for(i in seq_along(summaries)){
  summaries_metrics_df[[i]] <- lapply(summaries_test[[i]], function(summy_df){
    summy_df %>%
      summarise_at(vars(perc_error_y1:variogram_score), .funs = list(mean = mean,
                                                                     last = last))
  }) %>% 
    
    do.call(rbind, .) %>%
    mutate(model = models[i],
           model_index = 1:nrow(.)) %>%
    select(model, everything())
  
  
}

summaries_metrics_df <- do.call(rbind, summaries_metrics_df)

ggplot(summaries_metrics_df) +
  geom_line(aes(x= model_index, y = perc_error_y1_mean, col = model), size = 1) +
  scale_y_continuous(labels = scales::percent) +
  labs(y = "MAPE for Chikungunya",
       x = "Window",
       colour = "Model") + 
  theme_bw(base_size = 13) -> mape_chik
  

ggplot(summaries_metrics_df) +
  geom_line(aes(x= model_index, y = perc_error_y2_mean, col = model), size = 1) +
  scale_y_continuous(labels = scales::percent) +
  labs(y = "MAPE for Dengue",
       x = "Window",
       colour = "Model") + 
  theme_bw(base_size = 13) -> mape_dengue
  

ggplot(summaries_metrics_df) +
  geom_line(aes(x= model_index, 
                y = perc_error_y1_mean + perc_error_y2_mean, col = model), size = 1) +
  scale_y_continuous(labels = scales::percent) +
  labs(y = "MAPE sum",
       x = "Window",
       colour = "Model") + 
  theme_bw(base_size = 13)


ggplot(summaries_metrics_df) +
  geom_line(aes(x= model_index, y = relative_width_95_mean, col = model), size = 1) +
  labs(y = "Mean relative interval width",
       x = "Window",
       colour = "Model") + 
  theme_bw(base_size = 13) -> g_rel_width

ggplot(summaries_metrics_df) +
  geom_line(aes(x= model_index, y = relative_is_95_mean, col = model), size = 1) +
  labs(y = "Mean relative interval score",
       x = "Window",
       colour = "Model") + 
  theme_bw(base_size = 13) -> g_rel_IS

ggplot(summaries_metrics_df) +
  geom_line(aes(x= model_index, y = relative_width_95_last, col = model), size = 1) +
  labs(y = "Last relative interval width",
       x = "Window",
       colour = "Model") + 
  theme_bw(base_size = 13) -> g_rel_width_last


ggarrange(g_rel_width, g_rel_IS,
          mape_chik, mape_dengue,
          common.legend = T, legend = "bottom")


ggplot(summaries_metrics_df) +
  geom_line(aes(x= model_index, y = relative_is_95_last, col = model), size = 1) + 
  labs(y = "Last relative interval score",
       x = "Window",
       colour = "Model") + 
  theme_bw(base_size = 13)

ggplot(summaries_metrics_df) +
  geom_line(aes(x= model_index, y = relative_width_95_last, col = model), size = 1) +
  labs(y = "Last relative interval width",
       x = "Window",
       colour = "Model") + 
  theme_bw(base_size = 13)

ggplot(summaries_metrics_df) +
  geom_line(aes(x= model_index, y = coverage_95_y1_mean , col = model))

ggplot(summaries_metrics_df) +
  geom_line(aes(x= model_index, y = coverage_95_y2_mean , col = model))


ggplot(summaries_metrics_df) +
  geom_line(aes(x= model_index, y = coverage_95_both_mean, col = model))

ggplot(summaries_metrics_df) +
  geom_line(aes(x= model_index, y = coverage_95_atleastone_mean, col = model))

ggplot(summaries_metrics_df) +
  geom_line(aes(x= model_index, y = variogram_score_mean, col = model))

ggplot(summaries_metrics_df) +
  geom_line(aes(x= model_index, y = energy_score_mean, col = model))

ggplot(summaries_metrics_df) +
  geom_line(aes(x= model_index, y = variogram_score_last, col = model))

ggplot(summaries_metrics_df) +
  geom_line(aes(x= model_index, y = squared_error_y1_mean, col = model))  +
  labs(y = "MSE for Chikungunya",
       x = "Window",
       colour = "Model") + 
  theme_bw(base_size = 13)

  
ggplot(summaries_metrics_df) +
  geom_line(aes(x= model_index, y = squared_error_y2_mean, col = model), size = 1)  +
  labs(y = "MSE for Dengue",
       x = "Window",
       colour = "Model") + 
  theme_bw(base_size = 13)



summaries_metrics_df |> 
  select(model, MAPE_y1 = perc_error_y1_mean, MAPE_y2 = perc_error_y2_mean,
         coverage_95_y1_mean, coverage_95_y2_mean, coverage_95_both_mean,
         coverage_95_both_last,
         width_95_y1_mean, width_95_y2_mean,
         width_95_y1_last, width_95_y2_last,
         relative_width_95_mean, relative_is_95_mean,
         relative_width_95_last, relative_is_95_last,
         energy_score_mean, energy_score_last,
         variogram_score_mean, variogram_score_last
  ) |> 
  group_by(model) |> 
  summarise_all(mean) -> metrics_by_model


library(kableExtra)
metrics_by_model |> 
  mutate(MAPE_y1 = 100*MAPE_y1,
         MAPE_y2 = 100*MAPE_y2) |>
  #kableExtra::add_header_above("Window size") |> 
  select("MAPE (Chikungunya)" = MAPE_y1,
         "MAPE (Dengue)" = MAPE_y2, 
         "Coverage (Chikungunya)" = coverage_95_y1_mean,
         "Coverage (Dengue)" = coverage_95_y2_mean,
         "Coverage (both)" = coverage_95_both_mean,
         "Last coverage (both)" = coverage_95_both_last,
         "mean relWidth" = relative_width_95_mean,
         "mean relIS" = relative_is_95_mean,
         "mean last relWidth" = relative_width_95_last,
         "mean last relIS" = relative_is_95_last,
         "mean Energy Score" = energy_score_mean,
         "mean Variogram Score" = variogram_score_mean
  ) |> 
  t() |> 
  kableExtra::kbl(format = "latex", booktabs = T, digits = 3) |> 
  add_header_above(c(" " = 1,
                     "M0" = 1,
                     "M1" = 1,
                     "M2" = 1,
                     "M3" = 1)) |>
  add_header_above(c(" " = 1, "Model" = 4))

knitr::kable()  |> 
  
  xtable::xtable()




summary_list_60 <- summaries[[3]]


g_list_y1 <- list() -> g_list_y2
for (i in 1:length(summary_list_60)){
  data_no_delay <- data_model_list[[i]] |> 
    mutate(date = ymd("2017-01-01") + weeks(time - 1),
           y1_obs = as.numeric(y1_obs)) |> 
    group_by(date, time) |>
    summarise_at(vars(y1_obs:y2), sum, na.rm = T)
  
  data_real |> 
    filter(arbo == "Chik",
           date <= ymd("2019-10-01")) |>
    ggplot() +
    geom_line(aes(x = date, y = y, 
                  colour = "Real observed",
                  linetype = "Real observed"), size = 1) +
    geom_line(mapping = aes(x = date, y = y1_obs, 
                            colour = "Observed at the moment",
                            linetype = "Observed at the moment"), data = data_no_delay,
              size = 1) +
    geom_line(mapping = aes(x = date, y = Median_y1, 
                            colour = "Estimated", 
                            linetype = "Estimated"),
              size = 1, data = summary_list_60[[i]] |> inner_join(data_no_delay, by = "time")) +
    geom_ribbon(aes(ymin = LI_y1, ymax = LS_y1, x = date), fill = "black", alpha = 0.25,
                data = summary_list_60[[i]] |> inner_join(data_no_delay, by = "time")) +
    scale_x_date(limits = c(ymd("2019-01-01"), ymd("2019-10-01"))) +
    labs(y = ifelse(i %% 5 == 1,"Chikungunya cases",""), x= ifelse(i > 20, "Date","")) +
    scale_linetype_manual(values = c("dashed", "dotted", "solid")) +
    scale_color_manual(values = c("red", "black", "black")) +
    theme_bw() +
    guides(color = guide_legend(override.aes = list(linetype = c("dashed", "dotted", "solid")), 
                                title = "Legend"),
           linetype = "none")  -> g_list_y1[[i]]
  
  data_real |> 
    filter(arbo == "Den") |> 
    ggplot() +
    geom_line(aes(x = date, y = y, 
                  colour = "Real observed",
                  linetype = "Real observed"), size = 1) +
    geom_line(mapping = aes(x = date, y = y2_obs, 
                            colour = "Observed at the moment",
                            linetype = "Observed at the moment"), data = data_no_delay,
              size = 1) +
    geom_line(mapping = aes(x = date, y = Median_y2, 
                            colour = "Estimated", 
                            linetype = "Estimated"),
              size = 1, data = summary_list_60[[i]] |> inner_join(data_no_delay, by = "time")) +
    geom_ribbon(aes(ymin = LI_y2, ymax = LS_y2, x = date), fill = "black", alpha = 0.25,
                data = summary_list_60[[i]] |> inner_join(data_no_delay, by = "time")) +
    scale_x_date(limits = c(ymd("2019-01-01"), ymd("2019-10-01"))) +
    labs(y = ifelse(i %% 5 == 1,"Dengue cases",""), x= ifelse(i > 20, "Date","")) +
    scale_linetype_manual(values = c("dashed", "dotted", "solid")) +
    scale_color_manual(values = c("red", "black", "black")) +
    theme_bw() +
    guides(color = guide_legend(override.aes = list(linetype = c("dashed", "dotted", "solid")), 
                                title = "Legend"),
           linetype = "none") -> g_list_y2[[i]]
  
}

# leg1 <- cowplot::get_legend(g_list_y1[[1]])
# 
# g_list_y1 <- lapply(g_list_y1, \(x) x + theme(legend.position = "none"))


library(grid)
ggpubr::ggarrange(plotlist = lapply(g_list_y1, \(x) + rremove("ylab") + rremove("xlab")), 
                  common.legend = T, legend= "bottom",align = "hv") 



gridExtra::grid.arrange(grobs = g_list_y1, ncol = 5, bottom = leg1)
gridExtra::grid.arrange(grobs = g_list_y2, ncol = 5)

data_no_delay <- data_model_list[[1]] |> 
  mutate(date = ymd("2017-01-01") + weeks(time - 1),
         y1_obs = as.numeric(y1_obs)) |> 
  group_by(date, time) |>
  summarise_at(vars(y1_obs:y2), sum, na.rm = T)

data_real |> 
  filter(arbo == "Chik") |> 
  ggplot() +
  geom_line(aes(x = date, y = y)) +
  geom_line(mapping = aes(x = date, y = y1_obs), data = data_no_delay, linetype = "dashed") +
  geom_line(mapping = aes(x = date, y = Median_y1), colour = "red", linetype = "dashed",
            data = summary_list_60[[1]] |> inner_join(data_no_delay, by = "time")) +
  geom_ribbon(aes(ymin = LI_y1, ymax = LS_y1, x = date), fill = "black", alpha = 0.1,
              data = summary_list_60[[1]] |> inner_join(data_no_delay, by = "time")) +
  scale_x_continuous(limits = c(ymd("2019-01-01"), ymd("2019-12-31"))) + 
  theme_bw()


### fit
fit_files <- list.files(pattern = "bivariate_fit.*_25runs_regressora_fixo")
load(fit_files)

fit_files <- list.files(pattern = "bivariate_fit.*_25runs")
load("bivariate_fit_delay_15_window_70_25runs.RData")


output <- bivariate_fit[[1]]$inla_output



vectime = 1:139
effects_time <- data.frame(time = vectime,
                           #date = unique(data_real$date)[1:124],
                           time_chik =   output$summary.random$time$`0.5quant`[vectime] + output$summary.random$time1$`0.5quant`,
                           time_dengue = output$summary.random$time$`0.5quant`[vectime + length(vectime)] + output$summary.random$time2$`0.5quant`,
                           time_uni = output$summary.random$time1$`0.5quant`)

effects_time |> 
  ggplot() +
  geom_line(aes(x = time, y = time_chik, col = "Chikungunya"), size = 1) +
  geom_line(aes(x = time, y = time_dengue, col = "Dengue"), size = 1) +
  # geom_line(aes(x = time, y = time_uni, col = "common", linetype = "Common"),
  #           size = 1, linetype = 2) +
  labs(x = "Time", y = expression(paste("Effects ",(delta[t] + beta["i,t"]))), colour = "Disease") + 
  scale_x_continuous(breaks = seq(0, 140, by = 20), labels = seq(0, 140, by = 20)) + 
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11)  
  ) -> gtime

gtime


effects_delay <- data.frame(delay = 1:16,
                            delay_chik =    output$summary.fixed$`0.5quant`[3] + output$summary.random$delay$`0.5quant`[1:16] + output$summary.random$delay1$`0.5quant`,
                            delay_dengue =  output$summary.fixed$`0.5quant`[4] + output$summary.random$delay$`0.5quant`[1:16 + 16] + output$summary.random$delay2$`0.5quant`)

effects_delay |> 
  ggplot() +
  geom_line(aes(x = delay, y = delay_chik, col = "Chikungunya"), size = 1) +
  geom_line(aes(x = delay, y = delay_dengue, col = "Dengue"), size = 1) +
  scale_x_continuous(breaks = seq(1, 16, by = 2), labels = seq(0, 15, by = 2)) + 
  labs(x = "Delay", y = expression(paste("Effects ",(psi[d] + gamma["i,d"]))), colour = "Disease") + 
  theme_bw() +
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
#ggsave("effects_fulldata.pdf",height = 4.75, width = 9.75)



summaries <- bivar_nowcasting_summary(bivariate_samples[[1]])


i=1
data_no_delay <- data_model_list[[i]] |> 
  mutate(date = ymd("2017-01-01") + weeks(time - 1),
         y1_obs = as.numeric(y1_obs)) |> 
  group_by(date, time) |>
  summarise_at(vars(y1_obs:y2), sum, na.rm = T)

data_real |> 
  filter(arbo == "Chik",
         date <= ymd("2019-10-01")) |>
  ggplot() +
  geom_line(aes(x = date, y = y, 
                colour = "Real observed",
                linetype = "Real observed"), size = 1) +
  geom_line(mapping = aes(x = date, y = y1_obs, 
                          colour = "Observed at the moment",
                          linetype = "Observed at the moment"), data = data_no_delay,
            size = 1) +
  geom_line(mapping = aes(x = date, y = Median_y1, 
                          colour = "Estimated", 
                          linetype = "Estimated"),
            size = 1, data = summaries |> inner_join(data_no_delay, by = "time")) +
  geom_ribbon(aes(ymin = LI_y1, ymax = LS_y1, x = date), fill = "black", alpha = 0.25,
              data = summaries |> inner_join(data_no_delay, by = "time")) +
  scale_x_date(limits = c(ymd("2019-01-01"), ymd("2019-10-01"))) +
  labs(y = ifelse(i %% 5 == 1,"Chikungunya cases",""), x= ifelse(i > 20, "Date","")) +
  scale_linetype_manual(values = c("dashed", "dotted", "solid")) +
  scale_color_manual(values = c("red", "black", "black")) +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(linetype = c("dashed", "dotted", "solid")), 
                              title = "Legend"),
         linetype = "none")  -> g_y1

data_real |> 
  filter(arbo == "Den") |> 
  ggplot() +
  geom_line(aes(x = date, y = y, 
                colour = "Real observed",
                linetype = "Real observed"), size = 1) +
  geom_line(mapping = aes(x = date, y = y2_obs, 
                          colour = "Observed at the moment",
                          linetype = "Observed at the moment"), data = data_no_delay,
            size = 1) +
  geom_line(mapping = aes(x = date, y = Median_y2, 
                          colour = "Estimated", 
                          linetype = "Estimated"),
            size = 1, data =summaries |> inner_join(data_no_delay, by = "time")) +
  geom_ribbon(aes(ymin = LI_y2, ymax = LS_y2, x = date), fill = "black", alpha = 0.25,
              data = summaries |> inner_join(data_no_delay, by = "time")) +
  scale_x_date(limits = c(ymd("2019-01-01"), ymd("2019-10-01"))) +
  labs(y = ifelse(i %% 5 == 1,"Dengue cases",""), x= ifelse(i > 20, "Date","")) +
  scale_linetype_manual(values = c("dashed", "dotted", "solid")) +
  scale_color_manual(values = c("red", "black", "black")) +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(linetype = c("dashed", "dotted", "solid")), 
                              title = "Legend"),
         linetype = "none") -> g_y2


nb_list <- list(g_y1, g_y2)
# poisson_list <- list(g_y1, g_y2)

graphs_list <- c(nb_list, poisson_list)
ggarrange(plotlist = graphs_list, ncol=2, nrow=2,
          align = "hv",
          common.legend = TRUE, legend="bottom",
          labels = "AUTO") -> figure

ggarrange(plotlist = nb_list, ncol=2,
          align = "hv",
          common.legend = TRUE, legend = "bottom") 

ggsave("nb_compare_pois_lengeda_bottom.pdf") 

ggarrange(plotlist = poisson_list, ncol=2,
          align = "hv",
          common.legend = TRUE, legend="bottom") 

ggsave("pois_compare_legenda_bottom.pdf")


ggpubr::ggarrange(plotlist = g_list_y1,
                  common.legend = T, legend= "top") -> figure


library(grid)
annotate_figure(figure, 
                left = textGrob("Poisson                  Negative Binomial",
                                rot = 90, hjust = 0.33 , vjust = 0.9, gp = gpar(cex = 1.3)),
                )



output$summary.hyperpar |> 
  select(4,3,5) |> 
  xtable::xtable()
















