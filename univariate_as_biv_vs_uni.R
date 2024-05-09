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

M5 <- Y ~ -1 + alphai +
  + f(time,  model = "rw2", replicate = idrep,
      hyper = list(prec = list(param = c(0.001, 0.001)))) +
  + f(delay, model = "rw1", replicate = idrep,
      hyper = list(prec = list(param = c(0.001, 0.001)))) 

M1 <- Y ~ 1 +
  f(time,
    model = "rw2",
    hyper = list("prec" = list(prior = "loggamma",
                               param = c(0.001, 0.001))
    )) +
  f(delay, model = "rw1",
    hyper = list("prec" = list(prior = "loggamma",
                               param = c(0.001, 0.001)))
  )

fit_biv_uni <- run_inla_model(data = data_model, model = M5)
samples_biv_uni <- sample_predictive(inla_output = fit_biv_uni$inla_output,
                                     joint_data = fit_biv_uni$joint_data)

summary_biv_uni <- bivar_nowcasting_summary(samples_biv_uni)


dados_chik <- data_model %>% 
  select(time, delay, report_time, Y = y1_obs, count = y1)

dados_dengue <- data_model %>% 
  select(time, delay, report_time, Y = y2_obs, count = y2)


output_chik <- INLA::inla(M1, family = "nbinomial",
                          data = dados_chik,
                          control.predictor = list(link = 1, compute = T),
                          control.compute = list(config = T, waic=T, dic=T,
                                                 cpo = T),
                          control.family = list(
                            hyper = list("theta" = list(prior = "loggamma",
                                                        param = c(0.001,
                                                                  0.001))
                            )
                          )
)

output_dengue <- INLA::inla(M1, family = "nbinomial",
                            data = dados_dengue,
                            control.predictor = list(link = 1, compute = T),
                            control.compute = list(config = T, waic=T, dic=T,
                                                   cpo = T),
                            control.family = list(
                              hyper = list("theta" = list(prior = "loggamma",
                                                          param = c(0.001,
                                                                    0.001))
                              )
                            )
)

# fit_uni(dados, disease = "dengue", evol_time = "rw2", evol_delay = "rw1")


ts.plot(fit_biv_uni$inla_output$summary.random$time$`0.5quant`)
lines(c(output_chik$summary.random$time$`0.5quant`,output_dengue$summary.random$time$`0.5quant`), 
      col = 2)

cbind(fit_biv_uni$inla_output$summary.random$time$`0.5quant`,
      c(output_chik$summary.random$time$`0.5quant`,
        output_dengue$summary.random$time$`0.5quant`)) |> round(1)

ts.plot(fit_biv_uni$inla_output$summary.random$delay$`0.5quant`)
lines(c(output_chik$summary.random$delay$`0.5quant`,output_dengue$summary.random$delay$`0.5quant`), 
      col = 2)

cbind(fit_biv_uni$inla_output$summary.random$delay$`0.5quant`,
      c(output_chik$summary.random$delay$`0.5quant`,
        output_dengue$summary.random$delay$`0.5quant`)) |> round(2)


fit_biv_uni$inla_output$summary.hyperpar$mean
output_chik$summary.hyperpar$mean
output_dengue$summary.hyperpar$mean


fit_biv_uni$inla_output$summary.hyperpar$`0.5quant`
output_chik$summary.hyperpar$`0.5quant`
output_dengue$summary.hyperpar$`0.5quant`


#### AMOSTRA dengue ####
index.missing <- which(is.na(dados_dengue$Y))


## Step 1: Sampling from the approximate posterior distribution using INLA
srag.samples0.list <- INLA::inla.posterior.sample(n = 1000, output_dengue)

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
gg.age_uni <- function(x, dados.gg, idx){
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

## Step 4: Applying the age aggregation on each posterior
tibble.samples.0 <- lapply( X = vector.samples0,
                            FUN = gg.age_uni,
                            dados = dados_dengue,
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


summy_dengue <- total.summy


##### Amostra chik ####
index.missing <- which(is.na(dados_chik$Y))


## Step 1: Sampling from the approximate posterior distribution using INLA
srag.samples0.list <- INLA::inla.posterior.sample(n = 1000, output_chik)

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

## Step 4: Applying the age aggregation on each posterior
tibble.samples.0 <- lapply( X = vector.samples0,
                            FUN = gg.age_uni,
                            dados = dados_chik,
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


summy_chik <- total.summy


# save(summy_dengue,
#      summy_chik,
#      summary_biv_uni,
#      file = "Summaries_dengue_chik_uni.RData")

load("Summaries_dengue_chik_uni.RData")
# summary_biv_uni |> 
#   melt()
#   tidyr::pivot_longer(cols = starts_with("Median"),
#                       names_to = "Disease", values_to = "Count") 

uni_summy_join <- inner_join(summy_chik, summy_dengue, by = "time", suffix = c("_y1","_y2"))


pivot_longer(summary_biv_uni, cols=-1, names_pattern = "(.*)(..)$", names_to = c("limit", "name")) %>% 
    mutate(limit=ifelse(limit=="", "value", limit)) |> 
    mutate(model = "Univariado conjunto") -> long_biv_uni
  
pivot_longer(uni_summy_join, cols=-1, names_pattern = "(.*)(..)$", names_to = c("limit", "name")) %>% 
  mutate(limit=ifelse(limit=="", "value", limit)) |> 
  mutate(model = "Dois univariados") -> long_uni_join


bind_rows(long_biv_uni, long_uni_join) |> 
  mutate(limit = substr(limit, 1, nchar(limit) - 1)) |> 
  filter(limit %in% c("LI", "LS"),
         name == "y1") |> 
  pivot_wider(names_from = limit, values_from = value) |> 
  mutate(range_90 = LS - LI) |> 
  ggplot() +
  aes(x = time, y = range_90, colour = model) +
  geom_line() +
  scale_colour_discrete(labels = c("Two univariate models",
                                   "Bivariate model without common effects")) +
  labs(title = "Chikungunya", y = "95% interval width",
       x = "Time", col = "Model") +
  theme_bw(base_size = 13) -> chik_uni_compare


bind_rows(long_biv_uni, long_uni_join) |> 
  mutate(limit = substr(limit, 1, nchar(limit) - 1)) |> 
  filter(limit %in% c("LI", "LS"),
         name == "y2") |> 
  pivot_wider(names_from = limit, values_from = value) |> 
  mutate(range_90 = LS - LI) |> 
  ggplot() +
  aes(x = time, y = range_90, colour = model) +
  geom_line() +
  scale_colour_discrete(labels = c("Two univariate models",
                                   "Bivariate model without common effects")) +
  labs(title = "Dengue", y = "95% interval width",
       x = "Time", col = "Model") +
  theme_bw(base_size = 13) -> dengue_uni_compare


ggpubr::ggarrange(chik_uni_compare, dengue_uni_compare,
                  common.legend = TRUE, legend= "bottom")

ggsave("univariate_compare.pdf", height = 4.43, width = 9.3)





# #### RANGE 50
# bind_rows(long_biv_uni, long_uni_join) |> 
#   mutate(limit = substr(limit, 1, nchar(limit) - 1)) |> 
#   filter(limit %in% c("LIb", "LSb"),
#          name == "y1") |> 
#   pivot_wider(names_from = limit, values_from = value) |> 
#   mutate(range_50 = LSb - LIb) |> 
#   ggplot() +
#   aes(x = time, y = range_50, colour = model) +
#   geom_line() +
#   labs(title = "Chikungunya") +
#   theme_bw() 
# 
# 
# bind_rows(long_biv_uni, long_uni_join) |> 
#   mutate(limit = substr(limit, 1, nchar(limit) - 1)) |> 
#   filter(limit %in% c("LIb", "LSb"),
#          name == "y2") |> 
#   pivot_wider(names_from = limit, values_from = value) |> 
#   mutate(range_50 = LSb - LIb) |> 
#   ggplot() +
#   aes(x = time, y = range_50, colour = model) +
#   geom_line() +
#   labs(title = "Dengue") +
#   theme_bw() 


output <- fit_biv_uni$inla_output



vectime = unique(output$summary.random$time$ID)
effects_time <- data.frame(time = vectime,
                           #date = unique(data_real$date)[1:124],
                           time_chik =   output$summary.random$time$`0.5quant`[vectime],
                           time_dengue = output$summary.random$time$`0.5quant`[vectime + length(vectime)],
                           time_dengue_uni = output_dengue$summary.random$time$`0.5quant`,
                           time_chik_uni = output_chik$summary.random$time$`0.5quant`)

effects_time |> 
  ggplot() +
  geom_line(aes(x = time, y = time_chik, 
                col = "Chikungunya", lty = "Bivariate"), size = 1) +
  geom_line(aes(x = time, y = time_dengue, 
                col = "Dengue", lty = "Bivariate"), size = 1) +
  geom_line(aes(x = time, y = time_dengue_uni, lty = "Univariate"), 
            size = 1) +
  geom_line(aes(x = time, y = time_chik_uni, lty = "Univariate"), 
            size = 1) +
  scale_linetype_manual(labels = c("Univariate", 
                                   "Independent bivariate"),
                        values = c(1,2)) +
  labs(x = "Time", y = expression(paste("Effects ",(delta[t] + beta["i,t"]))),
       colour = "Disease",linetype = "Model") + 
  scale_x_continuous(breaks = seq(0, 140, by = 20), labels = seq(0, 140, by = 20)) + 
  theme_bw() +
  theme(legend.position = "bottom",
        legend.box="vertical",
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11)  
  ) -> gtime

gtime

effects_delay <- data.frame(delay = 1:16,
                            delay_chik =  output$summary.random$delay$`0.5quant`[1:16],
                            delay_dengue =  output$summary.random$delay$`0.5quant`[1:16 + 16],
                            delay_dengue_uni = output_dengue$summary.random$delay$`0.5quant`,
                            delay_chik_uni = output_chik$summary.random$delay$`0.5quant`)

effects_delay |> 
  ggplot() +
  geom_line(aes(x = delay, y = delay_chik, 
                col = "Chikungunya", lty = "Bivariate"), size = 1) +
  geom_line(aes(x = delay, y = delay_dengue, 
                col = "Dengue", lty = "Bivariate"), size = 1) +
  geom_line(aes(x = delay, y = delay_dengue_uni, lty = "Univariate"), 
            size = 1) +
  geom_line(aes(x = delay, y = delay_chik_uni, lty = "Univariate"), 
            size = 1) +
  scale_x_continuous(breaks = seq(1, 16, by = 2), labels = seq(0, 15, by = 2)) + 
  scale_linetype_manual(labels = c("Univariate", 
                                   "Independent bivariate"),
                        values = c(1,2)) + 
  labs(x = "Delay", 
       y = expression(paste("Effects ",(psi[d] + gamma["i,d"]))), 
       colour = "Disease",
       linetype = "Model") + 
  theme_bw() +
  theme(legend.position = "bottom",
        legend.box="vertical",
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11)  
  ) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE)) -> gdelay

gdelay

library(ggpubr)
ggarrange(gtime, gdelay, ncol=2, nrow=1, common.legend = TRUE, legend="bottom") 

ggsave("effects_univariate_compare.pdf", height = 4.43, width = 9.3)







