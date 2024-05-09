## Geração e tentativa de ajuste de modelo para os dados Bivariados 

# Pacotes
library(INLA)
library(ggplot2)
library(dplyr)
library(nowcaster)
library(truncnorm)

source("functions_multivar_nowcasting.R")



# Simulação ---------------------------------------------------------------


## Modelo: log(\lambda_{itd}) = \alpha_i + \beta_{i,t} + \gamma_{i,d} + \delta_{t} + \psi_{d}

max_t <- 200
D <- 15
y <- matrix(NA, max_t + D, D)


alpha_1 <- 2
alpha_2 <- 1
beta    <- matrix(nrow = max_t, ncol = 2)
gamma_d <- matrix(nrow = D + 1, ncol = 2)
delta   <- vector(length = max_t)
psi     <- vector(length = D + 1)
x <- vector(length = max_t)
beta_x <- vector(length = max_t)


set.seed(123)


beta[1, ] <- rnorm(2)
gamma_d[1,] <- rnorm(2)
delta[1] <- rnorm(1)
psi[1] <- rnorm(1)
x[1] <- rnorm(1)
beta_x[1] <- rnorm(1)


tau_time <- 50
tau_time1 <- 5
tau_delay <- 50
tau_delay1 <- 5
tau_beta_x <- 0.5
tau_x <- 10

phi1 <- 0.4
phi52 <- 0.3

# Geração dos efeitos/covariaveis
for(t in 2:max_t){
  # Efeitos do tempo (por doença e conjunto)
  beta[t, ] <- rnorm(2, mean = beta[t-1, ], sd = sqrt(1/tau_time))
  delta[t] <- rnorm(1, mean = delta[t-1], sd = sqrt(1/tau_time1))
  
  # Regressora e efeito da regressora
  x[t] <- x[t-1] + rnorm(1, sd = sqrt(1/tau_x))
  beta_x[t] <- beta_x[t-1] + rnorm(1, sd = sqrt(1/tau_beta_x))
  
  # if(t <= 52){
  #   X[t] <- phi1*X[t-1] + rnorm(1, 0, sqrt(1/tau_param))
  # }else{
  #   X[t] <- phi1*(X[t-1]) + phi52*(X[t-52])  + rnorm(1, 0, sqrt(.25))
  # }
}

for(d in 2:(D + 1)){
  gamma_d[d, ] <- rnorm(n = 2, mean = gamma_d[d-1, ], sd = sqrt(1/tau_delay))
  
  psi[d] <- rtruncnorm(n = 1, mean = psi[d - 1], sd = sqrt(1/tau_delay1), 
                       a = -Inf, b = psi[d - 1])
  
}



ts.plot(beta)
ts.plot(gamma_d)
ts.plot(delta)
ts.plot(psi)
ts.plot(x)
ts.plot(beta_x)

# beta <- apply(beta,2, \(x) x - mean(x))
# gamma_d <- apply(gamma_d,2, \(x) x - mean(x))
# delta <- delta - mean(delta)
# psi = psi - mean(psi)



x = x/10


input1 <- cbind(rep(alpha_1, max_t*(D+1)), 
                rep(gamma_d[,1], max_t),
                rep(beta[,1], each = D+1), 
                rep(delta, each = D+1),
                rep(psi, max_t), 
                rep(beta_x*x, each = D+1))

input2 <- cbind(rep(alpha_2, max_t*(D+1)), 
                rep(gamma_d[,2], max_t),
                rep(beta[,2], each = D+1), 
                rep(delta, each = D+1),
                rep(psi, max_t),
                rep(beta_x*x, each = D+1))

LP1 <- rowSums(input1)
LP2 <- rowSums(input2)

#y <- rnbinom(n = 2*(D+1)*max_t, size = 5, mu = exp(c(LP1, LP2)))

# Alternativa
y <- cbind(rnbinom(n = (D+1)*max_t, size = 5, mu = exp(LP1)),
           rnbinom(n = (D+1)*max_t, size = 5, mu = exp(LP2)))

cor(LP1, LP2)
cor(y)

colnames(y) <- c("y1", "y2")


y_df <- cbind.data.frame(y, 
                         time = rep(1:max_t, each = D+1), 
                         delay = rep(0:D, max_t)) %>% 
  mutate(report_time = time + delay)

y_df %>% 
  group_by(time) %>% 
  summarise(y1 = sum(y1),
            y2 = sum(y2)) %>% 
  tidyr::pivot_longer(y1:y2) %>%
  ggplot(aes(x = time, y = value, colour = name)) +
  geom_line() +
    scale_color_manual(values = c('y1' = '#F8766D',
                                  'y2' = '#00BFC4'),
                       labels = c(expression(N["1,t"]), expression(N["2,t"]))) +
  labs(x = "Time", y = "Count", colour = "") +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")
  
  
y_df %>% 
  mutate_at(1:2,scale) |> 
  group_by(time) %>% 
  summarise(y1 = sum(y1),
            y2 = sum(y2)) %>% 
  tidyr::pivot_longer(y1:y2) %>%
  ggplot(aes(x = time, y = value, colour = name)) +
  geom_line() +
  scale_color_manual(values = c('y1' = '#F8766D',
                                'y2' = '#00BFC4'),
                     labels = c(expression(N["1,t"]), expression(N["2,t"]))) +
  labs(x = "t", y = "Count", colour = "") +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")



y_new <- y_df %>% 
  mutate(y1_obs = ifelse(report_time > max_t, NA, y1),
         y2_obs = ifelse(report_time > max_t, NA, y2))


round(alpha_1 + mean(delta) + mean(psi) + mean(beta[,1]) + mean(gamma_d[,1]),2)
round(alpha_2 + mean(delta) + mean(psi) + mean(beta[,2]) + mean(gamma_d[,2]),2)
mean(beta_x)

# Joint_data --------------------------------------------------------

n <- nrow(y)
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
                             idrep = rep(1:2, each = n),
                             id_x = c(y_df$time, y_df$time),
                             id_x1 = c(y_df$time, rep(NA,n)),
                             id_x2 = c(rep(NA,n), y_df$time),
                             x1 = c(rep(x, each = D + 1), rep(NA,n)),
                             x2 = c(rep(NA,n)), rep(x, each = D + 1))


random_covariate <- list(time1  = c(y_df$time, rep(NA, n)),
                         delay1 = c(y_df$delay, rep(NA, n)),
                         time2  = c(rep(NA, n), y_df$time),
                         delay2 = c(rep(NA, n), y_df$delay),
                         time = rep(y_df$time, 2),
                         delay = rep(y_df$delay, 2),
                         x = rep(x, each = D + 1))




Y = y_bvar
joint_data <- cbind(Y, fixed_covariate, random_covariate)



# Modelos -----------------------------------------------------------------


# Modelo bivariado considerando dois efeitos diferentes para a covariável.
joint_model <- Y ~ -1 + alphai + 
  + f(time,  model = "rw1", replicate = idrep,
      hyper = list(prec = list(param = c(0.001, 0.001)))) +
  + f(delay, model = "rw1", replicate = idrep,
      hyper = list(prec = list(param = c(0.001, 0.001)))) +
  + f(time1, model = "rw1",
      hyper = list(prec = list(param = c(0.001, 0.001)))) +
  + f(time2, copy  = "time1", fixed = TRUE) + 
  + f(delay1, model = "rw1",
      hyper = list(prec = list(param = c(0.001, 0.001)))) +
  + f(delay2, copy  = "delay1", fixed = TRUE) +
  + x + 
  # + f(id_x, x, model = "rw1", replicate = idrep,
  #     hyper = list(prec = list(param = c(0.001, 0.001))), initial = 0) 
  + f(id_x1, x, model = "rw1",
      hyper = list(prec = list(param = c(0.001, 0.001))), initial = 0) +
  + f(id_x2, x, copy = "id_x1", fixed = T)
  

teste_regressora <- INLA::inla(joint_model, family = rep("nbinomial", 2),
                      data = joint_data,
                      control.predictor = list(link = 1, compute = T),
                      control.compute = list( config = T, waic=T, dic=T),
                      control.family = list(
                        list(hyper = list("theta" = list(prior = "loggamma", param = c(0.001, 0.001)))),
                        list(hyper = list("theta" = list(prior = "loggamma", param = c(0.001, 0.001))))
                      ),
                      verbose = F
                      
)



modelo_uni <- Y ~ -1 + alphai +
  + f(time,  model = "rw1", replicate = idrep,
      hyper = list("theta" = list(prior = "loggamma"))) +
  + f(delay, model = "rw1", replicate = idrep,
      hyper = list("theta" = list(prior = "loggamma"))) +
  + x +
  + f(id_x, x, model = "rw1", replicate = idrep,
      hyper = list(prec = list(param = c(0.001, 0.001))), initial = 0) 

teste_uni <- INLA::inla(modelo_uni, family = rep("nbinomial", 2),
                        data = joint_data,
                        control.predictor = list(link = 1, compute = T),
                        control.compute = list( config = T, waic=T, dic=T),
                        control.family = list(
                          list(hyper = list("theta" = list(prior = "loggamma", param = c(0.001, 0.001)))),
                          list(hyper = list("theta" = list(prior = "loggamma", param = c(0.001, 0.001))))
                        ),
                        verbose = F
                        
)


teste_poisson <- run_inla_model(data = NULL, family = "poisson",
                                model = joint_model,
                                data_inla = joint_data)

teste_poisson <- teste_poisson$inla_output




summary(teste_regressora)
summary(teste_uni)
summary(teste_poisson)


# Gráficos Bivariado ------------------------------------------------------



library(purrr)
samples <- inla.posterior.sample(teste_regressora, n = 1000)

map(samples, \(sample_1){
  sample_1$latent[grepl("id_x",rownames(sample_1$latent))] +
    sample_1$latent[grepl("\\bx\\b",rownames(sample_1$latent))]
  }) |> reduce(rbind) |>
  apply(MARGIN = 2, function(x){quantile(x, c(0.025, 0.5, 0.975))}) -> quantiles


rownames(quantiles) <- c("lower", "median", "upper")
quantiles[,1:200] |> 
  t() |> 
  as.data.frame() |> 
  ggplot() +
  geom_line(aes(x = 1:200, y = median, colour = "Posterior median")) +
  geom_ribbon(aes(x = 1:200, ymax = upper, ymin = lower),
              alpha = 0.2) +
  geom_line(data = NULL, 
            aes(x = 1:200, y  = beta_x, colour = "Real")) +
  labs(x = "t", y = expression(nu[t]), colour = "") +
  theme_bw(base_size = 13) -> gregressor



# lines(teste_regressora$summary.random$id_x1$`0.5quant` + teste_regressora$summary.fixed$`0.5quant`[3], ylim = c(-25,5), lty = 2, lwd = 2)
# lines(teste_regressora$summary.random$id_x1$`0.975quant`+ teste_regressora$summary.fixed$`0.5quant`[3], col = 2, lty = 2, lwd = 2)
# lines(teste_regressora$summary.random$id_x1$`0.025quant`+ teste_regressora$summary.fixed$`0.5quant`[3], col = 2, lty = 2, lwd = 2)
# abline(h=0)
# lines(beta_x, col = "green")


ggplot(data = teste_regressora$summary.random$id_x1) +
  geom_line(aes(x = ID, y = `0.5quant`, colour = "Posterior median")) +
  geom_ribbon(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`),
              alpha = 0.2) +
  geom_line(data = NULL, 
            aes(x = 1:200, y  = beta_x - mean(beta_x), colour = "Real")) +
  labs(x = "t", y = expression(nu[t]), colour = "") +
  theme_bw(base_size = 13) -> gregressor


## Regressora
ggplot(data = teste_regressora$summary.random$id_x1) +
  geom_line(aes(x = ID, y = `0.5quant`, colour = "Posterior median")) +
  geom_ribbon(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`),
              alpha = 0.2) +
  geom_line(data = NULL, 
            aes(x = 1:200, y  = beta_x - mean(beta_x), colour = "Real")) +
  labs(x = "t", y = expression(nu[t]), colour = "") +
  theme_bw(base_size = 13) -> gregressor

# ggplot(data = teste_regressora$summary.random$id_x[201:400,]) +
#   geom_line(aes(x = ID, y = `0.5quant`, colour = "Posterior median")) +
#   geom_ribbon(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`),
#               alpha = 0.2) +
#   geom_line(data = NULL, 
#             aes(x = 1:200, y  = beta_x - mean(beta_x), colour = "Real")) +
#   labs(x = "t", y = expression(nu[t]), colour = "") +
#   theme_bw(base_size = 13) -> gregressor2

## Delta (comum do tempo)
ggplot(data = teste_regressora$summary.random$time1) +
  geom_line(aes(x = ID, y = `0.5quant`, colour = "Posterior median")) +
  geom_ribbon(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`),
              alpha = 0.2) +
  geom_line(data = NULL, 
            aes(x = 1:200, y  = delta - mean(delta), colour = "Real")) +
  labs(x = "t", y = expression(delta[t]), color = "") +
  theme_bw(base_size = 13) -> gtime_delta

## beta_t,1
ggplot(data = teste_regressora$summary.random$time[1:200,]) +
  geom_line(aes(x = ID, y = `0.5quant`, colour = "Posterior median")) +
  geom_ribbon(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`),
              alpha = 0.2) +
  geom_line(data = NULL, 
            aes(x = 1:200, y  = beta[,1] - mean(beta[,1]), colour = "Real")) +
  labs(x = "t", y = expression(beta["1,t"]), color = "") +
  theme_bw(base_size = 13) -> gtime_beta1

## beta_t,2
ggplot(data = teste_regressora$summary.random$time[201:400,]) +
  geom_line(aes(x = ID, y = `0.5quant`, colour = "Posterior median")) +
  geom_ribbon(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`),
              alpha = 0.2) +
  geom_line(data = NULL, 
            aes(x = 1:200, y  = beta[,2] - mean(beta[,2]), colour = "Real")) +
  labs(x = "t", y = expression(beta["2,t"]), color = "")  +
  theme_bw(base_size = 13) -> gtime_beta2

## psi (comum delay)
ggplot(data = teste_regressora$summary.random$delay1) +
  geom_line(aes(x = ID, y = `0.5quant`, colour = "Posterior median")) +
  geom_ribbon(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`),
              alpha = 0.2) +
  geom_line(data = NULL, 
            aes(x = 0:15, y  = psi - mean(psi), colour = "Real")) +
  labs(x = "d", y = expression(psi[d]), colour = "") +
  theme_bw(base_size = 13) -> gdelay_psi

## gamma_d,1
ggplot(data = teste_regressora$summary.random$delay[1:16,]) +
  geom_line(aes(x = ID, y = `0.5quant`, colour = "Posterior median")) +
  geom_ribbon(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`),
              alpha = 0.2) +
  geom_line(data = NULL, 
            aes(x = 0:15, y  = gamma_d[,1] - mean(gamma_d[,1]), colour = "Real")) +
  labs(x = "d", y = expression(gamma["1,d"]), colour = "") +
  theme_bw(base_size = 13) -> gdelay_gamma1

## gamma_d,2
ggplot(data = teste_regressora$summary.random$delay[1:16 + 16,]) +
  geom_line(aes(x = ID, y = `0.5quant`, colour = "Posterior median")) +
  geom_ribbon(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`),
              alpha = 0.2) +
  geom_line(data = NULL, 
            aes(x = 0:15, y  = gamma_d[,2] - mean(gamma_d[,2]), colour = "Real")) +
  labs(x = "d", y = expression(gamma["2,d"]), colour = "") +
  theme_bw(base_size = 13) -> gdelay_gamma2

library(ggpubr)
library(patchwork)
ggarrange(gtime_delta, gdelay_psi,
          gtime_beta1, gdelay_gamma1,
          gtime_beta2, gdelay_gamma2,
          common.legend = TRUE, 
          legend = "bottom",
          ncol = 2,
          nrow = 3) / gregressor


(gtime_delta + gdelay_psi)/
(gtime_beta1 + gdelay_gamma1)/
(gtime_beta2 + gdelay_gamma2 )/ gregressor +   
  plot_layout(guides = "collect") &  theme(legend.position = "bottom")

# ggsave("simulado_effects_estimate_regressora_varia.pdf", height = 8, width = 6.5)

# Gráficos Univariado -----------------------------------------------------
## Regressora
ggplot(data = teste_uni$summary.random$id_x[1:200,]) +
  geom_line(aes(x = ID, y = `0.5quant`, colour = "Posterior median")) +
  geom_ribbon(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`),
              alpha = 0.2) +
  geom_line(data = NULL,
            aes(x = 1:200, y  = beta_x - mean(beta_x), colour = "Real")) +
  labs(x = "t", y = expression(nu["1,t"]), colour = "") +
  theme_bw(base_size = 13) -> uni_gregressor1
# 
ggplot(data = teste_uni$summary.random$id_x[1:200 + 200,]) +
  geom_line(aes(x = ID, y = `0.5quant`, colour = "Posterior median")) +
  geom_ribbon(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`),
              alpha = 0.2) +
  geom_line(data = NULL,
            aes(x = 1:200, y  = beta_x - mean(beta_x), colour = "Real")) +
  labs(x = "t", y = expression(nu["2,t"]), colour = "") +
  theme_bw(base_size = 13) -> uni_gregressor2
# 
# 
# ## beta_t,1
# ggplot(data = teste_uni$summary.random$time[1:200,]) +
#   geom_line(aes(x = ID, y = `0.5quant`, colour = "Posterior median")) +
#   geom_ribbon(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`),
#               alpha = 0.2) +
#   geom_line(data = NULL, 
#             aes(x = 1:200, y  = beta[,1] - mean(beta[,1]), colour = "Real")) +
#   labs(x = "t", y = expression(beta["1,t"]), color = "") +
#   theme_bw(base_size = 13) -> uni_gtime_beta1
# 
# ## beta_t,2
# ggplot(data = teste_uni$summary.random$time[201:400,]) +
#   geom_line(aes(x = ID, y = `0.5quant`, colour = "Posterior median")) +
#   geom_ribbon(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`),
#               alpha = 0.2) +
#   geom_line(data = NULL, 
#             aes(x = 1:200, y  = beta[,2] - mean(beta[,2]), colour = "Real")) +
#   labs(x = "t", y = expression(beta["2,t"]), color = "")  +
#   theme_bw(base_size = 13) -> uni_gtime_beta2
# 
# ## psi (comum delay)
# ggplot(data = teste_uni$summary.random$delay1) +
#   geom_line(aes(x = ID, y = `0.5quant`, colour = "Posterior median")) +
#   geom_ribbon(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`),
#               alpha = 0.2) +
#   geom_line(data = NULL, 
#             aes(x = 0:15, y  = psi - mean(psi), colour = "Real")) +
#   labs(x = "d", y = expression(psi[d]), colour = "") +
#   theme_bw(base_size = 13) -> uni_gdelay_psi
# 
# ## gamma_d,1
# ggplot(data = teste_uni$summary.random$delay[1:16,]) +
#   geom_line(aes(x = ID, y = `0.5quant`, colour = "Posterior median")) +
#   geom_ribbon(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`),
#               alpha = 0.2) +
#   geom_line(data = NULL, 
#             aes(x = 0:15, y  = gamma_d[,1] - mean(gamma_d[,1]), colour = "Real")) +
#   labs(x = "d", y = expression(gamma["1,d"]), colour = "") +
#   theme_bw(base_size = 13) -> uni_gdelay_gamma1
# 
# ## gamma_d,2
# ggplot(data = teste_uni$summary.random$delay[1:16 + 16,]) +
#   geom_line(aes(x = ID, y = `0.5quant`, colour = "Posterior median")) +
#   geom_ribbon(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`),
#               alpha = 0.2) +
#   geom_line(data = NULL, 
#             aes(x = 0:15, y  = gamma_d[,2] - mean(gamma_d[,2]), colour = "Real")) +
#   labs(x = "d", y = expression(gamma["2,d"]), colour = "") +
#   theme_bw(base_size = 13) -> uni_gdelay_gamma2
# 
# 
# ggarrange(uni_gtime_delta, uni_gdelay_psi,
#           uni_gtime_beta1, uni_gdelay_gamma1,
#           uni_gtime_beta2, uni_gdelay_gamma2,
#           common.legend = TRUE, 
#           legend = "bottom",
#           ncol = 2,
#           nrow = 3) / uni_gregressor
# 
# 
# (gtime_delta + gdelay_psi)/
#   (gtime_beta1 + gdelay_gamma1)/
#   (gtime_beta2 + gdelay_gamma2 )/ gregressor +   
#   plot_layout(guides = "collect") &  theme(legend.position = "bottom")




# Gráficos poisson --------------------------------------------------------
## Regressora
ggplot(data = teste_poisson$summary.random$id_x1) +
  geom_line(aes(x = ID, y = `0.5quant`, colour = "Posterior median")) +
  geom_ribbon(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`),
              alpha = 0.2) +
  geom_line(data = NULL, 
            aes(x = 1:200, y  = beta_x - mean(beta_x), colour = "Real")) +
  labs(x = "t", y = expression(nu[t]), colour = "") +
  theme_bw(base_size = 13) -> gregressor

## Delta (comum do tempo)
ggplot(data = teste_poisson$summary.random$time1) +
  geom_line(aes(x = ID, y = `0.5quant`, colour = "Posterior median")) +
  geom_ribbon(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`),
              alpha = 0.2) +
  geom_line(data = NULL, 
            aes(x = 1:200, y  = delta - mean(delta), colour = "Real")) +
  labs(x = "t", y = expression(delta[t]), color = "") +
  theme_bw(base_size = 13) -> gtime_delta

## beta_t,1
ggplot(data = teste_poisson$summary.random$time[1:200,]) +
  geom_line(aes(x = ID, y = `0.5quant`, colour = "Posterior median")) +
  geom_ribbon(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`),
              alpha = 0.2) +
  geom_line(data = NULL, 
            aes(x = 1:200, y  = beta[,1] - mean(beta[,1]), colour = "Real")) +
  labs(x = "t", y = expression(beta["1,t"]), color = "") +
  theme_bw(base_size = 13) -> gtime_beta1

## beta_t,2
ggplot(data = teste_poisson$summary.random$time[201:400,]) +
  geom_line(aes(x = ID, y = `0.5quant`, colour = "Posterior median")) +
  geom_ribbon(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`),
              alpha = 0.2) +
  geom_line(data = NULL, 
            aes(x = 1:200, y  = beta[,2] - mean(beta[,2]), colour = "Real")) +
  labs(x = "t", y = expression(beta["2,t"]), color = "")  +
  theme_bw(base_size = 13) -> gtime_beta2

## psi (comum delay)
ggplot(data = teste_poisson$summary.random$delay1) +
  geom_line(aes(x = ID, y = `0.5quant`, colour = "Posterior median")) +
  geom_ribbon(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`),
              alpha = 0.2) +
  geom_line(data = NULL, 
            aes(x = 0:15, y  = psi - mean(psi), colour = "Real")) +
  labs(x = "d", y = expression(psi[d]), colour = "") +
  theme_bw(base_size = 13) -> gdelay_psi

## gamma_d,1
ggplot(data = teste_poisson$summary.random$delay[1:16,]) +
  geom_line(aes(x = ID, y = `0.5quant`, colour = "Posterior median")) +
  geom_ribbon(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`),
              alpha = 0.2) +
  geom_line(data = NULL, 
            aes(x = 0:15, y  = gamma_d[,1] - mean(gamma_d[,1]), colour = "Real")) +
  labs(x = "d", y = expression(gamma["1,d"]), colour = "") +
  theme_bw(base_size = 13) -> gdelay_gamma1

## gamma_d,2
ggplot(data = teste_poisson$summary.random$delay[1:16 + 16,]) +
  geom_line(aes(x = ID, y = `0.5quant`, colour = "Posterior median")) +
  geom_ribbon(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`),
              alpha = 0.2) +
  geom_line(data = NULL, 
            aes(x = 0:15, y  = gamma_d[,2] - mean(gamma_d[,2]), colour = "Real")) +
  labs(x = "d", y = expression(gamma["2,d"]), colour = "") +
  theme_bw(base_size = 13) -> gdelay_gamma2


ggarrange(gtime_delta, gdelay_psi,
          gtime_beta1, gdelay_gamma1,
          gtime_beta2, gdelay_gamma2,
          common.legend = TRUE, 
          legend = "bottom",
          ncol = 2,
          nrow = 3) / gregressor


(gtime_delta + gdelay_psi)/
  (gtime_beta1 + gdelay_gamma1)/
  (gtime_beta2 + gdelay_gamma2 )/ gregressor +   
  plot_layout(guides = "collect") &  theme(legend.position = "bottom")

# ggsave("simulado_effects_estimate_poisson.pdf", height = 8, width = 6.5)

# Nowcasting Bivariado ----------------------------------------------------
set.seed(267)

load("samples_regressora_simulado_3models.RData")
# pred_samples <- sample_predictive(inla_output = teste_regressora, 
#                                   joint_data = joint_data)

summy <- bivar_nowcasting_summary(pred_samples)

dados_by_week <- y_new |>
  group_by(time) |>
  summarise_at(vars(starts_with("y")), .funs = list(~sum(., na.rm = T)))


summy |>
  ggplot() +
  geom_line(aes(x = time, y = Median_y1, colour = 'Nowcasting', 
                linetype = "nowcasting")) +
  geom_ribbon(aes(x = time, ymin = LI_y1, ymax = LS_y1), alpha = 0.2, show.legend = F) + 
  geom_line(data = dados_by_week %>% filter(time > 185),
            aes(x = time, y = y1_obs, colour = 'Observed at the moment',
                linetype = 'Observed at the moment')) +
  geom_line(data = dados_by_week %>% filter(time > 185), 
            aes(time, y = y1, colour = 'True', linetype = "True")) +
  labs(y = expression(N["1,t"]), x = "t") +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90)) +
  scale_linetype_manual(values = c("dashed", "dotted", "solid")) +
  scale_color_manual(values = c("red", "black", "black")) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(linetype = c("dashed", "dotted", "solid")),
                              title = ""),
         linetype = "none") -> g1_simul



summy |>
  ggplot() +
  geom_line(aes(x = time, y = Median_y2, colour = 'Nowcasting', linetype = "nowcasting")) +
  geom_ribbon(aes(x = time, ymin = LI_y2, ymax = LS_y2), alpha = 0.2, show.legend = F) + 
  geom_line(data = dados_by_week %>% filter(time > 185),
            aes(x = time, y = y2_obs, colour = 'Observed at the moment',
                linetype = 'Observed at the moment')) +
  geom_line(data = dados_by_week %>% filter(time > 185), 
            aes(time, y = y2, colour = 'True', linetype = "True")) +
  labs(y = expression(N["2,t"]), x = "t") +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90)) +
  scale_linetype_manual(values = c("dashed", "dotted", "solid")) +
  scale_color_manual(values = c("red", "black", "black")) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(linetype = c("dashed", "dotted", "solid")),
                              title = ""),
         linetype = "none") -> g2_simul



ggpubr::ggarrange(g1_simul, g2_simul, common.legend = TRUE,
                  legend = "bottom")

# ggsave("nowcasts_simulado_versao2.pdf")


  # scale_color_manual(values = c('grey50', 'black'), name = '') +
  # #scale_x_date(date_breaks = '2 weeks', date_labels = '%V/%y', name = 'Date in Weeks')+
  # labs(x = '', y = 'Nº Cases') 


# Nowcating Univariado ----------------------------------------------------

pred_samples_uni <- sample_predictive(inla_output = teste_uni,
                                      joint_data = joint_data)

summy_uni <- bivar_nowcasting_summary(pred_samples_uni)

dados_by_week <- y_new |>
  group_by(time) |>
  summarise_at(vars(starts_with("y")), .funs = list(~sum(., na.rm = T)))


summy_uni |>
  ggplot() +
  geom_line(aes(x = time, y = Median_y1, colour = 'Nowcasting', linetype = "nowcasting")) +
  geom_ribbon(aes(x = time, ymin = LI_y1, ymax = LS_y1), alpha = 0.2, show.legend = F) + 
  geom_line(data = dados_by_week %>% filter(time > 185),
            aes(x = time, y = y1_obs, colour = 'Observed at the moment',
                linetype = 'Observed at the moment')) +
  geom_line(data = dados_by_week %>% filter(time > 185), 
            aes(time, y = y1, colour = 'True', linetype = "True")) +
  labs(y = expression(N["1,t"]), x = "t") +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90)) +
  scale_linetype_manual(values = c("dashed", "dotted", "solid")) +
  scale_color_manual(values = c("red", "black", "black")) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(linetype = c("dashed", "dotted", "solid")),
                              title = ""),
         linetype = "none") -> g1_simul



summy_uni |>
  ggplot() +
  geom_line(aes(x = time, y = Median_y2, colour = 'Nowcasting', linetype = "nowcasting")) +
  geom_ribbon(aes(x = time, ymin = LI_y2, ymax = LS_y2), alpha = 0.2, show.legend = F) + 
  geom_line(data = dados_by_week %>% filter(time > 185),
            aes(x = time, y = y2_obs, colour = 'Observed at the moment',
                linetype = 'Observed at the moment')) +
  geom_line(data = dados_by_week %>% filter(time > 185), 
            aes(time, y = y2, colour = 'True', linetype = "True")) +
  labs(y = expression(N["2,t"]), x = "t") +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90)) +
  scale_linetype_manual(values = c("dashed", "dotted", "solid")) +
  scale_color_manual(values = c("red", "black", "black")) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(linetype = c("dashed", "dotted", "solid")),
                              title = ""),
         linetype = "none") -> g2_simul



ggpubr::ggarrange(g1_simul, g2_simul, common.legend = TRUE,
                  legend = "bottom")

# ggsave("nowcasts_simulado_versao2_univariate.pdf")
  

# Nowcasting poisson ------------------------------------------------------


# pred_samples_poisson <- sample_predictive(inla_output = teste_poisson, 
#                                   joint_data = joint_data)

summy_poisson <- bivar_nowcasting_summary(pred_samples_poisson)

dados_by_week <- y_new |>
  group_by(time) |>
  summarise_at(vars(starts_with("y")), .funs = list(~sum(., na.rm = T)))


summy_poisson |>
  ggplot() +
  geom_line(aes(x = time, y = Median_y1, colour = 'Nowcasting', 
                linetype = "nowcasting")) +
  geom_ribbon(aes(x = time, ymin = LI_y1, ymax = LS_y1), alpha = 0.2, show.legend = F) + 
  geom_line(data = dados_by_week %>% filter(time > 185),
            aes(x = time, y = y1_obs, colour = 'Observed at the moment',
                linetype = 'Observed at the moment')) +
  geom_line(data = dados_by_week %>% filter(time > 185), 
            aes(time, y = y1, colour = 'True', linetype = "True")) +
  labs(y = expression(N["1,t"]), x = "t") +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90)) +
  scale_linetype_manual(values = c("dashed", "dotted", "solid")) +
  scale_color_manual(values = c("red", "black", "black")) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(linetype = c("dashed", "dotted", "solid")),
                              title = ""),
         linetype = "none") -> g1_simul



summy_poisson |>
  ggplot() +
  geom_line(aes(x = time, y = Median_y2, colour = 'Nowcasting', linetype = "nowcasting")) +
  geom_ribbon(aes(x = time, ymin = LI_y2, ymax = LS_y2), alpha = 0.2, show.legend = F) + 
  geom_line(data = dados_by_week %>% filter(time > 185),
            aes(x = time, y = y2_obs, colour = 'Observed at the moment',
                linetype = 'Observed at the moment')) +
  geom_line(data = dados_by_week %>% filter(time > 185), 
            aes(time, y = y2, colour = 'True', linetype = "True")) +
  labs(y = expression(N["2,t"]), x = "t") +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90)) +
  scale_linetype_manual(values = c("dashed", "dotted", "solid")) +
  scale_color_manual(values = c("red", "black", "black")) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(linetype = c("dashed", "dotted", "solid")),
                              title = ""),
         linetype = "none") -> g2_simul



ggpubr::ggarrange(g1_simul, g2_simul, common.legend = TRUE,
                  legend = "bottom")

# ggsave("nowcasts_simulado_versao2_poisson.pdf")



# Comparação dos modelos --------------------------------------------------

metrics_summary <- mget(ls(pattern = "summy")) |> 
  map(~{selection_criteria_summary(.x, dados_by_week)}) |> 
  reduce(rbind) |> 
  mutate(model = rep(c("Bivariate", "Poisson bivariate", "Univariate"), 
                     each = 15))


Sys.setlocale("LC_ALL", locale = "Portuguese")
#Sys.setlocale("LC_ALL", locale = "English")


mget(ls(pattern = "pred_samples")) |> 
  map(~scores_posterior(.x, dados_by_week)) |> 
  reduce(rbind) |> 
  cbind(metrics_summary)|>  
  group_by(model) |> 
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


  # select(model, MAPE_y1 = perc_error_y1_mean, MAPE_y2 = perc_error_y2_mean,
  #        coverage_95_y1_mean, coverage_95_y2_mean, coverage_95_both_mean,
  #        coverage_95_both_last,
  #        width_95_y1_mean, width_95_y2_mean,
  #        width_95_y1_last, width_95_y2_last,
  #        is_95_y1_mean,
  #        is_95_y2_mean,
  #        relative_width_95_mean, relative_is_95_mean,
  #        relative_width_95_last, relative_is_95_last,
  #        energy_score_mean, energy_score_last,
  #        variogram_score_mean, variogram_score_last
  # )  |> t() |> 
  # 
  # |> View()


# summy |> 
#     mutate(model = "Bivariate") |> 
#     rbind(summy_uni |> mutate(model = "Independent")) |> 
#     mutate(range_95_y1 = )

mget(ls(pattern = "pred_samples")) |> 
    map(~scores_posterior(.x, dados_by_week)) |> 
    reduce(rbind) |> 
    cbind(metrics_summary) |> 
    ggplot() +
    aes(x = time, y = width_95_y1, colour = model) +
    geom_line()
  
  
  
  
  mget(ls(pattern = "pred_samples")) |> 
    map(~scores_posterior(.x, dados_by_week)) |> 
    reduce(rbind) |> 
    cbind(metrics_summary)|>  
    group_by(model) |> 
    summarise_all(.funs = list(mean = mean,
                               last = last)) |> 
    mutate(perc_error_y1_mean = perc_error_y1_mean*100,
           perc_error_y2_mean = perc_error_y2_mean*100) |> 
    select(
      "mean Interval score (N_1)" = is_95_y1_mean,
      "mean Interval score (N_2)" = is_95_y2_mean,
      "last Interval score (N_1)" = is_95_y1_last,
      "last Interval score (N_2)" = is_95_y2_last
    ) |> 
    t() |> View()
  kableExtra::kbl(format = "latex", booktabs = T, digits = 3) |> 
    add_header_above(c(" " = 1,
                       "Bivariate (NB)" = 1,
                       "Bivariate (poisson)" = 1,
                       "Univariate" = 1)) |>
    add_header_above(c(" " = 1, "Model" = 3))
  
  
  
  



