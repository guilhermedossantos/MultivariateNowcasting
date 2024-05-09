## Geração e tentativa de ajuste de modelo para os dados Bivariados 

# Pacotes
library(INLA)
library(ggplot2)
library(dplyr)
library(nowcaster)
library(truncnorm)


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



set.seed(123)


beta[1, ] <- rnorm(2)
gamma_d[1,] <- rnorm(2, mean = 2)
delta[1] <- rnorm(1)
psi[1] <- rnorm(1)

tau_param <- 4
phi1 <- 0.4
phi52 <- 0.3

# Geração dos efeitos/covariaveis
for(t in 2:max_t){
  beta[t, ] <- rnorm(2, mean = beta[t-1, ], sd = sqrt(1/tau_param))
  # if(t == 2){
  #   beta[t, ] <- rnorm(2, mean = beta[t-1, ], sd = sqrt(1/tau_param))
  #   delta[t] <- rnorm(1, mean = delta[t-1], sd = sqrt(1/tau_param))
  # }else{
  #   beta[t, ] <- rnorm(2, mean = 2*beta[t-1, ] - beta[t-2,], sd = sqrt(1/tau_param))
  #   delta[t]  <- rnorm(1, mean = 2*delta[t-1] - delta[t-2], sd = sqrt(1/tau_param))
  #   }
  
  delta[t] <- rnorm(1, mean = delta[t-1], sd = sqrt(1/tau_param))
}

for(d in 2:(D + 1)){
  gamma_d[d, ] <- rtruncnorm(n = 1, mean = gamma_d[d-1, ], sd = sqrt(1/tau_param), 
                             a = -Inf, b = gamma_d[d-1, ])
  
  psi[d] <- rtruncnorm(n = 1, mean = psi[d - 1], sd = sqrt(1/tau_param), 
                       a = -Inf, b = psi[d - 1])
  
}



ts.plot(beta)
ts.plot(gamma_d)
ts.plot(delta)
ts.plot(psi)


input1 <- cbind(rep(alpha_1, max_t*(D+1)), 
                rep(gamma_d[,1], max_t),
                rep(beta[,1], each = D+1), 
                rep(delta, each = D+1),
                rep(psi, max_t))

input2 <- cbind(rep(alpha_2, max_t*(D+1)), 
                rep(gamma_d[,2], max_t),
                rep(beta[,2], each = D+1), 
                rep(delta, each = D+1),
                rep(psi, max_t))
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
  ts.plot()
  


y_new <- y_df %>% 
  mutate(y1_obs = ifelse(report_time > max_t, NA, y1),
         y2_obs = ifelse(report_time > max_t, NA, y2))





# Ajuste Bivariado --------------------------------------------------------

# n <- nrow(y)
# y_bvar <- matrix(NA, nrow = 2 * n, ncol = 2)
# y_bvar[1:n, 1] <- y_new$y1_obs
# y_bvar[n + 1:n, 2] <- y_new$y2_obs
# 
# I <- matrix(NA, nrow = 2 * n, ncol = 2)
# I[1:n, 1] <- 1
# I[n + 1:n, 2] <- 1
# 
# 
# time <- matrix(NA, nrow = 2 * n, ncol = 2)
# time[1:n, 1] <- y_df$time
# time[n + 1:n, 2] <- y_df$time
# 
# id_time = matrix(NA, nrow = 2 * n, ncol = 2)
# id_time[1:n, 1] <-1
# id_time[n + 1:n, 2] <- 1
# #id_delay = rep(1:2, each = n)
# 
# # id_time = rep(1:n, 2)
# # id_delay = rep(1:n, 2)
# 
# id_delay = matrix(NA, nrow = 2 * n, ncol = 2)
# id_delay[1:n, 1] <-1
# id_delay[n + 1:n, 2] <- 1
# 
# 
# delay <- matrix(NA, nrow = 2 * n, ncol = 2)
# delay[1:n, 1] <- y_df$delay
# delay[n + 1:n, 2] <- y_df$delay
# 
# joint_time <- matrix(NA, nrow = 2 * n, ncol = 2)
# #joint_time <- cbind(rep(y_df$time, 2), rep(y_df$time, 2))   
# joint_time <- time
# 
# joint_delay <- matrix(NA, nrow = 2 * n, ncol = 2)
# #joint_delay <- cbind(rep(y_df$delay, 2), rep(y_df$delay, 2))   
# joint_delay <- delay
# 
# joint_model <- y ~ -1 + I + 
#   # f(inla.group(time), model = "rw2", scale.model=TRUE,
#   #   hyper = list("prec" = list(prior = "loggamma",
#   #                              param = c(0.001, 0.001))
#   #                             )) +
#   f(id_time, time,
#     model = "rw2d", n = 2*n,
#     hyper = list("prec" = list(prior = "loggamma",
#                                param = c(0.001, 0.001))
#     )) +
#   f(id_delay, delay, model = "rw2d",n = 2*n,
#     hyper = list("prec" = list(prior = "loggamma",
#                                param = c(0.001, 0.001)))
#   ) #+
#   # f(joint_time, model = "rw2",
#   #   hyper = list("prec" = list(prior = "loggamma",
#   #                              param = c(0.001, 0.001))
#   #   )) +
#   # f(joint_delay, model = "rw1",
#   #   hyper = list("prec" = list(prior = "loggamma",
#   #                              param = c(0.001, 0.001))
#   #  )) 
# 
# 
# #f(idx.age, age, model = "iid2d", n = 2)
# 
# 
# 
# 
# 
# 
# # 
# # data_inla <- list(y = y_bvar, 
# #                   I = I, 
# #                   time = time, 
# #                   id_time = id_time, 
# #                   delay = delay,
# #                   id_delay = id_delay,
# #                   joint_time = joint_time,
# #                   joint_delay = joint_delay)
# 
# 
# 
# 
# 
# output0 <- INLA::inla(joint_model, family = rep("nbinomial", 2),
#                       data = data_inla,
#                       control.predictor = list(link = 1, compute = T),
#                       control.compute = list( config = T, waic=F, dic=F),
#                       control.family = list(
#                         list(hyper = list("theta" = list(prior = "loggamma", param = c(0.001, 0.001)))),
#                         list(hyper = list("theta" = list(prior = "loggamma", param = c(0.001, 0.0001))))
#                         )
#                       
# )



###################
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
                       idrep = rep(1:2, each = n))


random_covariate <- list(time1  = c(y_df$time, rep(NA, n)),
                         delay1 = c(y_df$delay, rep(NA, n)),
                         time2  = c(rep(NA, n), y_df$time),
                         delay2 = c(rep(NA, n), y_df$delay),
                         time = rep(y_df$time, 2),
                         delay = rep(y_df$delay, 2))




Y = y_bvar
joint_data <- cbind(Y, fixed_covariate, random_covariate)


# joint_model <- Y ~ -1 + alphai +
#   + f(id_time1, time1, model = "rw2") +
#   + f(id_delay1, delay1, model = "rw1") +
#   + f(id_time2, time2, replicate  = "id_time1", fixed = FALSE) +
#   + f(id_delay2, delay2, replicate  = "id_delay1", fixed = FALSE)

joint_model <- Y ~ -1 + alphai +
  + f(time,  model = "rw2", replicate = idrep,
      hyper = list(prec = list(param = c(0.001, 0.001)))) +
  + f(delay, model = "rw1", replicate = idrep) +
  + f(time1, model = "rw2") +
  + f(time2, copy  = "time1", fixed = FALSE) + 
  + f(delay1, model = "rw1") +
  + f(delay2, copy  = "delay1", fixed = FALSE) 
  # + f(id_time1, time1, model = "rw2") +
  # + f(id_time2, time2, model = "rw2") +
  # + f(id_delay1, delay1, , model = "rw1") +
  # + f(id_delay2, delay2, , model = "rw1") 


# jointmodel <- y ~ -1 + alphai + f(id_time, time,
#     model = "rw2d", n = 2*n,
#     hyper = list("prec" = list(prior = "loggamma",
#                                param = c(0.001, 0.001))
#     )) +
#   f(id_delay, delay, model = "rw2d",n = 2*n,
#     hyper = list("prec" = list(prior = "loggamma",
#                                param = c(0.001, 0.001)))
  # )

output0 <- INLA::inla(joint_model, family = rep("nbinomial", 2),
                      data = joint_data,
                      control.predictor = list(link = 1, compute = T),
                      control.compute = list( config = T, waic=F, dic=F),
                      control.family = list(
                        list(hyper = list("theta" = list(prior = "loggamma", param = c(0.001, 0.001)))),
                        list(hyper = list("theta" = list(prior = "loggamma", param = c(0.001, 0.0001))))
                      )
                      
)


ts.plot(output0$summary.random$time$`0.5quant`[1:200])
ts.plot(output0$summary.random$time$`0.5quant`[201:400])
ts.plot(output0$summary.random$time1$`0.5quant`)
ts.plot(output0$summary.random$time2$`0.5quant`)

ts.plot(output0$summary.random$delay$`0.5quant`[1:16])
ts.plot(output0$summary.random$delay$`0.5quant`[17:32])
ts.plot(output0$summary.random$delay1$`0.5quant`)
ts.plot(output0$summary.random$delay2$`0.5quant`)

ts.plot(output0$summary.random$delay1$`0.5quant` + output0$summary.random$delay$`0.5quant`[1:16])
lines(gamma_d[,1] + psi, col = 2)
ts.plot(output0$summary.random$delay2$`0.5quant` + output0$summary.random$delay$`0.5quant`[17:32])
lines(gamma_d[,2] + psi, col = 2)


ts.plot(output0$summary.random$time1$`0.5quant` + output0$summary.random$time$`0.5quant`[1:200], ylim = c(-6, 9))
lines(beta[,1] + delta, col = 2)

ts.plot(output0$summary.random$time2$`0.5quant` + output0$summary.random$time$`0.5quant`[201:400])
lines(beta[,2] + delta, col = 2)

ts.plot(c(LP1, LP2))
lines(output0$summary.linear.predictor$`0.5quant`, col = 2)


## Amostra da posteriori aproximada
samples <- INLA::inla.posterior.sample(n = 1000, output0)



index_missing <- which(is.na(joint_data$`1`) & is.na(joint_data$`2`))


vector_samples_bvar <- lapply(X = samples,
                          FUN = function(x, idx = index_missing, n = 3200){
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


gg.age <- function(x, dados.gg, idx, n = 3200){
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


dados_by_week |> 
  ggplot()+
  geom_line(data = dados_by_week, aes(time, y = y1, col = 'Observed'))+
  theme_bw()+
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90)) +
  scale_color_manual(values = c('grey50', 'black'), name = '')+
  labs(x = '', y = 'Nº Cases')


total.summy |> 
  ggplot(aes(x = time, y = Median_y1, col = 'Nowcasting')) +
  geom_line(data = dados_by_week %>% filter(time > 185), aes(time, y = y1, col = 'Observed')) +
  geom_ribbon(aes(ymin = LI_y1, ymax = LS_y1, col = NA), alpha = 0.2, show.legend = F) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90)) +
  scale_color_manual(values = c('grey50', 'black'), name = '')+
  #scale_x_date(date_breaks = '2 weeks', date_labels = '%V/%y', name = 'Date in Weeks')+
  labs(x = '', y = 'Nº Cases')


total.summy |> 
  ggplot(aes(x = time, y = Median_y2, col = 'Nowcasting')) +
  geom_line(data = dados_by_week %>% filter(time > 185), aes(time, y = y2, col = 'Observed')) +
  geom_ribbon(aes(ymin = LI_y2, ymax = LS_y2, col = NA), alpha = 0.2, show.legend = F) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90)) +
  # scale_color_manual(values = c('grey50', 'black'), name = '')+
  #scale_x_date(date_breaks = '2 weeks', date_labels = '%V/%y', name = 'Date in Weeks')+
  labs(x = '', y = 'Nº Cases')



# Ajuste univariado? ------------------------------------------------------



model <- y1_obs ~ 1 +
  f(time,
    model = "rw1",
    hyper = list("prec" = list(prior = "loggamma",
                               param = c(0.001, 0.001))
    )) +
  f(delay, model = "rw1",
    hyper = list("prec" = list(prior = "loggamma",
                               param = c(0.001, 0.001)))
  )


output0 <- INLA::inla(model, family = "nbinomial",
                      data = y_new,
                      control.predictor = list(link = 1, compute = T),
                      control.compute = list( config = T, waic=F, dic=F),
                      control.family = list(
                        hyper = list("theta" = list(prior = "loggamma",
                                                    param = c(0.001, 0.001))
                        )
                      )
)


ts.plot(output0$summary.random$delay$`0.5quant`, ylim = c(-12,6))
lines(gamma_d[,1] + psi, col = 2)


ts.plot(output0$summary.random$time$`0.5quant`, ylim = c(-7,12))
lines(beta[,1] + delta, col = 2)



model <- y2_obs ~ 1 +
  f(time,
    model = "rw1",
    hyper = list("prec" = list(prior = "loggamma",
                               param = c(0.001, 0.001))
    )) +
  f(delay, model = "rw1",
    hyper = list("prec" = list(prior = "loggamma",
                               param = c(0.001, 0.001)))
  )


output0 <- INLA::inla(model, family = "nbinomial",
                      data = y_new,
                      control.predictor = list(link = 1, compute = T),
                      control.compute = list( config = T, waic=F, dic=F),
                      control.family = list(
                        hyper = list("theta" = list(prior = "loggamma",
                                                    param = c(0.001, 0.001))
                        )
                      )
)


ts.plot(output0$summary.random$delay$`0.5quant`, ylim = c(-12,6))
lines(gamma_d[,1] + psi, col = 2)


ts.plot(output0$summary.random$time$`0.5quant`, ylim = c(-7,12))
lines(beta[,1] + delta, col = 2)








index.missing <- which(is.na(y_new$y1_obs))
y_new$Y <- y_new$y2_obs


## Step 1: Sampling from the approximate posterior distribution using INLA
srag.samples0.list <- INLA::inla.posterior.sample(n = 1000, output0)

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

gg.age <- function(x, dados.gg, idx){
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
                            FUN = gg.age,
                            dados = y_new,
                            idx = index.missing)

srag.pred.0 <- dplyr::bind_rows(tibble.samples.0, .id = "sample")

#return(srag.pred.0)



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


output<- list()
output$total <- total.summy



dados_by_week <- y_new |>
  group_by(time) |>
  summarise(n = sum(Y, na.rm = T))


dados_by_week |>
  ggplot()+
  geom_line(data = dados_by_week, aes(time, y = n, col = 'Observed'))+
  theme_bw()+
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90)) +
  scale_color_manual(values = c('grey50', 'black'), name = '')+
  labs(x = '', y = 'Nº Cases')


output$total |>
  ggplot(aes(x = time, y = Median, col = 'Nowcasting')) +
  geom_line(data = dados_by_week %>% filter(time > 185), aes(time, y = n, col = 'Observed')) +
  geom_ribbon(aes(ymin = LI, ymax = LS, col = NA), alpha = 0.2, show.legend = F) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90)) +
  scale_color_manual(values = c('grey50', 'black'), name = '')+
  #scale_x_date(date_breaks = '2 weeks', date_labels = '%V/%y', name = 'Date in Weeks')+
  labs(x = '', y = 'Nº Cases')

# 
# 
# 
# 
# 
# 
# 
# 
