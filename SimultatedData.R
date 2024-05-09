library(INLA)
library(ggplot2)
library(dplyr)
library(nowcaster)

max_t <- 200
D <- 15
y <- matrix(NA, max_t + D, D)



### Geração de dados para o modelo de Morrison et al. (2016)

max_t <- 200
z <- vector(length = max_t)

set.seed(214)

z[1] <- rnorm(1)
tau_z <- 3 -> tau_x
rho <- 0.4
for(t in 2:max_t){z[t] <- rho*z[t-1] + rnorm(1, 0, sd = sqrt(1/tau_z))}
x <- rnorm(n = max_t, mean = z, sd = sqrt(1/tau_x))


alpha1 <- 0.5
beta1 <- 2
alpha2 <- -0.5
beta2 <- 3

lambda_1 <- exp(alpha1 + beta1*x)
lambda_2 <- exp(alpha2 + beta2*x)

dim_y <- 2
y <- matrix(nrow = max_t, ncol = dim_y)
y[,1] <- rpois(n = max_t, lambda = lambda_1)
y[,2] <- rpois(n = max_t, lambda = lambda_2)

ts.plot(y[,1])
lines(y[,2], col = 2)
cor(y)



### Tentando gerar dados adequados para Bastos et al. (2019)
library(truncnorm)

D <- 15 
max_t <- 200

beta <- vector(length = max_t)
gamma_d <- vector(length = D + 1)
alpha = 0

set.seed(456)

beta[1] <- rnorm(1, mean = 5)
gamma_d[1] <- rnorm(1)
for(t in 2:max_t){beta[t] <- rnorm(1, mean = beta[t-1], sd = 0.5)}
for(d in 2:(D + 1)){
  gamma_d[d] <- rtruncnorm(n = 1, mean = gamma_d[d-1], sd = 0.5, 
                           a = -Inf, b = gamma_d[d-1])
}

ts.plot(beta)
ts.plot(gamma_d)

y <- t(sapply(beta, 
            function(beta_t){
              rnbinom(n = D+1, size = 5, mu = exp(alpha + beta_t + gamma_d))
              }))


y_df <- as.data.frame(y)
colnames(y_df) <- paste0("D_", 0:15)

y_df %>% 
  mutate(time = 1:max_t) %>% 
  tidyr::pivot_longer(!time, names_to = "delay", values_to = "count") %>% 
  mutate(delay = readr::parse_number(delay),
         report_time = time + delay) -> y_df


y_new <- y_df %>% 
  mutate(Y = ifelse(report_time > 200, NA, count))




model <- Y ~ 1 +
  f(time,
    model = "rw2",
    hyper = list("prec" = list(prior = "loggamma",
                               param = c(0.001, 0.001))
    )) +
  f(delay, model = "rw1",
    hyper = list("prec" = list(prior = "loggamma",
                               param = c(0.001, 0.001)))
  )

## Age-Delay effects
## f(delay, model = "rw1",
## hyper = list("prec" = list(prior = "loggamma", param = c(0.001, 0.001))),
## group = fx_etaria.num, control.group = list(model = "iid"))

# if(zeroinflated){
#   ## Running the Zero inflated Negative Binomial model in INLA
#   output0 <- inla(model, family = "zeroinflatednbinomial0", data = dados.age,
#                   control.predictor = list(link = 1, compute = T),
#                   control.compute = list( config = T, waic=F, dic=F),
#                   control.family = list(
#                     hyper = list("theta1" = list(prior = "loggamma",
#                                                  param = c(0.001, 0.001))
#                                  ## "theta2" = list(prior = "gaussian",
#                                  ##                 param = c(-1, .1))
#                     )
#                   )
#   )
# }else{
## Running the Negative Binomial model in INLA
output0 <- INLA::inla(model, family = "nbinomial",
                      data = y_new,
                      control.predictor = list(link = 1, compute = T),
                      control.compute = list(config = T, waic=F, dic=F),
                      control.family = list(
                        hyper = list("theta" = list(prior = "loggamma",
                                                    param = c(0.001, 0.001))
                        )
                      )
)

index.missing <- which(is.na(y_new$Y))


## Step 1: Sampling from the approximate posterior distribution using INLA
srag.samples0.list <- INLA::inla.posterior.sample(n = 1000, output0)

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




## summary
ts.plot(output0$summary.random$delay$`0.5quant`, ylim = c(-5,2.5))
lines(gamma_d, col = 2)


ts.plot(output0$summary.random$time$`0.5quant`, ylim = c(-5,10))
lines(beta, col = 2)





## Nowcasting estimate
sample.now <- nowcasting_no_age(dados.age = dados.inla)

## Summary on the posteriors of nowcasting
now_summary<-nowcaster::nowcasting.summary(srag.pred.0,
                                age = F)
l<-1









#### exemplo pacote
library(nowcaster)
nowcasting_bh_no_age <- nowcasting_inla(dataset = sragBH, 
                                        date_onset = DT_SIN_PRI, 
                                        date_report = DT_DIGITA, 
                                        data.by.week = T)
head(nowcasting_bh_no_age$total)
















