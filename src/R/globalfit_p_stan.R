noAv = TRUE

if (noAv) {
  datTemp <-
    read.csv(file = "C:\\Users\\owner\\Documents\\Dorsal-Synthetics-Analysis\\dat\\datNoAv.csv")
} else{
  datTemp <-
    read.csv(file = "C:\\Users\\owner\\Documents\\Dorsal-Synthetics-Analysis\\dat\\dat.csv")
}

paramsTemp <-
  read.csv(file = "C:\\Users\\owner\\Documents\\Dorsal-Synthetics-Analysis\\dat\\params.csv")


x = datTemp$x

Y = datTemp$Y

nx = length(x)

nSets = 7


kdub = 1E4

kdp0 = 2E3

rub = 2E3

# w_inits = c(1E5, 1E3, 1E7)
# 
# paramsTemp$p0[1] = w_inits[1]
# paramsTemp$lb[1] = w_inits[2]
# paramsTemp$ub[1] = w_inits[3]

paramsTemp$ub[2:8] = kdub

paramsTemp$p0[2:8] = kdp0

paramsTemp$ub[9] = rub

#this is the p/kp parameter. p_inits is [p0, lb, ub]
p_inits = c(1E-5, 1E-7, 1E-3);

paramsTemp <- rbind(paramsTemp, c("p", p_inits))
paramsTemp$p0 = as.numeric(paramsTemp$p0)
paramsTemp$lb = as.numeric(paramsTemp$lb)
paramsTemp$ub = as.numeric(paramsTemp$ub)

dsid = datTemp$dsid


dat = list(
  "N" = nx,
  "x" = x,
  "Y" = Y,
  "dsid" = dsid,
  "p0" = paramsTemp$p0,
  "lb" = paramsTemp$lb,
  "ub" = paramsTemp$ub,
  "K" = length(paramsTemp$p0)
)

#
# n = 1;
# dsid = rep(0, times=nx)
# dsid[1] = 1;
#
# for (i in 2:nx){
#   if (x[i] < x[i-1]){
#     n = n + 1;
#   }
#   dsid[i] = n;
# }

df <- data.frame(x = x, y = Y, dsid = dsid)

parNames = c("w", "KD1", "KD2",
             "KD3", "KD4", "KD5",
             "KD6", "KD7", "R","p", "sigma");

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
fit <-
  stan(file = "C:\\Users\\owner\\Documents\\Dorsal-Synthetics-Analysis\\src\\globalSimpleFluo.stan",
       model_name = "simplefluo",
       data = dat,
       iter = 1E4,
       control = list(adapt_delta = 0.99))
print(fit, parNames)


Y_mean <- extract(fit, "Y_mean")
Y_mean_cred <- apply(Y_mean$Y_mean, 2, quantile, c(0.05, 0.95))
Y_mean_mean <- apply(Y_mean$Y_mean, 2, mean)

Y_pred <- extract(fit, "Y_pred")
Y_pred_cred <- apply(Y_pred$Y_pred, 2, quantile, c(0.05, 0.95))
Y_pred_mean <- apply(Y_pred$Y_pred, 2, mean)


predicted_df <- data.frame(
  Y_mean_mean = Y_mean_mean,
  Y_pred_mean = Y_pred_mean,
  Y_mean_cred_lower = Y_mean_cred[1,],
  Y_mean_cred_upper = Y_mean_cred[2,],
  Y_pred_cred_lower = Y_pred_cred[1,],
  Y_pred_cred_upper = Y_pred_cred[2,],
  dsid = dsid,
  x = x
)
library(extrafont)
library(viridis)
library(hrbrthemes)
library(ggthemes)
#font_import()
#loadfonts(device = "win")
#windowsFonts("Arial" = windowsFont("Arial"))

dfnames = c()

namesList = c("1Dg", "1DgS", "1DgW", "1DgAW3", "1DgSVW2", "1DgVVW", "1DgVW")

for (j in 1:nx) {
  dfnames[j] = namesList[dsid[j]]
}

df$dfnames = dfnames

predicted_df$dfnames = dfnames


windows()

df %>%
  ggplot(aes(x = x, y = Y, group = dfnames)) +
  # geom_point() +
  theme(
    legend.position = "bottom",
    panel.spacing = unit(0, "lines"),
    strip.text.x = element_text(size = 8, family = "ArialMT"),
    plot.title = element_text(size = 13, family = "ArialMT"),
    text = element_text(family = "ArialMT")
  ) +
  xlab("[Dorsal] (au)") +
  ylab("max fluorescence (au)") +
  facet_wrap( ~ dfnames) +
  geom_line(data = predicted_df,
            aes(x = x, y = Y_mean_mean, group = dfnames)) +
  geom_line(
    data = predicted_df,
    aes(x = x, y = Y_pred_cred_lower, group = dfnames),
    linetype = "dashed",
    color = "red"
  ) +
  geom_line(
    data = predicted_df,
    aes(x = x, y = Y_pred_cred_upper, group = dfnames),
    linetype = "dashed",
    color = "red"
  ) +
  geom_ribbon(
    data = predicted_df,
    aes(ymin = Y_mean_cred_lower, ymax = Y_mean_cred_upper, group = dfnames),
    fill = "blue",
    alpha = 0.3
  )
#legend=c("observation", "prediction", "mean prediction",
#        "90% mean cred. interval", "90% pred. cred. interval")

ggsave('C:\\Users\\owner\\Dropbox\\DorsalSyntheticsDropbox\\manuscript\\ggplot.pdf')


library(bayesplot)
posteriormat <- as.matrix(fit)
posteriorarr <- as.array(fit)

windows()

mcmc_pairs(posteriorarr,
           pars = c("w", "R", "KD1", "p", "sigma"))




plot(
  dat$Y ~ dat$x,
  xlab = "x",
  ylab = "Y",
  main = "simple weak fraction active",
  col = "green",
  bg = "green"
)
lines(dat$x, Y_mean_mean)
points(dat$x, Y_pred_mean, pch = 19)
lines(dat$x, Y_mean_cred[1, ], col = 4)
lines(dat$x, Y_mean_cred[2, ], col = 4)
lines(dat$x, Y_pred_cred[1, ], col = 2)
lines(dat$x, Y_pred_cred[2, ], col = 2)
legend(
  x = "bottomright",
  bty = "n",
  lwd = 2,
  lty = c(NA, NA, 1, 1, 1),
  legend = c(
    "observation",
    "prediction",
    "mean prediction",
    "90% mean cred. interval",
    "90% pred. cred. interval"
  ),
  col = c(1, 1, 1, 4, 2),
  pch = c(1, 19, NA, NA, NA)
)
#
# library(bayesplot)
# library(tidyverse)
#
# fit%>%
#   mcmc_trace()
#
# fit %>%
#   rhat() %>%
#   mcmc_rhat() +
#   yaxis_text()

np <- nuts_params(fit)
color_scheme_set("red")
mcmc_intervals(posteriormat,
               pars = parNames)
windows()
# 
# (p <- mcmc_pairs(
#   posteriorarr,
#   pars = parNames
# ))




mcmc_scatter(posteriorarr, pars = c("KD7", "R"))
mcmc_scatter(posteriorarr, pars = c("w", "p"))


windows()
mcmc_trace(posteriorarr,
           pars = parNames)

# 
# xaxis_title(size = 13, family = "sans") +
#   xaxis_ticks(on = FALSE) +
#   yaxis_ticks(on = FALSE)
# plot_title <- ggtitle("Posterior distributions",
#                       "with medians and 80% intervals")
# windows()
# mcmc_parcoord(
#   posteriormat,
#   transformations = "log",
#   pars = parNames
# )


# 
# windows()
# color_scheme_set("darkgray")
# mcmc_parcoord(
#   posteriormat,
#   transformations = "log",
#   pars = c("w", "R"),
#   np = np
# )
# 
# mcmc_areas(
#   posterior,
#   pars = parNames,
#   prob = 0.8,
#   transformations = "log"
# ) + plot_title

# mcmc_areas(posterior,
#            prob = 0.8) + plot_title
