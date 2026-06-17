library(dplyr)
library(ggplot2)
library(np)
library(mgcv)

freq <- read.csv("B:/AD/freMTPL2freq.csv")
sev  <- read.csv("B:/AD/freMTPL2sev.csv")

sev_agg <- sev %>% group_by(IDpol) %>% summarise(TotalClaimAmount = sum(ClaimAmount))
df_full <- freq %>% left_join(sev_agg, by = "IDpol") %>%
  mutate(TotalClaimAmount = ifelse(is.na(TotalClaimAmount), 0, TotalClaimAmount))

df_crashes <- df_full %>%
  filter(TotalClaimAmount > 0) %>%
  mutate(
    CappedClaim = pmin(TotalClaimAmount, quantile(TotalClaimAmount, 0.995, na.rm = TRUE)),
    log_claim = log(CappedClaim),
    logD_c = log(Density) - mean(log(Density), na.rm = TRUE),
    logBM = log(BonusMalus),
    VehPowerNum = as.numeric(as.character(VehPower)),
    VehPower = as.factor(VehPower),
    VehBrand = as.factor(VehBrand),
    LargeClaim = as.numeric(TotalClaimAmount > 3000)
  ) %>%
  filter(is.finite(log_claim), is.finite(logD_c), is.finite(logBM), is.finite(DrivAge), is.finite(VehAge))

lim <- df_crashes %>%
  summarise(
    logDc_low = quantile(logD_c, 0.01), logDc_high = quantile(logD_c, 0.99),
    logBM_low = quantile(logBM, 0.01),  logBM_high = quantile(logBM, 0.99),
    Age_low   = quantile(DrivAge, 0.01),  Age_high   = quantile(DrivAge, 0.99)
  )

df_crashes <- df_crashes %>%
  filter(logD_c >= lim$logDc_low, logD_c <= lim$logDc_high,
         logBM >= lim$logBM_low, logBM <= lim$logBM_high,
         DrivAge >= lim$Age_low, DrivAge <= lim$Age_high)

set.seed(67)
train_idx <- sample(nrow(df_crashes), 0.8 * nrow(df_crashes))
train_data <- df_crashes[train_idx, ]
test_data <- df_crashes[-train_idx, ]

bw_train_1d <- train_data %>% slice_sample(n = min(5000, nrow(train_data)))
x_grid_D <- seq(min(train_data$logD_c), max(train_data$logD_c), length.out = 300)

bw_nw_1d <- npregbw(log_claim ~ logD_c, data = bw_train_1d, regtype = "lc")
bw_ll_1d <- npregbw(log_claim ~ logD_c, data = bw_train_1d, regtype = "ll")
model_nw_1d <- npreg(bw_nw_1d); model_ll_1d <- npreg(bw_ll_1d)

pred_nw_1d <- predict(model_nw_1d, newdata = test_data)
pred_ll_1d <- predict(model_ll_1d, newdata = test_data)
mse_nw_1d <- mean((test_data$log_claim - pred_nw_1d)^2)
mse_ll_1d <- mean((test_data$log_claim - pred_ll_1d)^2)
cat(sprintf("MSE на тестовій вибірці: NW = %.4f | LL = %.4f\n", mse_nw_1d, mse_ll_1d))

fit_nw_1d <- predict(model_nw_1d, newdata = data.frame(logD_c = x_grid_D), se.fit = TRUE)
fit_ll_1d <- predict(model_ll_1d, newdata = data.frame(logD_c = x_grid_D), se.fit = TRUE)

plot_1d <- bind_rows(
  tibble(logD_c = x_grid_D, method = "Nadaraya-Watson", fit = fit_nw_1d$fit, se = fit_nw_1d$se.fit),
  tibble(logD_c = x_grid_D, method = "Local linear", fit = fit_ll_1d$fit, se = fit_ll_1d$se.fit)
) %>% mutate(lower = fit - 1.96 * se, upper = fit + 1.96 * se)

g1 <- ggplot(plot_1d, aes(logD_c, fit, color = method, fill = method)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.9) +
  labs(title = "1.1 Kernel regression: log(claim) on logD_c", x = "logD_c", y = "log(claim)") +
  theme_minimal()
print(g1)

cat("\nKernel regression для logD_c та logBM\n")
bw_train_2d <- train_data %>% slice_sample(n = min(3000, nrow(train_data)))

bw_nw_2d <- npregbw(log_claim ~ logD_c + logBM, data = bw_train_2d, regtype = "lc")
bw_ll_2d <- npregbw(log_claim ~ logD_c + logBM, data = bw_train_2d, regtype = "ll")
model_nw_2d <- npreg(bw_nw_2d); model_ll_2d <- npreg(bw_ll_2d)

pred_nw_2d <- predict(model_nw_2d, newdata = test_data)
pred_ll_2d <- predict(model_ll_2d, newdata = test_data)
mse_nw_2d <- mean((test_data$log_claim - pred_nw_2d)^2, na.rm=TRUE)
mse_ll_2d <- mean((test_data$log_claim - pred_ll_2d)^2, na.rm=TRUE)

bm_levels <- quantile(train_data$logBM, c(0.25, 0.50, 0.75))
grid_2d <- expand.grid(logD_c = x_grid_D, logBM = as.numeric(bm_levels))

add_fit <- function(model, method_name) {
  p <- predict(model, newdata = grid_2d, se.fit = TRUE)
  grid_2d %>% mutate(method = method_name, fit = drop(p$fit), se = drop(p$se.fit))
}

plot_df_2d <- bind_rows(add_fit(model_nw_2d, "Nadaraya-Watson"), add_fit(model_ll_2d, "Local linear")) %>% 
  mutate(logBM_factor = factor(logBM, labels = c("logBM Q25", "logBM median", "logBM Q75")))

g2 <- ggplot(plot_df_2d, aes(logD_c, fit, color = logBM_factor, fill = logBM_factor)) +
  geom_ribbon(aes(ymin = fit - 1.96 * se, ymax = fit + 1.96 * se), alpha = 0.16, color = NA) +
  geom_line(linewidth = 0.9) + facet_wrap(~ method) +
  labs(title = "1.2 Kernel regression: log(claim) on logD_c and logBM", x = "logD_c", y = "log(claim)") +
  theme_minimal()
print(g2)

cat("\nPartially linear regression\n")
bw_train_pl <- train_data %>% slice_sample(n = min(2500, nrow(train_data)))
bw_pl <- npplregbw(log_claim ~ DrivAge + VehAge + VehPower | logD_c + logBM, data = bw_train_pl, regtype = "ll")
pl_model <- npplreg(bw_pl)
cat("Лінійні коефіцієнти:\n"); print(summary(pl_model))

pred_pl <- predict(pl_model, newdata = test_data)
mse_pl <- mean((test_data$log_claim - pred_pl)^2, na.rm=TRUE)
par(mfrow = c(1, 2))
plot(pl_model, plot.errors.method = "asymptotic", plot.errors.style = "band", 
     col = "darkblue", main = "PLR Непараметричний вплив")
par(mfrow = c(1, 1))

# ==============================================================================
cat("\nGAM\n")
gam_model <- bam(log_claim ~ s(logD_c) + s(logBM) + s(DrivAge) + s(VehAge) + VehPower,
                 data = train_data, family = gaussian(), method = "fREML", discrete = TRUE)
cat("Оцінка:\n"); print(summary(gam_model))

pred_gam <- predict(gam_model, newdata = test_data)
mse_gam <- mean((test_data$log_claim - pred_gam)^2)

par(mfrow = c(2, 2)) 
plot(gam_model, shade = TRUE, shade.col = "lightblue", seWithMean = TRUE, scale = 0, main = "Частковий ефект")
par(mfrow = c(1, 1))

res_table_1 <- data.frame(
  Модель = c("Kernel: m(logD_c) NW", "Kernel: m(logD_c) LL", "Kernel: 2D NW", "Kernel: 2D LL", "PLR", "GAM"),
  Test_MSE = c(mse_nw_1d, mse_ll_1d, mse_nw_2d, mse_ll_2d, mse_pl, mse_gam)
)
print(res_table_1)