library(dplyr)
library(ggplot2)
library(np)
library(mgcv)
freq <- read.csv("B:/AD/freMTPL2freq.csv")
sev  <- read.csv("B:/AD/freMTPL2sev.csv")

sev_agg <- sev %>% group_by(IDpol) %>% summarise(TotalClaimAmount = sum(ClaimAmount))
df_full <- freq %>% left_join(sev_agg, by = "IDpol") %>% mutate(TotalClaimAmount = ifelse(is.na(TotalClaimAmount), 0, TotalClaimAmount))

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

lim <- df_crashes %>% summarise(logDc_low = quantile(logD_c, 0.01), logDc_high = quantile(logD_c, 0.99), logBM_low = quantile(logBM, 0.01), logBM_high = quantile(logBM, 0.99), Age_low = quantile(DrivAge, 0.01), Age_high = quantile(DrivAge, 0.99))
df_crashes <- df_crashes %>% filter(logD_c >= lim$logDc_low, logD_c <= lim$logDc_high, logBM >= lim$logBM_low, logBM <= lim$logBM_high, DrivAge >= lim$Age_low, DrivAge <= lim$Age_high)

set.seed(67)
train_idx <- sample(nrow(df_crashes), 0.8 * nrow(df_crashes))
train_data <- df_crashes[train_idx, ]
test_data <- df_crashes[-train_idx, ]

cat("\nKernel regression для DrivAge\n")
bw_train_1d_2 <- train_data %>% slice_sample(n = min(5000, nrow(train_data)))
x_grid_Age <- seq(min(train_data$DrivAge), max(train_data$DrivAge), length.out = 300)

bw_nw_1d_2 <- npregbw(log_claim ~ DrivAge, data = bw_train_1d_2, regtype = "lc")
bw_ll_1d_2 <- npregbw(log_claim ~ DrivAge, data = bw_train_1d_2, regtype = "ll")
model_nw_1d_2 <- npreg(bw_nw_1d_2); model_ll_1d_2 <- npreg(bw_ll_1d_2)

pred_nw_1d_2 <- predict(model_nw_1d_2, newdata = test_data)
pred_ll_1d_2 <- predict(model_ll_1d_2, newdata = test_data)
mse_nw_1d_2 <- mean((test_data$log_claim - pred_nw_1d_2)^2)
mse_ll_1d_2 <- mean((test_data$log_claim - pred_ll_1d_2)^2)

fit_nw_1d_2 <- predict(model_nw_1d_2, newdata = data.frame(DrivAge = x_grid_Age), se.fit = TRUE)
fit_ll_1d_2 <- predict(model_ll_1d_2, newdata = data.frame(DrivAge = x_grid_Age), se.fit = TRUE)

plot_1d_2 <- bind_rows(
  tibble(DrivAge = x_grid_Age, method = "Nadaraya-Watson", fit = fit_nw_1d_2$fit, se = fit_nw_1d_2$se.fit),
  tibble(DrivAge = x_grid_Age, method = "Local linear", fit = fit_ll_1d_2$fit, se = fit_ll_1d_2$se.fit)
) %>% mutate(lower = fit - 1.96 * se, upper = fit + 1.96 * se)

g3 <- ggplot(plot_1d_2, aes(DrivAge, fit, color = method, fill = method)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.9) +
  labs(title = "2.1 Kernel regression: log(claim) on DrivAge", x = "DrivAge", y = "log(claim)") + theme_minimal()
print(g3)

# ==============================================================================
cat("\nKernel regression для DrivAge та logBM\n")
bw_train_2d_2 <- train_data %>% slice_sample(n = min(3000, nrow(train_data)))

bw_nw_2d_2 <- npregbw(log_claim ~ DrivAge + logBM, data = bw_train_2d_2, regtype = "lc")
bw_ll_2d_2 <- npregbw(log_claim ~ DrivAge + logBM, data = bw_train_2d_2, regtype = "ll")
model_nw_2d_2 <- npreg(bw_nw_2d_2); model_ll_2d_2 <- npreg(bw_ll_2d_2)

pred_nw_2d_2 <- predict(model_nw_2d_2, newdata = test_data)
pred_ll_2d_2 <- predict(model_ll_2d_2, newdata = test_data)
mse_nw_2d_2 <- mean((test_data$log_claim - pred_nw_2d_2)^2, na.rm=TRUE)
mse_ll_2d_2 <- mean((test_data$log_claim - pred_ll_2d_2)^2, na.rm=TRUE)

bm_levels <- quantile(train_data$logBM, c(0.25, 0.50, 0.75))
grid_2d_2 <- expand.grid(DrivAge = x_grid_Age, logBM = as.numeric(bm_levels))

add_fit_2 <- function(model, method_name) {
  p <- predict(model, newdata = grid_2d_2, se.fit = TRUE)
  grid_2d_2 %>% mutate(method = method_name, fit = drop(p$fit), se = drop(p$se.fit))
}
plot_df_2d_2 <- bind_rows(add_fit_2(model_nw_2d_2, "Nadaraya-Watson"), add_fit_2(model_ll_2d_2, "Local linear")) %>% 
  mutate(logBM_factor = factor(logBM, labels = c("logBM Q25", "logBM median", "logBM Q75")))

g4 <- ggplot(plot_df_2d_2, aes(DrivAge, fit, color = logBM_factor, fill = logBM_factor)) +
  geom_ribbon(aes(ymin = fit - 1.96 * se, ymax = fit + 1.96 * se), alpha = 0.16, color = NA) +
  geom_line(linewidth = 0.9) + facet_wrap(~ method) +
  labs(title = "Kernel regression: log(claim) on DrivAge and logBM", x = "DrivAge", y = "log(claim)") + theme_minimal()
print(g4)

# ==============================================================================
cat("\nPartially linear regression\n")
bw_train_pl_2 <- train_data %>% slice_sample(n = min(2500, nrow(train_data)))
bw_pl_2 <- npplregbw(log_claim ~ VehAge + VehPower + logD_c | DrivAge + logBM, data = bw_train_pl_2, regtype = "ll")
pl_model_2 <- npplreg(bw_pl_2)
cat("Лінійні коефіцієнти:\n"); print(summary(pl_model_2))

pred_pl_2 <- predict(pl_model_2, newdata = test_data)
mse_pl_2 <- mean((test_data$log_claim - pred_pl_2)^2, na.rm=TRUE)

par(mfrow = c(1, 2))
plot(pl_model_2, plot.errors.method = "asymptotic", plot.errors.style = "band", 
     col = "darkorange", main = "PLR Непараметричний вплив")
par(mfrow = c(1, 1))

# ==============================================================================
cat("\nGAM\n")
gam_model_2 <- bam(log_claim ~ s(DrivAge) + s(logBM) + s(logD_c) + s(VehAge) + VehPower,
                   data = train_data, family = gaussian(), method = "fREML", discrete = TRUE)
pred_gam_2 <- predict(gam_model_2, newdata = test_data)
mse_gam_2 <- mean((test_data$log_claim - pred_gam_2)^2)

par(mfrow = c(2, 2)) 
plot(gam_model_2, shade = TRUE, shade.col = "bisque", seWithMean = TRUE, scale = 0, main = "Частковий ефект")
par(mfrow = c(1, 1))

res_table_2 <- data.frame(
  Модель = c("Kernel: m(DrivAge) NW", "Kernel: m(DrivAge) LL", "Kernel: 2D NW", "Kernel: 2D LL", "PLR", "GAM"),
  Test_MSE = c(mse_nw_1d_2, mse_ll_1d_2, mse_nw_2d_2, mse_ll_2d_2, mse_pl_2, mse_gam_2)
)
print(res_table_2)