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

lim <- df_crashes %>% summarise(
  logDc_low = quantile(logD_c, 0.01), logDc_high = quantile(logD_c, 0.99), 
  logBM_low = quantile(logBM, 0.01), logBM_high = quantile(logBM, 0.99), 
  Age_low = quantile(DrivAge, 0.01), Age_high = quantile(DrivAge, 0.99),
  VehAge_high = quantile(VehAge, 0.99)
)

df_crashes <- df_crashes %>% filter(
  logD_c >= lim$logDc_low, logD_c <= lim$logDc_high, 
  logBM >= lim$logBM_low, logBM <= lim$logBM_high, 
  DrivAge >= lim$Age_low, DrivAge <= lim$Age_high,
  VehAge <= lim$VehAge_high # Застосовуємо фільтр
)

set.seed(67)
train_idx <- sample(nrow(df_crashes), 0.8 * nrow(df_crashes))
train_data <- df_crashes[train_idx, ]
test_data <- df_crashes[-train_idx, ]

cat("\nKernel regression для VehAge\n")
bw_train_1d_3 <- train_data %>% slice_sample(n = min(5000, nrow(train_data)))

x_grid_VehAge <- seq(min(train_data$VehAge), max(train_data$VehAge), length.out = 200)

bw_nw_1d_3 <- npregbw(log_claim ~ VehAge, data = bw_train_1d_3, regtype = "lc")
bw_ll_1d_3 <- npregbw(log_claim ~ VehAge, data = bw_train_1d_3, regtype = "ll")
model_nw_1d_3 <- npreg(bw_nw_1d_3); model_ll_1d_3 <- npreg(bw_ll_1d_3)

pred_nw_1d_3 <- predict(model_nw_1d_3, newdata = test_data)
pred_ll_1d_3 <- predict(model_ll_1d_3, newdata = test_data)
mse_nw_1d_3 <- mean((test_data$log_claim - pred_nw_1d_3)^2)
mse_ll_1d_3 <- mean((test_data$log_claim - pred_ll_1d_3)^2)

fit_nw_1d_3 <- predict(model_nw_1d_3, newdata = data.frame(VehAge = x_grid_VehAge), se.fit = TRUE)
fit_ll_1d_3 <- predict(model_ll_1d_3, newdata = data.frame(VehAge = x_grid_VehAge), se.fit = TRUE)

plot_1d_3 <- bind_rows(
  data.frame(VehAge = x_grid_VehAge, fit = fit_nw_1d_3$fit, se = fit_nw_1d_3$se.fit, method = "Nadaraya-Watson"),
  data.frame(VehAge = x_grid_VehAge, fit = fit_ll_1d_3$fit, se = fit_ll_1d_3$se.fit, method = "Local Linear")
) %>% mutate(lower = fit - 1.96 * se, upper = fit + 1.96 * se)

g3_1 <- ggplot(plot_1d_3, aes(x = VehAge, y = fit, color = method, fill = method)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  labs(title = "3.1 Kernel regression: Вплив віку авто на тяжкість ДТП", x = "Вік автомобіля", y = "log(Збитку)") + theme_minimal()
print(g3_1)

# ==============================================================================
cat("\nKernel для VehAge та VehPowerNum\n")
bw_train_2d_3 <- train_data %>% slice_sample(n = min(3000, nrow(train_data)))

bw_nw_2d_3 <- npregbw(log_claim ~ VehAge + VehPowerNum, data = bw_train_2d_3, regtype = "lc")
bw_ll_2d_3 <- npregbw(log_claim ~ VehAge + VehPowerNum, data = bw_train_2d_3, regtype = "ll")
model_nw_2d_3 <- npreg(bw_nw_2d_3); model_ll_2d_3 <- npreg(bw_ll_2d_3)

pred_nw_2d_3 <- predict(model_nw_2d_3, newdata = test_data)
pred_ll_2d_3 <- predict(model_ll_2d_3, newdata = test_data)
mse_nw_2d_3 <- mean((test_data$log_claim - pred_nw_2d_3)^2, na.rm=TRUE)
mse_ll_2d_3 <- mean((test_data$log_claim - pred_ll_2d_3)^2, na.rm=TRUE)

power_levels <- quantile(train_data$VehPowerNum, c(0.25, 0.50, 0.75))
grid_2d_3 <- expand.grid(VehAge = x_grid_VehAge, VehPowerNum = as.numeric(power_levels))
fit_2d_3 <- predict(model_ll_2d_3, newdata = grid_2d_3, se.fit = TRUE)
plot_2d_3 <- grid_2d_3 %>% mutate(fit = fit_2d_3$fit, se = fit_2d_3$se.fit, Power_level = factor(VehPowerNum, labels = c("Малолітражки", "Середні", "Потужні")))

g3_2 <- ggplot(plot_2d_3, aes(x = VehAge, y = fit, color = Power_level, fill = Power_level)) +
  geom_ribbon(aes(ymin = fit - 1.96*se, ymax = fit + 1.96*se), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  labs(title = "3.2 Kernel 2D (LL): Старіння авто при різній потужності", x = "Вік авто", y = "log(Збитку)") + theme_minimal()
print(g3_2)

# ==============================================================================
cat("\nЧастково лінійна регресія\n")
bw_train_plr_3 <- train_data %>% slice_sample(n = min(2500, nrow(train_data)))
bw_pl_3 <- npplregbw(formula = log_claim ~ logD_c + logBM + DrivAge | VehAge + VehPowerNum, data = bw_train_plr_3, regtype = "ll")
pl_model_3 <- npplreg(bw_pl_3)
pred_pl_3 <- predict(pl_model_3, newdata = test_data)
mse_pl_3 <- mean((test_data$log_claim - pred_pl_3)^2, na.rm=TRUE)

par(mfrow = c(1, 2))
plot(pl_model_3, plot.errors.method = "asymptotic", plot.errors.style = "band", 
     col = "purple", main = "PLR Непараметричний вплив")
par(mfrow = c(1, 1))

# ==============================================================================
cat("\nУзагальнена адитивна модель\n")
gam_model_3 <- bam(log_claim ~ s(VehAge) + s(VehPowerNum) + s(logD_c) + s(logBM) + s(DrivAge), data = train_data, family = gaussian(), method = "fREML", discrete = TRUE)
pred_gam_3 <- predict(gam_model_3, newdata = test_data)
mse_gam_3 <- mean((test_data$log_claim - pred_gam_3)^2)

par(mfrow = c(2, 3))             
plot(gam_model_3, shade = TRUE, shade.col = "thistle", seWithMean = TRUE, scale = 0, main = "Частковий ефект")
par(mfrow = c(1, 1))

res_table_3 <- data.frame(
  Модель = c("1D Kernel NW", "1D Kernel LL", "2D Kernel NW", "2D Kernel LL", "PLR", "GAM"),
  MSE_Test = c(mse_nw_1d_3, mse_ll_1d_3, mse_nw_2d_3, mse_ll_2d_3, mse_pl_3, mse_gam_3)
)
print(res_table_3)