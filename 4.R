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

# ==============================================================================
cat("\nKernel regression\n")
bw_train_1d_4 <- train_data %>% slice_sample(n = min(5000, nrow(train_data)))

bw_nw_1d_4 <- npregbw(LargeClaim ~ logD_c, data = bw_train_1d_4, regtype = "lc")
bw_ll_1d_4 <- npregbw(LargeClaim ~ logD_c, data = bw_train_1d_4, regtype = "ll")
model_nw_1d_4 <- npreg(bw_nw_1d_4); model_ll_1d_4 <- npreg(bw_ll_1d_4)

pred_nw_1d_4 <- predict(model_nw_1d_4, newdata = test_data)
pred_ll_1d_4 <- predict(model_ll_1d_4, newdata = test_data)
mse_nw_1d_4 <- mean((test_data$LargeClaim - pred_nw_1d_4)^2)
mse_ll_1d_4 <- mean((test_data$LargeClaim - pred_ll_1d_4)^2)

x_grid_D <- seq(min(train_data$logD_c), max(train_data$logD_c), length.out = 200)
fit_nw_1d_4 <- predict(model_nw_1d_4, newdata = data.frame(logD_c = x_grid_D), se.fit = TRUE)
fit_ll_1d_4 <- predict(model_ll_1d_4, newdata = data.frame(logD_c = x_grid_D), se.fit = TRUE)

plot_1d_4 <- bind_rows(
  data.frame(logD_c = x_grid_D, fit = fit_nw_1d_4$fit, se = fit_nw_1d_4$se.fit, method = "Nadaraya-Watson"),
  data.frame(logD_c = x_grid_D, fit = fit_ll_1d_4$fit, se = fit_ll_1d_4$se.fit, method = "Local Linear")
) %>% mutate(lower = pmax(0, fit - 1.96 * se), upper = pmin(1, fit + 1.96 * se))

g4_1 <- ggplot(plot_1d_4, aes(x = logD_c, y = fit, color = method, fill = method)) +
  geom_hline(yintercept = mean(train_data$LargeClaim), linetype = "dashed", color = "gray50") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  labs(title = "4.1 Ймовірність катастрофічного збитку (>3000€)", x = "logD_c", y = "P(LargeClaim=1)") + theme_minimal()
print(g4_1)

# ==============================================================================
cat("\nДвовимірна ядрова регресія\n")
bw_train_2d_4 <- train_data %>% slice_sample(n = min(3000, nrow(train_data)))

bw_nw_2d_4 <- npregbw(LargeClaim ~ logD_c + VehPowerNum, data = bw_train_2d_4, regtype = "lc")
bw_ll_2d_4 <- npregbw(LargeClaim ~ logD_c + VehPowerNum, data = bw_train_2d_4, regtype = "ll")
model_nw_2d_4 <- npreg(bw_nw_2d_4); model_ll_2d_4 <- npreg(bw_ll_2d_4)

pred_nw_2d_4 <- predict(model_nw_2d_4, newdata = test_data)
pred_ll_2d_4 <- predict(model_ll_2d_4, newdata = test_data)
mse_nw_2d_4 <- mean((test_data$LargeClaim - pred_nw_2d_4)^2, na.rm=TRUE)
mse_ll_2d_4 <- mean((test_data$LargeClaim - pred_ll_2d_4)^2, na.rm=TRUE)

# ==============================================================================
cat("\nGAM - Логіт\n")
gam_model_4 <- bam(LargeClaim ~ s(logD_c) + s(VehPowerNum) + s(logBM) + s(DrivAge) + s(VehAge),
                   data = train_data, family = binomial(link = "logit"), method = "fREML", discrete = TRUE)

pred_gam_4 <- predict(gam_model_4, newdata = test_data, type = "response")
mse_gam_4 <- mean((test_data$LargeClaim - pred_gam_4)^2)
par(mfrow = c(2, 3)) 
plot(gam_model_4, shade = TRUE, shade.col = "lightpink", seWithMean = TRUE, scale = 0, main = "Частковий ефект (Логіт)")
par(mfrow = c(1, 1))

res_table_4 <- data.frame(
  Модель = c("1D Kernel NW", "1D Kernel LL", "2D Kernel NW", "2D Kernel LL", "GAM Logit"),
  Brier_Score_MSE = c(mse_nw_1d_4, mse_ll_1d_4, mse_nw_2d_4, mse_ll_2d_4, mse_gam_4)
)
print(res_table_4)