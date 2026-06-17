# Для цієї частини встановіть: install.packages(c("FactoMineR", "factoextra"))
library(dplyr)
library(ggplot2)
library(mgcv)
library(FactoMineR)
library(factoextra)

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
cat("\nФінальна узагальнена адитивна модель\n")
gam_model_full <- bam(log_claim ~ s(logD_c) + s(logBM) + s(DrivAge) + s(VehAge) + s(VehPowerNum),
                      data = train_data, family = gaussian(), method = "fREML", discrete = TRUE)
pred_gam_full <- predict(gam_model_full, newdata = test_data)
mse_gam_full <- mean((test_data$log_claim - pred_gam_full)^2)
cat(sprintf("MSE для Фінальної GAM моделі: %.4f\n", mse_gam_full))

# ==============================================================================
cat("\nPCA\n")
df_pca <- train_data %>%
  select(log_claim, logD_c, logBM, DrivAge, VehAge, VehPowerNum, LargeClaim) %>%
  tidyr::drop_na()

res.pca <- PCA(df_pca %>% select(-LargeClaim), scale.unit = TRUE, ncp = 5, graph = FALSE)

g5_1 <- fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50),
                 title = "5.2 Scree Plot: Відсоток поясненої дисперсії", ylab = "Дисперсія (%)", xlab = "Компоненти") + theme_minimal()
print(g5_1)

# ==============================================================================
cat("\nКореляційне коло змінних\n")
g5_2 <- fviz_pca_var(res.pca, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE, title = "5.3 Кореляційне коло ризиків") + theme_minimal()
print(g5_2)

# ==============================================================================
set.seed(42)
pca_sample_idx <- sample(nrow(df_pca), min(5000, nrow(df_pca)))
df_pca_sample <- df_pca[pca_sample_idx, ]

scores <- as_tibble(res.pca$ind$coord[pca_sample_idx, 1:2]) %>%
  rename(Dim1 = Dim.1, Dim2 = Dim.2) %>% mutate(LargeClaim = factor(df_pca_sample$LargeClaim, labels = c("Звичайний", "Катастрофа")))

var_coords <- as_tibble(res.pca$var$coord[, 1:2], rownames = "variable") %>% rename(Dim1 = Dim.1, Dim2 = Dim.2)
arrow_scale <- 4 
var_plot <- var_coords %>% mutate(xend = Dim1 * arrow_scale, yend = Dim2 * arrow_scale)

g5_3 <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5, color = "gray50") + geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5, color = "gray50") +
  geom_point(data = filter(scores, LargeClaim == "Звичайний"), aes(x = Dim1, y = Dim2), color = "gray40", alpha = 0.15, size = 1) +
  geom_point(data = filter(scores, LargeClaim == "Катастрофа"), aes(x = Dim1, y = Dim2), color = "red", alpha = 0.7, size = 1.5) +
  geom_segment(data = var_plot, aes(x = 0, y = 0, xend = xend, yend = yend), arrow = arrow(length = unit(0.2, "cm")), color = "black", linewidth = 0.8) +
  geom_text(data = var_plot, aes(x = xend * 1.1, y = yend * 1.1, label = variable), color = "black", fontface = "bold", size = 4) +
  labs(title = "5.4 Biplot: Проєкція водіїв (червоні - катастрофи)") + theme_minimal()
print(g5_3)