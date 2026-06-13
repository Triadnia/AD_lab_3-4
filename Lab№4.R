# ЛАБОРАТОРНА РОБОТА №4. НЕПАРАМЕТРИЧНА РЕГРЕСІЯ ТА PCA
# ДОСЛІДНИЦЬКЕ ПИТАННЯ: Чи впливає вік водія (DrivAge) на вартість 
# страхового збитку (Total_Claim_Amount)?

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  cat("Робочу директорію встановлено на папку зі скриптом:\n", getwd(), "\n\n")
} else {
  cat("Збереження відбудеться у поточну робочу директорію:\n", getwd(), "\n\n")
}

# 1. БІБЛІОТЕКИ
library(tidyverse)
library(factoextra)
library(mgcv)
library(np)

freq_data <- read_csv("C:/Users/Acer/Desktop/АД/freMTPL2freq.csv", show_col_types = FALSE)
sev_data  <- read_csv("C:/Users/Acer/Desktop/АД/freMTPL2sev.csv",  show_col_types = FALSE)

sev_agg <- sev_data %>%
  distinct() %>%
  group_by(IDpol) %>%
  summarise(Total_Claim_Amount = sum(ClaimAmount), .groups = "drop")

full_data <- freq_data %>%
  left_join(sev_agg, by = "IDpol") %>%
  mutate(
    Total_Claim_Amount = replace_na(Total_Claim_Amount, 0),
    VehBrand  = as.factor(VehBrand),
    VehGas    = as.factor(VehGas)
  ) %>%
  filter(!(ClaimNb == 0 & Total_Claim_Amount > 0)) %>%
  filter(!(ClaimNb > 0  & Total_Claim_Amount == 0))

# Відсікання топ-0.5% катастрофічних виплат
threshold_995 <- quantile(full_data$Total_Claim_Amount[full_data$Total_Claim_Amount > 0], 0.995)

df_claims <- full_data %>%
  filter(Total_Claim_Amount > 0, Total_Claim_Amount <= threshold_995) %>%
  mutate(
    logY = log(Total_Claim_Amount),
    logBM = log(BonusMalus),
    logDensity = log(Density + 1),
    Diesel = as.integer(VehGas == "Diesel")
  ) %>% drop_na()

cat(sprintf("N після підготовки: %d спостережень\n", nrow(df_claims)))

# 3. АНАЛІЗ ГОЛОВНИХ КОМПОНЕНТ (PCA)
cat("\n=== Виконання PCA ===\n")
df_pca <- df_claims %>%
  select(DrivAge, logBM, VehAge, VehPower, logDensity)

pca_result <- prcomp(df_pca, center = TRUE, scale. = TRUE)

# Графік 1: Власні числа (Scree Plot)
p_scree <- fviz_eig(pca_result, addlabels = TRUE, ylim = c(0, 50),
                    title = "Scree Plot (Графік власних чисел)",
                    xlab = "Головні компоненти",
                    ylab = "Відсоток поясненої дисперсії")
print(p_scree)
ggsave("plot_1_pca_scree.png", plot = p_scree, width = 8, height = 6, dpi = 300)

# Графік 2: Суміщений біграфік (Максимальна чіткість)
p_biplot <- fviz_pca_biplot(pca_result,
                            geom.ind = "point", 
                            pointshape = 16, 
                            pointsize = 0.8,              
                            col.ind = "gray60",
                            alpha.ind = 0.15,
                            
                            col.var = "#D32F2F",
                            arrowsize = 1.2,
                            labelsize = 5.5,
                            font.var = c(14, "bold"),
                            repel = TRUE,
                            
                            title = "PCA Biplot: Проєкція змінних та об'єктів",
                            ggtheme = theme_minimal(base_size = 14)) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    panel.grid.major = element_line(color = "gray90") # Легка сітка на тлі
  )

print(p_biplot)
ggsave("plot_2_pca_biplot.png", plot = p_biplot, width = 10, height = 8, dpi = 300)


# 4. ВІЗУАЛЬНЕ ПОРІВНЯННЯ: OLS (Парабола) vs LOESS (Непараметрична)
cat("\n=== Візуальне порівняння OLS та LOESS ===\n")
df_plot <- df_claims %>%
  group_by(DrivAge) %>%
  summarise(mean_logY = mean(logY), n = n(), .groups = "drop") %>%
  filter(n > 10)

p_compare <- ggplot(df_plot, aes(x = DrivAge, y = mean_logY)) +
  geom_point(aes(size = n), alpha = 0.5, color = "gray50") +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE, 
              aes(color = "OLS (Парабола)"), linewidth = 1.2) +
  geom_smooth(method = "loess", span = 0.5, se = TRUE, 
              aes(color = "LOESS"), fill = "#1F77B4", alpha = 0.2, linewidth = 1.2) +
  scale_color_manual(name = "Метод", 
                     values = c("OLS (Парабола)" = "#2CA02C", "LOESS" = "#1F77B4")) +
  labs(title = "Непараметрична регресія: Форма впливу DrivAge на logY",
       subtitle = "Порівняння жорсткої OLS-параболи з гнучким локальним згладжуванням (LOESS)",
       x = "Вік водія (DrivAge)",
       y = "Середній логарифм збитку") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

print(p_compare)
ggsave("plot_3_ols_vs_loess.png", plot = p_compare, width = 10, height = 6, dpi = 300)


# 5. ЯДРОВА РЕГРЕСІЯ ТА ЧАСТКОВО ЛІНІЙНА МОДЕЛЬ (ПАКЕТ `np`)
set.seed(42)
df_sample <- df_claims %>% sample_n(2000)

cat("\n=== Оцінка Ядрової регресії (Надараї-Вотсона та Локально-лінійна) ===\n")
bw_nw <- npregbw(formula = logY ~ DrivAge, data = df_sample, regtype = "lc")
model_nw <- npreg(bw_nw)

bw_ll <- npregbw(formula = logY ~ DrivAge, data = df_sample, regtype = "ll")
model_ll <- npreg(bw_ll)

graphics.off()

png("plot_4_NW.png", width = 6, height = 5, units = "in", res = 300)
plot(model_nw, plot.errors.method = "bootstrap", plot.errors.style = "band", 
     main = "Надараї-Вотсона", ylab = "logY", xlab = "DrivAge", col="darkblue")
dev.off()

png("plot_5_LL.png", width = 6, height = 5, units = "in", res = 300)
plot(model_ll, plot.errors.method = "bootstrap", plot.errors.style = "band", 
     main = "Локальна лінійна", ylab = "logY", xlab = "DrivAge", col="darkblue")
dev.off()
cat("Графіки ядрової регресії збережено як plot_4_NW.png та plot_5_LL.png\n")

cat("\n=== Оцінка PLM (Partially Linear Model) через пакет np ===\n")
bw_plm <- npplregbw(formula = logY ~ logBM + VehPower + VehAge + Diesel | DrivAge, 
                    data = df_sample, regtype = "ll")

plm_model_np <- npplreg(bw_plm)
print(summary(plm_model_np))

cat("\nКоефіцієнти лінійної частини PLM:\n")
print(coef(plm_model_np))


# 6. УЗАГАЛЬНЕНА АДИТИВНА МОДЕЛЬ (GAM)
cat("\n=== Оцінка повної GAM моделі ===\n")
gam_model <- gam(logY ~ s(DrivAge) + s(logBM) + s(VehPower) + s(VehAge) + s(logDensity) + Diesel + VehBrand,
                 data = df_claims, method = "REML")

print(summary(gam_model))

cat("\nЯкість GAM моделі:\n")
cat("GAM AIC:", AIC(gam_model), "\n")


# 7. ІЗОЛЬОВАНИЙ ВПЛИВ ВІКУ ВОДІЯ НА ЗБИТОК (ПРОГНОЗ)
cat("\n=== Побудова ізольованого прогнозу GAM ===\n")
typ_BM       <- median(df_claims$BonusMalus)
typ_VehAge   <- median(df_claims$VehAge)
typ_VehPower <- median(df_claims$VehPower)
typ_Density  <- median(df_claims$Density)
typ_Diesel   <- 0
typ_Brand    <- "B1"

new_data <- data.frame(
  DrivAge    = seq(18, 90, by = 1),
  logBM      = log(typ_BM),
  VehAge     = typ_VehAge,
  VehPower   = typ_VehPower,
  logDensity = log(typ_Density + 1),
  Diesel     = typ_Diesel,
  VehBrand   = factor(typ_Brand, levels = levels(df_claims$VehBrand))
)

predictions <- predict(gam_model, newdata = new_data, type = "response", se.fit = TRUE)

new_data$fit   <- predictions$fit
new_data$upper <- predictions$fit + 1.96 * predictions$se.fit
new_data$lower <- predictions$fit - 1.96 * predictions$se.fit

p_final <- ggplot(new_data, aes(x = DrivAge, y = fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "purple", alpha = 0.2) +
  geom_line(color = "darkblue", linewidth = 1.5) +
  labs(title = "Ізольований вплив віку водія на тяжкість збитку (Прогноз GAM)",
       subtitle = sprintf("Усі інші фактори зафіксовані (BM=%.0f, Вік авто=%.0f, Бренд=%s)", 
                          typ_BM, typ_VehAge, typ_Brand),
       x = "Вік водія (DrivAge)",
       y = "Прогнозований логарифм збитку (logY)") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold"))

print(p_final)
ggsave("plot_6_gam_forecast.png", plot = p_final, width = 10, height = 6, dpi = 300)

cat("\n=== РОБОТУ ЗАВЕРШЕНО! Усі графіки успішно збережено в папку зі скриптом. ===\n")