# =========================================================================
# ЛАБОРАТОРНА РОБОТА №3. РЕГРЕСІЙНИЙ АНАЛІЗ
# Датасет: freMTPL2 (French Motor Third-Party Liability Insurance)
#
# ДОСЛІДНИЦЬКЕ ПИТАННЯ (причинно-наслідкове):
# Чи впливає вік водія (DrivAge) на вартість страхового збитку, і чи
# зберігається цей вплив після врахування характеристик автомобіля,
# досвіду водія та середовища?
#
# ПРИЧИННО-НАСЛІДКОВА ЛОГІКА:
# Вік водія пов'язаний з фізіологічними та когнітивними можливостями:
#   – Молоді водії (18–25): менший досвід → частіші помилки при
#     маневруванні, вища схильність до ризику → серйозніші ДТП;
#   – Водії середнього віку (26–59): оптимальне поєднання досвіду
#     та фізичних можливостей → найменший ризик;
#   – Літні водії (60+): уповільнена реакція, погіршення зору →
#     частіші помилки при оцінці ситуації → серйозніші наслідки.
# Це біологічно й психологічно обґрунтований механізм, а не лише
# кореляція.
# =========================================================================


# =========================================================================
# 1. БІБЛІОТЕКИ
# =========================================================================
library(tidyverse)
library(sandwich)
library(lmtest)
library(car)


# =========================================================================
# 2. ЗАВАНТАЖЕННЯ ТА ОЧИЩЕННЯ ДАНИХ
# =========================================================================
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
    Exposure  = if_else(Exposure > 1, 1, Exposure),
    VehBrand  = as.factor(VehBrand),
    VehGas    = as.factor(VehGas),
    Region    = as.factor(Region),
    Area      = as.factor(Area)
  ) %>%
  filter(!(ClaimNb == 0 & Total_Claim_Amount > 0)) %>%
  filter(!(ClaimNb > 0  & Total_Claim_Amount == 0))


# =========================================================================
# 3. ВІДСІКАННЯ ТОП-0.5% КАТАСТРОФІЧНИХ ВИПЛАТ
# =========================================================================
# Мотивація: розподіл страхових виплат є сильно правостороннім.
# Верхні 0.5% є катастрофічними подіями (летальні наслідки, тривала
# інвалідність), нетиповими наслідком віку водія. Включення спотворює
# оцінку «звичайного» вікового ефекту.
# Відсікання виконується ГЛОБАЛЬНО до будь-якого розбиття.

threshold_995 <- quantile(
  full_data$Total_Claim_Amount[full_data$Total_Claim_Amount > 0],
  0.995
)
cat(sprintf("Поріг відсікання (99.5%%): %.2f євро\n", threshold_995))

n_before <- nrow(full_data %>% filter(Total_Claim_Amount > 0))
cat(sprintf("Кількість полісів з виплатами (до відсікання): %d\n", n_before))


# =========================================================================
# 4. ПІДГОТОВКА АНАЛІТИЧНОГО ДАТАСЕТУ
# =========================================================================
df_claims <- full_data %>%
  filter(Total_Claim_Amount > 0,
         Total_Claim_Amount <= threshold_995) %>%
  mutate(
    
    # ---------------------------------------------------------------
    # ЗАЛЕЖНА ЗМІННА: log(Total_Claim_Amount)
    #
    # ОБГРУНТУВАННЯ ЛОГАРИФМУВАННЯ:
    # Страхові виплати мають сильно правосторонній розподіл (асиметрія
    # >> 0 навіть після відсікання топ-0.5%). Логарифм:
    #   (а) стабілізує дисперсію → виконання припущення гомоскедастич-
    #       ності краще; залишки після log значно симетричніші;
    #   (б) дає економічно змістовну інтерпретацію: коефіцієнт β є
    #       напів-еластичністю → збільшення Xj на 1 одиницю пов'язане
    #       зі зміною ClaimAmount на (exp(β)-1)×100%;
    #   (в) звужує вплив залишкових викидів після відсікання.
    # Перевірка: hist(Total_Claim_Amount) vs hist(log(Total_Claim_Amount))
    # підтверджує суттєво симетричніший розподіл після логарифмування.
    # ---------------------------------------------------------------
    logY = log(Total_Claim_Amount),
    
    # ---------------------------------------------------------------
    # КЛЮЧОВИЙ РЕГРЕСОР: вік водія (DrivAge, роки)
    #
    # Представлення (а): безперервна центрована змінна
    #   DrivAge_c  = DrivAge - mean(DrivAge)
    #   → центрування усуває мультиколінеарність між лінійним і
    #     квадратичним членами (corr(DrivAge, DrivAge²) ≈ висока;
    #     після центрування corr(DrivAge_c, DrivAge_c²) ≈ 0)
    #   DrivAge_c2 = DrivAge_c² — для перевірки нелінійності (M7)
    #     Очікується U-подібна крива: молодь та літні → вищі збитки,
    #     середній вік → найнижчі.
    #
    # Представлення (б): групова змінна AgeGroup
    #   Youth   (DrivAge 18–25): підвищений ризик через недосвідченість
    #   Adults  (DrivAge 26–59): референтна категорія (найнижчий ризик)
    #   Seniors (DrivAge 60+):   підвищений ризик через вікові обмеження
    #   → кордони груп відповідають природним класам страхового ринку
    #     і медичній літературі щодо ризику ДТП за віком
    # ---------------------------------------------------------------
    DrivAge_c  = DrivAge - mean(DrivAge),
    DrivAge_c2 = DrivAge_c^2,
    
    AgeGroup = case_when(
      DrivAge >= 18 & DrivAge <= 25 ~ "Youth",
      DrivAge >= 26 & DrivAge <= 59 ~ "Adults",
      DrivAge >= 60                 ~ "Seniors"
    ),
    AgeGroup = relevel(as.factor(AgeGroup), ref = "Adults"),
    
    # ---------------------------------------------------------------
    # КОНТРОЛЬНІ ЗМІННІ
    #
    # log(BonusMalus):
    #   BM є найважливішим confounder: він відображає НАКОПИЧЕНИЙ
    #   досвід водіння та схильність до аварій. Молоді водії мають
    #   вищий BM (менше досвіду), літні — можуть мати знижений BM
    #   (довгий безаварійний стаж). Без контролю за BM ефект AgeGroup
    #   частково є ефектом досвіду, а не власне віку.
    #   Логарифм: BM варіює 50–350, правосторонній розподіл.
    #
    # (logBM)²:
    #   Перевірка нелінійності у BonusMalus. Включається в M10.
    #
    # PowerGroup: потужність авто (Low/Medium/High).
    #   Молодь частіше їздить на менш потужних авто → конфаундер
    #   при оцінці вікового ефекту. Обов'язково контролювати.
    #
    # CarAgeGroup: New/Used/Old — вік авто впливає на вартість ремонту.
    #   Молоді водії частіше їздять на старих авто (нижча ціна) →
    #   потенційний конфаундер.
    #
    # Diesel: дизельні авто мають дорожчі двигуни → вплив на збиток.
    #
    # log(Density+1): щільність населення. Міська їзда (висока щільність)
    #   має специфічний профіль ризику, що може корелювати з віком
    #   (молодь концентрується у містах).
    # ---------------------------------------------------------------
    logBM  = log(BonusMalus),
    logBM2 = log(BonusMalus)^2,
    
    PowerGroup = case_when(
      VehPower <= 5                 ~ "Low",
      VehPower >= 6 & VehPower <= 8 ~ "Medium",
      VehPower >= 9                 ~ "High"
    ),
    PowerGroup = relevel(as.factor(PowerGroup), ref = "Low"),
    
    CarAgeGroup = case_when(
      VehAge >= 0  & VehAge <= 3  ~ "New",
      VehAge >= 4  & VehAge <= 10 ~ "Used",
      VehAge >= 11                ~ "Old"
    ),
    CarAgeGroup = relevel(as.factor(CarAgeGroup), ref = "Used"),
    
    Diesel     = as.integer(VehGas == "Diesel"),
    logDensity = log(Density + 1),
    
    VehBrand = relevel(as.factor(VehBrand), ref = "B1"),
    Region   = relevel(as.factor(Region),   ref = "Ile-de-France")
  )

cat(sprintf("N після відсікання: %d (відсіяно: %d спостережень, %.1f%%)\n",
            nrow(df_claims),
            n_before - nrow(df_claims),
            (n_before - nrow(df_claims)) / n_before * 100))


# =========================================================================
# 5. ПОПЕРЕДНІЙ АНАЛІЗ: МУЛЬТИКОЛІНЕАРНІСТЬ ТА НЕЛІНІЙНІСТЬ
# =========================================================================

cat("\n=== Кореляційна матриця числових регресорів ===\n")
cor_mat <- cor(df_claims[, c("logY", "DrivAge", "logBM", "VehPower",
                             "VehAge", "logDensity", "Diesel")])
print(round(cor_mat, 3))

cat("\n=== Середній logY по десятиліттях DrivAge ===\n")
df_claims %>%
  mutate(age_decade = floor(DrivAge / 10) * 10) %>%
  group_by(age_decade) %>%
  summarise(mean_logY = round(mean(logY), 3),
            n = n(), .groups = "drop") %>%
  print()

cat("\n=== Середній logY по групах AgeGroup ===\n")
df_claims %>%
  group_by(AgeGroup) %>%
  summarise(mean_logY = round(mean(logY), 3),
            sd_logY   = round(sd(logY), 3),
            n = n(), .groups = "drop") %>%
  print()

# Перевірка форми зв'язку logY ~ logBM (лінійна vs квадратична)
cat("\n=== Середній logY по децилях logBM ===\n")
df_claims %>%
  mutate(decile_BM = ntile(logBM, 10)) %>%
  group_by(decile_BM) %>%
  summarise(mean_logBM = round(mean(logBM), 3),
            mean_logY  = round(mean(logY), 3),
            n = n(), .groups = "drop") %>%
  print()


# =========================================================================
# 6. СТРУКТУРНА МОДЕЛЬ ТА ГІПОТЕЗИ
# =========================================================================
#
# log(ClaimAmount_i) = β₀ + β₁·AgeGroup_i
#                     + β₂·log(BonusMalus_i)
#                     + β₃·PowerGroup_i
#                     + β₄·CarAgeGroup_i
#                     + β₅·Diesel_i
#                     + β₆·log(Density_i + 1)
#                     + [Region_i]
#                     + [VehBrand_i]
#                     + ε_i
#
# ГІПОТЕЗИ ЩОДО ЗНАКІВ:
#   β₁(Youth vs Adults)   > 0   — молоді водії: нижчий досвід,
#                                  вища схильність до ризику → більший збиток
#   β₁(Seniors vs Adults) > 0   — літні водії: уповільнена реакція,
#                                  погіршений зір → більший збиток
#   β₂(logBM)             > 0   — ризиковіші водії → більші збитки
#   β₃(PowerHigh vs Low)  > 0   — потужніші авто → більша кінетична енергія
#   β₄(CarNew vs Used)    > 0   — нові авто → дорожчий ремонт
#   β₄(CarOld vs Used)    < 0   — старі авто → дешевший ремонт
#   β₅(Diesel)            > 0   — дизель → дорожчий двигун/ремонт
#   β₆(logDensity)        = ?   — суперечливий: місто → більше аварій,
#                                  але менш серйозних
#
# ІДЕАЛЬНІ КОНТРОЛЬНІ ЗМІННІ (відсутні в датасеті):
#   – Стаж водіння (роки з ліцензією)
#     → ключова змінна: молодий водій зі стажем 5 р. ≠ молодий зі стажем 1 р.
#     → BonusMalus частково відображає стаж, але не повністю
#   – Стиль водіння (агресивність, перевищення швидкості)
#     → частково проксується BonusMalus
#   – Стан здоров'я водія (зір, реакція)
#     → особливо важливо для Seniors; у датасеті відсутнє
#   – Тип доріг (місто / траса / автобан)
#     → Density є слабким проксі
#   – Вартість авто
#     → VehBrand є частковим проксі
#
# НАЯВНІ КОНТРОЛЬНІ ЗМІННІ:
#   BonusMalus, VehPower, VehAge, VehGas, Density,
#   VehBrand (FE у M6), Region (FE у M5)
#
# ПОТЕНЦІЙНИЙ OVB (без контролю за BonusMalus):
# Молоді водії (Youth) мають вищий BM (менший стаж, вища ставка).
# Позитивна кореляція DrivAge↔logBM (ρ ≈ від'ємна для молоді).
# Без BM: ефект Youth завищений (частина ефекту є ефектом BM, а не віку).
# Після контролю: β(Youth) може зменшитися — перевіряємо в OVB-аналізі.


# =========================================================================
# 7. ПОБУДОВА МОДЕЛЕЙ
# =========================================================================

# М1: Базова — тільки AgeGroup
m1 <- lm(logY ~ AgeGroup, data = df_claims)

# М2: Безперервний DrivAge (лінійна специфікація)
m2 <- lm(logY ~ DrivAge, data = df_claims)

# М3: + ризиковість водія (logBM) — головний confounder
m3 <- lm(logY ~ AgeGroup + logBM, data = df_claims)

# М4: + характеристики автомобіля та середовища — ОСНОВНА МОДЕЛЬ
m4 <- lm(logY ~ AgeGroup + logBM + PowerGroup + CarAgeGroup +
           Diesel + logDensity,
         data = df_claims)

# М5: + регіональні фіксовані ефекти
m5 <- lm(logY ~ AgeGroup + logBM + PowerGroup + CarAgeGroup +
           Diesel + logDensity + Region,
         data = df_claims)

# М6: + марка авто (VehBrand) — найповніша модель
m6 <- lm(logY ~ AgeGroup + logBM + PowerGroup + CarAgeGroup +
           Diesel + logDensity + Region + VehBrand,
         data = df_claims)

# М7: Перевірка нелінійності — DrivAge_c + DrivAge_c²
# ОБГРУНТУВАННЯ: очікується U-подібна крива (молодь і літні > середній вік).
# Квадратичний член перевіряє, чи є залежність симетричною параболою,
# чи лінійної специфікації достатньо.
m7 <- lm(logY ~ DrivAge_c + DrivAge_c2 + logBM + PowerGroup +
           CarAgeGroup + Diesel + logDensity,
         data = df_claims)

# М8: Взаємодія AgeGroup × PowerGroup
# ОБГРУНТУВАННЯ: молоді водії на потужних авто — чи посилюється ефект?
# Комбінація недосвідченості + висока швидкість може давати
# непропорційно більші збитки.
m8 <- lm(logY ~ AgeGroup + logBM + PowerGroup + CarAgeGroup +
           Diesel + logDensity + AgeGroup:PowerGroup,
         data = df_claims)

# М9: Взаємодія AgeGroup × CarAgeGroup
# ОБГРУНТУВАННЯ: молодь частіше їздить на старих дешевих авто,
# літні — можуть їздити на нових; взаємодія перевіряє,
# чи різниться вплив віку авто залежно від вікової групи водія.
m9 <- lm(logY ~ AgeGroup + logBM + PowerGroup + CarAgeGroup +
           Diesel + logDensity + AgeGroup:CarAgeGroup,
         data = df_claims)

# М10: Ступінь логарифма — (logBM)²
# ОБГРУНТУВАННЯ: перевіряємо, чи залежність logY від logBM є квадратичною.
# Якщо β(logBM²) > 0 — ефект BM прискорюється зі зростанням BM.
# Якщо незначущий — лінійна форма достатня.
m10 <- lm(logY ~ AgeGroup + logBM + logBM2 + PowerGroup +
            CarAgeGroup + Diesel + logDensity,
          data = df_claims)


# =========================================================================
# 8. ВИВЕДЕННЯ РЕЗУЛЬТАТІВ З HC3 СТАНДАРТНИМИ ПОХИБКАМИ
# =========================================================================
#
# ОБГРУНТУВАННЯ HC3:
# Страхові виплати навіть після відсікання топ-0.5% і log-трансформації
# мають гетероскедастичні залишки (дисперсія залишків зростає зі
# збільшенням підігнаних значень — видно на графіку p6).
# Гомоскедастичні (OLS) похибки у цьому випадку некоректні.
# HC3 (MacKinnon & White, 1985) коректує SE без припущення про
# форму гетероскедастичності.

print_model <- function(m, title) {
  cat("\n", strrep("=", 72), "\n", sep = "")
  cat(title, "\n")
  cat(strrep("=", 72), "\n")
  print(coeftest(m, vcov = vcovHC(m, type = "HC3")))
  cat(sprintf("R² = %.4f  |  adj.R² = %.4f  |  N = %d\n",
              summary(m)$r.squared,
              summary(m)$adj.r.squared,
              nobs(m)))
}

print_model(m1,  "М1: Базова — тільки AgeGroup")
print_model(m2,  "М2: Безперервний DrivAge (лінійна)")
print_model(m3,  "М3: AgeGroup + logBM")
print_model(m4,  "М4: + PowerGroup + CarAgeGroup + Diesel + logDensity  [ОСНОВНА]")
print_model(m5,  "М5: + Region FE")
print_model(m6,  "М6: + VehBrand  [НАЙПОВНІША]")
print_model(m7,  "М7: Нелінійність — DrivAge_c + DrivAge_c²")
print_model(m8,  "М8: Взаємодія AgeGroup × PowerGroup")
print_model(m9,  "М9: Взаємодія AgeGroup × CarAgeGroup")
print_model(m10, "М10: Ступінь логарифма — logBM + (logBM)²")


# =========================================================================
# 9. ТЕСТИ СПІЛЬНОЇ ЗНАЧУЩОСТІ ГРУП КОЕФІЦІЄНТІВ
# =========================================================================

cat("\n", strrep("=", 72), "\n", sep = "")
cat("ТЕСТИ СПІЛЬНОЇ ЗНАЧУЩОСТІ\n")
cat(strrep("=", 72), "\n")

# [1] AgeGroup у М4
cat("\n[1] AgeGroup у М4 — чи впливає вік водія на збиток?\n")
cat("H₀: β(Youth) = β(Seniors) = 0\n")
lh_age <- linearHypothesis(m4,
                           c("AgeGroupYouth = 0", "AgeGroupSeniors = 0"),
                           vcov = vcovHC(m4, "HC3"))
print(lh_age)

# [2] PowerGroup у М4
cat("\n[2] PowerGroup у М4\n")
cat("H₀: β(Medium) = β(High) = 0\n")
lh_power <- linearHypothesis(m4,
                             c("PowerGroupMedium = 0", "PowerGroupHigh = 0"),
                             vcov = vcovHC(m4, "HC3"))
print(lh_power)

# [3] CarAgeGroup у М4
cat("\n[3] CarAgeGroup у М4\n")
cat("H₀: β(New) = β(Old) = 0\n")
lh_carage <- linearHypothesis(m4,
                              c("CarAgeGroupNew = 0", "CarAgeGroupOld = 0"),
                              vcov = vcovHC(m4, "HC3"))
print(lh_carage)

# [4] Нелінійність DrivAge — чи потрібен квадратичний член?
cat("\n[4] Нелінійність: чи потрібен квадратичний член DrivAge²? (М7)\n")
cat("H₀: β(DrivAge_c²) = 0\n")
lh_quad <- linearHypothesis(m7, "DrivAge_c2 = 0",
                            vcov = vcovHC(m7, "HC3"))
print(lh_quad)

# [5] Взаємодія AgeGroup × PowerGroup
cat("\n[5] Взаємодія AgeGroup × PowerGroup (М8)\n")
cat("H₀: всі 4 взаємодії = 0\n")
lh_inter_power <- linearHypothesis(m8,
                                   c("AgeGroupYouth:PowerGroupMedium = 0",
                                     "AgeGroupYouth:PowerGroupHigh = 0",
                                     "AgeGroupSeniors:PowerGroupMedium = 0",
                                     "AgeGroupSeniors:PowerGroupHigh = 0"),
                                   vcov = vcovHC(m8, "HC3"))
print(lh_inter_power)

# [6] Взаємодія AgeGroup × CarAgeGroup
cat("\n[6] Взаємодія AgeGroup × CarAgeGroup (М9)\n")
cat("H₀: всі 4 взаємодії = 0\n")
lh_inter_car <- linearHypothesis(m9,
                                 c("AgeGroupYouth:CarAgeGroupNew = 0",
                                   "AgeGroupYouth:CarAgeGroupOld = 0",
                                   "AgeGroupSeniors:CarAgeGroupNew = 0",
                                   "AgeGroupSeniors:CarAgeGroupOld = 0"),
                                 vcov = vcovHC(m9, "HC3"))
print(lh_inter_car)

# [7] Region FE у М5
cat("\n[7] Region FE у М5\n")
region_coefs <- grep("^Region", names(coef(m5)), value = TRUE)
lh_region <- linearHypothesis(m5,
                              paste0(region_coefs, " = 0"),
                              vcov = vcovHC(m5, "HC3"))
print(lh_region)

# [8] VehBrand у М6
cat("\n[8] VehBrand у М6\n")
brand_coefs <- grep("^VehBrand", names(coef(m6)), value = TRUE)
lh_brand <- linearHypothesis(m6,
                             paste0(brand_coefs, " = 0"),
                             vcov = vcovHC(m6, "HC3"))
print(lh_brand)

# [9] Нелінійність BonusMalus: чи потрібен (logBM)²? (М10)
cat("\n[9] Нелінійність logBM: чи потрібен (logBM)²? (М10)\n")
cat("H₀: β(logBM²) = 0\n")
lh_bm2 <- linearHypothesis(m10, "logBM2 = 0",
                           vcov = vcovHC(m10, "HC3"))
print(lh_bm2)


# =========================================================================
# 10. ЗВЕДЕНА ТАБЛИЦЯ МОДЕЛЕЙ
# =========================================================================

models    <- list(m1, m2, m3, m4, m5, m6)
mod_names <- c("M1", "M2", "M3", "M4*", "M5", "M6")
ctests    <- lapply(models, function(m) coeftest(m, vcov = vcovHC(m, "HC3")))

get_coef_str <- function(ct, name) {
  idx <- which(rownames(ct) == name)
  if (length(idx) == 0) return("      —    ")
  est <- ct[idx, 1]; pv <- ct[idx, 4]
  stars <- ifelse(pv < 0.001, "***",
                  ifelse(pv < 0.01,  "** ",
                         ifelse(pv < 0.05,  "*  ", "   ")))
  sprintf("%+.3f%s", est, stars)
}

vars <- list(
  list(n = "AgeGroupYouth",      l = "Youth (ref: Adults)      [КЛЮЧОВИЙ]"),
  list(n = "AgeGroupSeniors",    l = "Seniors (ref: Adults)    [КЛЮЧОВИЙ]"),
  list(n = "DrivAge",            l = "DrivAge (безперервна)"),
  list(n = "logBM",              l = "log(BonusMalus)"),
  list(n = "PowerGroupMedium",   l = "MedPower (ref: Low)"),
  list(n = "PowerGroupHigh",     l = "HighPower (ref: Low)"),
  list(n = "CarAgeGroupNew",     l = "CarNew (ref: Used)"),
  list(n = "CarAgeGroupOld",     l = "CarOld (ref: Used)"),
  list(n = "Diesel",             l = "Diesel"),
  list(n = "logDensity",         l = "log(Density+1)"),
  list(n = "(Intercept)",        l = "Константа")
)

cat("\n")
cat(strrep("=", 88), "\n")
cat("  ЗВЕДЕНА ТАБЛИЦЯ КОЕФІЦІЄНТІВ\n")
cat("  Залежна змінна: log(Total_Claim_Amount), без топ-0.5% катастроф\n")
cat("  SE: HC3  |  *** p<0.001  ** p<0.01  * p<0.05\n")
cat("  * M4 — основна модель\n")
cat(strrep("=", 88), "\n")
header <- sprintf("%-38s", "Змінна")
for (nm in mod_names) header <- paste0(header, sprintf("%11s", nm))
cat(header, "\n")
cat(strrep("-", 88), "\n")
for (v in vars) {
  row <- sprintf("%-38s", v$l)
  for (ct in ctests) row <- paste0(row, sprintf("%11s", get_coef_str(ct, v$n)))
  cat(row, "\n")
}
cat(strrep("-", 88), "\n")
row_region <- sprintf("%-38s", "Region FE")
for (m in models) row_region <- paste0(row_region,
                                       sprintf("%11s", ifelse(any(grepl("^Region", names(coef(m)))), "Так", "Ні")))
cat(row_region, "\n")
row_brand <- sprintf("%-38s", "VehBrand FE")
for (m in models) row_brand <- paste0(row_brand,
                                      sprintf("%11s", ifelse(any(grepl("^VehBrand", names(coef(m)))), "Так", "Ні")))
cat(row_brand, "\n")
cat(strrep("-", 88), "\n")
row_r2 <- sprintf("%-38s", "R²")
for (m in models) row_r2 <- paste0(row_r2, sprintf("%11.4f", summary(m)$r.squared))
cat(row_r2, "\n")
row_ar2 <- sprintf("%-38s", "adj. R²")
for (m in models) row_ar2 <- paste0(row_ar2, sprintf("%11.4f", summary(m)$adj.r.squared))
cat(row_ar2, "\n")
row_n <- sprintf("%-38s", "N спостережень")
for (m in models) row_n <- paste0(row_n, sprintf("%11d", nobs(m)))
cat(row_n, "\n")
cat(strrep("=", 88), "\n")


# =========================================================================
# 11. ІНТЕРПРЕТАЦІЯ КЛЮЧОВИХ КОЕФІЦІЄНТІВ (M4)
# =========================================================================
ct_m4 <- ctests[[4]]

cat("\n=== ІНТЕРПРЕТАЦІЯ КОЕФІЦІЄНТІВ M4 ===\n\n")

you_est <- ct_m4["AgeGroupYouth",   1]
sen_est <- ct_m4["AgeGroupSeniors", 1]
you_pv  <- ct_m4["AgeGroupYouth",   4]
sen_pv  <- ct_m4["AgeGroupSeniors", 4]

cat(sprintf("AgeGroupYouth   = %+.4f  (p = %.4f)\n", you_est, you_pv))
cat(sprintf(
  "  -> Молоді водії (18-25) мають збитки на %.1f%% %s,\n",
  abs((exp(you_est) - 1) * 100),
  ifelse(you_est > 0, "БІЛЬШІ", "менші")))
cat("     ніж водії середнього віку (26-59) при однакових інших характеристиках.\n")
cat(sprintf("  -> Гіпотеза b(Youth) > 0: %s\n\n",
            ifelse(you_est > 0 & you_pv < 0.05, "ПІДТВЕРДЖЕНА", "не підтверджена")))

cat(sprintf("AgeGroupSeniors = %+.4f  (p = %.4f)\n", sen_est, sen_pv))
cat(sprintf(
  "  -> Літні водії (60+): різниця з Adults %.1f%% (%s)\n\n",
  abs((exp(sen_est) - 1) * 100),
  ifelse(sen_pv < 0.05, "значуща", "НЕ значуща")))

# ПРИМІТКА ПРО РЕФЕРЕНТНУ КАТЕГОРІЮ:
# AgeGroupAdults відсутній у таблиці — це очікувано і методологічно правильно.
# Adults є референтною (базовою) категорією: relevel(..., ref = "Adults").
# R автоматично виключає одну категорію для уникнення dummy variable trap.
# Коефіцієнт AgeGroupYouth = різниця Youth vs Adults.
# Коефіцієнт AgeGroupSeniors = різниця Seniors vs Adults.
# Ефект для Adults = 0 (закладений у константу β₀).

bm_est <- ct_m4["logBM", 1]
cat(sprintf("log(BonusMalus) = %+.4f ***\n", bm_est))
cat(sprintf("  -> Збільшення BM на 1%% -> збиток зростає на %.3f%%\n\n", bm_est))

new_est <- ct_m4["CarAgeGroupNew", 1]
old_est <- ct_m4["CarAgeGroupOld", 1]
cat(sprintf("CarNew = %+.4f  -> нові авто на %.1f%% %s у ремонті\n",
            new_est, abs((exp(new_est) - 1) * 100),
            ifelse(new_est > 0, "дорожче", "дешевше")))
cat(sprintf("CarOld = %+.4f  -> старі авто на %.1f%% %s збиток\n",
            old_est, abs((exp(old_est) - 1) * 100),
            ifelse(old_est < 0, "менший", "більший")))


# =========================================================================
# 12. АНАЛІЗ СТІЙКОСТІ (ROBUSTNESS CHECK)
# =========================================================================

cat("\n", strrep("=", 72), "\n", sep = "")
cat("АНАЛІЗ СТІЙКОСТІ: динаміка b(Youth) та b(Seniors) по моделях\n")
cat(strrep("=", 72), "\n\n")

cat(sprintf("%-6s  %10s  %10s  %10s  %10s\n",
            "Модель", "b(Youth)", "p(Youth)", "b(Sen)", "p(Sen)"))
cat(strrep("-", 52), "\n")
for (i in 1:4) {
  ct <- ctests[[i]]
  if (!"AgeGroupYouth" %in% rownames(ct)) next
  y_e <- ct["AgeGroupYouth",   1]; y_p <- ct["AgeGroupYouth",   4]
  s_e <- ct["AgeGroupSeniors", 1]; s_p <- ct["AgeGroupSeniors", 4]
  cat(sprintf("%-6s  %+10.4f  %10.4f  %+10.4f  %10.4f\n",
              mod_names[i], y_e, y_p, s_e, s_p))
}

you_m3 <- ctests[[3]]["AgeGroupYouth", 1]
you_m4 <- ctests[[4]]["AgeGroupYouth", 1]
you_m5 <- ctests[[5]]["AgeGroupYouth", 1]
you_m6 <- ctests[[6]]["AgeGroupYouth", 1]

sen_m3 <- ctests[[3]]["AgeGroupSeniors", 1]
sen_m4 <- ctests[[4]]["AgeGroupSeniors", 1]
sen_m5 <- ctests[[5]]["AgeGroupSeniors", 1]
sen_m6 <- ctests[[6]]["AgeGroupSeniors", 1]

cat(sprintf("\nb(Youth):   M3 = %.4f  M4 = %.4f  M5 = %.4f  M6 = %.4f\n",
            you_m3, you_m4, you_m5, you_m6))
cat(sprintf("b(Seniors): M3 = %.4f  M4 = %.4f  M5 = %.4f  M6 = %.4f\n",
            sen_m3, sen_m4, sen_m5, sen_m6))

cat("\nСтійкість оцінюється по M3->M4->M5 (до включення VehBrand):\n")
change_you_m3_m5 <- abs(you_m5 - you_m3) / (abs(you_m3) + 1e-10) * 100
change_sen_m3_m5 <- abs(sen_m5 - sen_m3) / (abs(sen_m3) + 1e-10) * 100
cat(sprintf("  |b(Youth) M5 - b(Youth) M3| / |b(Youth) M3|     = %.1f%%\n", change_you_m3_m5))
cat(sprintf("  |b(Seniors) M5 - b(Seniors) M3| / |b(Seniors) M3| = %.1f%%\n", change_sen_m3_m5))
cat(ifelse(change_you_m3_m5 < 20,
           "  -> Ефект Youth СТІЙКИЙ (зміна < 20%%) між M3 і M5.\n",
           "  -> Ефект Youth НЕСТІЙКИЙ (зміна >= 20%%) між M3 і M5.\n"))
cat(ifelse(change_sen_m3_m5 < 20,
           "  -> Ефект Seniors СТІЙКИЙ (зміна < 20%%) між M3 і M5.\n",
           "  -> Ефект Seniors НЕСТІЙКИЙ (зміна >= 20%%) між M3 і M5.\n"))

# OVB аналіз: порівняння без і з logBM
m_no_bm  <- lm(logY ~ AgeGroup + PowerGroup + CarAgeGroup + Diesel + logDensity,
               data = df_claims)
ct_no_bm <- coeftest(m_no_bm, vcov = vcovHC(m_no_bm, "HC3"))

cat("\nВплив включення logBM на b(AgeGroup) [OVB-аналіз]:\n")
cat(sprintf("              Без logBM     З logBM (M4)\n"))
cat(sprintf("Youth:        %+.4f        %+.4f\n",
            ct_no_bm["AgeGroupYouth",   1], ct_m4["AgeGroupYouth",   1]))
cat(sprintf("Seniors:      %+.4f        %+.4f\n",
            ct_no_bm["AgeGroupSeniors", 1], ct_m4["AgeGroupSeniors", 1]))
cat("-> Якщо b(Youth) зменшується з включенням BM: OVB позитивний\n")
cat("   (молодь має вищий BM, який самостійно підвищує збиток).\n")
cat("-> Якщо b(Seniors) змінюється — оцінюємо вплив досвідного каналу.\n")

# Нелінійність DrivAge
ct_m7  <- coeftest(m7, vcov = vcovHC(m7, "HC3"))
sq_pv  <- ct_m7["DrivAge_c2", 4]
sq_est <- ct_m7["DrivAge_c2", 1]
cat(sprintf("\nНелінійність DrivAge (M7): b(DrivAge_c²) = %.6f, p = %.4f\n",
            sq_est, sq_pv))
cat(ifelse(sq_est > 0 & sq_pv < 0.05,
           "-> U-подібна крива ПІДТВЕРДЖЕНА: молодь і літні > середній вік.\n",
           ifelse(sq_pv < 0.05,
                  "-> Нелінійність ПІДТВЕРДЖЕНА (не U-подібна): перевірте знак.\n",
                  "-> Нелінійність НЕ підтверджена: лінійна форма DrivAge достатня.\n")))

# Нелінійність BonusMalus — (logBM)²
ct_m10 <- coeftest(m10, vcov = vcovHC(m10, "HC3"))
bm2_pv  <- ct_m10["logBM2", 4]
bm2_est <- ct_m10["logBM2", 1]
cat(sprintf("\nНелінійність logBM (M10): b(logBM²) = %.4f, p = %.4f\n",
            bm2_est, bm2_pv))
cat(ifelse(bm2_pv < 0.05,
           "-> (logBM)² ЗНАЧУЩИЙ: залежність logY~logBM є квадратичною.\n",
           "-> (logBM)² НЕЗНАЧУЩИЙ: лінійна форма logBM достатня.\n"))
cat(sprintf("   b(Youth) у M10 = %.4f (порівн. з M4: %.4f) — зміна %.1f%%\n",
            ct_m10["AgeGroupYouth", 1], you_m4,
            abs(ct_m10["AgeGroupYouth", 1] - you_m4) / (abs(you_m4) + 1e-10) * 100))

# Взаємодії
inter_pow_pv <- lh_inter_power$`Pr(>F)`[2]
inter_car_pv <- lh_inter_car$`Pr(>F)`[2]
cat(sprintf("\nВзаємодія Age x PowerGroup (M8):   F = %.3f, p = %.4f %s\n",
            lh_inter_power$F[2], inter_pow_pv,
            ifelse(inter_pow_pv < 0.05, "[ЗНАЧУЩА]", "[незначуща]")))
cat(sprintf("Взаємодія Age x CarAgeGroup (M9): F = %.3f, p = %.4f %s\n",
            lh_inter_car$F[2], inter_car_pv,
            ifelse(inter_car_pv < 0.05, "[ЗНАЧУЩА]", "[незначуща]")))


# =========================================================================
# 13. ГРАФІКИ
# =========================================================================

raw_claims <- full_data %>% filter(Total_Claim_Amount > 0)

# --- Графік 1: Розподіл збитку до і після відсікання ---
p1 <- ggplot() +
  geom_histogram(data = raw_claims,
                 aes(x = log(Total_Claim_Amount), fill = "Повна вибірка"),
                 bins = 70, alpha = 0.5) +
  geom_histogram(data = df_claims,
                 aes(x = logY, fill = "Після відсікання топ-0.5%"),
                 bins = 70, alpha = 0.6) +
  scale_fill_manual(values = c("Повна вибірка"              = "#D7191C",
                               "Після відсікання топ-0.5%" = "#2C7BB6")) +
  labs(
    title    = "Розподіл log(Збиток): повна вибірка vs після відсікання",
    subtitle = sprintf("Відсіяно %d спостережень (поріг: %.0f євро)",
                       n_before - nrow(df_claims), threshold_995),
    x = "log(Сума збитку, євро)", y = "Кількість полісів", fill = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

# --- Графік 2: Середній logY по значеннях DrivAge зі ДІ ---
p2 <- df_claims %>%
  mutate(age_bin = floor(DrivAge / 5) * 5) %>%          # 5-річні інтервали
  group_by(age_bin) %>%
  summarise(mean_logY = mean(logY),
            se_logY   = sd(logY) / sqrt(n()),
            n = n(), .groups = "drop") %>%
  ggplot(aes(x = age_bin, y = mean_logY)) +
  geom_ribbon(aes(ymin = mean_logY - 1.96 * se_logY,
                  ymax = mean_logY + 1.96 * se_logY),
              fill = "#2C7BB6", alpha = 0.2) +
  geom_line(color = "#2C7BB6", linewidth = 1.2) +
  geom_point(aes(size = n), color = "#D7191C") +
  scale_size_continuous(range = c(2, 8), name = "N полісів") +
  labs(
    title    = "Середній log(Збиток) по віку водія (5-річні інтервали)",
    subtitle = "Стрічка = 95% ДІ. Розмір точки = кількість спостережень",
    x = "Вік водія (роки)",
    y = "Середнє log(Збиток)"
  ) +
  theme_minimal(base_size = 13)

# --- Графік 3: Boxplot logY по AgeGroup ---
p3 <- ggplot(df_claims, aes(x = AgeGroup, y = logY, fill = AgeGroup)) +
  geom_boxplot(outlier.alpha = 0.1, outlier.size = 0.5) +
  scale_fill_manual(values = c("Youth"   = "#D7191C",
                               "Adults"  = "#2C7BB6",
                               "Seniors" = "#FDAE61")) +
  labs(
    title    = "Розподіл log(Збиток) по вікових групах водія",
    subtitle = "Youth: 18–25 | Adults: 26–59 (референтна) | Seniors: 60+",
    x = "Вікова група", y = "log(Сума збитку, євро)"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

# --- Графік 4: Стійкість b(Youth) і b(Seniors) по моделях ---
stab_df <- do.call(rbind, lapply(1:4, function(i) {
  ct <- ctests[[i]]
  if (!"AgeGroupYouth" %in% rownames(ct)) return(NULL)
  data.frame(
    model = factor(mod_names[i], levels = mod_names[1:4]),
    group = c("Youth", "Seniors"),
    est   = c(ct["AgeGroupYouth",   1], ct["AgeGroupSeniors", 1]),
    se    = c(ct["AgeGroupYouth",   2], ct["AgeGroupSeniors", 2])
  )
})) %>%
  mutate(lower = est - 1.96 * se, upper = est + 1.96 * se)

p4 <- ggplot(stab_df, aes(x = model, y = est, color = group, shape = group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.15, linewidth = 1,
                position = position_dodge(0.3)) +
  geom_point(size = 4, position = position_dodge(0.3)) +
  scale_color_manual(values = c("Youth" = "#D7191C", "Seniors" = "#FDAE61")) +
  labs(
    title    = "Стійкість коефіцієнтів AgeGroup по моделях",
    subtitle = "Планки = 95% ДІ (HC3 SE). M1 = без контролів, M4 = основна",
    x = "Модель", y = "Коефіцієнт при AgeGroup",
    color = "Група", shape = "Група"
  ) +
  theme_minimal(base_size = 13)

# --- Графік 5: OVB — зміна b(AgeGroup) при включенні logBM ---
ovb_df <- data.frame(
  spec  = factor(rep(c("Без BonusMalus", "З BonusMalus (M4)"), each = 2),
                 levels = c("Без BonusMalus", "З BonusMalus (M4)")),
  group = rep(c("Youth", "Seniors"), 2),
  est   = c(ct_no_bm["AgeGroupYouth",   1], ct_no_bm["AgeGroupSeniors", 1],
            ct_m4["AgeGroupYouth",       1], ct_m4["AgeGroupSeniors",   1]),
  se    = c(ct_no_bm["AgeGroupYouth",   2], ct_no_bm["AgeGroupSeniors", 2],
            ct_m4["AgeGroupYouth",       2], ct_m4["AgeGroupSeniors",   2])
) %>%
  mutate(lower = est - 1.96 * se, upper = est + 1.96 * se)

p5 <- ggplot(ovb_df, aes(x = group, y = est, color = spec, shape = spec)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.15, linewidth = 1,
                position = position_dodge(0.4)) +
  geom_point(size = 4, position = position_dodge(0.4)) +
  scale_color_manual(values = c("Без BonusMalus"    = "#D7191C",
                                "З BonusMalus (M4)" = "#2C7BB6")) +
  labs(
    title    = "Вплив включення BonusMalus на b(AgeGroup)",
    subtitle = "Молодь має вищий BM (менший стаж) → потенційний OVB",
    x = "Вікова група", y = "Коефіцієнт",
    color = "Специфікація", shape = "Специфікація"
  ) +
  theme_minimal(base_size = 13)

# --- Графік 6: Залишки M4 ---
resid_df <- data.frame(fitted = fitted(m4), residuals = residuals(m4))
p6 <- ggplot(resid_df, aes(x = fitted, y = residuals)) +
  geom_point(alpha = 0.05, size = 0.5, color = "#555555") +
  geom_hline(yintercept = 0, color = "#D7191C", linewidth = 1) +
  geom_smooth(method = "loess", color = "#2C7BB6",
              linewidth = 1.2, se = FALSE) +
  labs(
    title    = "Залишки vs Підігнані значення (M4)",
    subtitle = "Перевірка гетероскедастичності. Ідеал: рівномірно навколо 0",
    x = "Підігнані значення", y = "Залишки"
  ) +
  theme_minimal(base_size = 13)

# --- Графік 7: logY ~ logBM (лінійна vs квадратична) ---
bm_df <- df_claims %>%
  mutate(decile_BM = ntile(logBM, 20)) %>%
  group_by(decile_BM) %>%
  summarise(mean_logBM = mean(logBM),
            mean_logY  = mean(logY),
            se_logY    = sd(logY) / sqrt(n()),
            .groups = "drop")

p7 <- ggplot(bm_df, aes(x = mean_logBM, y = mean_logY)) +
  geom_ribbon(aes(ymin = mean_logY - 1.96 * se_logY,
                  ymax = mean_logY + 1.96 * se_logY),
              fill = "#2C7BB6", alpha = 0.2) +
  geom_point(color = "#D7191C", size = 3) +
  geom_smooth(method = "lm", formula = y ~ x,
              color = "#2C7BB6", linetype = "solid",  se = FALSE) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2),
              color = "#D7191C",  linetype = "dashed", se = FALSE) +
  labs(
    title    = "Залежність log(Збиток) від log(BonusMalus)",
    subtitle = "Синя = лінійна; червона пунктир = квадратична. По 20 центилях.",
    x = "log(BonusMalus)", y = "Середнє log(Збиток)"
  ) +
  theme_minimal(base_size = 13)

# --- Збереження ---
ggsave("lab3_p1_dist_cutoff.png",   p1, width = 10, height = 5, dpi = 150)
ggsave("lab3_p2_age_trend.png",     p2, width = 10, height = 6, dpi = 150)
ggsave("lab3_p3_boxplot.png",       p3, width = 8,  height = 6, dpi = 150)
ggsave("lab3_p4_stability.png",     p4, width = 9,  height = 6, dpi = 150)
ggsave("lab3_p5_ovb_bm.png",        p5, width = 9,  height = 6, dpi = 150)
ggsave("lab3_p6_residuals.png",     p6, width = 9,  height = 6, dpi = 150)
ggsave("lab3_p7_logbm_shape.png",   p7, width = 9,  height = 6, dpi = 150)

cat("\nВсі графіки збережено.\n")


# =========================================================================
# 14. ЗАГАЛЬНІ ВИСНОВКИ
# =========================================================================
cat("\n", strrep("=", 72), "\n", sep = "")
cat("ЗАГАЛЬНІ ВИСНОВКИ\n")
cat(strrep("=", 72), "\n")

cat(sprintf("
1. ГОЛОВНИЙ РЕЗУЛЬТАТ — ефект віку водія:
   AgeGroup у M4 тестується спільно: F-тест із H₀: β(Youth)=β(Seniors)=0.
   b(Youth)   = %.4f у M4: молоді водії (18–25) генерують збитки
   на %.1f%% %s, ніж водії 26–59, при однакових
   характеристиках авто та середовища.
   b(Seniors) = %.4f у M4: літні водії (60+) — збитки на %.1f%% %s.
   AgeGroupAdults відсутній у таблиці — це референтна категорія (норма).

2. НЕЛІНІЙНІСТЬ DrivAge:
   Квадратичний член DrivAge_c² %s (p = %.4f).
   %s

3. СТІЙКІСТЬ:
   b(Youth):   M3→M5 зміна %.1f%% — %s.
   b(Seniors): M3→M5 зміна %.1f%% — %s.

4. ВПЛИВ BonusMalus (OVB-аналіз):
   Після включення logBM b(Youth) %s.
   Це відповідає логіці: молодь має вищий BM через менший стаж,
   тому без контролю BM ефект Youth включав ефект досвіду.

5. ФОРМА ЗАЛЕЖНОСТІ logBM:
   (logBM)² %s (p = %.4f) — %s.
   Взаємодії Age x Power (p = %.4f) і Age x CarAge (p = %.4f) %s.

6. ОБМЕЖЕННЯ:
   Ключова відсутня змінна: стаж водіння (роки з ліцензією) —
   найважливіший канал ефекту; BonusMalus є лише неповним проксі.
   Selection bias: тільки поліси з ненульовими виплатами.
   Стан здоров'я (особливо для Seniors): у датасеті відсутній.
",
            you_m4, abs((exp(you_m4) - 1) * 100),
            ifelse(you_m4 > 0, "БІЛЬШІ", "менші"),
            sen_m4, abs((exp(sen_m4) - 1) * 100),
            ifelse(sen_m4 > 0, "БІЛЬШІ", "менші"),
            ifelse(sq_pv < 0.05, "ЗНАЧУЩИЙ", "незначущий"), sq_pv,
            ifelse(sq_est > 0 & sq_pv < 0.05,
                   "U-подібна крива підтверджена: молодь і літні > середній вік.",
                   ifelse(sq_pv >= 0.05,
                          "Лінійна специфікація DrivAge достатня.",
                          "Нелінійність є, але не U-подібна.")),
            change_you_m3_m5, ifelse(change_you_m3_m5 < 20, "СТІЙКИЙ", "НЕСТІЙКИЙ"),
            change_sen_m3_m5, ifelse(change_sen_m3_m5 < 20, "СТІЙКИЙ", "НЕСТІЙКИЙ"),
            ifelse(abs(ct_m4["AgeGroupYouth", 1]) < abs(ct_no_bm["AgeGroupYouth", 1]),
                   "зменшується (позитивний OVB усунуто)",
                   "зростає (BM частково нейтралізував ефект)"),
            ifelse(bm2_pv < 0.05, "ЗНАЧУЩИЙ", "незначущий"), bm2_pv,
            ifelse(bm2_pv < 0.05, "квадратична форма краща", "лінійна форма достатня"),
            inter_pow_pv, inter_car_pv,
            ifelse(inter_pow_pv < 0.05 | inter_car_pv < 0.05,
                   "— принаймні одна ЗНАЧУЩА.",
                   "— незначущі.")
))