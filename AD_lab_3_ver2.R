
packages <- c(
  "tidyverse",
  "sandwich",
  "lmtest",
  "broom",
  "scales",
  "pROC"
)

new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) {
  install.packages(new_packages)
}

library(tidyverse)
library(sandwich)
library(lmtest)
library(broom)
library(scales)
library(pROC)

freq_path <- "freMTPL2freq.csv"
sev_path  <- "freMTPL2sev.csv"

# Якщо файли лежать у папці data:
# freq_path <- "data/freMTPL2freq.csv"
# sev_path  <- "data/freMTPL2sev.csv"

out_dir <- "lab3_outputs"
plot_dir <- file.path(out_dir, "plots")
table_dir <- file.path(out_dir, "tables")

dir.create(out_dir, showWarnings = FALSE)
dir.create(plot_dir, showWarnings = FALSE)
dir.create(table_dir, showWarnings = FALSE)

set.seed(123)

stars <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    p < 0.1   ~ ".",
    TRUE ~ ""
  )
}

vcov_ols <- function(model) {
  sandwich::vcovHC(model, type = "HC3")
}

vcov_glm <- function(model) {
  sandwich::vcovHC(model, type = "HC1")
}

coef_table <- function(model, vcov_mat, model_name, model_type = c("ols_log", "logit", "linear")) {
  model_type <- match.arg(model_type)
  
  ct <- lmtest::coeftest(model, vcov. = vcov_mat)
  mat <- as.matrix(ct)
  
  df <- tibble(
    model = model_name,
    term = rownames(mat),
    estimate = as.numeric(mat[, 1]),
    std_error = as.numeric(mat[, 2]),
    statistic = as.numeric(mat[, 3]),
    p_value = as.numeric(mat[, 4])
  ) %>%
    mutate(
      ci_low = estimate - 1.96 * std_error,
      ci_high = estimate + 1.96 * std_error,
      stars = stars(p_value)
    )
  
  if (model_type == "ols_log") {
    df <- df %>%
      mutate(
        percent_effect = 100 * (exp(estimate) - 1),
        percent_ci_low = 100 * (exp(ci_low) - 1),
        percent_ci_high = 100 * (exp(ci_high) - 1)
      )
  }
  
  if (model_type == "logit") {
    df <- df %>%
      mutate(
        odds_ratio = exp(estimate),
        or_ci_low = exp(ci_low),
        or_ci_high = exp(ci_high)
      )
  }
  
  return(df)
}

model_summary_lm <- function(model, model_name = NULL) {
  # якщо випадково передали model_name першим, а model другим
  if (is.character(model) && inherits(model_name, "lm")) {
    tmp <- model
    model <- model_name
    model_name <- tmp
  }
  
  if (!inherits(model, "lm")) {
    stop("model_summary_lm отримала не lm-модель, а: ", paste(class(model), collapse = ", "))
  }
  
  s <- summary(model)
  
  tibble(
    model = as.character(model_name),
    n = nobs(model),
    r2 = s$r.squared,
    adj_r2 = s$adj.r.squared,
    aic = AIC(model),
    bic = BIC(model)
  )
}
#model_summary_glm <- function(model, model_name = NULL) {
#  
# if (!inherits(model, "glm")) {
#   stop(
#     "model_summary_glm отримала не glm-модель. ",
#     "model_name = ", model_name,
#     "; class(model) = ", paste(class(model), collapse = ", "),
#     "; value(model) = ", paste(model, collapse = " ")
#   )
# }
# 
# ll_model <- -0.5 * model$deviance
# ll_null  <- -0.5 * model$null.deviance
# 
# tibble(
#   model = as.character(model_name),
#   n = length(model$y),
#   logLik_approx = ll_model,
#   AIC = model$aic,
#   BIC = model$aic + (log(length(model$y)) - 2) * length(coef(model)),
#   McFadden_R2 = 1 - ll_model / ll_null
# )
#

wald_test_terms <- function(model, term_pattern, vcov_mat, test_name) {
  beta <- coef(model)
  idx <- grep(term_pattern, names(beta))
  
  if (length(idx) == 0) {
    return(tibble(
      test = test_name,
      terms_pattern = term_pattern,
      df = 0,
      wald_chi2 = NA_real_,
      p_value = NA_real_,
      note = "No matching terms"
    ))
  }
  
  b <- beta[idx]
  V <- vcov_mat[idx, idx, drop = FALSE]
  
  W <- as.numeric(t(b) %*% solve(V) %*% b)
  df <- length(b)
  p <- pchisq(W, df = df, lower.tail = FALSE)
  
  tibble(
    test = test_name,
    terms_pattern = term_pattern,
    df = df,
    wald_chi2 = W,
    p_value = p,
    note = ""
  )
}

auc_binary <- function(y, p) {
  # приводимо y до нормальних 0/1
  if (is.factor(y)) {
    y <- as.character(y)
  }
  
  if (is.character(y)) {
    y <- ifelse(y %in% c("1", "TRUE", "Yes", "yes"), 1, 0)
  }
  
  y <- as.numeric(y)
  p <- as.numeric(p)
  
  # якщо після factor вийшло 1/2, перетворюємо на 0/1
  if (all(na.omit(unique(y)) %in% c(1, 2))) {
    y <- y - 1
  }
  
  ok <- !is.na(y) & !is.na(p)
  y <- y[ok]
  p <- p[ok]
  
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  
  if (n1 == 0 || n0 == 0) {
    return(NA_real_)
  }
  
  r <- rank(p, ties.method = "average")
  auc <- (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  
  return(auc)
}
auc_from_scores <- function(y, p) {
  y_raw <- as.character(y)
  
  y01 <- dplyr::case_when(
    y_raw %in% c("1", "TRUE", "True", "true") ~ 1,
    y_raw %in% c("0", "FALSE", "False", "false") ~ 0,
    TRUE ~ NA_real_
  )
  
  p <- as.numeric(p)
  
  df <- tibble(y = y01, p = p) %>%
    filter(!is.na(y), !is.na(p))
  
  print(table(df$y))
  print(summary(df$p))
  
  n1 <- sum(df$y == 1)
  n0 <- sum(df$y == 0)
  
  if (n1 == 0 || n0 == 0) return(NA_real_)
  
  r <- rank(df$p, ties.method = "average")
  
  auc <- (sum(r[df$y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  return(auc)
}
auc_manual <- function(y, p) {
  if (is.factor(y)) y <- as.character(y)
  if (is.logical(y)) y <- as.integer(y)
  
  y <- as.numeric(y)
  p <- as.numeric(p)
  
  df <- tibble(y = y, p = p) %>%
    filter(!is.na(y), !is.na(p)) %>%
    mutate(
      y = case_when(
        y == 1 ~ 1,
        y == 0 ~ 0,
        TRUE ~ NA_real_
      )
    ) %>%
    filter(!is.na(y))
  
  print(table(df$y))
  print(summary(df$p))
  
  n_pos <- sum(df$y == 1)
  n_neg <- sum(df$y == 0)
  
  if (n_pos == 0 || n_neg == 0) {
    stop("AUC неможливо порахувати: є тільки один клас.")
  }
  
  ranks <- rank(df$p, ties.method = "average")
  
  auc <- (sum(ranks[df$y == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
  
  return(as.numeric(auc))
}
calibration_table <- function(data, y_col, p_col, bins = 10) {
  data %>%
    filter(!is.na(.data[[y_col]]), !is.na(.data[[p_col]])) %>%
    mutate(bin = ntile(.data[[p_col]], bins)) %>%
    group_by(bin) %>%
    summarise(
      n = n(),
      mean_pred = mean(.data[[p_col]], na.rm = TRUE),
      observed_rate = mean(.data[[y_col]], na.rm = TRUE),
      se = sqrt(observed_rate * (1 - observed_rate) / n),
      ci_low = observed_rate - 1.96 * se,
      ci_high = observed_rate + 1.96 * se,
      .groups = "drop"
    )
}

save_plot <- function(plot, filename, width = 8, height = 5) {
  ggsave(
    filename = file.path(plot_dir, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = 300
  )
}

freq <- read_csv(freq_path, show_col_types = FALSE)
sev <- read_csv(sev_path, show_col_types = FALSE)

claims_by_policy <- sev %>%
  group_by(IDpol) %>%
  summarise(
    TotalClaim = sum(ClaimAmount, na.rm = TRUE),
    ClaimCount_sev = n(),
    .groups = "drop"
  )

policy_data <- freq %>%
  left_join(claims_by_policy, by = "IDpol") %>%
  mutate(
    TotalClaim = replace_na(TotalClaim, 0),
    ClaimCount_sev = replace_na(ClaimCount_sev, 0),
    HasClaim = ifelse(ClaimNb > 0, 1, 0)
  ) %>%
  filter(
    !is.na(Exposure),
    Exposure > 0,
    !is.na(BonusMalus),
    BonusMalus >= 50,
    BonusMalus <= 200,
    !is.na(DrivAge),
    !is.na(VehAge),
    !is.na(VehPower),
    !is.na(Density),
    !is.na(Area),
    !is.na(VehGas),
    !is.na(VehBrand),
    !is.na(Region)
  ) %>%
  mutate(
    logExposure = log(Exposure),
    logDensity = log(Density + 1),
    logBonusMalus = log(BonusMalus),
    
    BonusGroup = case_when(
      BonusMalus >= 50 & BonusMalus <= 80  ~ "50_80_Super",
      BonusMalus >= 81 & BonusMalus <= 100 ~ "81_100_Normal",
      BonusMalus >= 101 & BonusMalus <= 120 ~ "101_120_Risk",
      BonusMalus >= 121 & BonusMalus <= 200 ~ "121_200_Bad",
      TRUE ~ NA_character_
    ),
    
    AgeGroup = case_when(
      DrivAge >= 18 & DrivAge <= 25 ~ "Youth_18_25",
      DrivAge >= 26 & DrivAge <= 59 ~ "Adults_26_59",
      DrivAge >= 60 ~ "Seniors_60plus",
      TRUE ~ NA_character_
    ),
    
    PowerGroup = case_when(
      VehPower <= 5 ~ "Low",
      VehPower >= 6 & VehPower <= 8 ~ "Medium",
      VehPower >= 9 ~ "High",
      TRUE ~ NA_character_
    ),
    
    CarAgeGroup = case_when(
      VehAge <= 3 ~ "New_0_3",
      VehAge >= 4 & VehAge <= 10 ~ "Used_4_10",
      VehAge >= 11 ~ "Old_11plus",
      TRUE ~ NA_character_
    ),
    
    Area = as.factor(Area),
    VehGas = as.factor(VehGas),
    VehBrand = as.factor(VehBrand),
    Region = as.factor(Region),
    
    BonusGroup = factor(
      BonusGroup,
      levels = c("50_80_Super", "81_100_Normal", "101_120_Risk", "121_200_Bad")
    ),
    
    AgeGroup = factor(
      AgeGroup,
      levels = c("Adults_26_59", "Youth_18_25", "Seniors_60plus")
    ),
    
    PowerGroup = factor(
      PowerGroup,
      levels = c("Low", "Medium", "High")
    ),
    
    CarAgeGroup = factor(
      CarAgeGroup,
      levels = c("Used_4_10", "New_0_3", "Old_11plus")
    )
  ) %>%
  filter(
    !is.na(BonusGroup),
    !is.na(AgeGroup),
    !is.na(PowerGroup),
    !is.na(CarAgeGroup)
  )

# Severity sample
claims_data_raw <- policy_data %>%
  filter(TotalClaim > 0)

q995 <- quantile(claims_data_raw$TotalClaim, 0.995, na.rm = TRUE)
q95 <- quantile(claims_data_raw$TotalClaim, 0.95, na.rm = TRUE)

claims_data <- claims_data_raw %>%
  mutate(
    CappedClaim = pmin(TotalClaim, q995),
    logY = log(CappedClaim),
    LargeClaim = ifelse(TotalClaim > q95, 1, 0),
    
    DrivAge_c = as.numeric(scale(DrivAge, center = TRUE, scale = FALSE)),
    VehAge_c = as.numeric(scale(VehAge, center = TRUE, scale = FALSE)),
    VehPower_c = as.numeric(scale(VehPower, center = TRUE, scale = FALSE)),
    logDensity_c = as.numeric(scale(logDensity, center = TRUE, scale = FALSE)),
    logBonusMalus_c = as.numeric(scale(logBonusMalus, center = TRUE, scale = FALSE))
  )

policy_data <- policy_data %>%
  mutate(
    DrivAge_c = as.numeric(scale(DrivAge, center = TRUE, scale = FALSE)),
    VehAge_c = as.numeric(scale(VehAge, center = TRUE, scale = FALSE)),
    VehPower_c = as.numeric(scale(VehPower, center = TRUE, scale = FALSE)),
    logDensity_c = as.numeric(scale(logDensity, center = TRUE, scale = FALSE)),
    logBonusMalus_c = as.numeric(scale(logBonusMalus, center = TRUE, scale = FALSE))
  )

# Загальна інформація
data_info <- tibble(
  metric = c(
    "N policies",
    "N policies with claims",
    "Share HasClaim",
    "Q95 TotalClaim",
    "Q99.5 TotalClaim"
  ),
  value = c(
    nrow(policy_data),
    nrow(claims_data),
    mean(policy_data$HasClaim),
    as.numeric(q95),
    as.numeric(q995)
  )
)

write_csv(data_info, file.path(table_dir, "00_data_info.csv"))

# Описова статистика
bonus_desc <- claims_data %>%
  group_by(BonusGroup) %>%
  summarise(
    n = n(),
    LargeClaimRate = mean(LargeClaim),
    MeanClaim = mean(TotalClaim),
    MedianClaim = median(TotalClaim),
    MeanCappedClaim = mean(CappedClaim),
    MeanLogY = mean(logY),
    .groups = "drop"
  )

write_csv(bonus_desc, file.path(table_dir, "00_bonusgroup_descriptive.csv"))


# БЛОК 1 HasClaim

claim_m1 <- glm(
  HasClaim ~ logBonusMalus + logExposure,
  family = binomial(link = "logit"),
  data = policy_data
)

claim_m2 <- glm(
  HasClaim ~ logBonusMalus + AgeGroup + logExposure,
  family = binomial(link = "logit"),
  data = policy_data
)

claim_m3 <- glm(
  HasClaim ~ logBonusMalus + AgeGroup + VehAge + VehPower + VehGas + logExposure,
  family = binomial(link = "logit"),
  data = policy_data
)

claim_m4 <- glm(
  HasClaim ~ logBonusMalus + AgeGroup + VehAge + VehPower + VehGas +
    logDensity + Area + logExposure,
  family = binomial(link = "logit"),
  data = policy_data
)

claim_m5 <- glm(
  HasClaim ~ logBonusMalus + AgeGroup + VehAge + VehPower + VehGas +
    logDensity + Area + Region + logExposure,
  family = binomial(link = "logit"),
  data = policy_data
)

claim_m6 <- glm(
  HasClaim ~ logBonusMalus + AgeGroup + VehAge + VehPower + VehGas +
    logDensity + Area + Region + VehBrand + logExposure,
  family = binomial(link = "logit"),
  data = policy_data
)

claim_models <- list(
  M1_logBM_exposure = claim_m1,
  M2_driver = claim_m2,
  M3_vehicle = claim_m3,
  M4_geo_main = claim_m4,
  M5_regionFE = claim_m5,
  M6_brandFE = claim_m6
)

claim_coef <- purrr::imap_dfr(
  claim_models,
  function(model_obj, model_nm) {
    coef_table(
      model = model_obj,
      vcov_mat = vcov_glm(model_obj),
      model_name = model_nm,
      model_type = "logit"
    )
  }
)

claim_summary <- bind_rows(
  lapply(names(claim_models), function(nm) {
    m <- claim_models[[nm]]
    
    tibble(
      model = nm,
      n = length(stats::fitted(m)),
      logLik = as.numeric(stats::logLik(m)),
      AIC = stats::AIC(m),
      BIC = stats::BIC(m),
      McFadden_R2 = 1 - as.numeric(stats::logLik(m)) / as.numeric(stats::logLik(update(m, . ~ 1)))
    )
  })
)
write_csv(claim_coef, file.path(table_dir, "01_claim_probability_logit_coefficients.csv"))
write_csv(claim_summary, file.path(table_dir, "01_claim_probability_model_summary.csv"))

# Wald-тести для HasClaim
claim_wald <- bind_rows(
  wald_test_terms(claim_m6, "^AgeGroup", vcov_glm(claim_m6), "HasClaim: AgeGroup in M6"),
  wald_test_terms(claim_m6, "^Area", vcov_glm(claim_m6), "HasClaim: Area in M6"),
  wald_test_terms(claim_m6, "^Region", vcov_glm(claim_m6), "HasClaim: Region FE in M6"),
  wald_test_terms(claim_m6, "^VehBrand", vcov_glm(claim_m6), "HasClaim: VehBrand FE in M6")
)

write_csv(claim_wald, file.path(table_dir, "01_claim_probability_wald_tests.csv"))

# Прогнозована ймовірність claim за BonusGroup
policy_data <- policy_data %>%
  mutate(pred_claim_m6 = predict(claim_m6, type = "response"))

claim_bonus_plot_data <- policy_data %>%
  group_by(BonusGroup) %>%
  summarise(
    n = n(),
    observed_rate = mean(HasClaim),
    predicted_rate = mean(pred_claim_m6),
    .groups = "drop"
  )

write_csv(claim_bonus_plot_data, file.path(table_dir, "01_claim_probability_by_bonusgroup.csv"))

p_claim_bonus <- claim_bonus_plot_data %>%
  pivot_longer(cols = c(observed_rate, predicted_rate),
               names_to = "type", values_to = "rate") %>%
  ggplot(aes(x = BonusGroup, y = rate, fill = type)) +
  geom_col(position = position_dodge(width = 0.8)) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(
    title = "Ймовірність claim за групами Bonus-Malus",
    subtitle = "Observed vs predicted from logit M6",
    x = "BonusGroup",
    y = "Ймовірність claim",
    fill = ""
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

save_plot(p_claim_bonus, "01_claim_probability_by_bonusgroup.png")
table(policy_data$HasClaim, useNA = "ifany")
summary(policy_data$pred_claim_m6)
# Calibration for HasClaim
claim_calib <- calibration_table(policy_data, "HasClaim", "pred_claim_m6", bins = 10)
claim_roc <- pROC::roc(
  response = policy_data$HasClaim,
  predictor = policy_data$pred_claim_m6,
  levels = c(0, 1),
  direction = "<",
  quiet = TRUE
)

claim_auc <- as.numeric(pROC::auc(claim_roc))

print(paste("HasClaim AUC =", round(claim_auc, 4)))

write_csv(claim_calib, file.path(table_dir, "01_claim_probability_calibration.csv"))

p_claim_calib <- ggplot(claim_calib, aes(x = mean_pred, y = observed_rate)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.002) +
  geom_point(size = 2.5) +
  scale_x_continuous(labels = percent_format(accuracy = 0.1)) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(
    title = paste0("Calibration plot for HasClaim logit M6; AUC = ", round(claim_auc, 3)),
    x = "Середня прогнозована ймовірність",
    y = "Фактична частка claims"
  ) +
  theme_minimal()

save_plot(p_claim_calib, "01_claim_probability_calibration.png")

# БЛОК 2. Типова severity: log(CappedClaim)
sev_m1 <- lm(
  logY ~ BonusGroup,
  data = claims_data
)

sev_m2 <- lm(
  logY ~ logBonusMalus,
  data = claims_data
)

sev_m3 <- lm(
  logY ~ BonusGroup + DrivAge,
  data = claims_data
)

sev_m4 <- lm(
  logY ~ BonusGroup + DrivAge + VehAge + VehPower +
    logDensity + Area + VehGas,
  data = claims_data
)

sev_m5 <- lm(
  logY ~ BonusGroup + DrivAge + VehAge + VehPower +
    logDensity + Area + VehGas + Region,
  data = claims_data
)

sev_m6 <- lm(
  logY ~ BonusGroup + DrivAge + VehAge + VehPower +
    logDensity + Area + VehGas + Region + VehBrand,
  data = claims_data
)

sev_models <- list(
  M1_BonusGroup_only = sev_m1,
  M2_logBonusMalus = sev_m2,
  M3_BonusGroup_DrivAge = sev_m3,
  M4_Main = sev_m4,
  M5_RegionFE = sev_m5,
  M6_RegionFE_BrandFE = sev_m6
)

sev_coef <- purrr::imap_dfr(
  sev_models,
  function(model_obj, model_nm) {
    coef_table(
      model = model_obj,
      vcov_mat = vcov_ols(model_obj),
      model_name = model_nm,
      model_type = "ols_log"
    )
  }
)

sev_summary <- bind_rows(
  lapply(names(sev_models), function(nm) {
    m <- sev_models[[nm]]
    s <- summary(m)
    
    tibble(
      model = nm,
      n = length(stats::fitted(m)),
      r2 = s$r.squared,
      adj_r2 = s$adj.r.squared,
      AIC = stats::AIC(m),
      BIC = stats::BIC(m)
    )
  })
)

write_csv(sev_coef, file.path(table_dir, "02_typical_severity_ols_coefficients.csv"))
write_csv(sev_summary, file.path(table_dir, "02_typical_severity_model_summary.csv"))

# Wald-тести для severity
sev_wald <- bind_rows(
  wald_test_terms(sev_m4, "^BonusGroup", vcov_ols(sev_m4), "Severity M4: BonusGroup"),
  wald_test_terms(sev_m5, "^BonusGroup", vcov_ols(sev_m5), "Severity M5: BonusGroup"),
  wald_test_terms(sev_m6, "^BonusGroup", vcov_ols(sev_m6), "Severity M6: BonusGroup"),
  wald_test_terms(sev_m6, "^Area", vcov_ols(sev_m6), "Severity M6: Area"),
  wald_test_terms(sev_m6, "^Region", vcov_ols(sev_m6), "Severity M6: Region FE"),
  wald_test_terms(sev_m6, "^VehBrand", vcov_ols(sev_m6), "Severity M6: VehBrand FE")
)

write_csv(sev_wald, file.path(table_dir, "02_typical_severity_wald_tests.csv"))

# Скориговані прогнозні виплати за BonusGroup на M4
adjusted_predictions_bonus <- map_dfr(levels(claims_data$BonusGroup), function(bg) {
  nd <- claims_data
  nd$BonusGroup <- factor(bg, levels = levels(claims_data$BonusGroup))
  
  pred_log <- predict(sev_m4, newdata = nd)
  
  tibble(
    BonusGroup = bg,
    adjusted_log_prediction = mean(pred_log, na.rm = TRUE),
    adjusted_typical_claim = exp(mean(pred_log, na.rm = TRUE))
  )
})

write_csv(adjusted_predictions_bonus, file.path(table_dir, "02_adjusted_predictions_bonusgroup.csv"))
adjusted_predictions_bonus <- adjusted_predictions_bonus %>%
  mutate(
    BonusGroup = factor(
      BonusGroup,
      levels = c("50_80_Super", "81_100_Normal", "101_120_Risk", "121_200_Bad")
    )
  )
p_sev_adj <- ggplot(adjusted_predictions_bonus,
                    aes(x = BonusGroup, y = adjusted_typical_claim)) +
  geom_col() +
  labs(
    title = "Скоригована типова виплата за BonusGroup",
    subtitle = "На основі OLS M4 для log(CappedClaim)",
    x = "BonusGroup",
    y = "Скоригована типова виплата"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

save_plot(p_sev_adj, "02_adjusted_typical_claim_by_bonusgroup.png")

# Стабільність BonusGroup у severity M1/M3/M4/M5/M6
sev_stability <- sev_coef %>%
  filter(
    model %in% c(
      "M1_BonusGroup_only",
      "M3_BonusGroup_DrivAge",
      "M4_Main",
      "M5_RegionFE",
      "M6_RegionFE_BrandFE"
    ),
    str_detect(term, "^BonusGroup")
  ) %>%
  mutate(
    model_short = recode(
      model,
      "M1_BonusGroup_only" = "M1",
      "M3_BonusGroup_DrivAge" = "M3",
      "M4_Main" = "M4",
      "M5_RegionFE" = "M5",
      "M6_RegionFE_BrandFE" = "M6"
    ),
    model_short = factor(model_short, levels = c("M1", "M3", "M4", "M5", "M6"))
  )

write_csv(sev_stability, file.path(table_dir, "02_bonusgroup_stability.csv"))

p_sev_stability <- ggplot(
  sev_stability,
  aes(x = model_short, y = estimate, group = term)
) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 2.4) +
  geom_line() +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.12) +
  facet_wrap(~ term, nrow = 1) +
  labs(
    title = "Стійкість коефіцієнтів BonusGroup у severity-моделях",
    subtitle = "Планки = 95% robust HC3 CI",
    x = "Модель",
    y = "Коефіцієнт у log-моделі"
  ) +
  theme_minimal()

save_plot(p_sev_stability, "02_bonusgroup_stability.png", width = 10, height = 4.5)

# БЛОК 3. LargeClaim
large_m1 <- glm(
  LargeClaim ~ BonusGroup,
  family = binomial(link = "logit"),
  data = claims_data
)

large_m2 <- glm(
  LargeClaim ~ BonusGroup + DrivAge,
  family = binomial(link = "logit"),
  data = claims_data
)

large_m3 <- glm(
  LargeClaim ~ BonusGroup + DrivAge + VehAge + VehPower + VehGas,
  family = binomial(link = "logit"),
  data = claims_data
)

large_m4 <- glm(
  LargeClaim ~ BonusGroup + DrivAge + VehAge + VehPower + VehGas +
    logDensity + Area,
  family = binomial(link = "logit"),
  data = claims_data
)

large_m5 <- glm(
  LargeClaim ~ BonusGroup + DrivAge + VehAge + VehPower + VehGas +
    logDensity + Area + Region,
  family = binomial(link = "logit"),
  data = claims_data
)

large_m6 <- glm(
  LargeClaim ~ BonusGroup + DrivAge + VehAge + VehPower + VehGas +
    logDensity + Area + Region + VehBrand,
  family = binomial(link = "logit"),
  data = claims_data
)

large_models <- list(
  M1_BonusGroup_only = large_m1,
  M2_DrivAge = large_m2,
  M3_Vehicle = large_m3,
  M4_Geo_main = large_m4,
  M5_RegionFE = large_m5,
  M6_BrandFE = large_m6
)

large_coef <- purrr::imap_dfr(
  large_models,
  function(model_obj, model_nm) {
    coef_table(
      model = model_obj,
      vcov_mat = vcov_glm(model_obj),
      model_name = model_nm,
      model_type = "logit"
    )
  }
)

large_summary <- bind_rows(
  lapply(names(large_models), function(nm) {
    m <- large_models[[nm]]
    
    tibble(
      model = nm,
      n = length(stats::fitted(m)),
      logLik = as.numeric(stats::logLik(m)),
      AIC = stats::AIC(m),
      BIC = stats::BIC(m),
      McFadden_R2 = 1 - as.numeric(stats::logLik(m)) / as.numeric(stats::logLik(update(m, . ~ 1)))
    )
  })
)

write_csv(large_coef, file.path(table_dir, "03_largeclaim_logit_coefficients.csv"))
write_csv(large_summary, file.path(table_dir, "03_largeclaim_model_summary.csv"))

large_wald <- bind_rows(
  wald_test_terms(large_m4, "^BonusGroup", vcov_glm(large_m4), "LargeClaim M4: BonusGroup"),
  wald_test_terms(large_m5, "^BonusGroup", vcov_glm(large_m5), "LargeClaim M5: BonusGroup"),
  wald_test_terms(large_m6, "^BonusGroup", vcov_glm(large_m6), "LargeClaim M6: BonusGroup"),
  wald_test_terms(large_m6, "^Area", vcov_glm(large_m6), "LargeClaim M6: Area"),
  wald_test_terms(large_m6, "^Region", vcov_glm(large_m6), "LargeClaim M6: Region FE"),
  wald_test_terms(large_m6, "^VehBrand", vcov_glm(large_m6), "LargeClaim M6: VehBrand FE")
)

write_csv(large_wald, file.path(table_dir, "03_largeclaim_wald_tests.csv"))

# Прогнозована ймовірність LargeClaim за BonusGroup
claims_data <- claims_data %>%
  mutate(pred_large_m6 = predict(large_m6, type = "response"))

large_bonus_plot_data <- claims_data %>%
  group_by(BonusGroup) %>%
  summarise(
    n = n(),
    observed_large_rate = mean(LargeClaim),
    predicted_large_rate = mean(pred_large_m6),
    .groups = "drop"
  )

write_csv(large_bonus_plot_data, file.path(table_dir, "03_largeclaim_by_bonusgroup.csv"))

p_large_bonus <- large_bonus_plot_data %>%
  pivot_longer(cols = c(observed_large_rate, predicted_large_rate),
               names_to = "type", values_to = "rate") %>%
  ggplot(aes(x = BonusGroup, y = rate, fill = type)) +
  geom_col(position = position_dodge(width = 0.8)) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(
    title = "Ймовірність LargeClaim за BonusGroup",
    subtitle = "Observed vs predicted from logit M6",
    x = "BonusGroup",
    y = "Ймовірність великої виплати",
    fill = ""
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

save_plot(p_large_bonus, "03_largeclaim_probability_by_bonusgroup.png")

# Calibration for LargeClaim
large_calib <- calibration_table(claims_data, "LargeClaim", "pred_large_m6", bins = 10)
large_roc <- pROC::roc(
  response = claims_data$LargeClaim,
  predictor = claims_data$pred_large_m6,
  levels = c(0, 1),
  direction = "<",
  quiet = TRUE
)

large_auc <- as.numeric(pROC::auc(large_roc))

write_csv(large_calib, file.path(table_dir, "03_largeclaim_calibration.csv"))

p_large_calib <- ggplot(large_calib, aes(x = mean_pred, y = observed_rate)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.002) +
  geom_point(size = 2.5) +
  scale_x_continuous(labels = percent_format(accuracy = 0.1)) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(
    title = paste0("Calibration plot for LargeClaim logit M6; AUC = ", round(large_auc, 3)),
    x = "Середня прогнозована ймовірність",
    y = "Фактична частка LargeClaim"
  ) +
  theme_minimal()

save_plot(p_large_calib, "03_largeclaim_calibration.png")

# БЛОК 4.
# 4.1 Нелінійності для typical severity
sev_nonlin <- lm(
  logY ~ logBonusMalus_c + I(logBonusMalus_c^2) +
    DrivAge_c + I(DrivAge_c^2) +
    VehAge_c + I(VehAge_c^2) +
    VehPower_c + I(VehPower_c^2) +
    logDensity_c + I(logDensity_c^2) +
    Area + VehGas,
  data = claims_data
)

sev_nonlin_coef <- coef_table(
  sev_nonlin,
  vcov_ols(sev_nonlin),
  "Severity_nonlinear",
  model_type = "ols_log"
)

write_csv(sev_nonlin_coef, file.path(table_dir, "04_severity_nonlinear_terms.csv"))

sev_nonlin_wald <- bind_rows(
  wald_test_terms(sev_nonlin, "I\\(logBonusMalus_c\\^2\\)", vcov_ols(sev_nonlin), "Severity: logBonusMalus squared"),
  wald_test_terms(sev_nonlin, "I\\(DrivAge_c\\^2\\)", vcov_ols(sev_nonlin), "Severity: DrivAge squared"),
  wald_test_terms(sev_nonlin, "I\\(VehAge_c\\^2\\)", vcov_ols(sev_nonlin), "Severity: VehAge squared"),
  wald_test_terms(sev_nonlin, "I\\(VehPower_c\\^2\\)", vcov_ols(sev_nonlin), "Severity: VehPower squared"),
  wald_test_terms(sev_nonlin, "I\\(logDensity_c\\^2\\)", vcov_ols(sev_nonlin), "Severity: logDensity squared")
)

write_csv(sev_nonlin_wald, file.path(table_dir, "04_severity_nonlinear_wald_tests.csv"))

# 4.2 Взаємодії для typical severity
sev_int_bonus_age <- lm(
  logY ~ BonusGroup * AgeGroup + VehAge + VehPower + logDensity + Area + VehGas,
  data = claims_data
)

sev_int_power_agecar <- lm(
  logY ~ BonusGroup + DrivAge + VehAge * VehPower + logDensity + Area + VehGas,
  data = claims_data
)

sev_int_density_area <- lm(
  logY ~ BonusGroup + DrivAge + VehAge + VehPower + logDensity * Area + VehGas,
  data = claims_data
)

interaction_models <- list(
  Severity_BonusGroup_x_AgeGroup = sev_int_bonus_age,
  Severity_VehAge_x_VehPower = sev_int_power_agecar,
  Severity_logDensity_x_Area = sev_int_density_area
)

interaction_coef <- purrr::imap_dfr(
  interaction_models,
  function(model_obj, model_nm) {
    coef_table(
      model = model_obj,
      vcov_mat = vcov_ols(model_obj),
      model_name = model_nm,
      model_type = "ols_log"
    )
  }
)

write_csv(interaction_coef, file.path(table_dir, "04_severity_interaction_coefficients.csv"))

interaction_wald <- bind_rows(
  wald_test_terms(sev_int_bonus_age, ":", vcov_ols(sev_int_bonus_age), "Severity: BonusGroup x AgeGroup"),
  wald_test_terms(sev_int_power_agecar, ":", vcov_ols(sev_int_power_agecar), "Severity: VehAge x VehPower"),
  wald_test_terms(sev_int_density_area, ":", vcov_ols(sev_int_density_area), "Severity: logDensity x Area")
)

write_csv(interaction_wald, file.path(table_dir, "04_severity_interaction_wald_tests.csv"))

# 4.3 Діагностика OLS severity-моделі

claims_data <- claims_data %>%
  mutate(
    sev_m4_fitted = fitted(sev_m4),
    sev_m4_resid = resid(sev_m4)
  )

# Residuals vs fitted
p_resid_scatter <- ggplot(claims_data, aes(x = sev_m4_fitted, y = sev_m4_resid)) +
  geom_point(alpha = 0.08, size = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_smooth(method = "loess", se = FALSE) +
  labs(
    title = "Residuals vs fitted values для severity M4",
    subtitle = "Ознаки нерівномірної дисперсії обґрунтовують HC3",
    x = "Fitted log(CappedClaim)",
    y = "Residual"
  ) +
  theme_minimal()

save_plot(p_resid_scatter, "04_severity_residuals_vs_fitted.png")

# Binned residuals
binned_residuals <- claims_data %>%
  mutate(bin = ntile(sev_m4_fitted, 20)) %>%
  group_by(bin) %>%
  summarise(
    n = n(),
    fitted_mid = mean(sev_m4_fitted),
    mean_resid = mean(sev_m4_resid),
    se_resid = sd(sev_m4_resid) / sqrt(n),
    ci_low = mean_resid - 1.96 * se_resid,
    ci_high = mean_resid + 1.96 * se_resid,
    .groups = "drop"
  )

write_csv(binned_residuals, file.path(table_dir, "04_severity_binned_residuals.csv"))
claim_stability <- claim_coef %>%
  filter(term == "logBonusMalus") %>%
  select(model, term, estimate, std_error, p_value, odds_ratio, or_ci_low, or_ci_high)

write_csv(claim_stability, file.path(table_dir, "01_claim_logbonusmalus_stability.csv"))
large_bonus_table <- large_coef %>%
  filter(str_detect(term, "^BonusGroup")) %>%
  select(model, term, estimate, std_error, p_value, odds_ratio, or_ci_low, or_ci_high)

write_csv(large_bonus_table, file.path(table_dir, "03_largeclaim_bonusgroup_coefficients.csv"))
large_observed_table <- claims_data %>%
  group_by(BonusGroup) %>%
  summarise(
    n = n(),
    observed_large_rate = mean(LargeClaim),
    .groups = "drop"
  )

write_csv(large_observed_table, file.path(table_dir, "03_largeclaim_observed_by_bonusgroup.csv"))
p_binned_resid <- ggplot(binned_residuals, aes(x = fitted_mid, y = mean_resid)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.03) +
  geom_point(size = 2.5) +
  labs(
    title = "Binned residuals для severity M4",
    subtitle = "Середні залишки за 20 групами fitted values",
    x = "Середнє fitted log(CappedClaim)",
    y = "Середній residual"
  ) +
  theme_minimal()

save_plot(p_binned_resid, "04_severity_binned_residuals.png")

# Breusch-Pagan test як додатковий формальний тест
bp_test <- lmtest::bptest(sev_m4)

bp_table <- tibble(
  test = "Breusch-Pagan for severity M4",
  statistic = as.numeric(bp_test$statistic),
  df = as.numeric(bp_test$parameter),
  p_value = as.numeric(bp_test$p.value)
)

write_csv(bp_table, file.path(table_dir, "04_severity_breusch_pagan.csv"))

# 4.4 Порівняльні ключові результати

# Ключові коефіцієнти для frequency
claim_key <- claim_coef %>%
  filter(
    model == "M6_brandFE",
    term %in% c("logBonusMalus", "AgeGroupYouth_18_25", "AgeGroupSeniors_60plus") |
      str_detect(term, "^Area")
  ) %>%
  mutate(block = "Claim probability")

# Ключові коефіцієнти для typical severity
sev_key <- sev_coef %>%
  filter(
    model == "M4_Main",
    str_detect(term, "^BonusGroup") |
      term %in% c("DrivAge", "VehAge", "VehPower", "logDensity")
  ) %>%
  mutate(block = "Typical severity")

# Ключові коефіцієнти для LargeClaim
large_key <- large_coef %>%
  filter(
    model == "M6_BrandFE",
    str_detect(term, "^BonusGroup") |
      term %in% c("DrivAge", "VehAge", "VehPower", "logDensity")
  ) %>%
  mutate(block = "LargeClaim probability")

key_results <- bind_rows(
  claim_key,
  sev_key,
  large_key
)

write_csv(key_results, file.path(table_dir, "05_key_results_for_report.csv"))

# Порівняльна таблиця по головних блоках
comparison_table <- tibble(
  factor = c("Bonus-Malus", "Driver age", "Vehicle power / age", "Area / Density", "VehBrand / Region"),
  claim_probability = c(
    "Перевірити за logit HasClaim: logBM / BonusGroup",
    "Перевірити AgeGroup / DrivAge",
    "Перевірити VehAge, VehPower",
    "Перевірити Area, logDensity",
    "Перевірити Region FE, VehBrand FE"
  ),
  typical_severity = c(
    "Перевірити BonusGroup у OLS log(CappedClaim)",
    "Перевірити DrivAge",
    "Перевірити VehAge, VehPower",
    "Перевірити Area, logDensity",
    "Перевірити Region FE, VehBrand FE"
  ),
  largeclaim_probability = c(
    "Перевірити BonusGroup у logit LargeClaim",
    "Перевірити DrivAge",
    "Перевірити VehAge, VehPower",
    "Перевірити Area, logDensity",
    "Перевірити Region FE, VehBrand FE"
  ),
  final_comment_template = c(
    "Якщо сильний у HasClaim і severity, але слабкий у LargeClaim: працює для frequency / typical risk, не для хвоста.",
    "Якщо змінюється після Bonus-Malus: є OVB / досвід відділяється від віку.",
    "Якщо слабкий у frequency, але сильний у severity: технічні фактори більше про розмір збитку.",
    "Якщо сильний у frequency, але слабкий у severity: простір більше про ймовірність ДТП.",
    "Якщо FE змінюють коефіцієнти: частина ефекту пояснюється регіоном/брендом."
  )
)

write_csv(comparison_table, file.path(table_dir, "05_final_comparison_template.csv"))


cat("\n================ DATA INFO ================\n")
print(data_info)

cat("\n================ Q THRESHOLDS ================\n")
cat("Q95 TotalClaim:", as.numeric(q95), "\n")
cat("Q99.5 TotalClaim:", as.numeric(q995), "\n")

cat("\n================ AUC ================\n")
cat("HasClaim AUC:", round(claim_auc, 4), "\n")
cat("LargeClaim AUC:", round(large_auc, 4), "\n")

cat("\n================ OUTPUTS SAVED ================\n")
cat("Tables:", table_dir, "\n")
cat("Plots:", plot_dir, "\n")
