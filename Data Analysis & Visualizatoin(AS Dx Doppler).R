
source(file = "Data Upload and Cleaning.R")
# ── 6. ADDITIONAL LIBRARIES ─────────────────────────────────
library(gtsummary)   # publication-ready summary tables
library(patchwork)   # combine ggplots
library(rstatix)     # pipe-friendly statistical tests
library(janitor)     # clean_names() for snake_case column names
library(blandr)      # Bland-Altman agreement analysis
library(pROC)        # ROC curve analysis
library(flextable)   # export tables to Word/HTML

# ── 7. BASELINE CHARACTERISTICS TABLE ───────────────────────
# Demographics and comorbidities summary (Table 1)
DVI_for_AS_Data %>%
  select(Age, Sex, BMI, HTN, Dyslipidemia, DM,
         `Coronary Artery Disease`, Stroke, CKD) %>%
  tbl_summary(
    missing = "no",
    statistic = Age ~ "{mean} ({sd})"   # Age is continuous → mean (SD)
  ) %>%
  bold_labels() %>%
  as_flex_table()

# Echocardiographic parameters summary (Table 2)
DVI_for_AS_Data %>%
  select(
    `Left Ventricular Outflow Tract Diameter (cm)`,
    `Aortic Valve Mean Pressure Gradient (mmHg)`,
    `Aortic Valve Velocity Time Integral (cm)`,
    `LVOT Velocity Time Integral (cm)`,
    `Mitral Valve Continuous Wave Velocity Time Integral (cm)`,
    `Doppler Velocity Index using LVOT`,
    `Doppler Velocity Index using Mitral Valve`,
    `Aortic Valve Area (cm2)`,
    `Ejection Fraction`
  ) %>%
  tbl_summary(
    # Normally distributed variables → mean (SD); rest default to median (IQR)
    statistic = list(
      `Left Ventricular Outflow Tract Diameter (cm)` ~ "{mean} ({sd})",
      `LVOT Velocity Time Integral (cm)`             ~ "{mean} ({sd})",
      `Doppler Velocity Index using LVOT`            ~ "{mean} ({sd})"
    )
  ) %>%
  bold_labels() %>%
  as_flex_table()

# ── 8. NORMALITY TESTING ─────────────────────────────────────
# Convert to snake_case for easier piping in rstatix/ggplot
Clean_Dvi_For_As_Data <- DVI_for_AS_Data %>% clean_names()

# Helper: plots histogram, density, Q-Q plot, then runs Shapiro-Wilk test
Normality_Testing <- function(data, x) {
  p1 <- data %>% ggplot(aes(x = {{ x }})) + geom_histogram()
  p2 <- data %>% ggplot(aes(x = {{ x }})) + geom_density()
  p3 <- data %>% ggplot(aes(sample = {{ x }})) + geom_qq()
  print(p1 / p2 / p3)
  data %>% shapiro_test({{ x }})
}

# Run normality checks on all continuous variables
Clean_Dvi_For_As_Data %>% Normality_Testing(age)
Clean_Dvi_For_As_Data %>% Normality_Testing(bmi)
Clean_Dvi_For_As_Data %>% Normality_Testing(left_ventricular_outflow_tract_diameter_cm)
Clean_Dvi_For_As_Data %>% Normality_Testing(lvot_velocity_time_integral_cm)
Clean_Dvi_For_As_Data %>% Normality_Testing(aortic_valve_mean_pressure_gradient_mm_hg)
Clean_Dvi_For_As_Data %>% Normality_Testing(aortic_valve_velocity_time_integral_cm)
Clean_Dvi_For_As_Data %>% Normality_Testing(mitral_valve_continuous_wave_velocity_time_integral_cm)
Clean_Dvi_For_As_Data %>% Normality_Testing(doppler_velocity_index_using_lvot)
Clean_Dvi_For_As_Data %>% Normality_Testing(doppler_velocity_index_using_mitral_valve)
Clean_Dvi_For_As_Data %>% Normality_Testing(aortic_valve_area_cm2)
Clean_Dvi_For_As_Data %>% Normality_Testing(ejection_fraction)

# ── 9. AVA CATEGORIZATION ────────────────────────────────────
# Categorize AVA into severity groups per ACC/AHA guidelines:
#   Severe: AVA ≤ 1.0 cm²  |  Moderate: 1.0–1.5 cm²  |  Mild: 1.5–2.0 cm²
DVI_for_AS_Data <- DVI_for_AS_Data %>%
  mutate(AVA_cat = cut(
    `Aortic Valve Area (cm2)`,
    breaks = c(-Inf, 1.0, 1.5, 2.0),
    labels = c("≤1.0", ">1.0–≤1.5", ">1.5–≤2.0"),
    right  = TRUE
  ))

# ── 10. AVA DISTRIBUTION PLOT ────────────────────────────────
DVI_for_AS_Data %>%
  filter(!is.na(AVA_cat)) %>%
  mutate(AVA_cat = factor(AVA_cat,
                          levels = c("≤1.0", ">1.0–≤1.5", ">1.5–≤2.0"),
                          labels = c("Severe AS\n(AVA ≤ 1.0 cm²)",
                                     "Moderate AS\n(AVA 1.0–1.5 cm²)",
                                     "Mild AS\n(AVA 1.5–2.0 cm²)"))) %>%
  ggplot(aes(x = AVA_cat)) +
  geom_bar(aes(fill = AVA_cat), show.legend = FALSE, width = 0.55) +
  geom_text(stat = "count", aes(label = after_stat(count)),
            vjust = -0.6, fontface = "bold", size = 4.5) +
  scale_fill_manual(values = c("#1E3A8A", "#2563EB", "#4E9AF1")) +  # dark→light = severe→mild
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Distribution of Aortic Valve Area Categories",
    x     = "Aortic Valve Area Category",
    y     = "Number of Patients"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title          = element_text(face = "bold", hjust = 0.5, size = 14,
                                       margin = margin(b = 12)),
    axis.text           = element_text(color = "black", size = 11),
    axis.title.x        = element_text(face = "bold", margin = margin(t = 10)),
    axis.title.y        = element_text(face = "bold", margin = margin(r = 10)),
    axis.line           = element_line(color = "black"),
    panel.grid.major.y  = element_line(color = "grey90", linetype = "dashed")
  )

# ── 11. CORRELATION ANALYSIS ─────────────────────────────────
# Pearson correlation: LVOT VTI vs Mitral Valve CW VTI
cor.test(
  DVI_for_AS_Data$`LVOT Velocity Time Integral (cm)`,
  DVI_for_AS_Data$`Mitral Valve Continuous Wave Velocity Time Integral (cm)`,
  method = "pearson"
)

# Pearson correlation: DVI-LVOT vs DVI-MV (primary agreement test)
cor.test(
  DVI_for_AS_Data$`Doppler Velocity Index using LVOT`,
  DVI_for_AS_Data$`Doppler Velocity Index using Mitral Valve`,
  method = "pearson"
)

# ── 12. BLAND-ALTMAN AGREEMENT ANALYSIS ──────────────────────
# Assesses agreement between DVI-LVOT and DVI-MV (method comparison)
stats_ba <- blandr.statistics(
  DVI_for_AS_Data$`Doppler Velocity Index using LVOT`,
  DVI_for_AS_Data$`Doppler Velocity Index using Mitral Valve`
)

# Manually compute BA components for custom ggplot
ba <- DVI_for_AS_Data %>%
  select(`Doppler Velocity Index using LVOT`,
         `Doppler Velocity Index using Mitral Valve`) %>%
  filter(complete.cases(.)) %>%
  mutate(
    mean = (`Doppler Velocity Index using LVOT` +
              `Doppler Velocity Index using Mitral Valve`) / 2,
    diff =  `Doppler Velocity Index using LVOT` -
      `Doppler Velocity Index using Mitral Valve`
  )

# Compute bias and 95% limits of agreement (LoA = bias ± 1.96 SD)
bias      <- mean(ba$diff)
sd_d      <- sd(ba$diff)
loa_upper <- bias + 1.96 * sd_d
loa_lower <- bias - 1.96 * sd_d

# Plot Bland-Altman
ggplot(ba, aes(x = mean, y = diff)) +
  geom_point(color = "#0B2C6B", alpha = 0.65, size = 2) +
  geom_hline(yintercept = bias,      color = "#0B2C6B", linewidth = 1.1) +
  geom_hline(yintercept = loa_upper, color = "#C1121F", linetype = 2, linewidth = 0.9) +
  geom_hline(yintercept = loa_lower, color = "#C1121F", linetype = 2, linewidth = 0.9) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", color = "black", hjust = 0.5),
    axis.title       = element_text(face = "bold", color = "black"),
    axis.text        = element_text(color = "black"),
    panel.grid.major = element_line(color = "#E6ECF5", linewidth = 0.4),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x     = "Mean of DVI-LVOT and DVI-MV",
    y     = "Difference (DVI-LVOT − DVI-MV)",
    title = "Bland–Altman Plot: DVI-LVOT vs DVI-MV"
  )

# ── 13. ROC ANALYSIS ─────────────────────────────────────────
# Create binary outcome variables for two severity definitions:
#   severe_AS_AVA: TRUE if AVA > 1 cm² (non-severe by AVA = reference standard)
#   severe_AS_DVI: TRUE if DVI-LVOT < 0.25 (severe by LVOT criterion)
# NOTE: direction is set to "auto" — pROC will determine whether
#       higher or lower DVI-MV values correspond to the case group.
#       Confirmed directions after fitting:
#         roc_ava_for_as → controls < cases (higher DVI-MV = non-severe AVA)
#         roc_dvi_for_as → controls > cases (lower DVI-MV = severe by LVOT)
DVI_for_roc <- DVI_for_AS_Data %>%
  mutate(
    severe_AS_AVA = `Aortic Valve Area (cm2)` > 1,
    severe_AS_DVI = `Doppler Velocity Index using LVOT` < 0.25
  ) %>%
  select(`Doppler Velocity Index using LVOT`, `Aortic Valve Area (cm2)`,
         `Doppler Velocity Index using Mitral Valve`,
         severe_AS_AVA, severe_AS_DVI)

# ROC 1: DVI-MV predicting non-severe AS by AVA (AVA > 1 cm²)
roc_ava_for_as <- roc(
  response  = DVI_for_roc$severe_AS_AVA,
  predictor = DVI_for_roc$`Doppler Velocity Index using Mitral Valve`,
  levels    = c(FALSE, TRUE),
  direction = "auto"
)

# ROC 2: DVI-MV predicting severe AS by DVI-LVOT (< 0.25)
roc_dvi_for_as <- roc(
  response  = DVI_for_roc$severe_AS_DVI,
  predictor = DVI_for_roc$`Doppler Velocity Index using Mitral Valve`,
  levels    = c(FALSE, TRUE),
  direction = "auto"
)

# Area Under the Curve for each model
auc(roc_ava_for_as)
auc(roc_dvi_for_as)

# Optimal cut-off via Youden Index (maximizes sensitivity + specificity)
coords(roc_ava_for_as, x = "best",
       ret = c("threshold", "sensitivity", "specificity", "ppv", "npv"),
       best.method = "youden")

coords(roc_dvi_for_as, x = "best",
       ret = c("threshold", "sensitivity", "specificity", "ppv", "npv"),
       best.method = "youden")

# DeLong test: compares the two AUCs for statistical significance
roc.test(roc_ava_for_as, roc_dvi_for_as)

# ── 14. ROC CURVE PLOT ───────────────────────────────────────
ggroc(
  list(
    "Aortic Valve Area < 1 cm²"          = roc_ava_for_as,
    "Doppler Velocity Index LVOT < 0.25" = roc_dvi_for_as
  ),
  legacy.axes = TRUE,  # x-axis = 1 - Specificity (conventional orientation)
  linewidth   = 0.9
) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "grey60", linewidth = 0.5) +
  scale_color_manual(values = c("#2563EB", "#DC2626")) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    x     = "1 \u2212 Specificity",
    y     = "Sensitivity",
    color = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text            = element_text(color = "black", size = 11),
    axis.title           = element_text(face = "bold", size = 12),
    axis.title.x         = element_text(margin = margin(t = 10)),
    axis.title.y         = element_text(margin = margin(r = 10)),
    legend.position      = c(0.97, 0.08),
    legend.justification = c("right", "bottom"),
    legend.text          = element_text(size = 10.5),
    legend.key.width     = unit(1.4, "cm")
  )

# ── 15. DIAGNOSTIC PERFORMANCE TABLE ────────────────────────
# Compile cut-off, sensitivity, specificity, PPV, NPV for both models
tab <- bind_rows(
  coords(roc_ava_for_as, x = "best", best.method = "youden",
         ret = c("threshold", "sensitivity", "specificity", "ppv", "npv"),
         transpose = FALSE) %>%
    mutate(Outcome = "Severe AS by AVA < 1 cm²", .before = 1),
  
  coords(roc_dvi_for_as, x = "best", best.method = "youden",
         ret = c("threshold", "sensitivity", "specificity", "ppv", "npv"),
         transpose = FALSE) %>%
    mutate(Outcome = "Severe AS by DVI_LVOT < 0.25", .before = 1)
) %>%
  # Convert proportions to percentages for display
  mutate(across(c(sensitivity, specificity, ppv, npv), ~ . * 100))

# Format and export as publication-ready flextable
flextable(tab) %>%
  colformat_double(j = "threshold",
                   digits = 3) %>%
  colformat_double(j = c("sensitivity", "specificity", "ppv", "npv"),
                   digits = 1) %>%
  set_header_labels(
    Outcome     = "Outcome definition",
    threshold   = "Optimal cut-off",
    sensitivity = "Sensitivity (%)",
    specificity = "Specificity (%)",
    ppv         = "PPV (%)",
    npv         = "NPV (%)"
  ) %>%
  theme_vanilla() %>%
  autofit()





