# ============================================================
# Project: Doppler Velocity Index (DVI) for Aortic Stenosis
# Description: Evaluates DVI using Mitral Valve as a predictor
#              of severe AS, compared against AVA and LVOT-DVI
# ============================================================

# ── 1. DEPENDENCIES ─────────────────────────────────────────
library(tidyverse)
library(readxl)

# ── 2. DATA IMPORT ──────────────────────────────────────────
DVI_for_AS_Data <- read_excel(path = "Project 10 data.xlsx")

# ── 3. VARIABLE RENAMING ────────────────────────────────────
# Replace short/coded column names with full clinical labels
DVI_for_AS_Data <- DVI_for_AS_Data %>% rename(
  "Left Ventricular Outflow Tract Diameter (cm)"           = `LVOT diameter`,
  "Aortic Valve Mean Pressure Gradient (mmHg)"             = `AV mean pressure gradient`,
  "Aortic Valve Velocity Time Integral (cm)"               = `AV VTI`,
  "LVOT Velocity Time Integral (cm)"                       = `LVOT VTI`,
  "Mitral Valve Continuous Wave Velocity Time Integral (cm)" = `MV CW VTI`,
  "Doppler Velocity Index using LVOT"                      = `DVI LVOT`,
  "Doppler Velocity Index using Mitral Valve"              = `DVI MV`,
  "Aortic Valve Area (cm2)"                                = AVA,
  "Dyslipidemia"                                           = dyslipidemia,
  "Coronary Artery Disease"                                = CAD,
  "Ejection Fraction"                                      = EF,
  "Stroke"                                                 = stroke
)

# ── 4. REMOVE IDENTIFIERS ───────────────────────────────────
# Drop patient identifiers before any analysis (de-identification)
DVI_for_AS_Data <- DVI_for_AS_Data %>% select(-MRN, -name)

# ── 5. FACTOR ENCODING ──────────────────────────────────────
# Encode binary comorbidities as labeled factors (0 = No, 1 = Yes)
DVI_for_AS_Data <- DVI_for_AS_Data %>%
  mutate(across(
    c(HTN, Dyslipidemia, DM, `Coronary Artery Disease`, Stroke, CKD),
    ~ factor(.x, levels = c("0", "1"), labels = c("No", "Yes"))
  ))

# Encode sex as labeled factor (0 = Male, 1 = Female)
DVI_for_AS_Data <- DVI_for_AS_Data %>%
  mutate(Sex = factor(Sex, levels = c("0", "1"), labels = c("Male", "Female")))





