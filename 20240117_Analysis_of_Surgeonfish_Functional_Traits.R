#### Analysis of Surgeonfish Functional Traits ####

# ======================================
# ENVIRONMENT SETUP USING RENV
# ======================================

# This section sets up a reproducible environment using 'renv'.
# If you're setting up this project for the first time, uncomment and run the necessary sections.

# Load the renv package
library(renv)  

# --- Initialization (only run once) ---
# This initializes a project-specific library for package management.

# renv::init()  

# Activate the renv environment
renv::activate()

# --- Load Libraries ---
library(lme4)
library(lmerTest)
library(ggplot2)
library(dplyr)
library(scales)
library(forcats)
library(patchwork)
library(car)
library(Matrix)
library(knitr)

# --- Snapshot (only run once) ---
# This saves the package info for reproducibility. 
# Uncomment below if adding new packages or setting up for the first time.
#renv::snapshot()

### Data download and cleaning ###

# The dataset is fully citeable as Lilkendey, J et al.: Foraging traits, lengths and 3D movement trajectories of 
# coral reef fishes obtained via stereo-video in Eilat, Gulf of Aqaba, Red Sea. 
# PANGAEA, https://doi.pangaea.de/10.1594/PANGAEA.957634 

## The Gulf_of_Aqaba_Acanthuridae dataset is available at: https://doi.pangaea.de/10.1594/PANGAEA.957631

# If the dataset is already in the workspace
data <- Gulf_of_Aqaba_Acanthuridae

head(Gulf_of_Aqaba_Acanthuridae)

# Removing unnecessary column and rows (Only keep the two focus species and the sampling stations at the IUI Housereef)
data <- data %>%
  dplyr::filter(`Species UID` != "Ctenochaetus striatus" & !Event %in% c("Gulf_of_Aqaba_NR1", "Gulf_of_Aqaba_NR2"))

# Rename the columns 'Fish wm [g] (Derived from the measured len...)' to 'Fish_Mass'
# and 'Fish TL [mm] (Elaborated after Neuswanger e...)' to 'Fish_Length'
Gulf_of_Aqaba_Acanthuridae <- Gulf_of_Aqaba_Acanthuridae %>%
  rename(Fish_Mass = `Fish wm [g] (Derived from the measured len...)`,
         Fish_Length = `Fish TL [mm] (Elaborated after Neuswanger e...)`)

# Checking unique values in the "Event" column after filtering
unique_events_after_filter <- unique(data$Event)
print(unique_events_after_filter)

### Bite Distance ###

# Create the jitter plot for Bite Distance
plot_bite_distance <- ggplot(data, aes(x = factor(`Species UID`), y = `Bite dist [mm] (Elaborated after Neuswanger e...)`)) +
  geom_jitter(position = position_jitter(width = 0.2), color = "blue", size = 3) +
  labs(title = "Jitter Plot of Bite Distances between Species",
       x = "Species",
       y = "Bite distance (mm)") +
  theme_minimal()

print(plot_bite_distance)

# Outlier handling
IQR <- IQR(data$`Bite dist [mm] (Elaborated after Neuswanger e...)`, na.rm = TRUE)
Q1 <- quantile(data$`Bite dist [mm] (Elaborated after Neuswanger e...)`, 0.25, na.rm = TRUE)
Q3 <- quantile(data$`Bite dist [mm] (Elaborated after Neuswanger e...)`, 0.75, na.rm = TRUE)
lower_fence <- Q1 - 1.5 * IQR
upper_fence <- Q3 + 1.5 * IQR
data_no_outliers_bite_distance <- data[data$`Bite dist [mm] (Elaborated after Neuswanger e...)` > lower_fence & data$`Bite dist [mm] (Elaborated after Neuswanger e...)` < upper_fence, ]

# Create the jitter plot for Bite Distance
plot_bite_distance <- ggplot(data_no_outliers, aes(x = factor(`Species UID`), y = `Bite dist [mm] (Elaborated after Neuswanger e...)`)) +
  geom_jitter(position = position_jitter(width = 0.2), color = "blue", size = 3) +
  labs(title = "Jitter Plot of Bite Distances between Species",
       x = "Species",
       y = "Bite distance (mm)") +
  theme_minimal()

print(plot_bite_distance)

# Log transformation of Bite Distances
data_no_outliers_bite_distance$`Log_Bite_dist_mm` <- log(data_no_outliers_bite_distance$`Bite dist [mm] (Elaborated after Neuswanger e...)`)

# Density Plot for Log-Transformed Bite Distances
density_plot_dist_with_log <- ggplot(data_no_outliers_bite_distance, aes(x = `Log_Bite_dist_mm`, fill = `Species UID`)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density Plot of Log-Transformed Bite Distances between Species",
       x = "Log-Transformed Bite distance [mm]") +
  theme_minimal()

print(density_plot_dist_with_log)


# Fit models for Log-transformed Bite Distance with and without Fish Mass

# Models including Fish Mass
bite_dist_model_mass_nested <- lmer(Log_Bite_dist_mm ~ `Species UID` + Fish_Mass + 
                                      (1 | `ID (Fish ID)`/`ID (Quadrat ID)`), data = data_no_outliers_bite_distance)

bite_dist_model_mass_fishID <- lmer(Log_Bite_dist_mm ~ `Species UID` + Fish_Mass + 
                                      (1 | `ID (Fish ID)`), data = data_no_outliers_bite_distance)

bite_dist_model_mass_noRE <- lm(Log_Bite_dist_mm ~ `Species UID` + Fish_Mass, 
                                data = data_no_outliers_bite_distance)

# Models excluding Fish Mass
bite_dist_model_no_mass_nested <- lmer(Log_Bite_dist_mm ~ `Species UID` + 
                                         (1 | `ID (Fish ID)`/`ID (Quadrat ID)`), data = data_no_outliers_bite_distance)

bite_dist_model_no_mass_fishID <- lmer(Log_Bite_dist_mm ~ `Species UID` + 
                                         (1 | `ID (Fish ID)`), data = data_no_outliers_bite_distance)

bite_dist_model_no_mass_noRE <- lm(Log_Bite_dist_mm ~ `Species UID`, 
                                   data = data_no_outliers_bite_distance)


# Model comparison using AIC
AIC_table <- data.frame(
  Model = c("Bite Distance with Mass Nested", "Bite Distance without Mass Nested", 
            "Bite Distance with Mass FishID", "Bite Distance without Mass FishID", 
            "Bite Distance with Mass NoRE", "Bite Distance without Mass NoRE"),
  AIC = c(AIC(bite_dist_model_mass_nested), AIC(bite_dist_model_no_mass_nested), 
          AIC(bite_dist_model_mass_fishID), AIC(bite_dist_model_no_mass_fishID), 
          AIC(bite_dist_model_mass_noRE), AIC(bite_dist_model_no_mass_noRE))
)

print(AIC_table)

# Print model summaries including mass
print(summary(bite_dist_model_mass_nested))
print(summary(bite_dist_model_mass_fishID))
print(summary(bite_dist_model_mass_noRE))

# Summary of model with lowest AIC
print(summary(bite_dist_model_no_mass_fishID))
print(summary(bite_dist_model_no_mass_noRE))

# Test model assumptions
plot(fitted(bite_dist_model_no_mass_noRE), resid(bite_dist_model_no_mass_noRE))
abline(h=0, col="red")
title("Residuals vs Fitted for Bite Distance Model")

qqnorm(resid(bite_dist_model_no_mass_noRE))
qqline(resid(bite_dist_model_no_mass_noRE), col = "red")
title("Normal Q-Q Plot for Bite Distance Model")


### Bite Rates ###

# Outlier handling for Bite Rates
IQR_rate <- IQR(data$`Bite rate [#/min] (Elaborated after Neuswanger e...)`, na.rm = TRUE)
Q1_rate <- quantile(data$`Bite rate [#/min] (Elaborated after Neuswanger e...)`, 0.25, na.rm = TRUE)
Q3_rate <- quantile(data$`Bite rate [#/min] (Elaborated after Neuswanger e...)`, 0.75, na.rm = TRUE)
lower_fence_rate <- Q1_rate - 1.5 * IQR_rate
upper_fence_rate <- Q3_rate + 1.5 * IQR_rate
data_no_outliers_bite_rate <- data[data$`Bite rate [#/min] (Elaborated after Neuswanger e...)` > lower_fence_rate & data$`Bite rate [#/min] (Elaborated after Neuswanger e...)` < upper_fence_rate, ]


# Density Plot for Bite Rates
density_plot_rates_no_log <- ggplot(data_no_outliers_bite_rate, aes(x = `Bite rate [#/min] (Elaborated after Neuswanger e...)`, fill = `Species UID`)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density Plot of Bite Rates between Species (Without Log Transformation)",
       x = "Bite rate [#/min]") +
  theme_minimal()

print(density_plot_rates_no_log)

# Log transformation of Bite Rates
data_no_outliers_bite_rate$`Log_Bite_rate_per_min` <- log(data_no_outliers_bite_rate$`Bite rate [#/min] (Elaborated after Neuswanger e...)`)

# Density Plot for Log-Transformed Bite Rates
density_plot_rates_with_log <- ggplot(data_no_outliers_bite_rate, aes(x = `Log_Bite_rate_per_min`, fill = `Species UID`)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density Plot of Log-Transformed Bite Rates between Species",
       x = "Log-Transformed Bite rate [#/min]") +
  theme_minimal()

print(density_plot_rates_with_log)

# Fit models for Log-transformed Bite Rates with and without Fish Mass

# Models including Fish Mass
bite_rate_model_mass_nested <- lmer(Log_Bite_rate_per_min ~ `Species UID` + Fish_Mass + 
                                      (1 | `ID (Fish ID)`/`ID (Quadrat ID)`), data = data_no_outliers_bite_rate)

bite_rate_model_mass_fishID <- lmer(Log_Bite_rate_per_min ~ `Species UID` + Fish_Mass + 
                                      (1 | `ID (Fish ID)`), data = data_no_outliers_bite_rate)

bite_rate_model_mass_noRE <- lm(Log_Bite_rate_per_min ~ `Species UID` + Fish_Mass, 
                                data = data_no_outliers_bite_rate)

# Models excluding Fish Mass
bite_rate_model_no_mass_nested <- lmer(Log_Bite_rate_per_min ~ `Species UID` + 
                                         (1 | `ID (Fish ID)`/`ID (Quadrat ID)`), data = data_no_outliers_bite_rate)

bite_rate_model_no_mass_fishID <- lmer(Log_Bite_rate_per_min ~ `Species UID` + 
                                         (1 | `ID (Fish ID)`), data = data_no_outliers_bite_rate)

bite_rate_model_no_mass_noRE <- lm(Log_Bite_rate_per_min ~ `Species UID`, 
                                   data = data_no_outliers_bite_rate)

# Model comparison using AIC
AIC_table <- data.frame(
  Model = c("Bite Rate with Mass Nested", "Bite Rate without Mass Nested", 
            "Bite Rate with Mass FishID", "Bite Rate without Mass FishID", 
            "Bite Rate with Mass NoRE", "Bite Rate without Mass NoRE"),
  AIC = c(AIC(bite_rate_model_mass_nested), AIC(bite_rate_model_no_mass_nested), 
          AIC(bite_rate_model_mass_fishID), AIC(bite_rate_model_no_mass_fishID), 
          AIC(bite_rate_model_mass_noRE), AIC(bite_rate_model_no_mass_noRE))
)

print(AIC_table)

# Print model summaries including mass
print(summary(bite_rate_model_mass_nested))
print(summary(bite_rate_model_mass_fishID))
print(summary(bite_rate_model_mass_noRE))

# Model summary with the lowest AIC
print(summary(bite_rate_model_no_mass_noRE))

# Test assumptions
plot(fitted(bite_rate_model_no_mass_noRE), resid(bite_rate_model_no_mass_noRE))
abline(h=0, col="red")
title("Residuals vs Fitted")

qqnorm(resid(bite_rate_model_no_mass_noRE))
qqline(resid(bite_rate_model_no_mass_noRE), col = "red")
title("Normal Q-Q Plot")



### Visualization ###

# Specify desired order
desired_order <- c("Acanthurus nigrofuscus", "Zebrazoma xanthurum")  # Add more species in the desired order if needed


# Bite Distance Violin Plot 
violin_bite_distance <- ggplot(data_no_outliers_bite_distance %>% 
                                 filter(!is.na(`Bite dist [mm] (Elaborated after Neuswanger e...)`)), 
                               aes(x = factor(`Species UID`) %>% fct_relevel(. , desired_order))) +
  geom_violin(aes(y = `Bite dist [mm] (Elaborated after Neuswanger e...)`, fill = "Bite Distance"), 
              position = "identity",
              alpha = 0.5, color = "NA") +
  geom_jitter(aes(y = `Bite dist [mm] (Elaborated after Neuswanger e...)`, color = "Bite Distance"), 
              position = position_jitter(width = 0.2), 
              alpha = 0.5) +
  labs(title = "",
       x = "",
       y = "Bite Distance [mm]",
       fill = "Variable",
       color = "Variable") +
  theme_minimal() +
  scale_fill_manual(values = "#16574D") +
  scale_color_manual(values = "#16574D") + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 200), breaks = seq(0, 200, by = 40), minor_breaks = seq(0, 200, by = 20)) +  # Manually setting minor breaks
  theme(
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.line.y = element_line(color = "black", size = 0.5),
    axis.ticks.y = element_line(color = "black"),
    panel.grid.major.x = element_blank(),  # Remove vertical grid lines
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),  # Remove horizontal grid lines
    panel.grid.minor.y = element_blank()
  )

print(violin_bite_distance)

# Bite Rate Violin Plot
violin_bite_rate <- ggplot(data_no_outliers_bite_rate %>% 
                             filter(!is.na(`Bite rate [#/min] (Elaborated after Neuswanger e...)`)), 
                           aes(x = factor(`Species UID`) %>% fct_relevel(. , desired_order))) +
  geom_violin(aes(y = `Bite rate [#/min] (Elaborated after Neuswanger e...)`, fill = "Bite Rate"), 
              position = "identity",
              alpha = 0.5, color = "NA") +
  geom_jitter(aes(y = `Bite rate [#/min] (Elaborated after Neuswanger e...)`, color = "Bite Rate"), 
              position = position_jitter(width = 0.2), 
              alpha = 0.5) +
  labs(title = "",
       x = "",
       y = "Bite Rate [bites min^-1]",
       fill = "Variable",
       color = "Variable") +
  theme_minimal() +
  scale_fill_manual(values = "#54AF9C") +  
  scale_color_manual(values = "#54AF9C") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 120), breaks = seq(0, 120, by = 20), minor_breaks = seq(0, 120, by = 10)) +
  theme(
    aspect.ratio = 1,
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.line.y = element_line(color = "black", size = 0.5),
    axis.ticks.y = element_line(color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )

common_theme <- theme(
  axis.title = element_text(size = 6),
  axis.text.x = element_text(size = 6, hjust = 0.5),  # Centered x-axis labels
  axis.text.y = element_text(size = 6),
  legend.title = element_text(size = 6),
  legend.text = element_text(size = 6),
  aspect.ratio = 0.75,
  axis.line.x = element_line(color = "black", size = 0.5),
  axis.line.y = element_line(color = "black", size = 0.5),
  axis.ticks.y = element_line(color = "black"),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.major.y = element_blank(),
  panel.grid.minor.y = element_blank()
)

# Adjusted Bite Distance Plot with the common theme and without the legend
violin_bite_distance <- violin_bite_distance + common_theme + theme(legend.position = "none")

# Adjusted Bite Rate Plot with the common theme and without the legend
violin_bite_rate <- violin_bite_rate + common_theme + theme(legend.position = "none")

# Combine plots using patchwork
combined_plot <- violin_bite_distance / violin_bite_rate 
combined_plot <- combined_plot + plot_layout(heights = c(1, 1))

print(combined_plot)

### Sumarize Bite Distances and Bite Rates ###

# Perform the summarization for bite distance, omitting NAs
summary_data_no_outliers_distance <- data_no_outliers_bite_distance %>%
  filter(!is.na(`Bite dist [mm] (Elaborated after Neuswanger e...)`)) %>%
  group_by(`Species UID`) %>%
  summarise(
    number_of_individuals = n_distinct(`ID (Fish ID)`),
    number_of_bites = sum(`Bites [#] (Elaborated after Neuswanger e...)`, na.rm = TRUE),
    mean_bite_distance = mean(`Bite dist [mm] (Elaborated after Neuswanger e...)`, na.rm = TRUE),
    sd_bite_distance = sd(`Bite dist [mm] (Elaborated after Neuswanger e...)`, na.rm = TRUE)
  )

# Print the summary data for bite distance
print(summary_data_no_outliers_distance)

# Perform the summarization for bite rate, omitting NAs
summary_data_no_outliers_rate <- data_no_outliers_bite_rate %>%
  filter(!is.na(`Bite rate [#/min] (Elaborated after Neuswanger e...)`)) %>%
  group_by(`Species UID`) %>%
  summarise(
    number_of_individuals = n_distinct(`ID (Fish ID)`),
    number_of_bites = sum(`Bites [#] (Elaborated after Neuswanger e...)`, na.rm = TRUE),
    mean_bite_rate = mean(`Bite rate [#/min] (Elaborated after Neuswanger e...)`, na.rm = TRUE),
    sd_bite_rate = sd(`Bite rate [#/min] (Elaborated after Neuswanger e...)`, na.rm = TRUE)
  )

# Print the summary data for bite rate
print(summary_data_no_outliers_rate)

# Export tables to file
write.csv(summary_data_no_outliers_distance, "surgeonfish_summary_bite_distance.csv")
write.csv(summary_data_no_outliers_rate, "surgeonfish_summary_bite_rate.csv")

### Correlation with Substrate Cover ###

# Fit and perform diagnostics for each benthos cover

## Cover rock [%] ##
lme_model_bite_distance_cover_rock <- lmer(`Log_Bite_dist_mm` ~ `Species UID` + `Cover rock [%]` +
                                                 (1 | `ID (Fish ID)`:`ID (Quadrat ID)`),
                                               data = data_no_outliers_bite_distance,
                                               control = lmerControl(optimizer = "Nelder_Mead", optCtrl = list(maxIter = 2000)))

qqnorm(resid(lme_model_bite_distance_cover_rock))
qqline(resid(lme_model_bite_distance_cover_rock))

plot(fitted(lme_model_bite_distance_cover_rock), resid(lme_model_bite_distance_cover_rock))
plot(data_no_outliers$`Cover rock [%]`, resid(lme_model_bite_distance_cover_rock))

acf(resid(lme_model_log_bite_distance_cover_rock))

print(summary(lme_model_log_bite_distance_cover_rock))


## Cover live cor [%] ##
lme_model_bite_distance_cover_live_cor <- lmer(`Log_Bite_dist_mm` ~ `Species UID` + `Cover live cor [%]` +
                                                 (1 | `ID (Fish ID)`:`ID (Quadrat ID)`),
                                               data = data_no_outliers_bite_distance,
                                               control = lmerControl(optimizer = "Nelder_Mead", optCtrl = list(maxIter = 2000)))

if (isSingular(lme_model_bite_distance_cover_live_cor)) {
  warning("Model did not converge; consider alternative optimizer or model specification.")
}

qqnorm(resid(lme_model_bite_distance_cover_live_cor))
qqline(resid(lme_model_bite_distance_cover_live_cor))

plot(fitted(lme_model_bite_distance_cover_live_cor), resid(lme_model_bite_distance_cover_live_cor))
plot(data_no_outliers$`Cover live cor [%]`, resid(lme_model_bite_distance_cover_live_cor))

acf(resid(lme_model_bite_distance_cover_live_cor))

print(summary(lme_model_bite_distance_cover_live_cor))

## Sand [#] ##
lme_model_bite_distance_sand <- lmer(`Log_Bite_dist_mm` ~ `Species UID` + `Sand [#]` +
                                       (1 | `ID (Fish ID)`:`ID (Quadrat ID)`),
                                     data = data_no_outliers_bite_distance,
                                     control = lmerControl(optimizer = "Nelder_Mead", optCtrl = list(maxIter = 2000)))

if (isSingular(lme_model_bite_distance_sand)) {
  warning("Model did not converge; consider alternative optimizer or model specification.")
}

qqnorm(resid(lme_model_bite_distance_sand))
qqline(resid(lme_model_bite_distance_sand))

plot(fitted(lme_model_bite_distance_sand), resid(lme_model_bite_distance_sand))
plot(data_no_outliers$`Sand [#]`, resid(lme_model_bite_distance_sand))

acf(resid(lme_model_bite_distance_sand))

print(summary(lme_model_bite_distance_sand))


## Cover dead cor [%] ##
lme_model_bite_distance_cover_dead_cor <- lmer(`Log_Bite_dist_mm` ~ `Species UID` + `Cover dead cor [%]` +
                                                 (1 | `ID (Fish ID)`:`ID (Quadrat ID)`),
                                               data = data_no_outliers_bite_distance,
                                               control = lmerControl(optimizer = "Nelder_Mead", optCtrl = list(maxIter = 2000)))

if (isSingular(lme_model_bite_distance_cover_dead_cor)) {
  warning("Model did not converge; consider alternative optimizer or model specification.")
}

qqnorm(resid(lme_model_bite_distance_cover_dead_cor))
qqline(resid(lme_model_bite_distance_cover_dead_cor))

plot(fitted(lme_model_bite_distance_cover_dead_cor), resid(lme_model_bite_distance_cover_dead_cor))
plot(data_no_outliers$`Cover dead cor [%]`, resid(lme_model_bite_distance_cover_dead_cor))

acf(resid(lme_model_bite_distance_cover_dead_cor))

print(summary(lme_model_bite_distance_cover_dead_cor))


## Cover rubble [%] ##
lme_model_bite_distance_cover_rubble <- lmer(`Log_Bite_dist_mm` ~ `Species UID` + `Cover rubble [%]` +
                                               (1 | `ID (Fish ID)`:`ID (Quadrat ID)`),
                                             data = data_no_outliers_bite_distance,
                                             control = lmerControl(optimizer = "Nelder_Mead", optCtrl = list(maxIter = 2000)))

if (isSingular(lme_model_bite_distance_cover_rubble)) {
  warning("Model did not converge; consider alternative optimizer or model specification.")
}

qqnorm(resid(lme_model_bite_distance_cover_rubble))
qqline(resid(lme_model_bite_distance_cover_rubble))

plot(fitted(lme_model_bite_distance_cover_rubble), resid(lme_model_bite_distance_cover_rubble))
plot(data_no_outliers$`Cover rubble [%]`, resid(lme_model_bite_distance_cover_rubble))

acf(resid(lme_model_bite_distance_cover_rubble))

print(summary(lme_model_bite_distance_cover_rubble))

## Bite Rate ##

## Cover rock [%] ##
lme_model_bite_rate_cover_rock <- lmer(`Log_Bite_rate_per_min` ~ `Species UID` + `Cover rock [%]` +
                                         (1 | `ID (Fish ID)`:`ID (Quadrat ID)`),
                                       data = data_no_outliers_bite_rate,
                                       control = lmerControl(optimizer = "Nelder_Mead", optCtrl = list(maxIter = 2000)))

if (isSingular(lme_model_bite_rate_cover_rock)) {
  warning("Model did not converge; consider alternative optimizer or model specification.")
}

qqnorm(resid(lme_model_bite_rate_cover_rock))
qqline(resid(lme_model_bite_rate_cover_rock))

plot(fitted(lme_model_bite_rate_cover_rock), resid(lme_model_bite_rate_cover_rock))
plot(data_no_outliers_rate$`Cover rock [%]`, resid(lme_model_bite_rate_cover_rock))

acf(resid(lme_model_bite_rate_cover_rock))

print(summary(lme_model_bite_rate_cover_rock))


## Cover live cor [%] ##
lme_model_bite_rate_cover_live_cor <- lmer(`Log_Bite_rate_per_min` ~ `Species UID` + `Cover live cor [%]` +
                                             (1 | `ID (Fish ID)`:`ID (Quadrat ID)`),
                                           data = data_no_outliers_bite_rate,
                                           control = lmerControl(optimizer = "Nelder_Mead", optCtrl = list(maxIter = 2000)))

if (isSingular(lme_model_bite_rate_cover_live_cor)) {
  warning("Model did not converge; consider alternative optimizer or model specification.")
}

qqnorm(resid(lme_model_bite_rate_cover_live_cor))
qqline(resid(lme_model_bite_rate_cover_live_cor))

plot(fitted(lme_model_bite_rate_cover_live_cor), resid(lme_model_bite_rate_cover_live_cor))
plot(data_no_outliers_rate$`Cover live cor [%]`, resid(lme_model_bite_rate_cover_live_cor))

acf(resid(lme_model_bite_rate_cover_live_cor))

print(summary(lme_model_bite_rate_cover_live_cor))


## Sand [#] ##
lme_model_bite_rate_sand <- lmer(`Log_Bite_rate_per_min` ~ `Species UID` + `Sand [#]` +
                                   (1 | `ID (Fish ID)`:`ID (Quadrat ID)`),
                                 data = data_no_outliers_bite_rate,
                                 control = lmerControl(optimizer = "Nelder_Mead", optCtrl = list(maxIter = 2000)))

if (isSingular(lme_model_bite_rate_sand)) {
  warning("Model did not converge; consider alternative optimizer or model specification.")
}

qqnorm(resid(lme_model_bite_rate_sand))
qqline(resid(lme_model_bite_rate_sand))

plot(fitted(lme_model_bite_rate_sand), resid(lme_model_bite_rate_sand))
plot(data_no_outliers_rate$`Sand [#]`, resid(lme_model_bite_rate_sand))

acf(resid(lme_model_bite_rate_sand))

print(summary(lme_model_bite_rate_sand))


## Cover dead cor [%] ##
lme_model_bite_rate_cover_dead_cor <- lmer(`Log_Bite_rate_per_min` ~ `Species UID` + `Cover dead cor [%]` +
                                             (1 | `ID (Fish ID)`:`ID (Quadrat ID)`),
                                           data = data_no_outliers_bite_rate,
                                           control = lmerControl(optimizer = "Nelder_Mead", optCtrl = list(maxIter = 2000)))

if (isSingular(lme_model_bite_rate_cover_dead_cor)) {
  warning("Model did not converge; consider alternative optimizer or model specification.")
}

qqnorm(resid(lme_model_bite_rate_cover_dead_cor))
qqline(resid(lme_model_bite_rate_cover_dead_cor))

plot(fitted(lme_model_bite_rate_cover_dead_cor), resid(lme_model_bite_rate_cover_dead_cor))
plot(data_no_outliers_rate$`Cover dead cor [%]`, resid(lme_model_bite_rate_cover_dead_cor))

acf(resid(lme_model_bite_rate_cover_dead_cor))

print(summary(lme_model_bite_rate_cover_dead_cor))


## Cover rubble [%] ##
lme_model_bite_rate_cover_rubble <- lmer(`Log_Bite_rate_per_min` ~ `Species UID` + `Cover rubble [%]` +
                                           (1 | `ID (Fish ID)`:`ID (Quadrat ID)`),
                                         data = data_no_outliers_bite_rate,
                                         control = lmerControl(optimizer = "Nelder_Mead", optCtrl = list(maxIter = 2000)))

if (isSingular(lme_model_bite_rate_cover_rubble)) {
  warning("Model did not converge; consider alternative optimizer or model specification.")
}

qqnorm(resid(lme_model_bite_rate_cover_rubble))
qqline(resid(lme_model_bite_rate_cover_rubble))

plot(fitted(lme_model_bite_rate_cover_rubble), resid(lme_model_bite_rate_cover_rubble))
plot(data_no_outliers_rate$`Cover rubble [%]`, resid(lme_model_bite_rate_cover_rubble))

acf(resid(lme_model_bite_rate_cover_rubble))

print(summary(lme_model_bite_rate_cover_rubble))

