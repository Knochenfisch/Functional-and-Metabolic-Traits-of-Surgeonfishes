#### Analysis of 3D Fish Trajectories to infer Surgeonfish Metabolic Traits ####

# ======================================
# ENVIRONMENT SETUP USING RENV
# ======================================

# This section sets up a reproducible environment using 'renv'.
# If you're setting up this project for the first time, uncomment and run the necessary sections.

# Load the renv package
library(renv)  

# --- Initialization (only run once) ---
# This initializes a project-specific library for package management.

#renv::init()  

# Activate the renv environment
renv::activate()


# --- Load Libraries ---
library(purrr)
library(codetools)
library(hms)
library(viridis)
library(plotly)
library(magrittr)
library(dplyr)
library(scales)
library(forcats)
library(tidyr)
library(ggplot2)
library(MASS)
library(patchwork)
library(stringr)
library(KFAS)
library(zoo)
library(car)
library(lmtest)
library(markdown)

# --- Snapshot (only run once) ---
# This saves the package info for reproducibility. 
# Uncomment below if adding new packages or setting up for the first time.

#renv::snapshot()

### Data download and cleaning ###

# The dataset is fully citeable as Lilkendey, J et al.: Foraging traits, lengths and 3D movement trajectories of 
# coral reef fishes obtained via stereo-video in Eilat, Gulf of Aqaba, Red Sea. 
# PANGAEA, https://doi.pangaea.de/10.1594/PANGAEA.957634 

## The Gulf_of_Aqaba_3D_track dataset is available at https://doi.pangaea.de/10.1594/PANGAEA.957632?format=textfile"

# If the dataset is already in the workspace
data_2 <- Gulf_of_Aqaba_3D_track

print(colnames(data_2))

# Rename columns and filter out specific data in one step
data_2 <- data_2 %>%
  rename(
    posX_no_interpolation = `Center pos x [cm]`,
    posY_no_interpolation = `Center pos y [cm]`,
    posZ_no_interpolation = `Center pos z [cm]`,
    `Fish ID` = 'ID'
  ) %>%
  filter(!(`Species UID` == "Acanthurus nigrofuscus" & `l [cm]` > 21))

### Initial outlier removal ###

# Calculate the 2.5% and 97.5% percentiles for each coordinate
x_range <- quantile(data_2$posX_no_interpolation, c(0.025, 0.975), na.rm = TRUE)
y_range <- quantile(data_2$posY_no_interpolation, c(0.025, 0.975), na.rm = TRUE)
z_range <- quantile(data_2$posZ_no_interpolation, c(0.025, 0.975), na.rm = TRUE)

# Filter the data for outliers
data_2 <- data_2 %>%
  mutate(
    filtered_posX = ifelse(posX_no_interpolation >= x_range[1] & posX_no_interpolation <= x_range[2], posX_no_interpolation, NA),
    filtered_posY = ifelse(posY_no_interpolation >= y_range[1] & posY_no_interpolation <= y_range[2], posY_no_interpolation, NA),
    filtered_posZ = ifelse(posZ_no_interpolation >= z_range[1] & posZ_no_interpolation <= z_range[2], posZ_no_interpolation, NA)
  )

### Median smoothing ###

# Applying a running median with a window of 5 data points for the filtered data
data_2 <- data_2 %>%
  group_by(`Fish ID`) %>%
  mutate(
    median_posX = rollapply(filtered_posX, width = 5, FUN = median, align = 'center', fill = NA),
    median_posY = rollapply(filtered_posY, width = 5, FUN = median, align = 'center', fill = NA),
    median_posZ = rollapply(filtered_posZ, width = 5, FUN = median, align = 'center', fill = NA)
  )

### Kalman filtering ###

# Kalman Filtering function
kalman_filter_data <- function(data_2){
  model <- SSModel(data_2 ~ SSMtrend(1, Q = list(free = TRUE, initial = "zero")), H = NA)
  fit <- fitSSM(inits = rep(0,2), model = model, method = 'BFGS', control = list(trace = 0))
  kfilter <- KFS(fit$model)
  return(kfilter$alphahat[,"level",drop=FALSE])
}

# Apply Kalman filter on the median smoothed data
data_2 <- data_2 %>%
  group_by(`Fish ID`) %>%
  mutate(
    kalman_median_posX = kalman_filter_data(median_posX),
    kalman_median_posY = kalman_filter_data(median_posY),
    kalman_median_posZ = kalman_filter_data(median_posZ)
  )

### Visualizing data transformation ###

# Filter for Fish ID "AN04"
fish_AN04_data <- dplyr::filter(data_2, `Fish ID` == "AN04")

# Plot X coordinates for Fish ID "AN04"
ggplot(fish_AN04_data, aes(x = seq_along(posX_no_interpolation), y = posX_no_interpolation)) +
  geom_line(aes(color = "Raw"), size = 1) +
  geom_line(aes(y = median_posX, color = "Running Median"), size = 1) +
  geom_line(aes(y = kalman_median_posX, color = "Kalman Filtered"), size = 1) +
  labs(title = "Comparison of Raw, Running Median, and Kalman Filtered X-coordinate (Fish AN5)",
       x = "Time Step",
       y = "X-coordinate",
       color = "Data") +
  scale_color_manual(values = c("Raw" = "blue", "Running Median" = "green", "Kalman Filtered" = "red")) +
  theme_minimal()

# Plot Y coordinates for Fish ID "AN04"
ggplot(fish_AN04_data, aes(x = seq_along(posY_no_interpolation), y = posY_no_interpolation)) +
  geom_line(aes(color = "Raw"), size = 1) +
  geom_line(aes(y = median_posY, color = "Running Median"), size = 1) +
  geom_line(aes(y = kalman_median_posY, color = "Kalman Filtered"), size = 1) +
  labs(title = "Comparison of Raw, Running Median, and Kalman Filtered Y-coordinate (Fish AN5)",
       x = "Time Step",
       y = "Y-coordinate",
       color = "Data") +
  scale_color_manual(values = c("Raw" = "blue", "Running Median" = "green", "Kalman Filtered" = "red")) +
  theme_minimal()

# Plot Z coordinates for Fish ID "AN04"
ggplot(fish_AN04_data, aes(x = seq_along(posZ_no_interpolation), y = posZ_no_interpolation)) +
  geom_line(aes(color = "Raw"), size = 1) +
  geom_line(aes(y = median_posZ, color = "Running Median"), size = 1) +
  geom_line(aes(y = kalman_median_posZ, color = "Kalman Filtered"), size = 1) +
  labs(title = "Comparison of Raw, Running Median, and Kalman Filtered Z-coordinate (Fish AN5)",
       x = "Time Step",
       y = "Z-coordinate",
       color = "Data") +
  scale_color_manual(values = c("Raw" = "blue", "Running Median" = "green", "Kalman Filtered" = "red")) +
  theme_minimal()


# Plotting Raw, Median smoothed, and Kalman filtered data for Fish AN04

# Filter data for Fish ID "AN04"
AN04_data <- dplyr::filter(data_2, `Fish ID` == "AN04")

# Create a 3D scatter plot for Fish ID "AN04"
plot_ly() %>%
  
  # Raw Data
  add_trace(data = AN04_data, 
            x = ~`filtered_posX`, 
            y = ~`filtered_posZ`, 
            z = ~`filtered_posY`,
            type = 'scatter3d', mode = 'markers', 
            name = 'Raw Data',
            marker = list(color = ~`Frame count, continuous`, colorscale = 'Viridis', size = 3, colorbar=list(title="Frame Number"))) %>%
  
  # Running Median Data
  add_trace(data = AN04_data, 
            x = ~median_posX, 
            y = ~median_posZ, 
            z = ~median_posY,
            type = 'scatter3d', mode = 'markers', 
            name = 'Running Median',
            marker = list(color = 'orange', size = 3)) %>%
  
  # Kalman Filtered Data
  add_trace(data = AN04_data, 
            x = ~kalman_median_posX, 
            y = ~kalman_median_posZ, 
            z = ~kalman_median_posY,
            type = 'scatter3d', mode = 'markers', 
            name = 'Kalman Filtered on Median',
            marker = list(color = 'red', size = 3)) %>%
  
  # Layout & Axes
  layout(scene = list(
    xaxis = list(title = "X"),
    yaxis = list(title = "Y"),
    zaxis = list(title = "Z")
  ))

# Plotting only Kalman filtered Data for AN01

# Filter data for fish AN01
AN01_data <- dplyr::filter(data_2, `Fish ID` == "AN01")

# Create a basic 3D scatter plot
scatter3d_plot <- plot_ly(
  AN01_data,
  x = ~kalman_median_posX,
  y = ~kalman_median_posZ,
  z = ~kalman_median_posY,
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 2, opacity = 0.6)
) %>%
  layout(
    scene = list(
      xaxis = list(title = "X", range = c(-70,30)),
      zaxis = list(title = "Y", range = c(-150,150)),
      yaxis = list(title = "Z", range = c(0,400)),
      aspectmode = "cube"
    ),
    title = "3D Scatter Plot of Fish AN5",
    margin = list(l = 0, r = 0, b = 0, t = 65)
  )

# Render the plot
scatter3d_plot


### Calculate Energy Expenditure [W] for each fish ###

data_2 <- data_2 %>%
  arrange(`Fish ID`) %>%  # Making sure data is ordered correctly
  group_by(`Fish ID`) %>% # Keeping calculations within each transect
  mutate(
    # Calculate velocity using filtered data
    velocity_x = c(0, diff(kalman_median_posX) * 60),
    velocity_y = c(0, diff(kalman_median_posY) * 60),
    velocity_z = c(0, diff(kalman_median_posZ) * 60),
    
    # Calculate acceleration using velocity from filtered data
    acceleration_x = c(0, diff(velocity_x) * 60),
    acceleration_y = c(0, diff(velocity_y) * 60),
    acceleration_z = c(0, diff(velocity_z) * 60),
    
    ODBA = abs(acceleration_x) + abs(acceleration_y) + abs(acceleration_z),
  )

### SMR - Mass relationships of surrogate species ###

# Dataset available at https://github.com/nschiett/activity_rate_fishes/blob/master/data/data_respirometry.csv
# Please cite: Schiettekatte NMD, Conte F, French B, Brandl SJ, Fulton CJ, Mercière A, Norin T, Villéger S, Parravicini V (2022) 
# Combining stereo‐video monitoring and physiological trials to estimate reef fish metabolic demands in the wild. 
# Ecology and Evolution 12. https://doi.org/10.1002/ece3.9084

data_respirometry

# Filter data for C. striatus
data_C_striatus <- subset(data_respirometry, Species == "Ctenochaetus striatus")

# Log-log regression for C. striatus
model_C_striatus <- lm(log(SMR_gd) ~ log(Weight_g), data = data_C_striatus)
summary_C_striatus <- summary(model_C_striatus)
coefs_C_striatus <- coef(model_C_striatus)
K_C_striatus <- exp(coefs_C_striatus[1])
E_C_striatus <- coefs_C_striatus[2]

# Filter data for Z. scopas
data_Z_scopas <- subset(data_respirometry, Species == "Zebrasoma scopas")

# Log-log regression for Z. scopas
model_Z_scopas <- lm(log(SMR_gd) ~ log(Weight_g), data = data_Z_scopas)
summary_Z_scopas <- summary(model_Z_scopas)
coefs_Z_scopas <- coef(model_Z_scopas)
K_Z_scopas <- exp(coefs_Z_scopas[1])
E_Z_scopas <- coefs_Z_scopas[2]

# Constants for surrogate species based on log-log regression
constant_data <- data.frame(
  `Species UID` = c("Acanthurus nigrofuscus", "Zebrasoma xanthurum"),
  K = c(K_C_striatus, K_Z_scopas), # Original K values in g
  E = c(E_C_striatus, E_Z_scopas) # Exponent E is unitless
)

# Adjust K for for mg O2
constant_data$K <- constant_data$K * 1000  # # Convert g O2/day to mg O2/day

# Check current column names in constant_data
constant_data

# Rename column in constant_data if necessary
constant_data <- rename(constant_data, `Species UID` = `Species.UID`)

# Join constant_data with data_2
data_2 <- left_join(data_2, constant_data, by = "Species UID")

# Check if the join worked correctly
View(data_2)


# Conversion factor from Joules to milligrams of oxygen
conversion_factor_X <- 14.1  # J per mg O2 taken from Elliott, J. M., and Wo Davison. "Energy equivalents of oxygen consumption in animal energetics." Oecologia 19 (1975): 195-201.

# Temperature correction using Q10 Value and Temperatures

# As we don't have a specific Q10 for our species, we used 1.92 as calculated for Zebrasoma scopas 
# Source: McFarlane, L. (2016). Thermal tolerance and metabolic responses of brushtail tang, Zebrasoma scopas (Master thesis, California State University, Northridge).

# Define Q10, reference temperature, and study temperature
Q10 <- 1.92  # Value for our surrogate species Z. scopas taken from McFarlane (2016)
ref_temp_C <- 28  # Reference temperature in Celsius (from Schiettekatte et al. (2022), recorded at Mo′orea,	French	Polynesia,	between	March	2018 and	February	 2019)
study_temp_C <- 21  # Study temperature in Celsius (from Eilat, Isreael, March 2018)

# Calculating the Temperature Adjustment Factor and apply to metabolic rate data

# Convert temperature difference to a factor of 10°C
temp_diff_factor <- (ref_temp_C - study_temp_C) / 10

# Calculate the temperature adjustment factor
temp_adjustment_factor <- Q10 ^ temp_diff_factor

# Calculate metabolic rate in Joules per day using the O2 consumption (in mg O2 per day) and conversion factor
data_2 <- data_2 %>%
  mutate(
    Metabolic_Rate_J_per_day = (K * (`Mass [g]` ^ E)) * ODBA * conversion_factor_X * temp_adjustment_factor, # Convert mg O2/day to J/day, Mass_g is mass in grams
    EE_W = Metabolic_Rate_J_per_day / (60 * 60 * 24),  # Convert J/day to W (J/s)
    velocity = sqrt(velocity_x^2 + velocity_y^2 + velocity_z^2)
  )

### Compare mean velocities and mean EE between species ###

individual_means <- data_2 %>%
  group_by(`Fish ID`, `Species UID`, `Station`, `Rack no`, `Part desc`) %>%
  summarise(
    mean_mass = mean(`Mass [g]`, na.rm = TRUE),
    mean_length = mean(`l [cm]`, na.rm = TRUE),
    mean_velocity = mean(velocity, na.rm = TRUE),
    sd_velocity = sd(velocity, na.rm = TRUE),
    mean_EE_W = mean(EE_W, na.rm = TRUE),
    sd_EE_W = sd(EE_W, na.rm = TRUE),
    mean_class_result = mean(`Class result`, na.rm = TRUE), # Calculate the mean classification result
    first_frame_count_continuous = first(`Frame count, continuous`), # Include the first entry for Frame count, continuous
    .groups = "drop" # This option ensures that the grouping is dropped after summarising
  )

print(individual_means)


# Export table to file
write.csv(individual_means, "individual_surgeonfish_summary_velocity_energy_expenditure.csv")


### Create Overview table to enhance comparison with literature ###

# Group by Species UID for species summary
species_summary_velocity_EE <- individual_means %>%
  group_by(`Species UID`) %>%
  summarise(
    number_of_individuals = n_distinct(`Fish ID`),
    species_mean_velocity = mean(mean_velocity, na.rm = TRUE),
    sd_velocity = sd(mean_velocity, na.rm = TRUE),
    species_mean_EE_W = mean(mean_EE_W, na.rm = TRUE),
    sd_EE_W = sd(mean_EE_W, na.rm = TRUE)
  )

# Convert mean EE from W to J min^-1
species_summary_velocity_EE <- species_summary_velocity_EE %>%
  mutate(
    mean_EE_J_min = species_mean_EE_W * 60  # Convert from W to J/min
  )

# Convert mean EE from J min^-1 to mg O2 kg^-1 min^-1 
species_summary_velocity_EE <- species_summary_velocity_EE %>%
  mutate(
    mean_EE_O2_kg_min = mean_EE_J_min / conversion_factor_X / 1000  # Convert from J/min to mg O2/min/kg
  )

# Convert mean EE from mg O2 kg^-1 min^-1 to mg O2 g^-1 h^-1
species_summary_velocity_EE <- species_summary_velocity_EE %>%
  mutate(
    mean_EE_O2_g_hr = mean_EE_O2_kg_min * 60 / 1000  # Convert from mg O2/kg/min to mg O2/g/hr
  )

# Convert mean EE from mg O2 g^-1 h^-1 to mg O2 g^-1 d^-1
species_summary_velocity_EE <- species_summary_velocity_EE %>%
  mutate(
    mean_EE_O2_g_day = mean_EE_O2_g_hr * 24 # Convert from mg O2/min/kg to mg O2/g/day (1440 min/day)
  )

# Print the species summary with new conversions
print(species_summary_velocity_EE) 

# Export table to file
write.csv(species_summary_velocity_EE, "surgeonfish_summary_velocity_energy_expenditure.csv")

## Violin plots

# Establish common theme

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

# Specify desired order
desired_order <- c("Acanthurus nigrofuscus", "Zebrasoma xanthurum")

# Velocity Plot 
violin_velocity <- ggplot(individual_means %>% 
                            dplyr::filter(!is.na(mean_velocity)), 
                          aes(x = factor(`Species UID`) %>% fct_relevel(. , desired_order))) +
  geom_violin(aes(y = mean_velocity, fill = "Velocity"), 
              position = "identity",
              alpha = 0.5, color = "NA") +
  geom_jitter(aes(y = mean_velocity, color = "Velocity"), 
              position = position_jitter(width = 0.2), 
              alpha = 0.5) +
  labs(title = "",
       x = "",
       y = "Velocity [cm s-1]") +
  theme_minimal() +
  scale_fill_manual(values = "#16574D") +
  scale_color_manual(values = "#16574D") + 
  common_theme +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(0, 50), expand = c(0,0))

print(violin_velocity)

# EE Plot
violin_EE <- ggplot(individual_means  %>% 
                      dplyr::filter(!is.na(mean_EE_W)), 
                    aes(x = factor(`Species UID`) %>% fct_relevel(. , desired_order))) +
  geom_violin(aes(y = mean_EE_W, fill = "Energy Expenditure"), 
              position = "identity",
              alpha = 0.5, color = "NA") +
  geom_jitter(aes(y = mean_EE_W, color = "Energy Expenditure"), 
              position = position_jitter(width = 0.2), 
              alpha = 0.5) +
  labs(title = "",
       x = "",
       y = "Energy Expenditure [W]") +
  theme_minimal() +
  scale_fill_manual(values = "#54AF9C") +
  scale_color_manual(values = "#54AF9C") + 
  common_theme +
  theme(legend.position = "none") +
  ylim(0, 40)

print(violin_EE)

# Combine Velocity and EE plots using patchwork
combined_plot_velocity_EE <- violin_velocity / violin_EE 
combined_plot_velocity_EE <- combined_plot_velocity_EE + plot_layout(heights = c(1, 1))

# Print the combined plot
print(combined_plot_velocity_EE)

## Statistical Copmarison

# Wilcoxon rank sum test using raw data (individual_means_filtered)
# For velocity
wilcox.test(mean_velocity ~ `Species UID`, data = individual_means)

# For EE
wilcox.test(mean_EE_W ~ `Species UID`, data = individual_means)

# Equal variances

# For velocity
group_velocity <- individual_means$`Species UID`
velocity_values <- individual_means$mean_velocity

leveneTest(velocity_values, group = group_velocity)

# For EE
group_EE <- individual_means$`Species UID`
EE_values <- individual_means$mean_EE_W

leveneTest(EE_values, group = group_EE)

str(individual_means)


# Check if the filtering worked
print(head(individual_means))

# Combined scatter plot with separate regression lines and confidence intervals for each species
# Using the filtered dataset
scatter_plot_mass_EE_combined <- ggplot(individual_means, aes(x = mean_mass, y = mean_EE_W, color = `Species UID`)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, aes(color = `Species UID`)) +  # Add separate colors and confidence intervals
  labs(title = "Correlation between Individual Mean Mass and Energy Expenditure",
       x = "Mean Mass [g]",
       y = "Mean Energy Expenditure [W]",
       color = "Species UID") +
  theme_minimal() +
  scale_color_manual(values = c("#16574D", "#54AF9C"))  # Specify colors for each species

# Print the combined scatter plot
print(scatter_plot_mass_EE_combined)

# Fit the ANCOVA model
ancova_model <- lm(mean_EE_W ~ mean_mass * `Species UID`, data = individual_means)

# Perform ANOVA on the ANCOVA model
ancova_anova <- anova(ancova_model)

# Print the ANOVA table
print(ancova_anova)




### Plotting all fish trajectories in the same Plot ###

is.ts(data_2$kalman_median_posX)
is.ts(data_2$kalman_median_posZ)
is.ts(data_2$kalman_median_posY)

# Convert time series into vector

data_2$kalman_median_posX <- as.vector(data_2$kalman_median_posX)
data_2$kalman_median_posZ <- as.vector(data_2$kalman_median_posZ)
data_2$kalman_median_posY <- as.vector(data_2$kalman_median_posY)

# Get unique fish IDs
unique_fish <- unique(data_2$`Fish ID`)

# Generate colors from the Viridis palette
num_fish <- length(unique_fish)
fish_colors <- viridis(num_fish)

# Create a named vector mapping fish IDs to colors
color_map <- setNames(fish_colors, unique_fish)

# Create the 3D scatter plot using the desired columns for positions
scatter3d_plot <- plot_ly(
  data = data_2,
  x = ~kalman_median_posX,
  y = ~kalman_median_posZ,
  z = ~kalman_median_posY,
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 2, opacity = 0.6, color = ~color_map[as.character(`Fish ID`)], showscale = TRUE),
  text = ~paste("Fish ID:", `Fish ID`)
) %>%
  layout(
    scene = list(
      xaxis = list(title = "X", range = c(-200,200)),
      yaxis = list(title = "Z", range = c(0,400)),
      zaxis = list(title = "Y", range = c(-100,100)),
      aspectmode = "cube"
    ),
    title = "3D Scatter Plot of All Fish",
    margin = list(l = 0, r = 0, b = 0, t = 65)
  )

# Display the plot
scatter3d_plot

### HTML display of all 3D fish trajectories with Energy Expenditure as coloration ###

data_per_fish <- split(data_2, data_2$`Fish ID`)

# Initialize the scatter3d_plot
scatter3d_plot <- plot_ly()

# Add data for all fishes to the plot
for (i in 1:length(unique_fish)) {
  
  # Raw traces
  scatter3d_plot <- scatter3d_plot %>%
    add_trace(
      data = data_per_fish[[i]],
      x = ~`filtered_posX`, 
      y = ~`filtered_posZ`, 
      z = ~`filtered_posY`,
      type = "scatter3d",
      mode = "markers",
      marker = list(color = '#16574D', size = 5, opacity = 0.6),
      text = ~paste("Trajectory ID:", `Fish ID`, "<br>Fish total length [cm]:", sprintf("%.1f", `l [cm]`)),
      visible = (i == 1),
      name = "Raw Data"
    )
  
  # Median traces
  scatter3d_plot <- scatter3d_plot %>%
    add_trace(
      data = data_per_fish[[i]],
      x = ~median_posX, 
      y = ~median_posZ, 
      z = ~median_posY,
      type = "scatter3d",
      mode = "markers",
      marker = list(color = '#54AF9C', size = 5, opacity = 0.6),
      text = ~paste("Trajectory ID:", `Fish ID`, "<br>Fish total length [cm]:", sprintf("%.1f", `l [cm]`)),
      visible = (i == 1),
      name = "Median Smoothed Data"
    )
  
  # Kalman traces
  scatter3d_plot <- scatter3d_plot %>%
    add_trace(
      data = data_per_fish[[i]],
      x = ~kalman_median_posX, 
      y = ~kalman_median_posZ, 
      z = ~kalman_median_posY,
      type = "scatter3d",
      mode = "markers",
      marker = list(
        color = ~log1p(EE_W),
        colorscale = 'Reds',
        size = 5,
        opacity = 0.6,
        colorbar = list(title = "log(1 + Energy Expenditure [W])")
      ),
      text = ~paste("Trajectory ID:", `Fish ID`, "<br>Fish total length [cm]:", sprintf("%.1f", `l [cm]`)),
      visible = (i == 1),
      name = "Kalman Filtered Data"
    )
}

# Create steps for the slider
steps <- list()
for(i in 1:length(unique_fish)) {
  step <- list(
    args = list('visible', rep(FALSE, length(unique_fish)*3)),
    label = as.character(unique_fish[i]),
    method = 'restyle'
  )
  step$args[[2]][(i*3-2):(i*3)] <- TRUE
  steps[[i]] <- step
}

# Finalize the plot layout
scatter3d_plot <- scatter3d_plot %>% 
  layout(
    scene = list(
      xaxis = list(title = "X [cm]", range = c(-200, 200)),
      yaxis = list(title = "Z [cm]", range = c(0, 400)),
      zaxis = list(title = "Y [cm]", range = c(-100, 100)),
      aspectmode = 'cube'
    ),
    sliders = list(
      list(
        active = 0,
        yanchor = "top",
        xanchor = "left",
        currentvalue = list(font = list(size = 20), prefix = "Trajectory ID: ", visible = TRUE),
        transition = list(duration = 300, easing = "cubic-in-out"),
        pad = list(b = 10, t = 50),
        len = 0.9,
        x = 0.1,
        y = 0,
        steps = steps
      )
    ),
    legend = list(
      x = 1,
      y = 0,
      xanchor = 'right',
      yanchor = 'bottom'
    )
  )

# Display the plot
scatter3d_plot
