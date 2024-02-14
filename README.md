# Herbivore Feeding Dynamics and Energy Expenditure on a Coral Reef: Insights from Remote Underwater Stereo-Video and AI-Driven 3D Tracking

This repository contains the code and data to reproduce species-level tables and figures presented in Lilkendey et al. "Herbivore feeding dynamics and energy expenditure on a coral reef: insights from remote underwater stereo-video and AI-driven 3D tracking".

**Authors:** Julian Lilkendey, Cyril Barrelet, Jingjing Zhang, Michael Meares, Houssam Larbi, Gérard Subsol, Marc Chaumont, Armagan Sabetian  
**DOI:** [10.1002/ece3.11070](https://doi.org/10.1002/ece3.11070)

**Short Abstract:** Our study unveils the complex interplay between animal movement ecology, feeding behaviour, and energy budgeting within coral reef ecosystems, focusing on the role of herbivorous fishes on coral reefs under anthropogenic stress. Employing advanced technologies such as remote underwater stereo-video and AI-driven 3D tracking, we investigate the feeding dynamics of Brown surgeonfish (*Acanthurus nigrofuscus*) and Yellowtail tang (*Zebrasoma xanthurum*), highlighting their contributions to reef health and resilience.

---

## R Script Documentation: Analysis of Surgeonfish Functional Traits

### Overview

This R script delves into the functional traits of surgeonfish, specifically examining bite distances and rates, and their associations with various substrate covers. The analysis employs linear mixed-effects models to account for nested data structures and potential random effects.

### Environment Setup

When working with an existing project that utilizes `renv` for environment management, you should first ensure that `renv` is installed in your R environment. Here's how to proceed:

1. **Install `renv`:** If you haven't already installed `renv`, you can do so by running the following command in your R console. This step is only necessary if `renv` is not already installed.

    ```R
    install.packages("renv")
    ```

2. **Open the Project:** Open the project in RStudio, which should automatically set your working directory to the root of the project. If you're not using RStudio, manually set your working directory to the root of the project in an R session.

3. **Restore the Environment:** Since the project is already set up with `renv`, it should include an `renv.lock` file, which specifies all the necessary R packages and their versions. Use the `renv::restore()` function to synchronize your local R environment with the project's `renv.lock` file.

    ```R
    renv::restore()
    ```

The `renv::restore()` command will check your current environment against the specifications in the `renv.lock` file and install any missing packages or required versions. This ensures that your environment matches the project's dependencies, facilitating reproducibility and consistency across different setups.

### Data Availability

The dataset is fully citeable as Lilkendey, J et al.: Foraging traits, lengths and 3D movement trajectories of coral reef fishes obtained via stereo-video in Eilat, Gulf of Aqaba, Red Sea. [PANGAEA.957634](https://doi.pangaea.de/10.1594/PANGAEA.957634)

The Gulf_of_Aqaba_Acanthuridae dataset is available at: [PANGAEA.957631](https://doi.pangaea.de/10.1594/PANGAEA.957631)

### Statistical Modeling

#### Bite Distance Analysis:

- Outliers in bite distances are identified and removed using interquartile range (IQR) methods.
- Linear mixed-effects models are fitted to log-transformed bite distances, considering species, fish mass, and nested random effects (Fish ID within Quadrat ID).
- Model comparisons are conducted using Akaike Information Criterion (AIC) to determine the best-fitting model.

#### Bite Rate Analysis:

- Similar to bite distances, outliers in bite rates are handled, and data is log-transformed to meet model assumptions.
- A series of linear mixed-effects models are fitted to the log-transformed bite rates, exploring the effects of species, fish mass, and nested random effects.
- Model diagnostics and comparisons are performed to ensure the validity and optimal fit of the models.

### Visualisation

- Jitter and violin plots are created to visualise the distributions of bite distances and rates among species.
- Density plots are used to illustrate the distribution of log-transformed bite distances and rates, providing insights into the variability within and between species.

### Correlation with Substrate Cover

- The script explores the relationship between bite metrics (distance and rate) and various types of substrate cover (e.g., rock, live coral, sand, dead coral, rubble) using linear mixed-effects models.
- Each substrate type is analysed separately, with models adjusted for species and nested random effects.
- Model diagnostics, including residual plots and Q-Q plots, are conducted to assess the assumptions and fit of each model.

---

## R Script Documentation: Analysis of 3D Fish Trajectories to Infer Surgeonfish Metabolic Traits

### Script Overview

This R script is dedicated to analysing 3D movement trajectories of surgeonfish to deduce their metabolic traits. The analysis encompasses data cleaning, outlier removal, median smoothing, Kalman filtering, and energy expenditure calculations, culminating in comprehensive data visualisations.

### Environment Setup

When working with an existing project that utilizes `renv` for environment management, you should first ensure that `renv` is installed in your R environment. Here's how to proceed:

1. **Install `renv`:** If you haven't already installed `renv`, you can do so by running the following command in your R console. This step is only necessary if `renv` is not already installed.

    ```R
    install.packages("renv")
    ```

2. **Open the Project:** Open the project in RStudio, which should automatically set your working directory to the root of the project. If you're not using RStudio, manually set your working directory to the root of the project in an R session.

3. **Restore the Environment:** Since the project is already set up with `renv`, it should include an `renv.lock` file, which specifies all the necessary R packages and their versions. Use the `renv::restore()` function to synchronize your local R environment with the project's `renv.lock` file.

    ```R
    renv::restore()
    ```

The `renv::restore()` command will check your current environment against the specifications in the `renv.lock` file and install any missing packages or required versions. This ensures that your environment matches the project's dependencies, facilitating reproducibility and consistency across different setups.

### Data Availability

The dataset is fully citeable as Lilkendey, J et al.: Foraging traits, lengths and 3D movement trajectories of coral reef fishes obtained via stereo-video in Eilat, Gulf of Aqaba, Red Sea. [PANGAEA.957634](https://doi.pangaea.de/10.1594/PANGAEA.957634)

The Gulf_of_Aqaba_3D_track dataset is available at [PANGAEA.957632](https://doi.pangaea.de/10.1594/PANGAEA.957632)

### Analytical Methods

- **Median Smoothing:** Applied to reduce noise in the trajectory data, enhancing the clarity of movement patterns.
- **Kalman Filtering:** Employed to further refine the trajectory data, providing a more accurate representation of fish movements.
- **Energy Expenditure Calculation:** Utilises the refined trajectory data to estimate the metabolic rates of the fish, incorporating factors such as velocity, acceleration, and overall dynamic body acceleration (ODBA).

### Visualizations

The script generates a series of plots to visualise the data transformation process, including comparisons of raw, median-smoothed, and Kalman-filtered data. Additionally, 3D scatter plots are created to depict the fish trajectories and their energy expenditures, offering insightful visual representations of the fish behaviour and metabolic traits.

### Statistical Analysis

The script conducts statistical comparisons to explore differences in velocity and energy expenditure between species, employing Wilcoxon rank sum tests and Levene's tests for equality of variances.

For further inquiries, Julian Lilkendey can be contacted at julian.lilkendey@icloud.com.


