rm(list = ls())

# Load required packages
# install the most recent version from GitHub
options(download.file.method = "wininet") 
#devtools::install_github("pet221/SSNbler", ref = "main")

library(SSNbler)
library(SSN2)
library(SSNbayes)
library(SSNdata)
library(sf)
library(tidyverse)
library(nhdplusTools)
library(elevatr)
library(raster)
library(exactextractr)

# Read shapefiles
CA_flow <- st_read("rmrs-flowline_ca18_nsi/Flowline_CA18_NSI.shp")

stan_flowline <- CA_flow |>
  filter(GNIS_NAME == "Stanislaus River")

CA_points <- st_read("rmrs-predictionpoints_ca18_nsi/PredictionPoints_CA18_NSI.shp")

stan_points <- CA_points |>
  filter(GNIS_NAME == "Stanislaus River")

stan_ref_points <- st_read("SSN_layers/StanRefPoints_final.shp")
stan_flow_line <- st_read("SSN_layers/StanFlowLine_final.shp")
snorkel <- read.csv("stanSnorkel 1.csv")

# Filter and process snorkel data
snorkel <- snorkel |>
  filter(sizeClass %in% c("Small (<150)", "Medium (150-300)")) |>
  filter(!is.na(Lat) & !is.na(Long)) |>
  group_by(Lat, Long, Date, Reach) |>             # Group by point and date
  summarise(
    counts = sum(count, na.rm = TRUE), # Sum counts
    Lngth_Yrd = mean(Lngth_Yrd, na.rm = TRUE), # Average Length of snorkel survey by grouping
    flow = mean(Flow, na.rm = TRUE),
    .groups = "drop"                    # Ungroup after summarizing
  ) |>
  mutate(meters = Lngth_Yrd*0.9144, trout_100m = round((counts/meters) * 100)) |>
  st_as_sf(coords = c("Long", "Lat"), crs = 4326) |> # Convert to sf object using longitude and latitude columns
  st_transform(geometry, crs = 3857)
  #filter(trout_100m > 1) #attempting this to see how big of a difference these higher desnities make in the model

# Summarize the snorkel dataframe
summary(snorkel)

# Check the structure of the dataframe
str(snorkel)

# Count missing values in each column
colSums(is.na(snorkel))

# Check unique values for key columns
unique(snorkel$Reach)
#unique(snorkel$sizeClass)

# Summary statistics for count
summary(snorkel$counts)

# Check the distribution of count
hist(snorkel$trout_100m, breaks = 20, main = "Distribution of Count", xlab = "Count", col = "blue")

# Boxplot of count
boxplot(snorkel$counts, main = "Boxplot of Count", ylab = "Count")

# Boxplot of count by Reach
ggplot(snorkel, aes(x = Reach, y = counts)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Count by Reach", x = "Reach", y = "Count")

# Scatterplot of count vs Flow
ggplot(snorkel, aes(x = flow, y = counts)) +
  geom_point() +
  geom_smooth(method = "lm", color = "red") +
  labs(title = "Count vs Flow", x = "Flow", y = "Count")

ggplot(snorkel, aes(x = flow, y = counts, color = Reach)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(title = "Count vs Flow by Reach", x = "Flow", y = "Count")

ggplot(snorkel, aes(x = Lngth_Yrd, y = counts, color = Reach)) +
  geom_point() +
  #geom_smooth(method = "lm") +
  labs(title = "Count vs Length by Reach", x = "Length of snorkel", y = "Count")

# Assuming stan_ref_points and snorkel are both sf objects with a geometry column
# Ensure they have the same CRS for comparison
stan_ref_points <- st_transform(stan_ref_points, crs = st_crs(snorkel))

# Step 1: Identify matching points
# Use st_equals or st_intersects to find matches
matching_indices <- st_equals(stan_ref_points, snorkel)

# Step 2: Filter out matching points
# st_equals returns a list, so check which indices have matches
non_matching_indices <- which(lengths(matching_indices) == 0)

# Create a new sf object with only non-matching points
filtered_stan_ref_points <- stan_ref_points[non_matching_indices, ]

# Step 3: Inspect the result
print(filtered_stan_ref_points)

# Optional: Save or export the filtered points
# st_write(filtered_stan_ref_points, "filtered_stan_ref_points.shp")

# Define a temporary directory for storing the Local Stream Network (LSN) files
lsn.path <- paste0(tempdir(), "/lsn")  
# Creates a temporary directory to store the LSN. 
# Using tempdir ensures that the directory is unique for this session.

# Convert the river network (streams) into an LSN object
edges <- lines_to_lsn(
  streams = stan_flow_line,   # Input river network (spatial lines for the Stanislaus River)
  lsn_path = lsn.path,        # Path to save the LSN object
  check_topology = TRUE,      # Ensure the network topology is valid (e.g., connections between edges are correct)
  snap_tolerance = 0.05,      # Maximum distance (in CRS units) to snap disconnected segments
  topo_tolerance = 20,        # Maximum tolerance for correcting topology errors
  overwrite = TRUE            # Overwrite existing files in the target directory
)

# Transform the CRS of the reference points (stan_ref_points) to match the CRS of the edges
stan_ref_points <- sf::st_transform(stan_ref_points, sf::st_crs(edges))
# The transformation ensures that the CRS (Coordinate Reference System) of the reference points aligns with the river network,
# which is essential for spatial operations like snapping, distance calculations, and upstream network analyses.


# Create prediction points by snapping reference points (stan_ref_points) to the river network (edges)
preds <- sites_to_lsn(
  sites = stan_ref_points,  # Reference points for predictions
  edges = edges,            # River network edges
  save_local = TRUE,        # Save output locally
  lsn_path = lsn.path,      # Path to local SSN directory
  file_name = "preds.gpkg", # Output file name for predictions
  snap_tolerance = 100,     # Snap points within 100 meters of the network
  overwrite = TRUE          # Overwrite existing files if they exist
)

# Transform snorkel data to match the CRS of the river network
snorkel <- sf::st_transform(snorkel, sf::st_crs(edges))

# Create observation points by snapping snorkel sampling points to the river network
obs <- sites_to_lsn(
  sites = snorkel,          # Observation points (snorkel data)
  edges = edges,            # River network edges
  lsn_path = lsn.path,      # Path to local SSN directory
  file_name = "obs",        # Output file name for observations
  snap_tolerance = 100,     # Snap points within 100 meters of the network
  save_local = TRUE,        # Save output locally
  overwrite = TRUE          # Overwrite existing files if they exist
)

# Compute upstream distances for each edge in the river network
edges <- updist_edges(
  edges = edges,            # River network edges
  save_local = TRUE,        # Save output locally
  lsn_path = lsn.path,      # Path to local SSN directory
  calc_length = TRUE        # Calculate segment lengths
)

# View the column names in the edges dataset
names(edges)  # Useful for inspecting available attributes

# Compute upstream distances for sites (obs and preds) based on river network
site.list <- updist_sites(
  sites = list(
    obs = obs,              # Observation points
    preds = preds           # Prediction points
  ),
  edges = edges,            # River network edges
  length_col = "Length",    # Column indicating edge lengths
  save_local = TRUE,        # Save output locally
  lsn_path = lsn.path       # Path to local SSN directory
)

# View the names of datasets in the site list
names(site.list)  # Useful for inspecting available datasets

# Check column names for observation points in the site list
names(site.list$obs)

# Plot the river network with upstream distances
ggplot() +
  geom_sf(data = edges, aes(color = upDist)) +             # Color edges by upstream distance
  geom_sf(data = site.list$obs, aes(color = upDist)) +     # Color observation points by upstream distance
  coord_sf(datum = st_crs(stan_flow_line)) +               # Coordinate reference system
  scale_color_viridis_c()                                  # Viridis color scale for continuous values

# Summarize AreaSqKM column in the edges dataset
summary(edges$AreaSqKM)  # Check for zeros or anomalies in the drainage area column

# Compute additive flow variance (AFV) for the edges
edges <- afv_edges(
  edges = edges,            # River network edges
  infl_col = "AreaSqKM",    # Column for drainage area (influencing flow variance)
  segpi_col = "areaPI",     # Column for proportional influence
  afv_col = "afvArea",      # Column name for output AFV
  lsn_path = lsn.path       # Path to local SSN directory
)

# View column names in the edges dataset
names(edges)

# Summarize the AFV column to check for zeros or anomalies
summary(edges$afvArea)

# Assign AFV to observation and prediction sites
site.list <- afv_sites(
  sites = site.list,        # Observation and prediction sites
  edges = edges,            # River network edges with AFV
  afv_col = "afvArea",      # AFV column name
  save_local = TRUE,        # Save output locally
  lsn_path = lsn.path       # Path to local SSN directory
)

# View column names in the observation dataset
names(site.list$obs)

# Summarize prediction dataset for AFV values
summary(site.list$preds)

# Assemble the Spatial Stream Network (SSN) object
stan_ssn <- ssn_assemble(
  edges = edges,            # River network edges
  lsn_path = lsn.path,      # Path to local SSN directory
  obs_sites = site.list$obs,# Observation sites
  preds_list = site.list["preds"], # List of prediction datasets
  ssn_path = "/stan.ssn",   # Path to save assembled SSN
  import = TRUE,            # Import the assembled SSN into R
  check = TRUE,             # Check for errors in the assembly process
  afv_col = "afvArea",      # AFV column name
  overwrite = TRUE          # Overwrite existing SSN if it exists
)


class(stan_ssn) ## Get class

names(stan_ssn)

names(stan_ssn$preds) ## print names of prediction datasets

# Step 1: Ensure preds is a data.frame
if (is.list(stan_ssn$preds)) {
  preds_df <- do.call(rbind, stan_ssn$preds) |> as.data.frame()
} else {
  preds_df <- stan_ssn$preds
}

# Step 2: Drop geometry if preds is an sf object
if ("geometry" %in% colnames(preds_df)) {
  preds_df <- st_drop_geometry(preds_df)
}

# Step 3: Select covariates
covariates <- preds_df |>
  dplyr::select(rid, AreaSqKM, TotDASqKM, DUP_Length, uniqueID, Position, slope_cm_p, slope_cm_1, Sep) 

# View covariates
print(head(covariates))

# Step 3: Join covariates to obs based on locID
stan_ssn$obs <- stan_ssn$obs |>
  left_join(covariates, by = "rid") # Perform the join

# Step 4: Validate the result
print(stan_ssn$obs)

ggplot() +
  # Add the river network (edges) to the plot
  geom_sf(
    data = stan_ssn$edges,
    color = "medium blue",
    aes(linewidth = AreaSqKM)
  )+
  # Add observation points to the plot
  geom_sf(
    data = stan_ssn$obs,
    size = 1.7,
    aes(color = trout_100m)
  ) +
  coord_sf(datum = st_crs(stan_flow_line)) +
  scale_color_viridis_c() +
  labs(color = "trout_100m", linewidth = "WS Area") +
  theme(
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10)
  ) +
  scale_linewidth(range = c(0.1, 2.5)) +
  # Add prediction points to the plot
  geom_sf(
    data = stan_ssn$preds$preds,
    size = 1.5,
    shape = 21,
    fill = "white",
    color = "dark grey"
  ) 

## Generate hydrologic distance matrices
ssn_create_distmat(stan_ssn, predpts = "preds", overwrite = TRUE)

# Create a data frame
data <- data.frame(
  TotDASqKM = scale(stan_ssn$obs$TotDASqKM),
  slope_cm_p = scale(stan_ssn$obs$slope_cm_p),
  AreaSqKM = scale(stan_ssn$obs$AreaSqKM),
  Sep = scale(stan_ssn$obs$Sep)
)

# Calculate the correlation matrix
cor_matrix <- cor(data)

# Print the correlation matrix
print(cor_matrix)

# Define the models-- covariance parameterizations were the same as Isaak et al. (2016)
ssn_1 <- ssn_glm(
  formula = trout_100m ~ slope_cm_p + Sep,
  ssn.object = stan_ssn,
  tailup_type = "none",
  taildown_type = "exponential",
  euclid_type = "exponential",
  family = "nbinomial",
  nugget_type = "none",
  random = ~ as.factor(Date)
)

ssn_2 <- ssn_glm(
  formula = trout_100m ~ slope_cm_p + Sep,
  ssn.object = stan_ssn,
  tailup_type = "none",
  taildown_type = "exponential",
  euclid_type = "exponential",
  family = "nbinomial",
  nugget_type = "nugget"
)

# Combine models into a list
models <- list(ssn_1 = ssn_1, ssn_2 = ssn_2)

# Extract AIC values
model_aic <- sapply(models, AIC)

# Create the summary table
model_summary <- data.frame(
  Model = names(models),
  AIC = model_aic
) %>%
  mutate(Delta_AIC = AIC - min(AIC))  # Calculate ΔAIC

# View the summary table
print(model_summary)

SSN2::loocv(ssn_1)
SSN2::loocv(ssn_2)

summary(ssn_1)

tidy(ssn_1, conf.int = TRUE)

# Extract the obs component from stan_ssn and ensure it's a data.frame
preds_df <- as.data.frame(stan_ssn$preds)

# Check the structure of the new data.frame
str(preds_df)

preds_mod <- predict(ssn_1, newdata = "preds", interval = "prediction")
head(preds_mod)

# Extract residuals
residuals_ssn_1 <- residuals(ssn_1, type = "response")

# Summary of residuals
summary(residuals_ssn_1)

# Standardized residuals
rstandard_ssn_1 <- rstandard(ssn_1)

# Summary of standardized residuals
summary(rstandard_ssn_1)

hist(residuals_ssn_1, breaks = 30, main = "Residual Histogram", xlab = "Residuals")

fitted_values <- fitted(ssn_1)
plot(fitted_values, residuals_ssn_1,
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted Values")
abline(h = 0, col = "red")

qqnorm(residuals_ssn_1, main = "Q-Q Plot of Residuals")
qqline(residuals_ssn_1, col = "red")

# Extract leverage values
hat_values <- hatvalues(ssn_1)

# Plot leverage values
plot(hat_values, main = "Leverage Values", xlab = "Index", ylab = "Leverage")
abline(h = 2 * mean(hat_values), col = "red", lty = 2)  # Rule of thumb for high leverage

# Influence measures
influence_measures <- influence(ssn_1)

# View influence diagnostics
head(influence_measures)

# Plot a Torgegram to examine spatial correlation in residuals
#torgegram_plot <- plot(Torgegram(stan_ssn))

# Combine observed and fitted values
comparison <- data.frame(
  observed = stan_ssn$obs$trout_100m,
  fitted = fitted_values
)

# Scatter plot
ggplot(comparison, aes(x = observed, y = fitted)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(title = "Observed vs Fitted Values", x = "Observed", y = "Fitted")

################################################
#Fitting zero-inflated model--much better fit
##############################################
library(pscl)

# Fit Zero-Inflated Negative Binomial model
zinb_model <- zeroinfl(trout_100m ~ slope_cm_p + Sep | 1, data = stan_ssn$obs, dist = "negbin")

# Summary of the zero-inflated model
summary(zinb_model)

# Extract fitted values
fitted_values <- fitted(zinb_model)

# Create a comparison data frame
comparison <- data.frame(
  observed = stan_ssn$obs$trout_100m,
  fitted = fitted_values
)

# Scatterplot: Observed vs Fitted
ggplot(comparison, aes(x = observed, y = fitted)) +
  geom_point(alpha = 0.7, color = "blue") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(
    title = "Observed vs Fitted Values (Negative Binomial)",
    x = "Observed Counts",
    y = "Fitted Counts"
  ) +
  theme_minimal()

##################
##################
#Bayesian SSN Example
##################
##################
options(mc.cores = parallel::detectCores())

# Extract the obs component from stan_ssn and ensure it's a data.frame
obs_df <- as.data.frame(stan_ssn$obs)

# Check the structure of the new data.frame
str(obs_df)

path <- system.file("extdata/clearwater.ssn", package = "SSNdata")

# Fit a Bayesian universal kriging model for steelhead counts as a function of discharge, CV of discharge, and temperature
bayesian_model <- ssnbayes(
  formula = trout_100m ~ AreaSqKM + slope_cm_p + TotDASqKM, 
  data = obs_df, 
  path = stan_ssn,
  time_method = list('ar', 'Date'),
  space_method = list('use_ssn', 'Exponential.taildown'),
  iter = 2000,
  warmup = 1000,
  chains = 4,
  addfunccol='afvArea')

# Check column names
colnames(stan_ssn$obs)

# Inspect the locID column
str(stan_ssn$obs$locID)


####GARAGE####
# Combine rows by rid and take the average of trout_100m
obs <- obs |>
  group_by(rid) |>
  summarise(
    trout_100m = mean(trout_100m, na.rm = TRUE), # Take average of trout_100m
    geometry = st_union(geometry),              # Combine geometries
    snapdist_mean = mean(snapdist),
    .groups = "drop"
  )

# View the resulting dataframe
print(head(obs))

obs <- sites_to_lsn(
  sites = obs,
  edges = edges,
  lsn_path = lsn.path,
  file_name = "obs",
  snap_tolerance = 100,
  save_local = TRUE,
  overwrite = TRUE
)

obs


# Fit models with varying configurations
for (formula in formulas) {
  for (taildown in taildown_types) {
    for (euclid in euclid_types) {
      # Generate a unique model name
      model_name <- paste(
        deparse(formula), taildown, euclid, sep = "_"
      )
      
      # Fit the model
      models[[model_name]] <- ssn_lm(
        formula = formula,
        ssn.object = stan_ssn,
        tailup_type = tailup_types,
        taildown_type = taildown,
        euclid_type = euclid,
        additive = "afvArea"
      )
    }
  }
}


# Define predictor sets
predictor_sets <- list(
  formula1 = trout_100m ~ 1,
  formula2 = trout_100m ~ slope_cm_p,
  formula3 = trout_100m ~ slope_cm_p + Sep
)

# Define taildown types
taildown_types <- c("none", "exponential", "spherical")

# Define Euclidean types
euclid_types <- c("none", "exponential", "gaussian")

# Create an empty list to store models
models <- list()

# Counter for model names
model_id <- 1

library(purrr)

# Define a grid of configurations based on your hypothesis
configurations <- expand.grid(
  predictors = predictor_sets,     # Predictor sets
  taildown = c("exponential", "spherical"),  # Taildown types
  euclid = c("none", "exponential"),        # Euclidean types
  nugget = c("none", "nugget"),             # Nugget types
  stringsAsFactors = FALSE
)

# Fit models for all configurations
models <- pmap(
  configurations,
  ~ ssn_glm(
    formula = ..1,                 # Predictor set
    ssn.object = stan_ssn,         # SSN object
    tailup_type = "none",          # Tailup omitted
    taildown_type = ..2,           # Taildown type
    euclid_type = ..3,             # Euclidean type
    nugget_type = ..4,             # Nugget type
    family = "nbinomial"           # Negative binomial family
  )
)

# Name models for easier identification
names(models) <- apply(configurations, 1, function(row) {
  paste("Predictor:", deparse(row$predictors), 
        "| Taildown:", row$taildown, 
        "| Euclid:", row$euclid, 
        "| Nugget:", row$nugget)
})
