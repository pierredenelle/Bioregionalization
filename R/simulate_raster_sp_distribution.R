
## Packages -------------------------------------------------------------------
library(gstat) # virtual rasters with autocorrelation
library(sp) # spatial plots
library(cowplot) # combining plots
library(ggplot2) # plot
library(raster)
library(virtualspecies) # simulate virtual species
library(dplyr)

## Virtual raster layer -------------------------------------------------------
# http://santiago.begueria.es/2010/10/generating-spatially-correlated-random-fields-with-r/

# XY coordinates
xy <- expand.grid(1:10, 1:10)
names(xy) <- c("x","y")
# Spatial model
g_ex <- gstat(formula = z ~ 1, locations  = ~x+y, dummy = TRUE,
              beta = 1, nmax = 20,
              model = vgm(psill = 0.025, range = 5, model = "Exp"))

# Four layers created
yy <- predict(g_ex, newdata = xy, nsim = 4)

# Plottind rasters
gridded(yy) = ~x+y

plot_grid(spplot(obj = yy[1]), spplot(obj = yy[2]),
          spplot(obj = yy[3]), spplot(obj = yy[4]),
          nrow = 2)

# Data.frame with pixel and environmental value
env_dat <- cbind(yy[1]@coords, yy[1]@data)
env_dat$site <- paste0("site", 1:nrow(env_dat))

## Virtual species ------------------------------------------------------------
# http://borisleroy.com/files/virtualspecies-tutorial.html

# Generating 100 species with mean all along Bio1 gradient
env_value <- env_dat$sim1
mean_gdt <- sort(rep(seq(min(env_value), max(env_value), length.out = 10), 10))
sd_gdt <- rep(seq(var(env_value), sd(env_value), length.out = 10), 10)

# Generating 100 species with 500 occurrences
sampled_points <- rep(500, length(env_dat$sim1))
env_raster <- stack(yy[1])

sp_df <- c()
for (i in 1:100){
  param_i <- formatFunctions(
    sim1 = c(fun = "dnorm", mean = mean_gdt[i], sd = sd_gdt[i]))
  # Generation of the virtual species
  sp_i <- generateSpFromFun(raster.stack = env_raster,
                            parameters = param_i, plot = TRUE)
  # plotResponse(sp_i)

  # Conversion to presence/absence
  pa_i <- convertToPA(sp_i, plot = FALSE, beta = 0.7)

  occ_i <- data.frame(sp = paste0("sp", i),
                      rasterToPoints(pa_i$suitab.raster),
                      pa = data.frame(rasterToPoints(pa_i$pa.raster))$layer)
  colnames(occ_i)[colnames(occ_i) == "layer"] <- "suitab"

  # Binding results
  sp_df <- rbind(sp_df, occ_i)
  print(i)
}

# Removing lines with null abundances
sp_df <- sp_df[which(sp_df$suitab != 0), ]

# Add site and environment columns
sp_df <- left_join(sp_df, env_dat, by = c("x", "y"))

# Renaming environmental column
colnames(sp_df)[colnames(sp_df) == "sim1"] <- "env"

# If suitability inferior to 0.15, suitab = 0 (make some species absent)
# sp_df$suitab <- ifelse(sp_df$suitab < 0.15, 0, sp_df$suitab)
# Remove absent species
# sp_df <- sp_df[which(sp_df$suitab != 0), ]

# Save sp_df as RData
save(sp_df, file = "data/virtual_sp.RData")

# Map of species richness

