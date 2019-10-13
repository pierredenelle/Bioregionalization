
## Packages -------------------------------------------------------------------
library(gstat) # virtual rasters with autocorrelation
library(sp) # spatial plots
library(cowplot) # combining plots
library(ggplot2) # plot
library(raster)
library(virtualspecies) # simulate virtual species

## Virtual raster layer -------------------------------------------------------
# http://santiago.begueria.es/2010/10/generating-spatially-correlated-random-fields-with-r/

# XY coordinates
xy <- expand.grid(1:100, 1:100)
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

# yy2 <- rasterToPoints(yy[1])

## Virtual species ------------------------------------------------------------
# http://borisleroy.com/files/virtualspecies-tutorial.html

# Generating 100 species with mean all along Bio1 gradient
env_value <- yy[1]@data$sim1
mean_gdt <- sort(rep(seq(min(env_value), max(env_value), length.out = 10), 10))
sd_gdt <- rep(seq(sd(env_value), 2*sd(env_value), length.out = 10), 10)

# Generating 100 species with 500 occurrences
sampled_points <- rep(500, length(yy[1]@data$sim1))
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
  # pa_i <- convertToPA(sp_i, plot = FALSE, beta = 0.7)

  occ_i <- data.frame(rasterToPoints(sp_i$suitab.raster))
  occ_i$sp <- paste0("sp", i)
  # Binding results
  sp_df <- rbind(sp_df, occ_i)
  print(i)
}

colnames(sp_df)[colnames(sp_df) == "layer"] <- "ab"

# Removing lines with null abundances
sp_df <- sp_df[which(sp_df$ab != 0), ]

# Add site column
xy$site <- paste0("site", 1:nrow(xy))
sp_df <- left_join(sp_df, xy, by = c("x", "y"))

# If suitability inferior to 0.15, ab = 0 (make some species absent)
sp_df$ab <- ifelse(sp_df$ab < 0.15, 0, sp_df$ab)

# Remove absent species
sp_df <- sp_df[which(sp_df$ab != 0), ]

# Save sp_df as RData
save(sp_df, file = "data/virtual_sp.RData")

# Map of species richness

