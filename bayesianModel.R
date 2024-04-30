# Build and run Bayesian GLMs for three stability metrics with ecosystem characteristics, dsiturbance legacy and climate legacy as explanaotry variables. 
# Author: Liezl M. Vermeulen

# Install and load required packages
library(spBayes) 
library (raster)
library(rgdal)
library (rts)
library(terra)
library(spdep)
library(brms)
library(ggplot2)
library(ggpubr)
library(MASS)
library(lme4)
library(gtools)
library(devtools)
library(DHARMa)
library(dplyr)
library(tidyr)
library(car)
library(nlme)
library(spatialreg)
library(glmmTMB)
library(rsample)
library(tidybayes)
library(effects)

# Set working directory
wd.eco <- 'path/to/ecosystem/characteristics'
wd.dist <- 'path/to/dsiturbance/legacy'
wd.climate <- 'path/to/climate/legacy'
wd.response <- 'path/to/response/metrics'

# Load your study area shapefile
knp <- shapefile('path/to/study/area')

carbon.factor <- 10
sand.factor <- 10
pH.factor <- 10

# Ecosystem characteristics
soil.carbon <- raster(file.path(wd.eco, "soil/SoilGrids/soilCarbon_250m_KNP_utm.tif")) / carbon.factor
soil.pH <- raster(file.path(wd.eco, "soil/SoilGrids/soilpH_250m_KNP_utm.tif")) / pH.factor
soil.sand <- raster(file.path(wd.eco, "soil/SoilGrids/soilSand_250m_KNP_utm.tif")) / sand.factor
geology <- shapefile(file.path(wd.eco, "soil/basalt_granite.shp"))
elev <- raster(file.path(wd.eco, "terrain/elev_knp_utm.tif"))
slope <- raster(file.path(wd.eco, "terrain/slope_knp_utm.tif"))
twi <- raster(file.path(wd.eco, "terrain/twi_knp_utm.tif"))

# Climate legacy
drought <- stack(file.path(wd.climate, "spi/spi12_drought_sev_stack.tif")) # 1981 - 2015
drought.sum.2015 <- abs(sum(drought[[19:34]])) # 1981 - 2015 (for resistance/resilience)
drought.sum.full <- abs(sum(drought[[19:36]])) # 1981- 2017 (for woody cover)

#drought.sum <- abs(sum(drought[[1:36]]))
wetness <- stack(file.path(wd.climate, "spi/spi12_wetness_sev_stack.tif")) # 1981 - 2015
wetness.sum.2015 <- abs(sum(wetness[[19:34]])) # 1981 - 2015 (for resistance/resilience)
wetness.sum.full <- abs(sum(wetness[[19:42]])) # 1981- 2017 (for woody cover)

# Disturbance legacy
elephant <- raster(file.path(wd.dist, "elephants/elephant_impact_heatmap.tif")) / 100
fire.freq.files <- list.files(path=paste0(wd.dist,"/fire/output/fire_freq"),pattern="*.tif", full.names=TRUE, recursive=FALSE) 
fire.freq <- stack(fire.freq.files) 
fire.freq.sum.2015 <- sum(fire.freq[[19:72]]) # 1981 - 2015 (for resistance/resilience)
fire.freq.sum.full <- sum(fire.freq[[19:74]]) # 1981- 2017 (for woody cover)

# Response metrics
resistance <- raster(file.path(wd.response,"resistance.tif"))
resilience <- raster(file.path(wd.response,"resilience.tif"))
woodyChange <- raster(file.path(wd.response,"woody_cover/Venter_2017/woodyChangeTrend_Venter2017_buff.tif"))

# Project everything to CHIRPS res, local projection
chirps.crs <- crs(drought.sum.2015)
utm.crs <- crs(fire.freq)

base.r <- mask(projectRaster(resistance, res=c(1000,1000), crs=utm.crs, method="bilinear"), knp)

resistance.rp <- mask(projectRaster(resistance, base.r, method="bilinear"), knp)
resilience.rp <- mask(projectRaster(resilience, base.r, method="bilinear"), knp)
woody.change.rp <- mask(projectRaster(woodyChange, base.r, method="bilinear"), knp)


fire.freq.2015.rp <- round(mask(projectRaster(fire.freq.sum.2015, base.r, method="bilinear"), knp))
fire.freq.full.rp <- round(mask(projectRaster(fire.freq.sum.full, base.r, method="bilinear"), knp))
elephant.rp <- mask(projectRaster(elephant, base.r, method="bilinear"), knp)

drought.2015.rp <- mask(projectRaster(drought.sum.2015, base.r, method="bilinear"), knp)
drought.full.rp <- mask(projectRaster(drought.sum.full, base.r, method="bilinear"), knp)
wetness.2015.rp <- mask(projectRaster(wetness.sum.2015, base.r, method="bilinear"), knp)
wetness.full.rp <- mask(projectRaster(wetness.sum.full, base.r, method="bilinear"), knp)

elev.rp <- mask(projectRaster(elev, base.r, method="bilinear"), knp)
slope.rp <- mask(projectRaster(slope, base.r, method="bilinear"), knp)
twi.rp <- mask(projectRaster(twi, base.r, method="bilinear"), knp)
soil.sand.rp <- mask(projectRaster(soil.sand, base.r, method="bilinear"), knp)
soil.pH.rp <- mask(projectRaster(soil.pH, base.r, method="bilinear"), knp)
soil.carbon.rp <- mask(projectRaster(soil.carbon, base.r, method="bilinear"), knp)
geology.rp <- mask(rasterize(as(geology, "SpatVector"), rast(base.r), field="GEOLOGY"), vect(knp))

r.names <- c("resist" , "resil", "woodyChange",
             "elev", "slope", "twi", "soilSand", "soilpH", "soilCarbon", "geology",
             "fireFreq2015", "fireFreqFull", "elephant",
             "drought2015","droughtFull", "wetness2015", "wetnessFull")

r.list <- c(rast(resistance.rp), rast(resilience.rp), rast(woody.change.rp),
            rast(elev.rp), rast(slope.rp), rast(twi.rp), rast(soil.sand.rp), rast(soil.pH.rp), rast(soil.carbon.rp),geology.rp,
            rast(fire.freq.2015.rp), rast(fire.freq.full.rp), rast(elephant.rp),
            rast(drought.2015.rp), rast(drought.full.rp), rast(wetness.2015.rp), rast(wetness.full.rp))

names(r.list) <- r.names

# Plot all layers
plotAll<- lapply(r.list, function(x){
  plot(x, main=names(x))
  plot(knp, add=TRUE)
})

rivers <- mtq <- st_read("path/to/rivers.shp") %>% 
  dplyr::summarise() 

# Build data matrix
model.stack <- r.list
# the masking object:
mask <- st_bbox(model.stack) %>% # take extent of your raster... 
  st_as_sfc() %>% # make it a sf object
  st_set_crs(st_crs(mtq)) %>% # in CRS of your polygon 
  st_difference(mtq) %>% # intersect with the polygon object
  st_as_sf() # interpret as sf (and not sfc) object

model.stack.masked <- model.stack %>% 
  mask(mask)
plot(model.stack.masked)
model.df <- as.data.frame(model.stack.masked,xy=TRUE)

model.df2 <- cbind(ID = 1:nrow(model.df), model.df)
model.df2 <- na.omit(model.df2)

model.df2 <- model.df2[!(model.df2$soilSand<40),]
model.df2 <- model.df2[!(model.df2$soilpH<6),]
model.df2 <- model.df2[!(model.df2$soilCarbon>20),]

# Scale and center predictor variables
model.df2$soilCarbon2 <- scale(model.df2$soilCarbon, center = TRUE, scale = TRUE)
model.df2$soilSand2 <- scale(model.df2$soilSand, center = TRUE, scale = TRUE)
model.df2$soilpH2 <- scale(model.df2$soilpH, center = TRUE, scale = TRUE)
model.df2$slope2 <- scale(model.df2$slope, center = TRUE, scale = TRUE)
model.df2$elev2 <- scale(model.df2$elev, center = TRUE, scale = TRUE)
model.df2$twi2 <- scale(model.df2$twi, center = TRUE, scale = TRUE)
model.df2$fireFreq20152 <- scale(model.df2$fireFreq2015, center = TRUE, scale = TRUE)
model.df2$fireFreqFull2 <- scale(model.df2$fireFreqFull, center = TRUE, scale = TRUE)
model.df2$elephant2 <- scale(model.df2$elephant, center = TRUE, scale = TRUE)
model.df2$drought20152 <- scale(model.df2$drought2015, center = TRUE, scale = TRUE)
model.df2$droughtFull2 <- scale(model.df2$droughtFull, center = TRUE, scale = TRUE)
model.df2$wetness20152 <- scale(model.df2$wetness2015, center = TRUE, scale = TRUE)
model.df2$wetnessFull2 <- scale(model.df2$wetnessFull, center = TRUE, scale = TRUE)

data <- model.df2[,c(43,7,8)]
colnames(data) <- c("Resistance", "Resilience", "Woody Change")
groups <- model.df2[,17]

# Function to add correlation coefficients
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  Cor <- abs(cor(x, y)) # Remove abs function if desired
  txt <- paste0(prefix, format(c(Cor, 0.123456789), digits = digits)[1])
  if(missing(cex.cor)) {
    cex.cor <- 0.4 / strwidth(txt)
  }
  text(0.5, 0.5, txt,
       cex = 1 + cex.cor * Cor) # Resize the text by level of correlation
}

# Plotting the correlation matrix
pairs(data,
      upper.panel = panel.cor,    # Correlation panel
      lower.panel = panel.smooth) # Smoothed regression lines

library(psych)
corPlot(data, cex = 1.2)

library(corrplot)

corrplot.mixed(cor(data),
               lower = "color", 
               upper = "number",
               tl.col = "black")


corrplot(cor(data), type="upper", tl.col="black")

scale_color_manual(name="Geology", values = c("Basalt" = "darkorange4","Granite"="darkgoldenrod2")) +
  scale_fill_manual(name="Geology", values = c("Basalt" = "darkorange4","Granite"="darkgoldenrod2"))

# Bayesian model -----------------------------------------------------
set.seed(2021)
split <- rsample::initial_split(model.df2, prop = 0.37, strata = geology)
sample.df <- rsample::training(split)
test.df <- sample.df

# Create neighbourhood structure for autocorrelation
xy.df2 <- distinct(as.data.frame(cbind(test.df$ID, test.df$x, test.df$y)))
names(xy.df2) <- c("ID", "x", "y")
coordinates(xy.df2) <- ~ x + y
knea <- knearneigh(coordinates(xy.df2), longlat = FALSE)
W.nb <- knn2nb(knea, sym=TRUE)
W <- nb2mat(W.nb, style="B", zero.policy=FALSE)
rownames(W) <- unique(test.df$ID)
listW <- nb2listw(W.nb, style="W")

# Test for spatial correlation
globalMoran <- moran.test(test.df$woodyChange, listW)
globalMoran

# Resistance model -----------------------------------------------------
resist.formula <- resist ~ soilSand2 + slope2 + twi2 + fireFreq20152 + elephant2 + drought20152 + wetness20152 + geology + car(W, gr=ID)
resist_model <- brm(formula = resist.formula,  
                    data = test.df,
                    data2 = list(W=W),
                    family= Beta(link="logit"),
                    warmup = 1500, 
                    iter = 5000, 
                    chains = 4, 
                    cores=4,
                    seed = 123)
summary(resist_model)
plot(resist_model)
pp_check(resist_model)
#saveRDS(resist_model, file = "resist_model.rda")

# Extract marginal effects
soil_me <- marginal_effects(resist_model, effects="soilSand2")
soil_plot <- plot(soil_me, plot = FALSE)[[1]] + 
                theme_bw() + 
                labs(x = "Soil sand content", y = "Estimate") +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
                theme(axis.title=element_text(size=12)) +
                theme(axis.text=element_text(size=10))
soil_plot

slope_me <- marginal_effects(resist_model, effects="slope2")
slope_plot <- plot(slope_me, plot = FALSE)[[1]] + 
  theme_bw() + 
  labs(x = "Slope", y = "Estimate") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(axis.title=element_text(size=12)) +
  theme(axis.text=element_text(size=10))
slope_plot

twi_me <- marginal_effects(resist_model, effects="twi2")
twi_plot <- plot(twi_me, plot = FALSE)[[1]] + 
  theme_bw() + 
  labs(x = "Topographic wetness index", y = "Estimate") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(axis.title=element_text(size=12)) +
  theme(axis.text=element_text(size=10))
twi_plot
  
fire_me <- marginal_effects(resist_model, effects="fireFreq20152")
fire_plot <- plot(fire_me, plot = FALSE)[[1]] + 
  theme_bw() + 
  labs(x = "Fire frequency", y = "Estimate") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(axis.title=element_text(size=12)) +
  theme(axis.text=element_text(size=10))
fire_plot
  
elephant_me <- marginal_effects(resist_model, effects="elephant2")
elephant_plot <- plot(elephant_me, plot = FALSE)[[1]] + 
  theme_bw() + 
  labs(x = "Elephant density", y = "Estimate") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(axis.title=element_text(size=12)) +
  theme(axis.text=element_text(size=10))
elephant_plot

wetness_me <- marginal_effects(resist_model, effects="wetness20152")
wetness_plot <- plot(wetness_me, plot = FALSE)[[1]] + 
  theme_bw() + 
  labs(x = "Cumulative annual wetness severity", y = "Estimate") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(axis.title=element_text(size=12)) +
  theme(axis.text=element_text(size=10))
wetness_plot
  
drought_me <- marginal_effects(resist_model, effects="drought20152")
drought_plot <- plot(drought_me, plot = FALSE)[[1]] + 
  theme_bw() + 
  labs(x = "Cumulative annual drought severity", y = "Estimate") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(axis.title=element_text(size=12)) +
  theme(axis.text=element_text(size=10))
drought_plot
  
geology_me <- marginal_effects(resist_model, effects="geology")
geology_plot <- plot(geology_me, plot = FALSE)[[1]] + 
  theme_bw() + 
  labs(x = "Geology", y = "Estimate") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(axis.title=element_text(size=12)) +
  theme(axis.text=element_text(size=10))
geology_plot

figure <- ggarrange(soil_plot, slope_plot, twi_plot, fire_plot, elephant_plot, wetness_plot, drought_plot, geology_plot,
                    labels = c("a)", "b)", "c)", "d)", "e)", "f)", "g)", "h)"),
                    ncol = 2, nrow = 4)
figure

# Resilence model -----------------------------------------------------
resil.formula <- resil ~ soilSand2 + slope2 + twi2 + fireFreq20152 + elephant2 + drought20152 + wetness20152 + geology + car(W, gr=ID)
resil_model <- brm(formula = resil.formula,  
                   data = test.df,
                   data2 = list(W=W),
                   family=Beta(link="logit"),
                   warmup = 1500, 
                   iter = 5000, 
                   chains = 4, 
                   cores= 4,
                   seed = 123)
summary(resil_model)
plot(resil_model)
pp_check(resil_model)
# saveRDS(resil_model, file = "resil_model.rda")

# Extract marginal effects
soil_me <- marginal_effects(resil_model, effects="soilSand2")
soil_plot <- plot(soil_me, plot = FALSE)[[1]] + 
  theme_bw() + 
  labs(x = "Soil sand content", y = "Estimate") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(axis.title=element_text(size=12)) +
  theme(axis.text=element_text(size=10))
soil_plot

slope_me <- marginal_effects(resil_model, effects="slope2")
slope_plot <- plot(slope_me, plot = FALSE)[[1]] + 
  theme_bw() + 
  labs(x = "Slope", y = "Estimate") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(axis.title=element_text(size=12)) +
  theme(axis.text=element_text(size=10))
slope_plot

twi_me <- marginal_effects(resil_model, effects="twi2")
twi_plot <- plot(twi_me, plot = FALSE)[[1]] + 
  theme_bw() + 
  labs(x = "Topographic wetness index", y = "Estimate") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(axis.title=element_text(size=12)) +
  theme(axis.text=element_text(size=10))
twi_plot

fire_me <- marginal_effects(resil_model, effects="fireFreq20152")
fire_plot <- plot(fire_me, plot = FALSE)[[1]] + 
  theme_bw() + 
  labs(x = "Fire frequency", y = "Estimate") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(axis.title=element_text(size=12)) +
  theme(axis.text=element_text(size=10))
fire_plot

elephant_me <- marginal_effects(resil_model, effects="elephant2")
elephant_plot <- plot(elephant_me, plot = FALSE)[[1]] + 
  theme_bw() + 
  labs(x = "Elephant density", y = "Estimate") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(axis.title=element_text(size=12)) +
  theme(axis.text=element_text(size=10))
elephant_plot

wetness_me <- marginal_effects(resil_model, effects="wetness20152")
wetness_plot <- plot(wetness_me, plot = FALSE)[[1]] + 
  theme_bw() + 
  labs(x = "Cumulative annual wetness severity", y = "Estimate") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(axis.title=element_text(size=12)) +
  theme(axis.text=element_text(size=10))
wetness_plot

drought_me <- marginal_effects(resil_model, effects="drought20152")
drought_plot <- plot(drought_me, plot = FALSE)[[1]] + 
  theme_bw() + 
  labs(x = "Cumulative annual drought severity", y = "Estimate") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(axis.title=element_text(size=12)) +
  theme(axis.text=element_text(size=10))
drought_plot

geology_me <- marginal_effects(resil_model, effects="geology")
geology_plot <- plot(geology_me, plot = FALSE)[[1]] + 
  theme_bw() + 
  labs(x = "Geology", y = "Estimate") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(axis.title=element_text(size=12)) +
  theme(axis.text=element_text(size=10))
geology_plot

figure <- ggarrange(soil_plot, slope_plot, twi_plot, fire_plot, elephant_plot, wetness_plot, drought_plot, geology_plot,
                    labels = c("a)", "b)", "c)", "d)", "e)", "f)", "g)", "h)"),
                    ncol = 2, nrow = 4)
figure

# Woody change model -----------------------------------------------------
woodyChange.formula <- woodyChange ~ soilSand2 + slope2 + twi2 + fireFreq20152 + elephant2 + drought20152 + wetness20152 + geology + car(W, gr=ID)
woodyChange_model <- brm(formula = woodyChange.formula,  
                         data = test.df,
                         data2 = list(W=W),
                         family=gaussian(),
                         warmup = 1500, 
                         iter = 5000, 
                         chains = 4, 
                         cores=4,
                         seed = 123)
woodyChange_model
summary(woodyChange_model)
plot(woodyChange_model)
pp_check(woodyChange_model)
saveRDS(resil_model, file = "C:/Users/u0142455/Documents/PhD/Processing/ch2/models/woody_model.rda")

# Extract marginal effects
soil_me <- marginal_effects(woodyChange_model, effects="soilSand2")
soil_plot <- plot(soil_me, plot = FALSE)[[1]] + 
  theme_bw() + 
  labs(x = "Soil sand content", y = "Estimate") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(axis.title=element_text(size=12)) +
  theme(axis.text=element_text(size=10))
soil_plot

slope_me <- marginal_effects(woodyChange_model, effects="slope2")
slope_plot <- plot(slope_me, plot = FALSE)[[1]] + 
  theme_bw() + 
  labs(x = "Slope", y = "Estimate") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(axis.title=element_text(size=12)) +
  theme(axis.text=element_text(size=10))
slope_plot

twi_me <- marginal_effects(woodyChange_model, effects="twi2")
twi_plot <- plot(twi_me, plot = FALSE)[[1]] + 
  theme_bw() + 
  labs(x = "Topographic wetness index", y = "Estimate") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(axis.title=element_text(size=12)) +
  theme(axis.text=element_text(size=10))
twi_plot

fire_me <- marginal_effects(woodyChange_model, effects="fireFreq20152")
fire_plot <- plot(fire_me, plot = FALSE)[[1]] + 
  theme_bw() + 
  labs(x = "Fire frequency", y = "Estimate") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(axis.title=element_text(size=12)) +
  theme(axis.text=element_text(size=10))
fire_plot

elephant_me <- marginal_effects(woodyChange_model, effects="elephant2")
elephant_plot <- plot(elephant_me, plot = FALSE)[[1]] + 
  theme_bw() + 
  labs(x = "Elephant density", y = "Estimate") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(axis.title=element_text(size=12)) +
  theme(axis.text=element_text(size=10))
elephant_plot

wetness_me <- marginal_effects(woodyChange_model, effects="wetness20152")
wetness_plot <- plot(wetness_me, plot = FALSE)[[1]] + 
  theme_bw() + 
  labs(x = "Cumulative annual wetness severity", y = "Estimate") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(axis.title=element_text(size=12)) +
  theme(axis.text=element_text(size=10))
wetness_plot

drought_me <- marginal_effects(woodyChange_model, effects="drought20152")
drought_plot <- plot(drought_me, plot = FALSE)[[1]] + 
  theme_bw() + 
  labs(x = "Cumulative annual drought severity", y = "Estimate") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(axis.title=element_text(size=12)) +
  theme(axis.text=element_text(size=10))
drought_plot

geology_me <- marginal_effects(woodyChange_model, effects="geology")
geology_plot <- plot(geology_me, plot = FALSE)[[1]] + 
  theme_bw() + 
  labs(x = "Geology", y = "Estimate") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(axis.title=element_text(size=12)) +
  theme(axis.text=element_text(size=10))
geology_plot

figure <- ggarrange(soil_plot, slope_plot, twi_plot, fire_plot, elephant_plot, wetness_plot, drought_plot, geology_plot,
                    labels = c("a)", "b)", "c)", "d)", "e)", "f)", "g)", "h)"),
                    ncol = 2, nrow = 4)
figure


