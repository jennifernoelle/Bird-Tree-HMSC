# Uses interactions to improve occurrene predictions
# Atlantic frugivory data
# Jennifer Kampe
# May 18, 2022


### To do:
  ## Modeling/Debugging
    # consider cleaning Rebeirovegtype if we're using it
    # Weird error when i specify longlat = TRUE in hmsc random level
    # Consider changing the spatial covariancce with nn or predicted GP (sData = xycoords, method =...), better for large number of locations
    # Add second order terms for all continuous covariates (use poly(x, degree = 2, raw = TRUE))
    # Add Otso batch code to run multiple at once

  ## Data
    # Use otso weather code to get weather at exact locations
    # Reduce habitat and ecosystem vars with vars which we can obtain on a grid to enable new locations: geobr?
    # Import Otso's landuse data X 
    # Get rainfall/temperature data on a grid, take 10 year averages or something
    # Get bird phylogeny  from birdtree.org
    # Or get posterior estimates of interaction prob from Georgia's model
    # See if Georgia's code has bird and tree phylogeny in it

 ## Questions
    # Do I want predictEtaMean = TRUE in predictions?
    # Is there a command to get posterior means or just arrays of predictions for each post sample?
    # Tree covariates: poserior means or binary?

# Load packages 
library(Hmsc)
library(tidyverse)
library(magrittr)
library(tibble)
library(igraph)
library(zoo)
library(geobr)
library(abind)
library(vioplot)
library(colorspace)

# Load data
birds <- read.csv("./Data/atlantic_birds/birds_wide_genus_sp_covars.csv", header = TRUE)
taxa <- read.csv("./Data/atlantic_birds/taxa_and_commonnames.csv", header = TRUE)

trees <- read.csv("./Data/Trees/trees_wide_covars.csv", header = TRUE)
trees_T <- read.csv("./Data/Trees/T_tree_co.csv", header = TRUE, sep = ";") # traits for 1016 tree sp 
trees_X <- read.csv("./Data/Trees/X_tree_co.csv", header = TRUE, sep = ";") # environmental covars at 574 sites

frug <- read.csv("./Data/Frugivory/frug_clean.csv", header = TRUE)

#################### STEP 0: VISUALIZATION ##################################


### Visualize locations and select appropriate subset of birds data

# Load Brazilian states
states <- read_state(year=2020)

# Plot all Brazilian states with sites colored by ecosystem type
birds_ecosys <- birds %>% group_by(Longitude_x, Latitude_y, Olsoneconame) %>%
  summarize(n = n()) %>%
  rename(Ecosystem = Olsoneconame)

ggplot() +
  geom_sf(data=states, fill="white", color="black", size=.15, show.legend = FALSE) +
  geom_point(data = birds_ecosys, aes(x = Longitude_x, y=Latitude_y, color = Ecosystem, size = n), alpha = 0.3)  + 
  labs(subtitle="Sampling Sites by Ecosystem and Observations", size=12) +
  theme_minimal() +
  scale_size(guide = "none") 


##################### STEP 1: USE HMSC TO PREDICT TREES AT BIRD LOCATIONS ##################

### DATA PREPARATION
apply(trees_X, 2, function(x) sum(is.na(x))) # all covars present

### SPECIES OCCURRENCE: select just the most common trees and convert to binary
trees_occ <- trees[, which(colnames(trees)=="Abarema.cochliacarpos"):which(colnames(trees) =="Zygia.latifolia")]
trees_occ <- ifelse(trees_occ > 0, 1, 0)
ntrees <- 25
commontrees <- colnames(trees_occ[, order(colSums(-trees_occ))])[1:ntrees]
Y.trees <- as.matrix(trees_occ[, which(colnames(trees_occ) %in% commontrees)])
Y.trees <- Y.trees[, order(colnames(Y.trees))] # make sure it's alphabetized

### STUDY DESIGN AND SPATIAL STRUCTURE: if sites are too close together it causes numerical instability
# Create small geographical blocks and assign the same xy-coords to each study in the block

## DIVIDE THE DATA INTO BLOCKS, ADD LOCATIONS FOR CENTERS
max.lat.trees <- ceiling(max(trees$Latitude_y))
min.lat.trees <- floor(min(trees$Latitude_y))
max.long.trees <- ceiling(max(trees$Longitude_x))
min.long.trees <- floor(min(trees$Longitude_x))
n_groups <- 100

# Groups go from left bottom to left top, next left bottom to next left top, ...
trees.df <- within(trees, {
  grp.lat = cut(Latitude_y, seq(min.lat.trees, max.lat.trees, 
                                (max.lat.trees - min.lat.trees)/n_groups), labels = FALSE)
  grp.long = cut(Longitude_x, seq(min.long.trees, max.long.trees, 
                                (max.long.trees - min.long.trees)/n_groups), labels = FALSE)}) %>%
  mutate(block.id = as.numeric(interaction(grp.lat, grp.long)))

## Compute centers and attach block ids to block centers
grid.centers <- data.frame(center_lat_y = rollsum(seq(min.lat.trees, max.lat.trees, 
                                  (max.lat.trees - min.lat.trees)/n_groups), 2)/2, 
                      center_long_x = rollsum(seq(min.long.trees, max.long.trees, 
                     (max.long.trees - min.long.trees)/n_groups), 2)/2)

# Blocks move from bottom left to top left, ..., to top right
# For each long, we visit all lats, then next long
centers_long <- rep(grid.centers$center_long_x, times = rep(n_groups, times = n_groups))
centers_lat <- rep(grid.centers$center_lat_y, times = n_groups)
block_centers <- data.frame(cbind(centers_long, centers_lat, 1:(n_groups^2)))
colnames(block_centers) <- c("centers_long", "centers_lat", "block.id")

## Plot to make sure this is going right: don't plot with full number of blocks (10,000 - too many for gg)
# ggplot() + 
#   geom_point(data = trees.df, aes(x = Longitude_x, y = Latitude_y, color = as.factor(block.id))) +
#   geom_point(data = block_centers, aes(x= centers_long, y = centers_lat, fill = as.factor(block.id) ), 
#              shape = 21, colour = "black") + 
#   scale_y_continuous(breaks = seq(min.lat.trees, max.lat.trees, 
#                                   (max.lat.trees - min.lat.trees)/n_groups)) +
#   scale_x_continuous(breaks = seq(min.long.trees, max.long.trees, 
#       (max.long.trees - min.long.trees)/n_groups))
#   

# Merge block centers with trees
trees.df <- left_join(trees.df, block_centers)

## STUDY DESIGN: n_obs rows assigning each obs to a study design unit
studyDesign.trees = data.frame(trees.df[, "block.id"])
colnames(studyDesign.trees) <- c("block.id")
studyDesign.trees = data.frame(lapply(studyDesign.trees,as.factor))
block.id = HmscRandomLevel(units = studyDesign.trees$block.id)
ranlevels = list(block.id = block.id) 
ranlevels

# SPATIAL RANDOM EFFECTS: length(levels(study design units)) rows giving location of each unit
# rownames are study design units
xycoords.trees <- trees.df %>% 
                  mutate(block.id = as.factor(block.id)) %>%
                  group_by(block.id) %>%
                  # mean just collapses becaus loc constant within groups
                  summarize(centers_x_long = mean(centers_long), 
                            centers_y_lat = mean(centers_lat)) %>%
                  tibble::column_to_rownames(var = "block.id") 

rL.spatial.trees = HmscRandomLevel(sData = xycoords.trees)
rL.spatial.trees = setPriors(rL.spatial.trees,nfMin=1,nfMax=1) 

# Environmental features: must exist at birds locations
X.env.trees <- trees_X[, c( 
                "Precipitation", "Temperature")]
                #"Declivity", "Northness", "Eastness",
                #"Bioclimatic_stress", "Human_influence", "Soil_quality"
            
  #mutate(Soil_quality = as.factor(Soil_quality))
XFormula.trees = ~  poly(Temperature,degree = 2, raw=TRUE) + poly(Precipitation,degree = 2, raw=TRUE)
                   # Declivity + Northness + Eastness + Bioclimatic_stress + Human_influence + Soil_quality


# Traits
Tree.species.full <- colnames(trees.df[,which(colnames(trees.df)=="Abarema.cochliacarpos"):which(colnames(trees.df) =="Zygia.latifolia")])
Tr.temp.trees <- trees_T %>%
  set_rownames(Tree.species.full) %>% 
  filter(rownames(.) %in% commontrees)  %>%
  select(Maximum_height, Leaf_area, Seed_length, Wood_density, Leaf_type) 
Tr.trees <- Tr.temp.trees[order(row.names(Tr.temp.trees)), ] %>% 
  mutate(Leaf_type = as.factor(Leaf_type))

TrFormula.trees = ~ Maximum_height + Leaf_area + Seed_length + Wood_density + Leaf_type


# Define model 
model1.trees <- Hmsc(Y=Y.trees, XData = X.env.trees, XFormula = XFormula.trees, 
               TrData = Tr.trees, TrFormula = TrFormula.trees, 
               studyDesign = studyDesign.trees,
               ranLevels = list("block.id" = rL.spatial.trees), 
               distr = "probit" )

# Fit model
## Run the models
nChains = 2
test.run = FALSE
if (test.run){
  #with this option, the vignette runs fast but results are not reliable
  thin = 1
  samples = 10
  transient = 5
  verbose = 1
} else {
  #with this option, the vignette evaluates slow but it reproduces the results of the
  #.pdf version
  thin = 10
  samples = 250
  transient = 250*thin
  verbose = 10
}

mod1trees_HMSC = sampleMcmc(model1.trees,
                       samples = samples,
                       thin = thin,
                       transient = transient,
                       nChains = nChains, 
                       nParallel = nChains
                       )

save(mod1trees_HMSC,file="Models/mod1trees_HMSC.Rdata")


##################### STEP 2: BIRDS HMSC USING TREES AS A COVARIATE ###############


## Data Preparation: clean vars and drop NA's
birds.df <- birds %>% 
  mutate(Habitat2 = ifelse(Olsoneconame == "Alto Parana Atlantic forests" 
                           | Olsoneconame == "Xingu-Tocantins-Araguaia moist forests", "Tropical Moist Forest", 
                           ifelse(str_detect(birds$Olsoneconame, "interior"), "Interior Forest", 
                                  ifelse(Olsoneconame == "Cerrado" | str_detect(birds$Olsoneconame, "savanna"), "Savanna", 
                                         ifelse(str_detect(birds$Olsoneconame, "coastal")
                                                | Olsoneconame == "Southern Atlantic mangroves", "Coastal Forest", 
                                                ifelse(Olsoneconame =="Araucaria moist forests", "Temperate Forest", NA
                                                       #ifelse(Olsoneconame == "Caatinga", "Shrubland", NA)
                                                )))))) %>%
  select(-c(Methods, Year_start, Year_finish, 
            Habitat, State, replicate, year, name_biome, code_biome)) %>% 
  relocate(Habitat2, .after = "ID_codref") %>%
  na.omit


## Species Occurrence: Y Select just the most common birds
birds_occ <- birds.df[, which(colnames(birds.df)=="Aburria.jacutinga"):which(colnames(birds.df) == "Zonotrichia.capensis")]
nbirds <- 10
commonbirds <- colnames(birds_occ[, order(colSums(-birds_occ))])[1:nbirds]
Y <- as.matrix(birds_occ[, which(colnames(birds_occ) %in% commonbirds)])
Y <- Y[, order(colnames(Y))] # make sure it's alphabetized


### STUDY DESIGN AND SPATIAL STRUCTURE: if sites are too close together it causes numerical instability

## CREATE BLOCKS AND LOCATIONS AS BLOCK CENTERS
max.lat.birds <- ceiling(max(birds.df$Latitude_y))
min.lat.birds <- floor(min(birds.df$Latitude_y))
max.long.birds <- ceiling(max(birds.df$Longitude_x))
min.long.birds <- floor(min(birds.df$Longitude_x))
n_groups <- 100

# Groups go from left bottom to left top, next left bottom to next left top, ...
birds.df <- within(birds.df, {
  grp.lat = cut(Latitude_y, seq(min.lat.birds, max.lat.birds, 
                                (max.lat.birds - min.lat.birds)/n_groups), labels = FALSE)
  grp.long = cut(Longitude_x, seq(min.long.birds, max.long.birds, 
                                  (max.long.birds - min.long.birds)/n_groups), labels = FALSE)}) %>%
  mutate(block.id = as.numeric(interaction(grp.lat, grp.long)))

## Compute centers and attach block ids to block centers
grid.centers <- data.frame(center_lat_y = rollsum(seq(min.lat.birds, max.lat.birds, 
                                                      (max.lat.birds - min.lat.birds)/n_groups), 2)/2, 
                           center_long_x = rollsum(seq(min.long.birds, max.long.birds, 
                                                       (max.long.birds - min.long.birds)/n_groups), 2)/2)

# Blocks move from bottom left to top left, ..., to top right
centers_long <- rep(grid.centers$center_long_x, times = rep(n_groups, times = n_groups))
centers_lat <- rep(grid.centers$center_lat_y, times = n_groups)
block_centers <- data.frame(cbind(centers_long, centers_lat, 1:(n_groups^2)))
colnames(block_centers) <- c("centers_long", "centers_lat", "block.id")

# Merge block centers with birds
birds.df <- left_join(birds.df, block_centers)


# Define Study Design: n_obs rows assigning each obs to a study design unit
studyDesign.birds = data.frame(birds.df[, "block.id"])
colnames(studyDesign.birds) <- c("block.id")
studyDesign.birds = data.frame(lapply(studyDesign.birds,as.factor))
block.id = HmscRandomLevel(units = studyDesign.birds$block.id)
ranlevels = list(block.id = block.id) 
ranlevels

# Spatial data: length(levels(stuyd design units)) rows giving location of each unit
# rownames are study design units?
xycoords.birds <- birds.df %>% 
  mutate(block.id = as.factor(block.id)) %>%
  group_by(block.id) %>%
  # mean just collapses becaus loc constant within groups
  summarize(centers_x_long = mean(centers_long), 
            centers_y_lat = mean(centers_lat)) %>%
  tibble::column_to_rownames(var = "block.id") 

rL.spatial.birds = HmscRandomLevel(sData = xycoords.birds)
rL.spatial.birds = setPriors(rL.spatial.birds,nfMin=1,nfMax=1) 


### FIXED EFFECTS

## TREES AS COVARIATES tree predictions at bird locations

# Predict trees at unique bird locations
xycoords <- birds.df[, c("centers_long", "centers_lat")]
unique.indices <- which(!duplicated(xycoords))
X.trees.new.unique <- birds.df[unique.indices, c("Annual_temperature", "Annual_rainfall")]
colnames(X.trees.new.unique) <- c("Temperature", "Precipitation")

xy.coords.new.unique <- birds.df[unique.indices, ]%>%
  select(centers_long, centers_lat)
Gradient.unique = prepareGradient(hM = mod1trees_HMSC, XDataNew = X.trees.new.unique, sDataNew = list(block.id = xy.coords.new.unique))
pred.trees <- predict(mod1trees_HMSC, Gradient = Gradient.unique, expected = TRUE)

Epred.trees <- apply(abind(pred.trees, along = 3), c(1,2), mean) # Get posterior mean predictions
Epred.trees <- cbind(Epred.trees, xy.coords.new.unique)
birds.df <- left_join(birds.df, Epred.trees) # Join by centers_long, centers_lat to expand preds to duplicated locations

# Get adjacency matrix
frug.common <- frug %>%
  mutate(Bird_species = gsub(" ", ".", Bird_species),
         Tree_species = gsub(" ", ".", Tree_species)) %>% 
  filter(Bird_species %in% commonbirds, Tree_species %in% commontrees)

# Create binary interaction network
el <- as.matrix(frug.common[, c(2,3)], ncol = 2)
g <- graph.data.frame(el)
A1 <- get.adjacency(g,sparse=FALSE)
A1 <- A1[rownames(A1) %in% commonbirds, colnames(A1) %in% commontrees]

A.full <- matrix(0, nrow = nbirds, ncol = ntrees, dimnames = list(commonbirds, commontrees))

present_bird_indices <- sapply(dimnames(A1)[[1]], function(x) which(x==commonbirds))
present_tree_indices <- sapply(dimnames(A1)[[2]], function(x) which(x==commontrees))
A.full[present_bird_indices, present_tree_indices] <- A1
A.full <- ifelse(A.full > 0, 1, 0) 
A.full <- A.full[order(row.names(A.full)), order(colnames(A.full))] # alphabetize, make sure correct
  

# Use adjacency matrix for variable selection: simple version for now
# Eventually use posterior interaction probablities from Georgia's model
# If no observed interaction, set prior prob for tree covar to 0.5, otherwise include tree covar with prob 1
A.prior <- ifelse(A.full == 0, 0.4, 1)
Xsel <- list()
for(j in 1:ntrees){
 Xsel[[j]] <- list(covGroup = c(j), spGroup = 1:nbirds, q = A.prior[, j] ) 
}


# Environmental features: set up for models with and without trees
X.env.trees <- birds.df[, c(which(colnames(birds.df) == "Alchornea.triplinervia"):
                    which(colnames(birds.df) == "Zanthoxylum.rhoifolium"))] %>%
         cbind(.,birds.df[, c("Habitat2","Annual_temperature", "Annual_rainfall", 
                              "Altitude")]) %>%
         mutate(Habitat2 = as.factor(Habitat2)) 

XFormula.trees = as.formula(paste(" ~ poly(Annual_temperature,degree = 2, raw=TRUE) +
                                  poly(Annual_rainfall, degree = 2, raw=TRUE) + 
                                  Altitude + Habitat2 + ", 
                                  paste(colnames(A.prior), collapse=' + ')))

X.env <- birds.df[, c("Habitat2","Annual_temperature", "Annual_rainfall", "Altitude")] %>%
  mutate(Habitat2 = as.factor(Habitat2)) 

XFormula = ~ poly(Annual_temperature,degree = 2, raw=TRUE) + 
             poly(Annual_rainfall, degree = 2, raw=TRUE) + Altitude + Habitat2 



# Traits
Tr.temp <- taxa %>%
      mutate(Species = gsub(" ", ".", taxa$Species)) %>%
      set_rownames(.$Species) %>% 
      filter(Species %in% commonbirds)  %>%
      select(family.common, Order) 
Tr <- Tr.temp[order(row.names(Tr.temp)), ] %>% 
  mutate(Order = as.factor(Order), family.common = as.factor(family.common))

TrFormula = ~ Order



## Define models

# 1 - uses trees and interaction information
model1 <- Hmsc(Y=Y, XData = X.env.trees, XFormula = XFormula.trees, XSelect = Xsel,
               TrData = NULL, TrFormula = NULL, 
               studyDesign = studyDesign.birds,
               ranLevels = list("block.id" = rL.spatial.birds), distr = "probit" )

# 2 - trees but no selection
model2 <- Hmsc(Y=Y, XData = X.env.trees, XFormula = XFormula.trees,
               TrData = NULL, TrFormula = NULL, 
               studyDesign = studyDesign.birds,
               ranLevels = list("block.id" = rL.spatial.birds), distr = "probit" )

# 3 - no trees
model3 <- Hmsc(Y=Y, XData = X.env, XFormula = XFormula,
               TrData = NULL, TrFormula = NULL, 
               studyDesign = studyDesign.birds,
               ranLevels = list("block.id" = rL.spatial.birds), distr = "probit" )

## Run the models
nChains = 2
test.run = FALSE
if (test.run){
  #with this option, the vignette runs fast but results are not reliable
  thin = 1
  samples = 10
  transient = 5
  verbose = 1
} else {
  #with this option, the vignette evaluates slow but it reproduces the results of the
  #.pdf version
  thin = 10
  samples = 250
  transient = 250*thin
  verbose = 10
}

# Model 1: trees and interactions
mod1_HMSC = sampleMcmc(model1,
                      samples = samples,
                      thin = thin,
                      transient = transient,
                      nChains = nChains, 
                      nParallel = nChains)
save(mod1_HMSC,file="Models/mod1_HMSC.Rdata")

# MCMC diagnostics
m1post <- convertToCodaObject(mod1_HMSC)
hist(effectiveSize(m1post$Beta), main="ess(beta)")
hist(gelman.diag(m1post$Beta, multivariate=FALSE)$psrf, main="psrf(beta)")

# Look at params
postBeta = getPostEstimate(mod1_HMSC, parName = "Beta")
plotBeta(mod1_HMSC, post = postBeta, param = "Support", supportLevel = 0.95)

# Look at fit
preds.spatial.1 = computePredictedValues(mod1_HMSC)
MF.spatial.1= evaluateModelFit(hM=mod1_HMSC, predY=preds.spatial.1)
MF.spatial.1


mod2_HMSC = sampleMcmc(model2,
                       samples = samples,
                       thin = thin,
                       transient = transient,
                       nChains = nChains, 
                       nParallel = nChains)
save(mod2_HMSC,file="Models/mod2_HMSC.Rdata")


mod3_HMSC = sampleMcmc(model3,
                       samples = samples,
                       thin = thin,
                       transient = transient,
                       nChains = nChains, 
                       nParallel = nChains)
save(mod3_HMSC,file="Models/mod3_HMSC.Rdata")


  

############################ DEBUGGING CODE ####################################
# Debugging: simulate site and plot data and see if that works better: this works
# Lesson learned: aim to have no fewer than 15 observations at the finest level

# Variable selection 
studyDesign <- data.frame("site" = rep(1:5, times = c(100,100,100,100, nrow(birds.df) - 400))) %>%
  group_by(site) %>%
  mutate(plot = rep(1:4, times = c(15,15,15,n() - 45)))
studyDesign = data.frame(lapply(studyDesign,as.factor))

# Variable structuring 
site = HmscRandomLevel(units = studyDesign$site)
plot = HmscRandomLevel(units = studyDesign$plot)
ranlevels = list(site = site) 
ranlevels

# Get block counts
empty_blocks <- data.frame(block.id = setdiff(1:(n_groups^2), unique(trees.df$block.id)), block.count = 0)
block_counts <- data.frame(table(trees.df$block.id))
colnames(block_counts) <- c("block.id", "block.count")
block_counts$block.id <- as.numeric(as.character(block_counts$block.id))
block_counts <- rbind(block_counts, empty_blocks) %>% 
  arrange(block.id)


########################### OLD CODE ###########################################

# Previously had to keep only unique obs
xycoords <- birds[, c("Longitude_x", "Latitude_y")]
birds.df <- birds[!duplicated(xycoords), ]

mutate(Olsoneconame = ifelse(Olsoneconame %in% c("Pernambuco interior forests", "Uruguayan savanna", 
                                                 "Caatinga", "Southern Atlantic mangroves", "Campos Rupestres montane savanna", 
                                                 "Pernambuco interior forests", "Unknown")
                             , NA, Olsoneconame))


# Attempt 1: Trees at all bird locations (including repeats)
X.trees.new <- birds.df[, c("Annual_temperature", "Annual_rainfall")]
colnames(X.trees.new) <- c("Temperature", "Precipitation")

xy.coords.new <- birds.df%>%
  select(centers_long, centers_lat) %>%
  rename(centers_x_long = centers_long, centers_y_lat = centers_lat)
Gradient = prepareGradient(hM = mod1trees_HMSC, XDataNew = X.trees.new, sDataNew = list(block.id = xy.coords.new))
pred.trees <- predict(mod1trees_HMSC, Gradient = Gradient, expected = TRUE)
