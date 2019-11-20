#Package prep
packages<-c("cowplot", "dplyr", "geosphere", "ggplot2", "ggExtra", "maps", "maptools", "readxl", "rgdal", "rgeos", "sf", "sp", "spatialreg", "spdep", "tidyr", "viridis")
sapply(packages, require, character.only=T)

#Designating format for every column
data <- read.csv('childpov18_southfull.csv', 
                 colClasses = c("character", "character", "character", 
                                "numeric", "numeric", "numeric", "numeric",
                                "numeric", "numeric", "numeric", "numeric",
                                "numeric", "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", "numeric",
                                "numeric", "numeric", "numeric", "numeric",
                                "numeric", "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", "numeric",
                                "numeric", "numeric", "numeric", "numeric"))
head(data)
names(data)[names(data)=="X2016.child.poverty"] <- "child.pov.2016"

#Now we pick a state, in this case I have chosen Colorado 
ga_pov <- data %>% subset(State == "GA")
summary(ga_pov)

#Ordinary least squares
equation <- child.pov.2016 ~ rural + urban + lnmanufacturing + lnag + 
  lnretail + lnhealthss + lnconstruction + lnlesshs + 
  lnunemployment + lnsinglemom + lnblack + lnhispanic + 
  lnuninsured + lnincome_ratio + lnteenbirth + lnunmarried

options(scipen = 5)

ols <- lm(equation, data=ga_pov)
summary(ols)

#Now to create a list of contiguity neighbors
#Obtain FIPS Codes by county 
fips <- county.fips

#Create county polygons
georgia <- map(database = "county", regions = "georgia", fill=T, plot=F)
IDs <- sub("^georgia,","",georgia$names)

#Add FIPS codes to the county polygons
fips.codes <- separate(data = fips, col = polyname, into = c("state", "county"), sep = ",")
ga_fips <- subset(fips.codes, state=="georgia", select=fips)
names <- fips.codes$county
ga_IDs <- unique(ga_fips$fips)

#Create spatial polygons
ga_sp = map2SpatialPolygons(georgia,ga_fips$fips,CRS("+proj=longlat"))
names(ga_sp@polygons) <- ga_IDs

#Create neighbor weights using the queens case
neighb.data <- poly2nb(ga_sp, queen=T)
names(neighb.data) <- names(ga_sp@polygons)

#Create list of neighbors
cont.neighb <- nb2listw(neighb.data,style="W", zero.policy = TRUE)

lm.morantest(ols, cont.neighb)

lm.LMtests(ols, cont.neighb, test="all")

#Spatially lagged x model
SLX.model <- lmSLX(equation, data=ga_pov, cont.neighb)
summary(SLX.model)

summary(impacts(SLX.model, cont.neighb), zstats = TRUE)[["pzmat"]]

#Spatial lag model

sp.lag.model <- spatialreg::lagsarlm(equation, data=ga_pov, cont.neighb)
summary(sp.lag.model, Nagelkerke = TRUE)

summary(impacts(sp.lag.model, listw = cont.neighb, R=100), zstats = TRUE)[["pzmat"]]

#Spatial error model
sp.err.model <- spatialreg::errorsarlm(equation, data=ga_pov, cont.neighb)
summary(sp.err.model, Nagelkerke = TRUE)

#We have selected the Spatial lag model based on the fact it had the best pvalue

spatialreg::Hausman.test(sp.lag.model)

sd.lag <- spatialreg::lagsarlm(equation, ga_pov, cont.neighb)
sdm <- spatialreg::lagsarlm(equation, ga_pov, cont.neighb, type = "mixed")

summary(sd.lag, Nagelkerke = TRUE)

summary(spatialreg::impacts(sd.lag, listw = cont.neighb, R = 100), zstats = TRUE)[["pzmat"]]
spatialreg::LR.sarlm(sd.lag,sp.lag.model)

#Now onto spatial regression looking at k-nearest neighbors
all.xy <- centroid(ga_sp)
#tx_IDs <- unique(tx_fips$fips) this value was created in the contiguity section but would be needed here if only using distance functions. See "creating list of contiguity neighbors" for details.
rownames(all.xy) <- ga_IDs
colnames(all.xy) <- cbind("x","y")

#Create neighbors
all.dist.k1 <- knn2nb(knearneigh(all.xy, k=1, longlat = TRUE))
all.dist.k3 <- knn2nb(knearneigh(all.xy, k=3, longlat = TRUE))
all.dist.k5 <- knn2nb(knearneigh(all.xy, k=5, longlat = TRUE))

#Determine max k distance value to neighbor

all.max.k1 <- max(unlist(nbdists(all.dist.k1, all.xy, longlat=TRUE)))
all.max.k3 <- max(unlist(nbdists(all.dist.k3, all.xy, longlat=TRUE)))
all.max.k5 <- max(unlist(nbdists(all.dist.k5, all.xy, longlat=TRUE)))

#Calculate neighbors based on distance
all.sp.dist.k1 <- dnearneigh(all.xy, d1=0, d2=1 * all.max.k1, longlat = TRUE)
all.sp.dist.k3 <- dnearneigh(all.xy, d1=0, d2=1 * all.max.k3, longlat = TRUE)
all.sp.dist.k5 <- dnearneigh(all.xy, d1=0, d2=1 * all.max.k5, longlat = TRUE)

#Create neighbor list
all.dist.neighb.k1 <- nb2listw(all.sp.dist.k1,style="W", zero.policy = TRUE)
all.dist.neighb.k3 <- nb2listw(all.sp.dist.k3,style="W", zero.policy = TRUE)
all.dist.neighb.k5 <- nb2listw(all.sp.dist.k5,style="W", zero.policy = TRUE)

all.dist.lag.k1 <- lagsarlm(equation, data = ga_pov, listw = all.dist.neighb.k1)
all.dist.lag.k3 <- lagsarlm(equation, data = ga_pov, listw = all.dist.neighb.k3)
all.dist.lag.k5 <- lagsarlm(equation, data = ga_pov, listw = all.dist.neighb.k5)

summary(all.dist.lag.k1, Nagelkerke = TRUE)
#Distance error modeling
all.dist.lag.k1 <- errorsarlm(equation, data = ga_pov, listw = all.dist.neighb.k1)
all.dist.lag.k3 <- errorsarlm(equation, data = ga_pov, listw = all.dist.neighb.k3)
all.dist.lag.k5 <- errorsarlm(equation, data = ga_pov, listw = all.dist.neighb.k5)

summary(all.dist.lag.k1, Nagelkerke = TRUE)

dist.lag.data <- summary(all.dist.lag.k1, correlation=TRUE, Nagelkerke = TRUE)

dist.lag.output <- cbind.data.frame(ga_pov$FIPS,
                                    dist.lag.data$fitted.values, 
                                    dist.lag.data$residual, 
                                    ga_pov$child.pov.2016, 
                                    ga_pov$lnsinglemom, 
                                    ga_pov$lnunmarried, 
                                    ga_pov$lnlesshs, 
                                    ga_pov$lnincome_ratio,
                                    stringsAsFactors = FALSE)

#Renaming columns
colnames(dist.lag.output) <- c("fips","fitted","resid","childpov",
                               "single_mom","unmarried","less_hs","income_ratio")

#Create quantiles
quantiles_sm <- dist.lag.output %>%
  pull(single_mom) %>%
  quantile(probs = seq(0, 1, length.out = 4), na.rm = TRUE)

quantiles_lh <- dist.lag.output %>%
  pull(less_hs) %>%
  quantile(probs = seq(0, 1, length.out = 4), na.rm = TRUE)

#Create ranks
sm_rank <- cut(dist.lag.output$single_mom, 
               breaks= quantiles_sm, 
               labels=c("1", "2", "3"), 
               na.rm = TRUE, 
               include.lowest = TRUE)

lh_rank <- cut(dist.lag.output$less_hs, 
               breaks= quantiles_pov, 
               labels=c("1", "2", "3"), 
               na.rm = TRUE,
               include.lowest = TRUE)

#Join ranks and combined column to dataset
dist.lag.output$sm_score <- as.numeric(sm_rank)
dist.lag.output$lh_score <- as.numeric(lh_rank)
dist.lag.output$sm_lh <- paste(as.numeric(dist.lag.output$lh_score), 
                               "-", 
                               as.numeric(dist.lag.output$sm_score))
#creating the legend for 
legend_colors <- tibble(
  x = c(3,2,1,3,2,1,3,2,1),
  y = c(3,3,3,2,2,2,1,1,1),
  z = c("#574249", "#627f8c", "#64acbe", "#985356", "#ad9ea5", "#b0d5df", "#c85a5a", "#e4acac", "#e8e8e8"))

xlabel <- "Less than HS edu,Low \u2192 High"
xlabel <- gsub(",", "\n", xlabel)
ylabel <- "Single Mother Household,Low \u2192 High"
ylabel <- gsub(",", "\n", ylabel)

legend <- ggplot(legend_colors, aes(x,y)) + 
  geom_tile(aes(fill=z)) + 
  theme_minimal() + theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  labs(x = xlabel, y = ylabel) + 
  scale_fill_identity() +
  ggtitle("Legend") +
  theme(axis.title.y = element_text(face = "italic", hjust = 0.5, size = 8)) +
  theme(axis.title.x = element_text(face = "italic", hjust = 0.5, size = 8)) +
  theme(plot.title = element_text(face="bold", hjust = 0.5, size = 10))

#We need the FIPS codes
world <- map_data("world")
states <- map_data("state")
counties <- map_data("county")

counties$polyname <- paste(counties$region, counties$subregion, sep = ",")
counties <- counties %>% left_join(fips, by = c("polyname" = "polyname"))
counties$fips <- as.character(counties$fips)
counties <- counties %>% left_join(ga_pov, by = c("fips" = "FIPS"))

southern_states <- subset(states, region %in% 
                            c("texas", "arkansas", "louisiana", "mississippi", 
                              "alabama", "georgia", "florida", "north carolina",
                              "south carolina", "tennessee", "oklahoma", 
                              "kentucky", "west virginia", "virginia", 
                              "maryland", "delaware", "district of columbia"))

southern_counties <- subset(counties, region %in% 
                              c("texas", "arkansas", "louisiana", "mississippi", 
                                "alabama", "georgia", "florida", "north carolina",
                                "south carolina", "tennessee", "oklahoma", 
                                "kentucky", "west virginia", "virginia", 
                                "maryland", "delaware", "district of columbia"))
ga_counties <- subset(southern_counties, region == "georgia")

#Attach the data via the FIPS column and fortify the polygon
ga_poly <- ga_counties %>% 
  left_join(dist.lag.output, by = c("fips" = "fips")) %>%
  fortify

#Add custom color scheme based on ranks
bivariate_color_scale <- tibble(
  "3 - 3" = "#574249", 
  "2 - 3" = "#627f8c",
  "1 - 3" = "#64acbe",
  "3 - 2" = "#985356",
  "2 - 2" = "#ad9ea5",
  "1 - 2" = "#b0d5df",
  "3 - 1" = "#c85a5a",
  "2 - 1" = "#e4acac",
  "1 - 1" = "#e8e8e8") %>%
  gather("group", "fill")

ga_poly <- ga_poly %>% 
  left_join(bivariate_color_scale, by = c("sm_lh" = "group"))

sm_lh_map <- ggplot() + 
  geom_polygon(data = world, aes(x=long,y=lat, group=group), fill = "gray95", color = "white") +
  geom_polygon(data = states, aes(x=long,y=lat, group=group), fill = "gray", color = "white") +
  geom_polygon(data = ga_poly, aes(x=long, y=lat, group=group, fill = fill)) + 
  geom_polygon(data = southern_states, aes(x=long,y=lat, group=group), fill = NA, color = "white") +
  geom_polygon(data = ga_counties, aes(x=long,y=lat, group=group), fill = NA, color = "black", size = 0.05) +
  coord_map("conic", lat0 = 30, xlim=c(-90,-75), ylim=c(30,36)) +
  scale_fill_identity() +
  theme_grey() + theme(legend.position="bottom") + theme(legend.title.align=0.5) +
  theme(panel.background = element_rect(fill = 'deepskyblue'),
        panel.grid.major = element_line(colour = NA)) +
  labs(x = "Longitude", y = "Latitude", fill = "Less than HS degree", 
       title = "Bivariate Map of no HS degree and Single Mother Households") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

#mom_pov_map use to preview the map
final_map <- ggdraw() +
  draw_plot(sm_lh_map, x = 0, y = 0, width = 1, height = 1) +
  draw_plot(legend, x = 1.0, y = 0.25, width = 0.2, height = 0.25) 

final_map
