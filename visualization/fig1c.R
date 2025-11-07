setwd("E:\\Genome\\Comparative_genomics\\Species_distribution") # set working directory

library(dismo)
library(raster)
library(geodata)

#bioclim_data <- worldclim_global(var = "bio",res = 10, path = "data/")  #download
                                
bioclim_data <-list.files(path = "D:\\Population_643\\all_761_info\\local_adaption\\wc2.1_2.5m_bio\\", pattern = "*.tif$", full.names=TRUE)
bioclim_data <- stack(bioclim_data)
predNames <- c('wc2.1_2.5m_bio_19','wc2.1_2.5m_bio_8','wc2.1_2.5m_bio_9','wc2.1_2.5m_bio_15','wc2.1_2.5m_bio_2','wc2.1_2.5m_bio_17','wc2.1_2.5m_bio_18') ##only use some bioclim data
bioclim_data=bioclim_data[[predNames]]
##test model use 397 data sample,397sample_info.txt
###real model use data download from GBIF.org,

# Determine geographic extent of our data (whole world)
max_lat = 90
min_lat = -60
max_lon = 180
min_lon = -170
#max_lat <- ceiling(max(obs_data$Latitude))
#min_lat <- floor(min(obs_data$Latitude))
#max_lon <- ceiling(max(obs_data$Longitude))
#min_lon <- floor(min(obs_data$Longitude))

geographic_extent <- extent(x = c(min_lon, max_lon, min_lat, max_lat))
# Crop the bioclim data to geographic extent of my data
bioclim_data <- raster::crop(x = bioclim_data, y = geographic_extent)

tif_files <- list.files(path = "D:\\Population_643\\all_761_info\\local_adaption\\wc2.1_2.5m_bio\\", 
                        pattern = "*.tif$", 
                        full.names = TRUE)

# We only need one file, so use the first one in the list of .tif files
mask <- raster(tif_files[1])
set.seed(1)

obs_data <- read.table(file = "use_m_caerulea.txt",header = T,sep = "\t")  ###alfalfa:use_medicago_sativa_USDA_ARS.txt
obs_data=obs_data[,c("Longitude","Latitude")] #only use latitude and longitude


# Randomly sample points, not overlap with our present data
background <- randomPoints(mask = mask,     # Provides resolution of sampling points
                           n = nrow(obs_data),      # Number of random points
                           ext = geographic_extent, # Spatially restricts sampling
                           extf = 1.25)             # Expands sampling a little bit

# Arbitrarily assign group 1 as the testing data group
testing_group <- 1

# Create vector of group memberships,test using random point, real data need use representative data
group_presence <- kfold(x = obs_data, k = 5) # kfold is in dismo package

# Separate observations into training and testing groups
presence_train <- obs_data[group_presence != testing_group, ]
presence_test <- obs_data[group_presence == testing_group, ]


# Repeat the process for pseudo-absence points
group_background <- kfold(x = background, k = 5)

background_train <- background[group_background != testing_group, ]
background_test <- background[group_background == testing_group, ]
colnames(background_train)=c("Longitude","Latitude")
train <- rbind(presence_train, background_train)
pb_train <- c(rep(1, nrow(presence_train)), rep(0, nrow(background_train)))
envtrain <- raster::extract(bioclim_data, train[,c("Longitude","Latitude")])
envtrain <- data.frame( cbind(pa=pb_train, envtrain) )
envtrain=na.omit(envtrain) ##remove na
testpres <- data.frame( extract(bioclim_data, presence_test) )
testbackg <- data.frame( extract(bioclim_data, background_test) )
testpres=na.omit(testpres)
testbackg=na.omit(testbackg)
##random forest
library(randomForest)
#rf1 <- randomForest(pa~wc2.1_2.5m_bio_1+wc2.1_2.5m_bio_2+wc2.1_2.5m_bio_3+wc2.1_2.5m_bio_4+wc2.1_2.5m_bio_5+wc2.1_2.5m_bio_6+wc2.1_2.5m_bio_7+wc2.1_2.5m_bio_8+wc2.1_2.5m_bio_9+wc2.1_2.5m_bio_10+wc2.1_2.5m_bio_11+wc2.1_2.5m_bio_12+wc2.1_2.5m_bio_13+wc2.1_2.5m_bio_14+wc2.1_2.5m_bio_15+wc2.1_2.5m_bio_16+wc2.1_2.5m_bio_17+wc2.1_2.5m_bio_18+wc2.1_2.5m_bio_19, data=envtrain)
rf1 <- randomForest(pa~wc2.1_2.5m_bio_19+wc2.1_2.5m_bio_8+wc2.1_2.5m_bio_9+wc2.1_2.5m_bio_15+wc2.1_2.5m_bio_2+wc2.1_2.5m_bio_17+wc2.1_2.5m_bio_18, data=envtrain)
#rf1 <- randomForest(pa~wc2.1_10m_bio_19+wc2.1_10m_bio_8+wc2.1_10m_bio_9+wc2.1_10m_bio_15+wc2.1_10m_bio_2+wc2.1_10m_bio_17+wc2.1_10m_bio_18, data=envtrain)

erf <- evaluate(testpres, testbackg, rf1)
erf
pr <- dismo::predict(bioclim_data, rf1, ext=geographic_extent) ##need some time
tr <- threshold(erf, 'spec_sens')
pr_caerulea=pr
tr_caerulea=tr
#saveRDS(pr_caerulea,file="pr_caerulea_info.rds") ##tr_caerulea=0.5036333
#saveRDS(pr_sativa,file="pr_sativa_info.rds") ##tr_sativa=0.4726333

##combine sativa and caerulea
pdf(file = "RF_model_present_distribution_combine_sativa_caerulea.pdf",   width = 6, height = 4)
plot(pr_sativa > (tr_sativa*0.8), main='',xlab="Longitude",ylab="Latitude",legend = FALSE,col=c("gray95","lightsteelblue1"),ylim=c(-60,90),xlim=c(-180,180)) ##这个分布阈值可以根据实际分布区域调整
plot(pr_caerulea > (tr_caerulea*1.7), main='',add = TRUE,legend = FALSE,col=c("#00000000","#00868B"),ylim=c(-60,90),xlim=c(-180,180)) ##这个分布阈值可以根据实际分布区域调整

dev.off()
