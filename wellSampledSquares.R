library (plyr)
library (dplyr)
library (rgdal)
library (sp)
library (compare)
library (rgeos)
library (maptools)
library (udunits2)
library (geojsonio)
library (ggplot2)
library (fossil)
library (reshape2)
# It has too many libraries. I have not removed ones that are not needed

#-------------------------------------------------------------------
#  LOAD EBIRD DATA AND PROCESS IT TO ADD GRIDS
#  DO NOT DO IT EVERYTIME AS ITS TIME CONSUMING
#-------------------------------------------------------------------

##### IMPORTANT: THIS NEEDS AN AWK SCRIPT TO BE RUN BEFORE HAND
##### LOOK AT ebdStrip.R

# Load eBird Data
g_ebd <- readRDS ("ebd.rds")
# Keep g_ebd intact, if some further processing is needed
ebd   <- g_ebd
# Configurable square resolution. 0.5 degree to start with
squareResolution <- 0.5

# Create a block column in the database
ebd <- within (ebd, 
                         BLOCK <-  paste( squareResolution * as.integer(LONGITUDE/squareResolution + 0.5), 
                                          squareResolution * as.integer(LATITUDE/squareResolution  + 0.5),
                                          sep=':'))
ebd$LATITUDE <- NULL
ebd$LONGITUDE<- NULL

# Create an era column in the database. Configuring this function provides multiple ERAs
ebd <- within (ebd, 
               ERA <- ifelse(as.numeric(format(as.Date(ebd$OBSERVATION.DATE),"%Y"))< 2010,
                             ifelse(as.numeric(format(as.Date(ebd$OBSERVATION.DATE),"%Y")) < 2000, 
                             ifelse(as.numeric(format(as.Date(ebd$OBSERVATION.DATE),"%Y")) < 1980, 1980,2000),
                             2010),
                             2018))

# Use only species, issf (subspecies group), domestic (Rock Pigeon)
ebd <- rbind (ebd [ebd$CATEGORY=="species ",],
              ebd [ebd$CATEGORY=="issf ",],
              ebd [ebd$CATEGORY=="domestic ",])

ebd$CATEGORY <- NULL

saveRDS(ebd, "min_ebd.rds")
---------------------------------------------------------------
  
  
  
#--------------------------------------------------------------  
# This function calculates the Chao2 value for a specified block/grid
chaoPerBlock <- function (Block, DataSet)
{
  print (Block)
  DataSet <- DataSet[DataSet$BLOCK == Block,]
  DataCast <- dcast (DataSet, 
                     OBSERVATION.DATE ~ COMMON.NAME,
                     value.var = "BLOCK",
                     fun.aggregate = function(u){length(u[u>0])} )
  return (chao2(DataCast, taxa.row=TRUE))  
}
#---------------------------------------------------------------


#---------------------------------------------------------------
# A VECTORED CHAO FUNCTION THAT CALCULATES PER ERA
#---------------------------------------------------------------

chaoPerEra <- function (Era, DataSet)
{
  print(paste("Processing Era:",Era))
  ebd <- DataSet[DataSet$ERA == Era,]
  
  # Create a species vs grid matrix - with values showing the number of days when it was found
  species_per_grid <- dcast(ebd, 
                            BLOCK ~ COMMON.NAME, 
                            value.var = "OBSERVATION.DATE",
                            fun.aggregate = length)

  # Number of species per grid
  n_species_per_grid <- as.data.frame(rowSums(species_per_grid!=0))
  colnames(n_species_per_grid) <- c("NO.OF.SPECIES")

  # Number of records per grid
  block_count <- dcast (ebd, BLOCK~BLOCK) #Diagonal will have the record counts

  # Extract the diagonal as numeric data frame
  n_records_per_grid <- as.data.frame (as.numeric ( diag(as.matrix(block_count)[,2:ncol(block_count)])))
  rownames(n_records_per_grid) <- rownames(n_species_per_grid)
  colnames(n_records_per_grid) <- c("NO.OF.RECORDS")

  # Chao2 per grid - this is a costly function. Can we optimize it?
  n_species_chao2 <- mapply (chaoPerBlock, 
                             Block = species_per_grid$BLOCK, 
                             MoreArgs = list(DataSet = ebd)) 

  # Save the data for intermediate processing
  saveRDS(n_species_per_grid, paste(Era,"n_species.rds"))
  saveRDS(n_records_per_grid, paste(Era,"n_records.rds"))
  saveRDS(n_species_chao2, paste(Era,"chao.rds"))
  
  n_species_per_grid  = readRDS(paste(Era,"n_species.rds"))
  n_records_per_grid  = readRDS(paste(Era,"n_records.rds"))
  n_species_chao2     = readRDS(paste(Era,"chao.rds"))

  # C value per grid
  c_value = n_species_chao2 / n_species_per_grid 
  colnames(c_value) <- c("C.VALUE")
  saveRDS(c_value, paste(Era,"c_value.rds"))

  well_sampled_cell = (c_value > 0.9) & (n_records_per_grid > 200)

  colnames(well_sampled_cell) = c("WELL.SAMPLED")

  chao <- as.data.frame (readRDS(paste(Era,"chao.rds")))

  chao <- cbind(chao, rownames(chao))

  colnames(chao) <- c("chao", "block")

  chao <- cbind(chao, as.numeric(sapply(strsplit(as.character(chao$block),":"), `[`, 1)))
  chao <- cbind(chao, as.numeric(sapply(strsplit(as.character(chao$block),":"), `[`, 2)))
  chao$block <- NULL

  chao <- chao[,c(2, 3, 1)]

  colnames(chao) <- c("LONGITUDE", "LATITUDE", "chao")

  chao <- cbind(chao, well_sampled_cell, c_value, n_species_per_grid, n_records_per_grid)

  nchao <- chao[chao$WELL.SAMPLED==1,]

  saveRDS(nchao, paste(Era,"nchao.rds"))
  return (chao)
}

#---------------------------------------------------------------
# Load the India map on to which we need to project. Need only for ggplot2
setwd('data')
india_map<-   readRDS('india_map.rds')
setwd('..')
#-------------------------------------------------------------

# Actual execution starts here
ebd <- readRDS("min_ebd.rds")

# Calculate the nchao_list for each ERA. You need to updated the Era list for new years
nchao_list <- mapply (chaoPerEra, 
                           Era = c(1980,2000, 2010, 2018), 
                           MoreArgs = list(DataSet = ebd)) 

# Next ggplot code is dirty - should be cleaned up, formatted well and vectorised
ggplot (aes(x=LONGITUDE,y=LATITUDE,fill=NO.OF.RECORDS),
        data=as.data.frame(nchao_list[,1])) + 
  scale_color_gradient(low="blue", high="red")  +
  geom_tile() +
  coord_fixed (ratio = 1.0) +
  geom_polygon(data=india_map,aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0)
ggsave("No_Records_1800_1980.png", width = 500, height = 500, units = "mm", dpi = 90)

ggplot (aes(x=LONGITUDE,y=LATITUDE,fill=NO.OF.RECORDS),
        data=as.data.frame(nchao_list[,2])) + 
  scale_color_gradient(low="blue", high="red")  +
  geom_tile() +
  coord_fixed (ratio = 1.0) +
  geom_polygon(data=india_map,aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0)
ggsave("No_Records_1980_2000.png", width = 500, height = 500, units = "mm", dpi = 90)

ggplot (aes(x=LONGITUDE,y=LATITUDE,fill=NO.OF.RECORDS),
        data=as.data.frame(nchao_list[,3])) + 
  scale_color_gradient(low="blue", high="red")  +
  geom_tile() +
  coord_fixed (ratio = 1.0) +
  geom_polygon(data=india_map,aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0)
ggsave("No_Records_2000_2010.png", width = 500, height = 500, units = "mm", dpi = 90)

ggplot (aes(x=LONGITUDE,y=LATITUDE,fill=NO.OF.RECORDS),
        data=as.data.frame(nchao_list[,4])) + 
  scale_color_gradient(low="blue", high="red")  +
  geom_tile() +
  coord_fixed (ratio = 1.0) +
  geom_polygon(data=india_map,aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0)
ggsave("No_Records_2010_2018.png", width = 500, height = 500, units = "mm", dpi = 90)

  
ggplot (aes(x=LONGITUDE,y=LATITUDE,fill=WELL.SAMPLED),
        data=as.data.frame(nchao_list[,1])) + 
  scale_color_gradient(low="blue", high="red")  +
  geom_tile() +
  coord_fixed (ratio = 1.0) +
  geom_polygon(data=india_map,aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0)
ggsave("Well_Sampled_1800_1980.png", width = 500, height = 500, units = "mm", dpi = 90)

ggplot (aes(x=LONGITUDE,y=LATITUDE,fill=WELL.SAMPLED),
        data=as.data.frame(nchao_list[,2])) + 
  scale_color_gradient(low="blue", high="red")  +
  geom_tile() +
  coord_fixed (ratio = 1.0) +
  geom_polygon(data=india_map,aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0)
ggsave("Well_Sampled_1980_2000.png", width = 500, height = 500, units = "mm", dpi = 90)

ggplot (aes(x=LONGITUDE,y=LATITUDE,fill=WELL.SAMPLED),
        data=as.data.frame(nchao_list[,3])) + 
  scale_color_gradient(low="blue", high="red")  +
  geom_tile() +
  coord_fixed (ratio = 1.0) +
  geom_polygon(data=india_map,aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0)
ggsave("Well_Sampled_2000_2010.png", width = 500, height = 500, units = "mm", dpi = 90)

ggplot (aes(x=LONGITUDE,y=LATITUDE,fill=WELL.SAMPLED),
        data=as.data.frame(nchao_list[,4])) + 
  scale_color_gradient(low="blue", high="red")  +
  geom_tile() +
  coord_fixed (ratio = 1.0) +
  geom_polygon(data=india_map,aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0)
ggsave("Well_Sampled_2010_2018.png", width = 500, height = 500, units = "mm", dpi = 90)


