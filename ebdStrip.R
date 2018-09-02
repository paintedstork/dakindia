##### Required to do once for an eBird set. After that use ebd.rds

#######################Configurations###########################
#Name of the eBird quarterly archive WITHOUT .zip
ebd_file_name <- 'ebd_IN_relMay-2018'

################################################################

#Unzip and read eBird records. This assumes the file is stored in a top level folder called 'data'
unzip(paste('..\\data\\',ebd_file_name,'.zip',sep=''))

# Strip eBird file using an awk script as the file cant load fully in R due to excessive memory requirements
# The field that interest us are indicated by the column numbers
awkscript <- '\"BEGIN {FS=\\"\\t\\"};{print $4,\\"\\t\\",$5,\\"\\t\\",$6,\\"\\t\\",$7,\\"\\t\\",$26,\\"\\t\\",$27,\\"\\t\\",$28,\\"\\t\\",$31,\\"\\t\\",$40,\\"\\t\\",$42,\\"\\t\\",$43;}\"'

system2 ("awk", 
         # args = c('-f', 'ebdStrip.awk', paste(ebd_file_name,".txt", sep='')),
         args = c(awkscript, paste(ebd_file_name,".txt", sep='')),
         wait = TRUE, # Wait till awk is completed
         stdout = paste(ebd_file_name,"_strip.txt", sep=''))


# Read the stripped file
ebd <- read.delim(paste(ebd_file_name,'_strip.txt',sep=''), na.strings = c("NA", "", "null"), as.is=TRUE, quote="")

#Add unique list identifier for removing duplicates
ebd <- within (ebd, UNIQUE_SAMPLING_ID <-  ifelse(is.na(GROUP.IDENTIFIER),SAMPLING.EVENT.IDENTIFIER,GROUP.IDENTIFIER))

#If subspecies, copy subspecies common name to roll up issfs to species level (Comment this line if you want issf explicitly)
ebd <- within (ebd, COMMON.NAME <-  ifelse(CATEGORY=='issf',SUBSPECIES.COMMON.NAME,COMMON.NAME))

#Remove entries from shared lists (Comment this line if you dont want all lists)
ebd_records   <- ebd[!duplicated(ebd[c("UNIQUE_SAMPLING_ID","COMMON.NAME")]),]

# Reduce dataset to only columns of interest
ebd <- subset(ebd, select = c("LATITUDE", 
                              "LONGITUDE", 
                              "COMMON.NAME", 
                              "CATEGORY", 
                              "OBSERVATION.DATE", 
                              "APPROVED"))

# Remove UnApproved records - which are of exotic species 
ebd <- ebd [ebd$APPROVED==1,]

ebd <- subset(ebd, select = c("LATITUDE", 
                              "LONGITUDE", 
                              "COMMON.NAME", 
                              "CATEGORY", 
                              "OBSERVATION.DATE"
))

saveRDS(ebd,"ebd.rds")
