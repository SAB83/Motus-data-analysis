# Motus-data-analysis
#Using motus package other related packages and other analyses to calculate movements of deployed motus tags on migratory birds
library(tidyverse)
library(lubridate)
library(motus)
library(motusData)
library(DBI)
library(RSQLite)
library(ggmap)
library(maps) 

install.packages("rnaturalearth")

getwd()
setwd("/home/baniasss/164")

Sys.setenv(TZ = "UTC")

proj.num <- 164
#we downloded the project but we want to update it and specify a different location to save
# the data by entering your preferred filepath

sql.motus <- tagme(projRecv = proj.num, new = TRUE, update = TRUE, dir = "C:/Users/name/OneDrive/Documents")

Username
Password

# specify the filepath where your .motus file is
# saved, and the file name.

file.name <- dbConnect(SQLite(), "C:/Users/name/OneDrive/Documents/project-164.motus")

# get a list of tables in the .motus file specified
# above.

dbListTables(file.name)

# get a list of variables in the 'species' table in
# the .motus file.
dbListFields(file.name, "alltags")

#The metadata() function can be used to add the complete Motus metadata to your .motus file.
metadata(sql.motus, projectIDs = 164)

# this retrieves the 'alltags' table from the
# 'sql.motus' SQLite file we read in earlier

tbl.alltags <- tbl(sql.motus, "alltags")  # virtual table


#To omit runs identified as dubious by motusFilter we can use an anti-join.
# Number of rows with runs 3 hits or less

filter(tbl(sql.motus, "alltags"), runLen <= 4) %>% collect() %>% nrow()
820733


# Identify runs to remove
to_remove <- tbl(sql.motus, "runs") %>%
  select(runID, motusFilter) %>%
  filter(motusFilter == 0)

# Use anti-join to remove those runs from the alltags

tbl_filtered <- anti_join(tbl(sql.motus, "alltags"), to_remove, by = "runID")

# Number of rows with runs 4 hits or less after being filtered

filter(tbl_filtered, runLen <= 4) %>% collect() %>% nrow() 


#For example, the following code, adds a probability column to the sample project data,
#which is identical to the motusFilter column (i.e., by default filterByActivity() uses the same conditions). 

tbl_motusFilter <- filterByActivity(sql.motus, return = "all")

#These parameters can also be less strict (i.e., exclude more detections). 
#This example excludes all runs of length 2 or less and will exclude any runs less than length 4 from 
#hourBins which have more than 500 runs and where at least 95% of those runs have a run length of 2. 

tbl_relaxed <- filterByActivity(sql.motus, minLen = 2, maxLen = 4, 
                                maxRuns = 500, ratio = 0.95, return = "all")


#Preparing the data
#First we will use the filterByActivity() function to label dubious detections. 
#This returns all the data in the alltags view with a new probability column.

tbl.alltags <- filterByActivity(sql.motus, return = "all")


#FILTER AND MAKE A FLAT DATAFRAME

year2018 <- as.numeric(as.POSIXct("2018-01-01 00:00:00 UTC"))

df.alltags <- tbl.alltags %>% filter(ts> year2018) %>%
  
  df.alltags <- tbl.alltags %>% 
  mutate(recvLat = if_else((is.na(gpsLat)|gpsLat == 0|gpsLat == 999), 
                           recvDeployLat, gpsLat),
         recvLon = if_else((is.na(gpsLon)|gpsLon == 0|gpsLon == 999), 
                           recvDeployLon, gpsLon),
         recvAlt = if_else(is.na(gpsAlt), recvDeployAlt, gpsAlt)) %>%
  select(-noise, -slop, -burstSlop, -done, -bootnum, 
         -codeSet, -mfg, -nomFreq, -markerNumber, -markerType, 
         -deviceID, -recvDeployAlt, -speciesGroup,- recvAlt) %>%
  collect() %>%
  as.data.frame() %>%
  mutate(ts = as_datetime(ts),  # work with dates AFTER transforming to flat file
         tagDeployStart = as_datetime(tagDeployStart),
         tagDeployEnd = as_datetime(tagDeployEnd), 
         recvLat = plyr::round_any(recvLat, 0.05), 
         recvLon = plyr::round_any(recvLon, 0.05),
         recvDeployName = if_else(is.na(recvDeployName), 
                                  paste(recvLat, recvLon, sep=":"), 
                                  recvDeployName))

#Now that we have a nice clean data frame, let’s filter it to the ‘good’ detections. 
#We will filter to 1 for detections to keep and 0 for dubious detections. column.

df.alltags.sub <- filter(df.alltags, probability == 1)

#Let us also save the excluded detections for later analysis.

df.block.0 <- filter(df.alltags, probability == 0) %>%
  select(motusTagID, runID) %>%
  distinct()

#_______________________________________________________________________________________________________________

# simplify the data by summarizing by the runID. 
# If you want to summarize at a finer/coarser scale, you can also create other groups.  
# The simplest alternative is a rounded timestamp variable; for example by using 
# mutate(ts.h = plyr::round_any(ts, 3600) function call. 
# Other options are to just use date (e.g date = as_date(ts))


fun.getpath <- function(df) {
  df %>%
    filter(tagProjID == proj.num, # keep only tags registered to the sample project
           !is.na(recvLat) | !(recvLat == 0)) %>% # drops data without lon/lat
    group_by(motusTagID, runID, recvDeployName, ambigID, 
             tagDepLon, tagDepLat, recvLat, recvLon) %>%
    #summarizing by runID to get max run length and mean time stamp:
    summarize(max.runLen = max(runLen), ts.h = mean(ts)) %>% 
    arrange(motusTagID, ts.h)
} # end of function call

df.alltags.path <- fun.getpath(df.alltags.sub)

ggplot(data = filter(df.alltags.path, 
                     motusTagID %in% c(26067)), 
       aes(x = ts.h, y = recvLat)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  geom_point() + 
  geom_path() +
  facet_wrap(~ motusTagID, scales = "free", ncol = 2)

#Alternatively, you can produce a list of sites where the maximum run length of detections
#was never greater than (say) 4, which may sometimes (but not always!) indicate they are simply false detections.

df.alltags.sub %>%
  mutate(month = month(ts)) %>%
  group_by(recvDeployName, month) %>%
  summarize(max.rl = max(runLen)) %>%
  filter(max.rl < 5) %>%
  spread(key = month, value = max.rl)

#check for ambiguous detections: _____________________________________________________________________________

clarify(sql.motus)

#Let's get a vector of these,

df.ambigTags <- select(df.alltags.sub, ambigID, motusTagID) %>% 
  filter(!is.na(ambigID)) %>% distinct()

#Using our getpath function, we'll create paths and then plot these detections
# create a boolean variable for ambiguous detections:

df.alltags.path <- fun.getpath(filter(df.alltags.sub, 
                                      motusTagID %in% df.ambigTags$motusTagID, tagProjID == 
                                        proj.num)) %>% # create a boolean variable for ambiguous
  # detections:
  mutate(Ambiguous = !(is.na(ambigID)))

# to put all ambiguous tags from the same project
# on the same plot together, we need to create a
# new 'ambig tag' variable we call 'newID.

ambigTags.2 <- filter(df.alltags.sub) %>% select(ambigID, 
                                                 motusTagID) %>% filter(!is.na(ambigID)) %>% distinct() %>% 
  group_by(ambigID) %>% summarize(newID = paste(unique(ambigID), 
                                                toString(motusTagID), sep = ": ")) %>% left_join(df.ambigTags, 
                                                                                                 by = "ambigID")

# and merge that with df.alltags.path
df.alltags.path <- left_join(df.alltags.path, ambigTags.2, 
                             by = "motusTagID") %>% arrange(ts.h)

p <- ggplot(data = df.alltags.path, aes(ts.h, recvLat, 
                                        group = Ambiguous, colour = Ambiguous))
p + geom_point() + geom_path() + theme_bw() + facet_wrap(~newID, 
                                                         scales = "free", ncol = 2) + theme(axis.text.x = element_text(angle = 45, 
                                                                                                                       vjust = 1, hjust = 1))
#_______________________________________________________________________________________________________________
#Blocking:

#The detections in Nov. Dec. and Jan. are false detections so we filter them

df.block.1 <- filter(df.alltags.sub, month(ts) %in% 
                       c(11, 12, 1), motusTagID %in% c(26032, 26033, 26040, 26041, 26044, 26050, 26053, 
                                                       26055, 26056, 26058, 26060, 26085, 26061, 
                                                       26062, 26067, 26068, 26070, 26072, 26075, 26079,
                                                       26080,30297, 30291, 30285, 30280,30281, 30289, 30295, 30297)) %>% 
  select(mfgID, motusTagID, runID) %>% distinct()



df.block.2 <- filter(df.alltags.sub, month(ts) == 7, day(ts) %in% c(1, 2, 3, 4, 5, 6), 
                     motusTagID %in% c(26035, 26085, 26047, 26031, 26060, 26061, 26079, 26055)) %>% 
  select(motusTagID, runID) %>% distinct()

df.block.3 <- filter(df.alltags.sub, month(ts) == 7, day(ts) == 5, 
                     motusTagID == 26059) %>% 
  select(motusTagID, runID) %>% distinct()

df.block.17 <- filter(df.alltags.sub, month(ts) == 7, day(ts) == 7, 
                      motusTagID == 26035) %>% 
  select(motusTagID, runID) %>% distinct()

#1. 26083 (472):it is for BARS

df.block.4 <- filter(df.alltags.sub, ambigID == -359) %>% 
  select(motusTagID, runID) %>% distinct()

#2 & 3. 26075 (463) and 26082 (471): for BAR

df.block.5 <- filter(df.alltags.sub, ambigID == -358) %>% 
  select(motusTagID, runID) %>% distinct()

df.block.6 <- filter(df.alltags.sub, ambigID == -350) %>% 
  select(motusTagID, runID) %>% distinct()


#4. 26068 (456):it is for PUMA

df.block.7 <- filter(df.alltags.sub, ambigID == -356, 
                     motusTagID == 25007) %>% select(motusTagID, runID) %>% 
  distinct()

df.block.8 <- filter(df.alltags.sub, month(ts) == 12, motusTagID == 26058) %>% 
  select(mfgID, motusTagID, runID) %>% distinct()

#5. 26076 (464):it is for BARS


df.block.9 <- filter(df.alltags.sub, ambigID == -355) %>% 
  select(motusTagID, runID) %>% distinct()


#6. 26073 (461):it is for BARS

df.block.10 <- filter(df.alltags.sub, ambigID == -353) %>% 
  select(motusTagID, runID) %>% distinct()

#7. 26071 (459): it is for PUMA

df.block.11 <- filter(df.alltags.sub, ambigID == -352) %>% 
  select(motusTagID, runID) %>% distinct()

#8. 26077 (465):it is for BAR

df.block.12 <- filter(df.alltags.sub, ambigID == -351) %>% 
  select(motusTagID, runID) %>% distinct()

#9. 26084 (473):it is for BAR


df.block.13 <- filter(df.alltags.sub, ambigID == -342) %>% 
  select(motusTagID, runID) %>% distinct()

#10. 26081(470): it is for BARS

df.block.14 <- filter(df.alltags.sub, ambigID == -340) %>% 
  select(motusTagID, runID) %>% distinct()

#11.26070 (458):it is for PUMA

df.block.15 <- filter(df.alltags.sub, ambigID == -343, 
                      motusTagID == 25009) %>% select(motusTagID, runID) %>% 
  distinct() 

df.block.16 <- filter(df.alltags.sub, month(ts) == 9, day(ts) %in% c(18, 19, 20, 21, 22, 23, 24, 25, 26, 27), motusTagID == 26058) %>% select(motusTagID, runID) %>% distinct()



#12.26079 (467): it is for PUMA

df.block.18 <- filter(df.alltags.sub, ambigID == -360) %>% 
  select(motusTagID, runID) %>% distinct() 

df.block.19 <- filter(df.alltags.sub, month(ts) == 7, day(ts) == 4, motusTagID == 26035) %>% 
  select(motusTagID, runID) %>% distinct()

df.block.20 <- filter(df.alltags.sub, ambigID == -367) %>% 
  select(motusTagID, runID) %>% distinct() 

df.block.21 <- filter(df.alltags.sub, ambigID == -376) %>% 
  select(motusTagID, runID) %>% distinct() 

df.block.22 <- filter(df.alltags.sub, ambigID == -377) %>% 
  select(motusTagID, runID) %>% distinct()
################################################

df.block.23 <- filter(df.alltags.sub, ambigID == -1055) %>% 
  select(motusTagID, runID) %>% distinct()


df.block.24 <- filter(df.alltags.sub, ambigID == -1272) %>% 
  select(motusTagID, runID) %>% distinct()


df.block.25 <- filter(df.alltags.sub, ambigID == -1276) %>% 
  select(motusTagID, runID) %>% distinct()

df.block.26 <- filter(df.alltags.sub, ambigID == -1277) %>% 
  select(motusTagID, runID) %>% distinct()

#_____________________________________________________________________________________________________________


# combine our df.block data frames into a single
# dataframe, and add probability = 0 for filtered
# records.

df.block.all <- bind_rows(df.block.0, df.block.1, df.block.2, 
                          df.block.3, df.block.5, df.block.6, 
                          df.block.7, df.block.8, df.block.9, df.block.10, df.block.11, df.block.12,
                          df.block.13, df.block.14, df.block.15, df.block.17, df.block.18,
                          df.block.19, df.block.21,df.block.22
                          ,df.block.23,df.block.24,df.block.25,df.block.26) %>% mutate(probability = 0)

df.alltags.sub <- anti_join(df.alltags, df.block.all, by = c("runID", "motusTagID"))


#SAVE

saveRDS(df.alltags.sub, file = "/home/baniasss/164/dfAlltagsSub.rds")

#Read the tbl

df.alltags.sub <- readRDS("/home/baniasss/164/dfAlltagsSub.rds")

#___________________________________________________________________________________________________________

#We can manipulate existing variables or create new ones with dplyr's mutate function,
#here we'll convert ts to a POSIXct format, then make a new variable for year and day of year (doy).

df.alltags.sub <- df.alltags.sub %>%
  mutate(ts = as_datetime(ts, tz = "UTC"), # convert ts to POSIXct format
         year = year(ts), # extract year from ts
         doy = yday(ts)) %>% # extract numeric day of year from ts
  filter(!is.na(recvLat))
head(df.alltags.sub)

#Mapping__________________________________________________________________________________

fun.getpath <- function(df) {
  df %>%
    filter(tagProjID == proj.num, # keep only tags registered to the sample project
           !is.na(recvLat) | !(recvLat == 0)) %>% 
    group_by(motusTagID, runID, recvDeployName, ambigID, 
             tagDepLon, tagDepLat, recvLat, recvLon) %>%
    summarize(max.runLen = max(runLen), ts.h = mean(ts)) %>%
    arrange(motusTagID, ts.h) %>%
    data.frame()
} # end of function call

df.alltags.path <- fun.getpath(df.alltags.sub)


df.alltags.sub.path <- df.alltags.sub %>%
  filter(tagProjID == proj.num) %>% # only tags registered to project
  arrange(motusTagID, ts) %>%       # order by time stamp for each tag
  mutate(date = as_date(ts)) %>%    # create date variable
  group_by(motusTagID, date, recvDeployName, ambigID, 
           tagDepLon, tagDepLat, recvLat, recvLon)

df.alltags.path <- fun.getpath(df.alltags.sub.path)


#ggmap_______________________________________________________________________________________

library(ggmap)

ggmap::register_google('googleAPIkey')


gmap <-  get_stamenmap(bbox = c(left = -85, right = -75, bottom = 40, top = 45),
                       maptype = "terrain-background", # select maptype
                       zoom = 6) # zoom, must be a whole number

#_____________________________________________________________________________________________

#This year2017
#26030, 26031, 26032, 26033, 26034, 26035, 26036, 26037, 26038, 26039, 26040, 26041, 26042,
#26043, 26044, 26045, 26046, 26047, 26048, 26049, 26050, 26051, 26052, 26053, 26054, 26055, 26056,
#26057, 26058, 26059, 26060, 26061, 26062, 26063, 26064, 26065, 26066, 26067, 26068, 26069, 26070,
#26071, 26072, 26073, 26074, 26075, 26076, 26077, 26078, 26079, 26080, 26081, 26082, 26083, 26084, 26085

#This year2018
#30245
#[58] 30246 30247 30248 30249 30250 30251 30252 30253 30254 30255 30256 30257 30258 30259 30260 30261 30262 30263 30264
#[77] 30265 30266 30267 30268 30269 30270 30271 30272 30273 30274 30275 30276 30277 30278 30279 30280 30281 30282 30283
#[96] 30284 30285 30286 30287 30288 30289 30290 30291 30292 30293 30294 30295 30296 30297 30298 30875


#30245, 30246, 30247, 30248, 30249, 30250, 30251, 30252,
#30253, 30254, 30255, 30256, 30257, 30258, 30259, 30260, 30261,
#30262, 30263, 30264,30265, 30266, 30267, 30268, 30269, 30270,
#30271, 30272, 30273, 30274, 30275, 30276, 30277, 30278, 30279,
#30280, 30281, 30282, 30283, 30284, 30285, 30286, 30287, 30288,
#30289, 30290, 30291, 30292, 30293, 30294, 30295, 30296, 30297,
#30298, 30875

# just use the tags that we have examined carefully and filtered (in the previous chapter)
df.tmp <- filter(df.alltags.path, 
                 motusTagID %in% c(71168))

#OR

#df.tmp <- filter(df.alltags.path, 
#   motusTagID == 30291)
#__________________________________________________________________

df.tmp <- arrange(df.tmp, ts.h) # arange by hour
df.tmp <- as.data.frame(df.tmp)


p <- ggmap(gmap)
p + geom_point(data=df.tmp, 
               aes(recvLon, recvLat), pch=21, colour = "black", fill = "yellow") +
  geom_path(data=df.tmp, 
            aes(recvLon, recvLat, group=motusTagID, col = as.factor(motusTagID))) +
  theme_bw() + 
  scale_color_discrete(name="MotusTagID")

#Creating simple outline maps

na.lakes <- map_data(map = "lakes")
na.lakes <- na.lakes %>% mutate(long = long - 360)

# Include all of the Americas to begin
na.map <- map_data(map = "world2")
na.map <- na.map %>% filter(region %in% c("Canada", 
                                          "USA")) %>% mutate(long = long - 360)

# set limits to map based on locations of
# detections, ensuring they include the deployment
# locations
xmin <- min(df.tmp$recvLon, na.rm = TRUE) - 2
xmax <- max(df.tmp$recvLon, na.rm = TRUE) + 2
ymin <- min(df.tmp$recvLat, na.rm = TRUE) - 2
ymax <- max(df.tmp$recvLat, na.rm = TRUE) + 2


# map
ggplot(na.lakes, aes(long, lat)) + geom_polygon(data = na.map, 
                                                aes(long, lat, group = group), colour = "grey", 
                                                fill = "grey98") + geom_polygon(aes(group = group), 
                                                                                colour = "grey", fill = "white") + coord_map(projection = "mercator", 
                                                                                                                             xlim = c(xmin, xmax), ylim = c(ymin, ymax)) + xlab("") + 
  ylab("") + theme_bw() + geom_path(data = df.tmp, 
                                    aes(recvLon, recvLat, group = as.factor(motusTagID), 
                                        colour = as.factor(motusTagID))) + geom_point(data = df.tmp, 
                                                                                      aes(tagDeployLon, tagDeployLat), colour = "black", 
                                                                                      shape = 4) + scale_colour_discrete("motusTagID")


#_______________________________________________________________

#transition:calculating the flight distances

transitions <- siteTrans(filter(df.alltags, motusTagID == 26060), 
                         latCoord = "recvDeployLat", lonCoord = "recvDeployLon")


transitions<- filter(transitions, rate<= 6.8)

transitions<- select(transitions, motusTagID,ts.x, ts.y,lat.x, lon.x, lat.y, lon.y, dist, rate, bearing) 



#transitions<-transitions[-c(11), ]

max(transitions$dist)/1000


sum(transitions$dist )/1000


write_csv(transitions, "/home/baniasss/164/26060.CSV")

head(transitions)


sum(transitions$dist)

#_________________________________________________________________________
#making some plots for better underestanding


#This year2017
#  [1] 26030 26031 26032 26033 26034 26035 26036 26037 26038 26039 26040 26041 26042 26043 26044 26045 26046 26047 26048
#[20] 26049 26050 26051 26052 26053 26054 26055 26056 26057 26058 26059 26060 26061 26062 26063 26064 26065 26066 26067
#[39] 26068 26069 26070 26071 26072 26073 26074 26075 26076 26077 26078 26079 26080 26081 26082 26083 26084 26085


#This year2018
#30245
#[58] 30246 30247 30248 30249 30250 30251 30252 30253 30254 30255 30256 30257 30258 30259 30260 30261 30262 30263 30264
#[77] 30265 30266 30267 30268 30269 30270 30271 30272 30273 30274 30275 30276 30277 30278 30279 30280 30281 30282 30283
#[96] 30284 30285 30286 30287 30288 30289 30290 30291 30292 30293 30294 30295 30296 30297 30298 30875 

#station and signals

df.site <- filter(df.alltags.sub, motusTagID == 71169)

p <- ggplot(data = df.site, aes(ts, sig, colour = paste(recvDeployLat, 
                                                        recvDeployLon,
                                                        recvDeployName, 
                                                        sep = ": ")))

p + geom_point() + scale_colour_discrete(name = "Lat/Lon and\nStation Name") + 
  theme_bw() + facet_wrap(~as_date(ts), scales = "free_x")

#signal strength and time at a particular month

filter(df.alltags.sub, motusTagID == 71168, 
       month(ts) == 7) %>% ggplot(aes(ts, sig, group = recvDeployName, 
                                      colour = recvDeployName)) + geom_point() +
  theme_bw() + 
  xlab("Time") + ylab("Signal strength") + facet_grid(recvDeployLon ~ .)

#signal, ts, location

df.30295.ts <- filter(df.alltags.sub, motusTagID == 71168) %>% 
  mutate(LatLonStationName = paste(recvDeployLat, recvDeployLon, 
                                   recvDeployName, sep = ": "))

p <- ggplot(data = df.30295.ts, aes(ts, sig, colour = LatLonStationName))

p + geom_point() + theme_bw()

#another plot with ts and longt

df.30295.lo <- filter(df.alltags.sub, motusTagID == 71168) %>% 
  mutate(sig = ifelse(sig > 0, sig * -1, sig))

p <- ggplot(data = df.30295.lo , aes(recvLon, ts, colour = paste(recvLat, 
                                                                 recvLon, recvDeployName, sep = ": ")))
p + geom_point() + theme_bw() + scale_colour_discrete(name = "Lat/Lon and\nStation Name")


plotAllTagsCoord(filter(tbl.alltags, motusTagID == 30295), tagsPerPanel = 1)

# plotRouteMap:

plotRouteMap(sql.motus, maptype = "terrain", latCentre = 44, 
             lonCentre = -70, zoom = 5, recvStart = "2018-07-01", 
             recvEnd = "2018-11-01")

#Plot signal strength vs time for all tags detected at a specified site, coloured by antenna
#Plot select tags for site : Piskwamish


df.alltags.sub.30295 <- mutate(df.alltags.sub)

p <- ggplot(filter(df.alltags.sub, motusTagID == 71168, 
                   ts > ymd("2022-01-01"), ts < ymd("2022-09-01")), 
            aes(ts, sig, col = as.factor(antBearing)))
p + theme_bw() + geom_point() + xlab("Time of day") + 
  ylab("Signal strength") + scale_color_discrete(name = "Antenna bearing") + 
  facet_grid(recvDeployName ~ .)

#Creates a data.frame consisting of only detections of tags that are detected
#at two or more receivers at the same time.

simSites <- simSiteDet(df.alltags.sub)
head(simSites)

#Creates a summary of the first and last detection at a site, 
#the length of time between first and last detection, the number of tags,
#and the total number of detections at a site. Plots total number of detections
#across all tags, and total number of tags detected at each site.

site_summary <- siteSum(filter(df.alltags.sub, recvDeployName %in% 
                                 c("Niapiskau", "Netitishi", "Old Cur", "Washkaugou")), 
                        units = "mins")

#Creates a summary of the first and last daily detection at a site, 
#the length of time between first and last detection, the number of tags,
#and the total number of detections at a site for each day. Same as siteSum, but daily by site.

daily_site_summary <- siteSumDaily(tbl.alltags, units = "mins")
head(daily_site_summary)

# Creates a summary for each tag of it's first and last detection time,
#first and last detection site, length of time between first and last detection,
#straight line distance between first and last detection site, rate of movement, and bearing

df.tag_summary <- tagSum(df.alltags.sub)
head(tag_summary)

#Creates a summary for each tag of it's first and last detection time at each site, 
#length of time between first and last detection of each site,
#and total number of detections at each site.

tag_site_summary <- tagSumSite(filter(df.alltags.sub, motusTagID %in% 
                                        c(71168)))
head(tag_site_summary)




#To update your existing packages:
update.packages()




#To see the tags:______________________________________________________________

tbl.tags <- tbl(sql.motus, "tags")
df.tags <- tbl.tags %>% filter(projectID == proj.num) %>% 
  collect() %>% as.data.frame()

unique(df.tags$tagID)

write_csv(df.tags, "C:/motus/164/164.CSV")



#36029,	36030,	36031,	36032,	36033,	36034,	36035,	36036,	36037,	36038,	36039,	36040,	36041,	36042,
#36043,	36044,	36045,	36046,	36047,	36048,	36049,	36050,	36051,	36052,	36053,	36054,	36055,	36056,
#36057,	36058,	36059,	36060,	36061,	36062,	36063,	36064,	36065,	36066,	36067,	36068,	36069,	36070,
#36071,	36072,	36073,	36074,	36075,	36076,	36077,	36078,	36128,	36079,	36080,	36081,	36082,	36083,
#36084,	36085,	36086,	36087,	36088,	36090,	36093,	36109,	35364,	35365,	35366,	35367,	35368,	35369,
#35370,	35371,	35372,	35373,	35374,	35375,	35376,	35377,	35378,	35379,	35380,	35381,	35382,	35383,
#35384,	35385,	35386,	35387,	35388,	35389,	35390,	35391,	35392,	35393,	35394,	35395,	35396,	35397,
#35398,	35399,	35400,	35401,	35402,	35403,	35404,	35405,	35406,	35407,	35408,	35409,	35410,	35411,
#35412,	35413,	35414,	35415,	35416,	35417,	35418,	35419,	35420,	35421,	35422,	35423,	35424,	35426,
#35427,	35428,	35429,	35430,	35431,	35432,	35433,	35434,	35435,	35436,	35437,	35438,	35439,	35440,
#35463,	35441,	35442,	35443,	35445,	35446,	35447,	35448,	35449,	35450,	35451,	35452,	35453,	35454,
#35455,	35456,	35457,	35458,	35459,	35460,	35461,	35462

#______________________________________________________________________________________________________
#______________________________________________________________________________________________________
#Simplify the data for plotting

fun.getpath <- function(df) 
{
  df %>%
    filter(tagProjID == proj.num, # keep only tags registered to the sample project
           !is.na(recvLat) | !(recvLat == 0)) %>% # drops data without lon/lat
    group_by(motusTagID, runID, recvDeployName, ambigID, 
             tagDeployLon, tagDeployLat, recvLat, recvLon) %>%
    #summarizing by runID to get max run length and mean time stamp:
    summarize(max.runLen = max(runLen), ts.h = mean(ts)) %>% 
    arrange(motusTagID, ts.h)
} # end of function call

df.alltags.path <- fun.getpath(df.alltags.sub)

p <- ggplot(data = filter(df.alltags.path, motusTagID %in% 
                            c(36029)), aes(ts.h, recvLat))
p + geom_point() + geom_path() + theme_bw() + facet_wrap(~motusTagID, 
                                                         scales = "free", ncol = 2) + theme(axis.text.x = element_text(angle = 45, 
                                                                                                                       vjust = 1, hjust = 1))
#__________________________________________________________________________________



t <- filter(df.alltags.sub, motusTagID == 26040)


write_csv(t, "/home/baniasss/164/point/26040.CSV")


transitions<- siteTrans(filter(df.alltags, motusTagID == 26040),
                        latCoord = "recvDeployLat", lonCoord = "recvDeployLon")


write_csv(transitions, "/home/baniasss/164/26040.CSV")

transitions <- siteTrans(filter(df.alltags, motusTagID == 16037), 
                         latCoord = "recvDeployLat", lonCoord = "recvDeployLon")

