library(dplyr)
library(readr)
library(tidyverse)
library(stringr)
library(geometry)
library(RANN)
library(reshape2)
library(zoo)
library(magrittr)
library(data.table)
library(purrr)

Vidsync_Data1 <- read.csv(file = "Downloads/R-Program/NRES_746/NRES746-FinalProject/rawData_processing/VidSync_Raw-csv/before/control/Porter_BACI_JesusToast_30June2018_Part1.csv",
                         skip = 2, 
                         col.names = c("objects", "event", "timecode", "time", "X", "Y", "Z", "pld_error", "projection_error", "nearest_camera_distance", "screen_coordinates"),
                         colClasses = c("character", "character", "character", "double", "double", "double", "double", "double", "double", "double", "double")) %>% 
  mutate(sample_event = "before") %>% #These are the columns to be added to each video prior to analysis
  mutate(site = "18.1")

Vidsync_Data2 <- read.csv(file = "Downloads/R-Program/NRES_746/NRES746-FinalProject/rawData_processing/VidSync_Raw-csv/before/control/Porter_BACI_GolfBall_30June2018_Part3.csv",
                          skip = 2, 
                          col.names = c("objects", "event", "timecode", "time", "X", "Y", "Z", "pld_error", "projection_error", "nearest_camera_distance", "screen_coordinates"),
                          colClasses = c("character", "character", "character", "double", "double", "double", "double", "double", "double", "double", "double")) %>% 
  mutate(sample_event = "before") %>% #These are the columns to be added to each video prior to analysis
  mutate(site = "18.2")

Vidsync_Data3 <- read.csv(file = "Downloads/R-Program/NRES_746/NRES746-FinalProject/rawData_processing/VidSync_Raw-csv/before/control/Porter_BACI_Waterfall_18.4_June_14_18.csv",
                          skip = 2, 
                          col.names = c("objects", "event", "timecode", "time", "X", "Y", "Z", "pld_error", "projection_error", "nearest_camera_distance", "screen_coordinates"),
                          colClasses = c("character", "character", "character", "double", "double", "double", "double", "double", "double", "double", "double")) %>% 
  mutate(sample_event = "before") %>% #These are the columns to be added to each video prior to analysis
  mutate(site = "18.3")

Vidsync_Data4 <- read.csv(file = "Downloads/R-Program/NRES_746/NRES746-FinalProject/rawData_processing/VidSync_Raw-csv/before/impact/Porter_BACI_RoachRun_29June2018_Part3.csv",
                          skip = 2, 
                          col.names = c("objects", "event", "timecode", "time", "X", "Y", "Z", "pld_error", "projection_error", "nearest_camera_distance", "screen_coordinates"),
                          colClasses = c("character", "character", "character", "double", "double", "double", "double", "double", "double", "double", "double")) %>% 
  mutate(sample_event = "before") %>% #These are the columns to be added to each video prior to analysis
  mutate(site = "18.4")

Vidsync_Data5 <- read.csv(file = "Downloads/R-Program/NRES_746/NRES746-FinalProject/rawData_processing/VidSync_Raw-csv/before/impact/Porter_BACI_BigBrego_29June2018_Part1.csv",
                          skip = 2, 
                          col.names = c("objects", "event", "timecode", "time", "X", "Y", "Z", "pld_error", "projection_error", "nearest_camera_distance", "screen_coordinates"),
                          colClasses = c("character", "character", "character", "double", "double", "double", "double", "double", "double", "double", "double")) %>% 
  mutate(sample_event = "before") %>% #These are the columns to be added to each video prior to analysis
  mutate(site = "18.5")

Vidsync_Data6 <- read.csv(file = "Downloads/R-Program/NRES_746/NRES746-FinalProject/rawData_processing/VidSync_Raw-csv/before/impact/Porter_BACI_HalfTire_29June2018_Part2.csv",
                          skip = 2, 
                          col.names = c("objects", "event", "timecode", "time", "X", "Y", "Z", "pld_error", "projection_error", "nearest_camera_distance", "screen_coordinates"),
                          colClasses = c("character", "character", "character", "double", "double", "double", "double", "double", "double", "double", "double")) %>% 
  mutate(sample_event = "before") %>% #These are the columns to be added to each video prior to analysis
  mutate(site = "18.6")

Vidsync_Data7 <- read.csv(file = "Downloads/R-Program/NRES_746/NRES746-FinalProject/rawData_processing/VidSync_Raw-csv/after/control/Porter_BACI_JesusToast_6July2018_Part1.csv",
                          skip = 2, 
                          col.names = c("objects", "event", "timecode", "time", "X", "Y", "Z", "pld_error", "projection_error", "nearest_camera_distance", "screen_coordinates"),
                          colClasses = c("character", "character", "character", "double", "double", "double", "double", "double", "double", "double", "double")) %>% 
  mutate(sample_event = "after") %>% #These are the columns to be added to each video prior to analysis
  mutate(site = "18.1")

Vidsync_Data8 <- read.csv(file = "Downloads/R-Program/NRES_746/NRES746-FinalProject/rawData_processing/VidSync_Raw-csv/after/control/Porter_BACI_Golfball_6July2018_Part1.csv",
                          skip = 2, 
                          col.names = c("objects", "event", "timecode", "time", "X", "Y", "Z", "pld_error", "projection_error", "nearest_camera_distance", "screen_coordinates"),
                          colClasses = c("character", "character", "character", "double", "double", "double", "double", "double", "double", "double", "double")) %>% 
  mutate(sample_event = "after") %>% #These are the columns to be added to each video prior to analysis
  mutate(site = "18.2")

Vidsync_Data9 <- read.csv(file = "Downloads/R-Program/NRES_746/NRES746-FinalProject/rawData_processing/VidSync_Raw-csv/after/control/Porter_BACI_Waterfall_18.4_July_6_18.csv",
                          skip = 2, 
                          col.names = c("objects", "event", "timecode", "time", "X", "Y", "Z", "pld_error", "projection_error", "nearest_camera_distance", "screen_coordinates"),
                          colClasses = c("character", "character", "character", "double", "double", "double", "double", "double", "double", "double", "double")) %>% 
  mutate(sample_event = "after") %>% #These are the columns to be added to each video prior to analysis
  mutate(site = "18.3")

Vidsync_Data10 <- read.csv(file = "Downloads/R-Program/NRES_746/NRES746-FinalProject/rawData_processing/VidSync_Raw-csv/after/impact/Porter_BACI_RoachRun_5July2018_Part1.csv",
                          skip = 2, 
                          col.names = c("objects", "event", "timecode", "time", "X", "Y", "Z", "pld_error", "projection_error", "nearest_camera_distance", "screen_coordinates"),
                          colClasses = c("character", "character", "character", "double", "double", "double", "double", "double", "double", "double", "double")) %>% 
  mutate(sample_event = "after") %>% #These are the columns to be added to each video prior to analysis
  mutate(site = "18.4")

Vidsync_Data11 <- read.csv(file = "Downloads/R-Program/NRES_746/NRES746-FinalProject/rawData_processing/VidSync_Raw-csv/after/impact/Porter_BACI_BigBrego_5July2018_Part1.csv",
                          skip = 2, 
                          col.names = c("objects", "event", "timecode", "time", "X", "Y", "Z", "pld_error", "projection_error", "nearest_camera_distance", "screen_coordinates"),
                          colClasses = c("character", "character", "character", "double", "double", "double", "double", "double", "double", "double", "double")) %>% 
  mutate(sample_event = "after") %>% #These are the columns to be added to each video prior to analysis
  mutate(site = "18.5")

Vidsync_Data12 <- read.csv(file = "Downloads/R-Program/NRES_746/NRES746-FinalProject/rawData_processing/VidSync_Raw-csv/after/impact/Porter_BACI_HalfTire_5July2018_Part3.csv",
                          skip = 2, 
                          col.names = c("objects", "event", "timecode", "time", "X", "Y", "Z", "pld_error", "projection_error", "nearest_camera_distance", "screen_coordinates"),
                          colClasses = c("character", "character", "character", "double", "double", "double", "double", "double", "double", "double", "double")) %>% 
  mutate(sample_event = "after") %>% #These are the columns to be added to each video prior to analysis
  mutate(site = "18.6")

master_df <- rbind(Vidsync_Data1, Vidsync_Data2, Vidsync_Data3, Vidsync_Data4, Vidsync_Data5, Vidsync_Data6,
                   Vidsync_Data7, Vidsync_Data8, Vidsync_Data9, Vidsync_Data10, Vidsync_Data11, Vidsync_Data12)

VidsyncFormat <- function(dataframe=MasterFormat){
    dataframe %>% 
    mutate(species = str_extract(objects, "Omykiss|Okisutch")) %>% 
    mutate(subsample = str_extract(objects, "\\d")) %>%
    mutate(index = str_extract(objects, "\\h\\d{1,2}")) %>%
    transform(index = as.numeric(index), subsample = as.numeric(subsample)) %>% 
    mutate(behavior = str_extract(event, pattern = "[A-z]+")) %>% 
    select(sample_event, site, subsample, index, species, behavior, time, X, Y, Z) %>% 
    na.omit(species) %>% 
    arrange(subsample, index, time)
} 
MasterFormat <- VidsyncFormat(master_df)

ExtractVolume <- function(dataframe){ #Needs fix for multiple sample events
    dataframe %>% 
    filter(!behavior == "Length") %>%
    arrange(subsample, index, time) %>%
    select(sample_event, site, subsample, index, behavior, X, Y, Z) %>%
    group_by(sample_event, site, index) %>% 
    add_count(index) %>% 
    filter(!n <= 3) %>% 
    split(.$index) %>% 
    map(~ as.matrix(.[ , c("X","Y", "Z")])) %>%
    map(~ convhulln(., option="FA")) %>% 
    map_dfr("vol") %>% 
    tidyr::gather(index, volume) %>% 
    transform(index = as.numeric(index)) %>% 
    left_join(dataframe, by = "index") %>% 
    distinct(volume, .keep_all = T) %>% 
    select(sample_event, site, subsample, index, volume) 
}
vol1 <- ExtractVolume(VidsyncFormat(Vidsync_Data1))
vol2 <- ExtractVolume(VidsyncFormat(Vidsync_Data2))
vol3 <- ExtractVolume(VidsyncFormat(Vidsync_Data3))
vol4 <- ExtractVolume(VidsyncFormat(Vidsync_Data4))
vol5 <- ExtractVolume(VidsyncFormat(Vidsync_Data5))
vol6 <- ExtractVolume(VidsyncFormat(Vidsync_Data6))
vol7 <- ExtractVolume(VidsyncFormat(Vidsync_Data7))
vol8 <- ExtractVolume(VidsyncFormat(Vidsync_Data8))
vos9 <- ExtractVolume(VidsyncFormat(Vidsync_Data9))
vol10 <- ExtractVolume(VidsyncFormat(Vidsync_Data10))
vol11 <- ExtractVolume(VidsyncFormat(Vidsync_Data11))
vol12 <- ExtractVolume(VidsyncFormat(Vidsync_Data12))
MasterVolume <- rbind(vol1, vol2, vol3, vol4, vol5, vol6, vol7, vol8, vol9, vol10, vol11, vol12) 

DistancePerTime <- function(dataframe=MasterFormat){
    dataframe %>% 
    filter(!behavior == "Length") %>%
    group_by(sample_event, site, subsample, index) %>% 
    mutate(distance_travelled_X_cm = X - lag(X, default = first(X))) %>%
    mutate(distance_travelled_Y_cm = Y - lag(Y, default = first(Y))) %>%
    mutate(distance_travelled_Z_cm = Z - lag(Z, default = first(Z))) %>%
    mutate(fish_distance_travelled_cm = sqrt((distance_travelled_X_cm)^2
                                             + (distance_travelled_Y_cm)^2
                                             + (distance_travelled_Z_cm)^2)) %>%
    mutate(distance_cm_per_sec = fish_distance_travelled_cm /(time - lag(time, default = first(time)))) %>%
    filter(!distance_cm_per_sec == Inf) %>% 
    group_by(sample_event, site, index) %>% 
    mutate(DistPerTime_Median = median(distance_cm_per_sec)) %>%
    mutate(DistPerTime_Mean = mean(distance_cm_per_sec)) %>% 
    select(sample_event, site, subsample, index, time, distance_cm_per_sec, DistPerTime_Mean, DistPerTime_Median) %>% 
    unique()
}
MasterDistance <- DistancePerTime(MasterFormat)

BehaviorTypes <- function(dataframe=MasterFormat) {
    dataframe %>% 
    select(sample_event, site, subsample, index, time, behavior) %>% 
    filter(!behavior == "Length") %>%
    filter(!behavior == "Point") %>% 
    mutate(drift_forage = if_else(behavior == "Drift_Forage", 1, 0)) %>% 
    mutate(benthic_forage = if_else(behavior == "Benthic_Forage", 1, 0)) %>% 
    mutate(search_forage = if_else(behavior == "Search_Forage", 1, 0)) %>% 
    mutate(movement = if_else(behavior == "Movement", 1, 0)) %>% 
    mutate(surface_strike = if_else(behavior == "Surface_Strike", 1, 0)) %>% 
    mutate(attack = if_else(behavior == c("Attack","Nip"), 1, 0)) #%>% 
    #mutate(nip = if_else(behavior == "Nip", 1, 0))
}
MasterBehaviors <- BehaviorTypes(MasterFormat)

GetLength <- function(dataframe=MasterFormat){
    dataframe %>% 
    filter(behavior == "Length") %>% 
    group_by(site, sample_event, index, time) %>% 
    mutate(length = 10*sqrt(((X-lag(X, default = first(X)))^2)+
                              ((Y-lag(Y, default = first(Y)))^2)+
                              ((Z-lag(Z, default = first(Z)))^2))) %>% 
    filter(!length == "0") %>% 
    ungroup() %>%  
    group_by(site, sample_event, index) %>% 
    mutate(length_mm = mean(length)) %>% 
    select(sample_event, site, subsample, index, length_mm) %>% 
    distinct(index, .keep_all = TRUE)
} #last join
MasterLengths <- GetLength(MasterFormat)

GetNND <- function(dataframe=MasterFormat){ #Needs fix for multiple sample events
    dataframe %>% 
    dplyr::filter(!behavior == "Length") %>% 
    dplyr::arrange(time) %>% 
    dplyr::select(time, X, Y, Z) %>% 
    dplyr::add_count(time) %>% 
    dplyr::filter(!n <= 1) %>% 
    split(.$time) %>% 
    purrr::map(~ as.matrix(.[ , c("X","Y","Z")])) %>%
    purrr::map(~ dist(., method = "euclidean")) %>% 
    purrr::map(~ melt(as.matrix(.), varnames = c("row", "col"))) %>%
    purrr::map(~ filter(., !value == "0")) %>% 
    purrr::map(~ filter(., value == min(value))) %>% 
    purrr::map(~ distinct(., value)) %>% 
    map_dfr(~ as.numeric(str_extract(., pattern = "\\d{1,3}.\\d{1,5}"))) %>% 
    gather(time, NND_cm) %>% 
    transform(time = as.numeric(time)) %>% 
    left_join(dataframe, by = "time") %>% 
    arrange(site, desc(sample_event), subsample, index, time) %>% 
    select(sample_event, site, subsample, index, time, NND_cm) 
}

NND1 <- GetNND(VidsyncFormat(Vidsync_Data1))
NND2 <- GetNND(VidsyncFormat(Vidsync_Data2))
#NND3 <- GetNND(VidsyncFormat(Vidsync_Data3))
NND4 <- GetNND(VidsyncFormat(Vidsync_Data4))
NND5 <- GetNND(VidsyncFormat(Vidsync_Data5))
NND6 <- GetNND(VidsyncFormat(Vidsync_Data6))
NND7 <- GetNND(VidsyncFormat(Vidsync_Data7))
NND8 <- GetNND(VidsyncFormat(Vidsync_Data8))
NND9 <- GetNND(VidsyncFormat(Vidsync_Data9))
NND10 <- GetNND(VidsyncFormat(Vidsync_Data10))
NND11 <- GetNND(VidsyncFormat(Vidsync_Data11))
NND12 <- GetNND(VidsyncFormat(Vidsync_Data12))
MasterNND <- rbind(NND1, NND2, NND4, NND5, NND6, NND7, NND8, NND9, NND10, NND11, NND12)

df_Final1 <-  left_join(MasterBehaviors, MasterDistance, by = c("sample_event" = "sample_event", 
                                                           "site" = "site", 
                                                           "subsample" = "subsample",
                                                           "index" = "index", 
                                                           "time" = "time"))
df_Final2 <- left_join(df_Final1, MasterNND, by = c("sample_event" = "sample_event", 
                                                    "site" = "site", 
                                                    "subsample" = "subsample",
                                                    "index" = "index", 
                                                    "time" = "time"))
df_Final3 <- left_join(df_Final2, MasterLengths, by = c("sample_event" = "sample_event", 
                                                     "site" = "site", 
                                                     "subsample" = "subsample",
                                                     "index" = "index"))
df_Final <- left_join(df_Final3, MasterVolume, by = c("sample_event" = "sample_event", 
                                                      "site" = "site", 
                                                      "subsample" = "subsample",
                                                      "index" = "index")) %>% 
  select(sample_event, site, subsample, index, time, length_mm, behavior, distance_cm_per_sec, DistPerTime_Mean, DistPerTime_Median, NND_cm, volume, drift_forage, search_forage, benthic_forage, movement, surface_strike, attack) %>% 
  arrange(desc(sample_event), site, subsample, index, time)
filename <- file.path("Downloads/R-Program/NRES_746/NRES746-FinalProject/rawData_processing/", "finalDataFrame.csv")
write.csv(x = df_Final, file = filename)