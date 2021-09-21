```{r setup, include=FALSE, message=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r}
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
```

```{r}    
Vidsync_Data <- read.csv(file = "VidSync_Raw-csv/after/impact/Porter_BACI_HalfTire_5July2018_Part3.csv",
                         skip = 2, 
                         col.names = c("objects", "event", "timecode", "time", "X", "Y", "Z", "pld_error", "projection_error", "nearest_camera_distance", "screen_coordinates"),
                         colClasses = c("character", "character", "character", "double", "double", "double", "double", "double", "double", "double", "double")) %>% 
  mutate(date = as.Date("2018-07-05")) %>% #These are the columns to be added to each video prior to analysis
  mutate(sample_event = "After") %>% #These are the columns to be added to each video prior to analysis
  mutate(site = "18.8")
Vidsync_Data
```
head(Vidsync_Data)

#Import VidSync Data into proper format
```{r}
VidsyncFormat <- function(dataframe){
  dataframe <- dataframe %>% 
    mutate(species = str_extract(objects, "Omykiss|Okisutch")) %>% 
    mutate(subsample = str_extract(objects, "\\d")) %>%
    mutate(index = str_extract(objects, "\\h\\d{1,2}")) %>%
    transform(index = as.numeric(index), subsample = as.numeric(subsample)) %>% 
    mutate(behavior = str_extract(event, pattern = "[A-z]+")) %>% 
    select(sample_event, date, site, subsample, index, species, behavior, time, X, Y, Z) %>% 
    na.omit(species) %>% 
    arrange(subsample, index, time)
} #We would run the function "VidsyncFormat(Vidsync_Data)" to create a clean dataset to work with for the rest of the functions

test <- VidsyncFormat(Vidsync_Data)
head(test)
```
#Volume Macro
```{r}
ExtractVolume <- function(dataframe){
  dataframe %>% 
  filter(!behavior == "Length") %>%
  filter(!behavior == "Surface_Shots") %>%
  filter(!behavior == "Attack") %>%
  filter(!behavior == "Nip") %>%
  arrange(subsample, index, time) %>%
  select(sample_event, date, site, subsample, index, behavior, X, Y, Z) %>% 
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
  select(site, sample_event, subsample, index, volume) %>% 
  arrange(site, desc(sample_event), subsample, index)
}
volTest <- ExtractVolume(test)
```

#Swimming Speed Macro
```{r}
DistancePerTime <- function(dataframe){
    dataframe %>% 
    filter(!behavior == "Length") %>%
    filter(!behavior == "Surface_Strike") %>% 
    filter(!behavior == "") %>% 
    group_by(subsample, index) %>% 
    mutate(distance_travelled_X_cm = X - lag(X, default = first(X))) %>%
    mutate(distance_travelled_Y_cm = Y - lag(Y, default = first(Y))) %>%
    mutate(distance_travelled_Z_cm = Z - lag(Z, default = first(Z))) %>%
    mutate(fish_distance_travelled_cm = sqrt((distance_travelled_X_cm)^2
                                        + (distance_travelled_Y_cm)^2
                                        + (distance_travelled_Z_cm)^2)) %>%
    mutate(distance_cm_per_sec = fish_distance_travelled_cm /(time - lag(time, default = first(time)))) %>%
    filter(!distance_cm_per_sec == Inf) %>% 
    mutate(DistPerTime_Median = median(fish_distance_travelled_cm)) %>%
    mutate(DistPerTime_Mean = mean(fish_distance_travelled_cm)) %>% 
    select(sample_event, date, site, subsample, index, time, behavior, DistPerTime_Mean, DistPerTime_Median) %>% 
    unique()
}
test2 <- DistancePerTime(test)
test2
```  

#Behavior Type Macro
```{r}
BehaviorTypes <- function(dataframe) {
    dataframe %>% 
    select(subsample, index, behavior) %>% 
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
test4 <- BehaviorTypes(VidsyncFormat(Vidsync_Data))
head(test4)
```

#Length Macro
```{r}
GetLength <- function(dataframe){
  length <- dataframe %>% 
    filter(behavior == "Length") %>% 
    group_by(index, time) %>% 
    mutate(length = 10*sqrt(((X-lag(X, default = first(X)))^2)+
                            ((Y-lag(Y, default = first(Y)))^2)+
                            ((Z-lag(Z, default = first(Z)))^2))) %>% 
    filter(!length == "0") %>% 
    ungroup() %>% 
    group_by(index) %>% 
    mutate(length_mm = mean(length)) %>% 
    select(sample_event, date, site, subsample, index, length_mm) %>% 
    distinct(index, .keep_all = TRUE)
}

test3 <- GetLength(VidsyncFormat(Vidsync_Data))
test3
```
#Modified NND
```{r}
GetNND1 <- function(dataframe){
for (i in 1:63) {
    dataframe %>% 
    filter(!behavior == "Length") %>% 
    arrange(time) %>% 
    mutate(time_count = rleid(time)) %>% 
    filter(i == time_count) %>% 
    select(X, Y, Z) %>% 
    dist(p = 2L) %>% 
    as.matrix() %>% 
    melt(varnames = c("row", "col")) %>% 
    filter(!value == 0) %>% 
    filter(value == min(value)) %>% 
    distinct(value) %>% 
    as.data.frame() %>% 
    print()
  }
}
GetNND1(test)

GetNND <- function(dataframe){
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
  select(site, sample_event, subsample, index, time, NND_cm) 
}

GetNND(test)

test %>% 
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
  left_join(test, by = "time") %>% 
  arrange(site, desc(sample_event), subsample, index, time) %>% 
  select(site, sample_event, subsample, index, time, NND_cm)
```
