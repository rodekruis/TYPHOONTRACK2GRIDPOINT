library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(tmap)
library(viridis)
library(maps)
library(ggmap)
library(httr)
library(sf)

# read data
TRACK_DATA <- read.csv("C:/Users/ATeklesadik/Rode Kruis/510 - Data preparedness - Documents/[PRJ] FbF - Philippines - (PMF)/Analysis/historical_best_track_metoc_navy/best_track.csv")

#wind_track_haiyan<- read.csv("C:/Users/ATeklesadik/Rode Kruis/510 - Data preparedness - Documents/[PRJ] FbF - Philippines - (PMF)/model_improvement/wind_track_haiyan.csv")
grid_points <- read.csv("C:/Users/ATeklesadik/Rode Kruis/510 - Data preparedness - Documents/[PRJ] FbF - Philippines - (PMF)/model_improvement/county_points_php.csv", sep=";")
grid_points_adm3<- read.csv("C:/Users/ATeklesadik/Rode Kruis/510 - Data preparedness - Documents/[PRJ] FbF - Philippines - (PMF)/model_improvement/grid_points_admin3.csv", sep=",")

TRACK_DATA<- TRACK_DATA %>%
  dplyr::mutate(typhoon = paste0(TRACK_DATA$STORMNAME,substr(TRACK_DATA$YYYYMMDDHH, 1, 4)))

landmask <- readr::read_csv("C:/SOFTWARES/stormwindmodel/data-raw/landseamask_ph1.csv",
                            col_names = c("longitude", "latitude", "land")) %>%
  dplyr::mutate(land = factor(land, levels = c(1, 0), labels = c("land", "water")))


typhoon_names <- list('Bopha','Goni','Hagupit','Haima','Haiyan','Kalmaegi','Koppu','Melor','Nock-Ten','Rammasun','Sarika','Utor')

typhoon_names<-tolower(typhoon_names)
 
TRACK_DATA<-TRACK_DATA[tolower(TRACK_DATA$STORMNAME) %in% typhoon_names, ]

#setwd("C:/Users/ATeklesadik/Rode Kruis/510 - Data preparedness - Documents/[PRJ] FbF - Philippines - (PMF)/Analysis/historical_best_track_metoc_navy")
#BASIN,CY,YYYYMMDDHH,TECHNUM,TECH,TAU,LatN/S,LonE/W,vmax,mslp,ty,RAD,WINDCODE,RAD1,RAD2,RAD3,RAD4,RADP,RRP,MRD,GUSTS,EYE,SUBREGION,MAXSEAS,INITITIALS,DIR,STORMNAME,DEPTH


typhoon_names<-unique(TRACK_DATA$typhoon)


radians_to_degrees<- function(radians){
  degrees  <-  radians *180/pi
  return(degrees)
}

degrees_to_radians<- function(degrees){
  radians <- degrees * pi / 180
  return(radians)
}

create_full_track <- function(hurr_track, tint = 0.5){  
  
  hurr_track <- hurr_track %>%
    mutate(
      index                          = 1:nrow(hurr_track),
      date                           = lubridate::ymd_hm(YYYYMMDDHH),
      tclat                           = abs(as.numeric(LAT)),
      tclon                           = as.numeric(LON),
      tclon                           = ifelse(tclon < 180, tclon, tclon - 180),
      latitude                       = as.numeric(LAT),#as.numeric(LAT),#sprintf("%03d", LAT)#,#as.numeric(as.character(TRACK_DATA$LAT)),
      longitude                      = as.numeric(LON),#as.numeric(LON),# as.numeric(as.character(TRACK_DATA$LON)),
      vmax                           = as.numeric(VMAX)*0.514444,#as.numeric(VMAX),# as.numeric(as.character(TRACK_DATA$VMAX)),
      typhoon_name                   = tolower(STORMNAME),
      wind                           = as.numeric(VMAX)*0.514444
    ) %>%
    select (date                                    ,
            latitude                                ,
            tclat                                   ,
            tclon                                   ,
            longitude                               ,
            typhoon_name                            ,
            wind                                    ,
            typhoon                                 ,
            vmax      ) 
  
  interp_df <- (floor(nrow(hurr_track)/3)-1)
  #interp_df<- 34 # 108 --30    
  interp_date <- seq(from = min(hurr_track$date),
                     to = max(hurr_track$date),
                     by = tint * 3600) # Date time sequence must use `by` in
  # seconds
  interp_date <- data.frame(date = interp_date)
  
  tclat_spline <- stats::glm(tclat ~ splines::ns(date, df = interp_df),  data = hurr_track)
  interp_tclat <- stats::predict.glm(tclat_spline, newdata = interp_date)  
  tclon_spline <- stats::glm(tclon ~ splines::ns(date, df = interp_df),  data = hurr_track)
  interp_tclon <- stats::predict.glm(tclon_spline, newdata = interp_date)
  
  vmax_spline <- stats::glm(vmax ~ splines::ns(date, df = interp_df),    data = hurr_track)
  interp_vmax <- stats::predict.glm(vmax_spline, newdata = interp_date)
  typhoon_name <- as.vector(TRACK_DATA1$STORMNAME[1])
  full_track <- data.frame(typhoon_name=typhoon_name,date = interp_date, tclat = interp_tclat, tclon = interp_tclon, vmax = interp_vmax)
  return(full_track)
}



latlon_to_km<- function(tclat_1, tclon_1, tclat_2, tclon_2, Rearth = 6378.14){
  tclat_1 <- degrees_to_radians(tclat_1)
  tclon_1 <- degrees_to_radians(tclon_1)
  tclat_2 <- degrees_to_radians(tclat_2)
  tclon_2 <- degrees_to_radians(tclon_2)
  
  delta_L <- tclon_1 - tclon_2
  delta_tclat <- tclat_1 - tclat_2
  
  hav_L <- sin(delta_L / 2) ^ 2
  hav_tclat <- sin(delta_tclat / 2) ^ 2
  
  hav_gamma <- hav_tclat + cos(tclat_1) * cos(tclat_2) * hav_L
  gamma <- 2 * asin(sqrt(hav_gamma))
  
  dist <- Rearth * gamma
  return(dist)
}

calc_forward_speed<- function(tclat_1, tclon_1, time_1, tclat_2, tclon_2, time_2){
  dist <- latlon_to_km(tclat_1, tclon_1, tclat_2, tclon_2) * 1000
  time <- as.numeric(difftime(time_2, time_1, units = "secs"))
  tcspd <- dist / time
  return(tcspd)
}

##Calculate the direction of gradient winds at each location

calc_bearing<- function(tclat_1, tclon_1, tclat_2, tclon_2){
  tclat_1 <- degrees_to_radians(tclat_1)
  tclon_1 <- degrees_to_radians(tclon_1)
  tclat_2 <- degrees_to_radians(tclat_2)
  tclon_2 <- degrees_to_radians(tclon_2)
  
 
  
  S <- cos(tclat_2) * sin(tclon_1 - tclon_2)
  C <- cos(tclat_1) * sin(tclat_2) - sin(tclat_1) * cos(tclat_2) * cos(tclon_1 - tclon_2)
  
  theta_rad <- atan2(S, C)
  theta <- radians_to_degrees(theta_rad) + 90
  theta <- theta %% 360 # restrict to be between 0 and 360 degrees
  return(theta)
}


remove_forward_speed<- function(vmax, tcspd){
  vmax_sfc_sym <- vmax - 0.5 * tcspd
  vmax_sfc_sym[vmax_sfc_sym < 0] <- 0
  return(vmax_sfc_sym)
}

calc_gradient_speed<- function(vmax_sfc_sym, over_land){
  reduction_factor <- 0.9
  if(over_land){
    reduction_factor <- reduction_factor * 0.8
  }
  vmax_gl <- vmax_sfc_sym / reduction_factor
  return(vmax_gl)
}


check_over_land<- function(tclat, tclon){
  lat_diffs <- abs(tclat - landmask$latitude)
  closest_grid_lat <- landmask$latitude[which(lat_diffs ==
                                                min(lat_diffs))][1]
  
  lon_diffs <- abs(tclon - (360 - landmask$longitude))
  closest_grid_lon <- landmask$longitude[which(lon_diffs ==
                                                 min(lon_diffs))][1]
  
  over_land <- landmask %>%
    dplyr::filter_(~ latitude == closest_grid_lat &
                     longitude == closest_grid_lon) %>%
    dplyr::mutate_(land = ~ land == "land") %>%
    dplyr::select_(~ land)
  over_land <- as.vector(over_land$land[1])
  
  return(over_land)
}


#Rmax: Radius from the storm center to the point at which the maximum wind occurs (km)

will7a<- function(vmax_gl, tclat){
  Rmax <- 46.4 * exp(-0.0155 * vmax_gl + 0.0169 * tclat)
  return(Rmax)
}

# X1, which is a parameter needed for the Willoughby model. This is done using equation 10a from Willoughby et al. (2006):
# X1=317.1a^2.026Vmax,G+1.915Ï

will10a<-function(vmax_gl, tclat){
  X1 <- 317.1 - 2.026 * vmax_gl + 1.915 * tclat
  return(X1)
}
#
#Next, the code calculates another Willoughby parameter, n. This is done with equation 10b from Willoughby et al. (2006):
#n=0.406+0.0144Vmax,Gâ 0.0038Ï 

will10b<- function(vmax_gl, tclat){
  n <- 0.4067 + 0.0144 * vmax_gl - 0.0038 * tclat
  return(n)
}
#Next, the code calculates another Willoughby parameter, A, with equation 10c from Willoughby et al. (2006)
will10c<- function(vmax_gl, tclat){
  A <- 0.0696 + 0.0049 * vmax_gl - 0.0064 * tclat
  A[A < 0 & !is.na(A)] <- 0
  return(A)
}

will3_right<- function(n, A, X1, Rmax){
  eq3_right <- (n * ((1 - A) * X1 + 25 * A)) /
    (n * ((1 - A) * X1 + 25 * A) + Rmax)
  return(eq3_right)
}


will3_deriv_func<- function(xi, eq3_right){
  deriv <- 70 * 9 * xi ^ 8 - 315 * 8 * xi ^ 7 + 540 * 7 * xi ^ 6 -
    420 * 6 * xi ^ 5 + 126 * 5 * xi ^ 4
  func <- 70 * xi ^ 9 - 315 * xi ^ 8 + 540 * xi ^ 7 - 420 * xi ^ 6 +
    126 * xi ^ 5 - eq3_right
  deriv_func <-c(deriv, func)
  return(deriv_func)
}

solve_for_xi<- function(xi0 = 0.5, eq3_right, eps = 10e-4, itmax = 100){
  if(is.na(eq3_right)){
    return(NA)
  } else{
    i <- 1
    xi <- xi0
    while(i <= itmax){
      deriv_func <- will3_deriv_func(xi, eq3_right)
      if(abs(deriv_func[2]) <= eps){ break }
      xi <- xi - deriv_func[2] / deriv_func[1]
    }
    if(i < itmax){
      return(xi)
    } else{
      warning("Newton-Raphson did not converge.")
      return(NA)
    }
  }
}


##While the Newton-Raphson method can sometimes perform poorly in finding global maxima, in this case the function for which we are trying ##to find the root is well-behaved. Across tropical storms from 1988 to 2015, the method never failed to converge, and identified roots ##were consistent across storms (typically roots are for Î¾ of 0.6--0.65). We also tested using the optim function in the stats R package ##and found similar estimated roots but slower convergence times than when using the Newton-Raphson method.

##Now an equation from the Willoughby et al. 2006 paper can be used to calculate R1 (Willoughby, Darling, and Rahn 2006):

#R1â=â Rmaxâ â â¾(R2â â â R1)

#For this function, the package code assumes that R2âââR1 (the width of the transition region) is 25 kilometers when Rmax is larger than 20 #kilometers and 15 kilometers otherwise.

calc_R1<- function(Rmax, xi){
  R2_minus_R1 <- ifelse(Rmax > 20, 25, 15)
  R1 <- Rmax - xi * R2_minus_R1
  return(R1)
}


##Determine radius for end of transition region
##Next, the code estimates the radius for the end of the transition period, R2. We assume that smaller storms (Rmaxââ¤â20 km) have a ##transition region of 15 km while larger storms (Rmaxâ>â20 km) have a transition region of 25 km:

##$$ R_2 = \begin{cases} R_1 + 25 & \text{ if } R_{max} > 20\mbox{ km}\\ R_1 + 15 & \text{ if } R_{max} \le 20\mbox{ km} \end{cases} $$

##where:

##R1: Radius to the start of the transition region (km)
##R2: Radius to the end of the transition region (km)
##R2 = ifelse(Rmax > 20, R1 + 25, R1 + 15)
##Calculate wind speed at each grid point
##Next, the code models the wind speed at a location (e.g., a county center). As a note, this function calculates wind characteristics at a ##single location; a later function applies this function across many grid points):

##After calculating the grid wind time series for a grid point, you can input the time series for a grid point into summarize_grid_wind to ##generate overall storm summaries for the grid point. This functions calculate wind characteristics at each grid point (or county center ##location) for every storm observation. These characteristics are:

##vmax_gust: Maximum value of surface-level (10 meters) sustained winds, in meters per second, over the length of the storm at the given ##location
##vmax_sust: Maximum value of surface-level (10 meters) gust winds, in meters per second, over the length of the storm at the given location
##gust_dur: Length of time, in minutes, that surface-level sustained winds were above a certain wind speed cutoff (e.g., 20 meters per ##second)
##sust_dur: Length of time, in minutes, that surface-level gust winds were above a certain wind speed cutoff (e.g., 20 meters per second)

##Determine gradient wind speed at each location
##Next, the package calculates VG(r), the gradient level 1-minute sustained wind at the grid point, which is at radius r from the tropical ##cyclone center (Cdis**t for the grid point). Note there are different equations for VG(r) for (1) the eye to the start of the transition ##region; (2) outside the transition region; and (3) within the transition region.


will1<- function(cdist, Rmax, R1, R2, vmax_gl, n, A, X1, X2 = 25){
  
  if(is.na(Rmax) || is.na(vmax_gl) ||
     is.na(n) || is.na(A) || is.na(X1)){
    return(NA)
  } else {
    
    Vi <- vmax_gl * (cdist / Rmax) ^ n
    Vo <- vmax_gl * ((1 - A) * exp((Rmax - cdist)/X1) + A * exp((Rmax - cdist) / X2))
    
    if(cdist < R1){
      wind_gl_aa <- Vi
    } else if (cdist > R2){
      wind_gl_aa <- Vo
    } else {
      eps <- (cdist - R1) / (R2 - R1)
      w <- 126 * eps ^ 5 - 420 * eps ^ 6 + 540 * eps ^ 7 - 315 *
        eps ^ 8 + 70 * eps ^ 9
      wind_gl_aa <- Vi * (1 - w) + Vo * w
    }
    
    wind_gl_aa[wind_gl_aa < 0 & !is.na(wind_gl_aa)] <- 0
    
    return(wind_gl_aa)
  }
}

##Estimate surface level winds from gradient winds

gradient_to_surface<- function(wind_gl_aa, cdist){
  if(cdist <= 100){
    reduction_factor <- 0.9
  } else if(cdist >= 700){
    reduction_factor <- 0.75
  } else {
    reduction_factor <- 0.90 - (cdist - 100) * (0.15/ 600)
  }
  # Since all counties are over land, reduction factor should
  # be 20% lower than if it were over water
  reduction_factor <- reduction_factor * 0.8
  wind_sfc_sym <- wind_gl_aa * reduction_factor
  return(wind_sfc_sym)
}


##Next, the function calculates the gradient wind direction based on the bearing of a location from the storm. This gradient wind direction is##calculated by adding 90 degrees to the bearing of the grid point from the storm center.

#gwd = (90 + chead) %% 360

##Calculate the surface wind direction
##The next step is to change from the gradient wind direction to the surface wind direction. To do this, the function adds an inflow angle ##to the gradient wind direction (making sure the final answer is between 0 and 360 degrees). This step is necessary because surface ##friction changes the wind direction near the surface compared to higher above the surface.

##The inflow angle is calculated as 

add_inflow<- function(gwd, cdist, Rmax){
  if(is.na(gwd) | is.na(cdist) | is.na(Rmax)){
    return(NA)
  }
  
  # Calculate inflow angle over water based on radius of location from storm
  # center in comparison to radius of maximum winds (Phadke et al. 2003)
  if(cdist < Rmax){
    inflow_angle <- 10 + (1 + (cdist / Rmax))
  } else if(Rmax <= cdist & cdist < 1.2 * Rmax){
    inflow_angle <- 20 + 25 * ((cdist / Rmax) - 1)
  } else {
    inflow_angle <- 25
  }
  
  # Add 20 degrees to inflow angle since location is over land, not water
  overland_inflow_angle <- inflow_angle + 20
  
  # Add inflow angle to gradient wind direction
  gwd_with_inflow <- (gwd + overland_inflow_angle) %% 360
  
  return(gwd_with_inflow)
}


#Add back in wind component from forward speed of storm Next, to add back in the storm's forward motion at each grid point, the code #reverses the earlier step that used the Phadke correction factor (equation 12, Phadke et al. 2003). The package calculates a constant #correction factor (correction_factor), as a function of r, radius from the storm center to the grid point, and Rmax, radius from storm #center to maximum winds.

#$$ U(r) = \frac{R_{max}r}{R_{max}^2 + r^2}F $$ where:


add_forward_speed <- function(wind_sfc_sym, tcspd_u, tcspd_v, swd, cdist, Rmax){
  # Calculate u- and v-components of surface wind speed
  swd<- swd * pi / 180
  wind_sfc_sym_u <- wind_sfc_sym * cos(swd)
  wind_sfc_sym_v <-  wind_sfc_sym * sin(swd)
  
  # Add back in component from forward motion of the storm
  correction_factor <- (Rmax * cdist) / (Rmax ^ 2 + cdist ^ 2)
  
  # Add tangential and forward speed components and calculate
  # magnitude of this total wind
  wind_sfc_u <- wind_sfc_sym_u + correction_factor * tcspd_u
  wind_sfc_v <- wind_sfc_sym_v + correction_factor * tcspd_v
  wind_sfc <- sqrt(wind_sfc_u ^ 2 + wind_sfc_v ^ 2)
  
  # Reset any negative values to 0
  wind_sfc <- ifelse(wind_sfc > 0, wind_sfc, 0)
  
  return(wind_sfc)
}


#Calculate 3-second gust wind speed from sustained wind speed

#Here is a table with gust factors based on location (Harper, Kepert, and Ginger 2010):

#Location	Gust factor (G3,â 60)
#In-land	1.49
#Just offshore	1.36
#Just onshore	1.23
#At sea	1.11
#The stormwindmodel package uses the "in-land" gust factor value throughout.
 
get_grid_winds<- function(hurr_track , grid_df ,tint = 0.5,gust_duration_cut = 20,sust_duration_cut = 20){
  full_track <- create_full_track(hurr_track = hurr_track, tint = tint)
  with_wind_radii <- add_wind_radii(full_track = full_track)
  
  grid_winds <- plyr::adply(grid_df, 1, calc_and_summarize_grid_wind,
                            with_wind_radii = with_wind_radii,
                            tint = tint,
                            gust_duration_cut = gust_duration_cut,
                            sust_duration_cut = sust_duration_cut)
  
  return(grid_winds)
}


add_wind_radii <-function(full_track = create_full_track()){
  
  with_wind_radii <- full_track %>%
    dplyr::mutate_(tcspd = ~ calc_forward_speed(tclat, tclon, date,
                                                dplyr::lead(tclat),
                                                dplyr::lead(tclon),
                                                dplyr::lead(date)),
                   tcdir = ~ calc_bearing(tclat, tclon,
                                          dplyr::lead(tclat),
                                          dplyr::lead(tclon)),
                   tcspd_u = ~ tcspd * cos(tcdir* pi / 180),
                   tcspd_v = ~ tcspd * sin(tcdir* pi / 180),
                   vmax_sfc_sym = ~ remove_forward_speed(vmax, tcspd),
                   over_land = ~ mapply(check_over_land, tclat, tclon),
                   vmax_gl = ~ mapply(calc_gradient_speed,
                                      vmax_sfc_sym = vmax_sfc_sym,
                                      over_land = over_land),
                   Rmax = ~ will7a(vmax_gl, tclat),
                   X1 = ~ will10a(vmax_gl, tclat),
                   n = ~ will10b(vmax_gl, tclat),
                   A = ~ will10c(vmax_gl, tclat),
                   eq3_right = ~ will3_right(n, A, X1, Rmax),
                   xi = ~ mapply(solve_for_xi, eq3_right = eq3_right),
                   R1 = ~ calc_R1(Rmax, xi),
                   R2 = ~ ifelse(Rmax > 20, R1 + 25, R1 + 15)
    ) %>%
    dplyr::select_(quote(-vmax), quote(-tcspd), quote(-vmax_sfc_sym),
                   quote(-over_land), quote(-eq3_right), quote(-xi))
  return(with_wind_radii)
}

calc_grid_wind <- function(grid_point,
                           with_wind_radii = add_wind_radii()){
  
  grid_wind <- dplyr::mutate_(with_wind_radii,
                              # Calculated distance from storm center to location
                              cdist = ~ latlon_to_km(tclat, tclon,
                                                     grid_point$glat, grid_point$glon),
                              # Calculate gradient winds at the point
                              wind_gl_aa = ~ mapply(will1, cdist = cdist, Rmax = Rmax,
                                                    R1 = R1, R2 = R2, vmax_gl = vmax_gl,
                                                    n = n, A = A, X1 = X1),
                              # calculate the gradient wind direction (gwd) at this
                              # grid point
                              chead = ~ calc_bearing(tclat, tclon,
                                                     grid_point$glat, - grid_point$glon),
                              gwd = ~ (90 + chead) %% 360,
                              # Bring back to surface level (surface wind reduction factor)
                              wind_sfc_sym = ~ mapply(gradient_to_surface,
                                                      wind_gl_aa = wind_gl_aa,
                                                      cdist = cdist),
                              # Get surface wind direction
                              swd = ~ mapply(add_inflow, gwd = gwd, cdist = cdist,
                                             Rmax = Rmax),
                              # Add back in storm forward motion component
                              windspeed = ~ add_forward_speed(wind_sfc_sym,
                                                              tcspd_u, tcspd_v,
                                                              swd, cdist, Rmax)) %>%
    dplyr::select_(~ date, ~ windspeed,~ cdist)
  return(grid_wind)
}

get_grid_winds <- function(hurr_track,
                           grid_df,
                           tint = 0.5,
                           gust_duration_cut = 20,
                           sust_duration_cut = 20){
  full_track <- create_full_track(hurr_track = hurr_track, tint = tint)
  with_wind_radii <- add_wind_radii(full_track = full_track)
  
  grid_winds <- plyr::adply(grid_df, 1, calc_and_summarize_grid_wind,
                            with_wind_radii = with_wind_radii,
                            tint = tint,
                            gust_duration_cut = gust_duration_cut,
                            sust_duration_cut = sust_duration_cut)
  
  return(grid_winds)
}


summarize_grid_wind <- function(grid_wind, tint = 0.5, gust_duration_cut = 20,
                                sust_duration_cut = 20){
  grid_wind_summary <- grid_wind %>%
    dplyr::mutate_(gustspeed = ~ windspeed * 1.49) %>%
    # Determine max of windspeed and duration of wind over 20
    dplyr::summarize_(vmax_gust = ~ max(gustspeed, na.rm = TRUE),
                      vmax_sust = ~ max(windspeed, na.rm = TRUE),
                      dist_track = ~ min(cdist, na.rm = TRUE),
                      gust_dur = ~ 60 * sum(gustspeed > gust_duration_cut,
                                            na.rm = TRUE),
                      sust_dur = ~ 60 * sum(windspeed > sust_duration_cut,
                                            na.rm = TRUE)) %>%
    dplyr::mutate_(gust_dur = ~ gust_dur * tint,
                   sust_dur = ~ sust_dur * tint)
  grid_wind_summary <- as.matrix(grid_wind_summary)
  return(grid_wind_summary)
}

calc_and_summarize_grid_wind <- function(grid_point,
                                         with_wind_radii = add_wind_radii(),
                                         tint = 0.5, gust_duration_cut = 20,
                                         sust_duration_cut = 20){
  grid_wind <- calc_grid_wind(grid_point = grid_point,
                              with_wind_radii = with_wind_radii)
  grid_wind_summary <- summarize_grid_wind(grid_wind = grid_wind, tint = tint,
                                           gust_duration_cut = gust_duration_cut,
                                           sust_duration_cut = sust_duration_cut)
  return(grid_wind_summary)
  
}


map_wind <- function(grid_winds, value = "vmax_sust", break_point = NULL,
                     wind_metric = "mps"){
  
  grid_winds$value <- grid_winds[ , value]
  if(wind_metric != "mps"){
    grid_winds$value <- weathermetrics::convert_wind_speed(grid_winds$value,
                                                           old_metric = "mps",
                                                           new_metric = wind_metric)
  }
  
  if(!is.null(break_point)){
    cut_values <- cut(grid_winds$value,
                      breaks = c(0, break_point, max(grid_winds$value)),
                      include.lowest = TRUE)
    grid_winds$value <- cut_values
    num_colors <- 2
  } else {
    if(wind_metric == "mps"){
      breaks <- c(0, seq(15, 45, 5))
      exposure_palette <- c("#FEE5D9", "#FCBBA1", "#FC9272", "#FB6A4A",
                            "#DE2D26", "#A50F15")
    } else if(wind_metric == "knots"){
      breaks <- c(0, 34, 50, 64, 100)
      exposure_palette <- c("#FEE0D2", "#FC9272", "#DE2D26")
    }
    palette_name <- "Reds"
    
    # Adjust for right outliers
    if(max(grid_winds$value) > max(breaks)){
      breaks <- c(breaks, max(grid_winds$value))
    }
    
    exposure_palette <- c("#ffffff", exposure_palette, "#1a1a1a")
    
    grid_winds <- grid_winds %>%
      dplyr::mutate_(value = ~ cut(value, breaks = breaks,
                                   include.lowest = TRUE))
  }
  
  map_data <- dplyr::mutate_(grid_winds,
                             fips = ~ fips) %>%
    dplyr::select_(~ fips, ~ value)
  us_counties<- php_admin3 %>%
    dplyr::left_join(map_data, by = "fips")
  
  us_counties2<- php_admin3[php_admin3$adm3_pcode %in% list("PH175302000","PH175307000","PH175309000","PH175322000","PH175310000","PH175313000","PH060401000",
                                                            "PH060402000","PH060403000","PH060404000","PH060405000","PH060406000","PH060407000","PH060408000",
                                                            "PH060409000","PH060410000","PH060412000","PH060411000","PH060413000","PH060414000","PH060415000","PH060416000",
                                                            "PH060417000","PH060602000","PH060604000","PH060605000","PH060606000","PH060609000","PH060610000","PH060611000",
                                                            "PH060612000","PH060614000","PH060615000","PH060617000","PH060618000","PH061901000","PH061902000","PH061903000",
                                                            "PH061904000","PH061905000","PH061906000","PH061907000","PH061908000","PH061909000","PH061910000","PH061911000",
                                                            "PH061912000","PH061913000","PH061914000","PH061915000","PH061916000","PH061917000","PH063001000","PH063002000",
                                                            "PH063003000","PH063004000","PH063005000","PH063006000","PH063007000","PH063008000","PH063009000","PH063010000",
                                                            "PH063012000","PH063013000","PH063014000","PH063015000","PH063016000","PH063017000","PH063018000","PH063019000",
                                                            "PH063023000","PH063025000","PH063027000","PH063029000","PH063031000","PH063032000","PH063035000","PH063037000",
                                                            "PH063038000","PH063039000","PH063042000","PH063044000","PH063047000","PH064504000","PH064508000","PH064509000",
                                                            "PH064518000","PH064523000","PH064526000","PH064531000","PH072209000","PH072211000","PH072213000","PH072221000",
                                                            "PH072228000","PH072231000","PH072236000","PH072238000","PH072242000","PH072243000","PH072247000","PH072244000",
                                                            "PH072248000","PH072249000","PH072252000","PH072253000","PH083701000","PH083702000","PH083703000","PH083705000",
                                                            "PH083706000","PH083710000","PH083713000","PH083714000","PH083715000","PH083717000","PH083718000","PH083722000",
                                                            "PH083723000","PH083724000","PH083725000","PH083726000","PH083728000","PH083729000","PH083730000","PH083731000",
                                                            "PH083733000","PH083735000","PH083736000","PH083739000","PH083740000","PH083741000","PH083742000","PH083743000",
                                                            "PH083744000","PH083745000","PH083746000","PH083748000","PH083749000","PH083750000","PH083751000","PH083708000",
                                                            "PH083738000","PH083747000","PH082602000","PH082603000","PH082607000","PH082608000","PH082609000","PH082610000",
                                                            "PH082612000","PH082613000","PH082615000","PH082616000","PH082618000","PH082619000","PH086002000","PH086006000",
                                                            "PH086010000","PH086017000","PH086019000","PH086021000","PH086416000","PH087801000","PH087802000","PH087803000",
                                                            "PH087804000","PH087805000","PH087808000","PH168505000","PH054103000"),]%>%
    dplyr::left_join(map_data, by = "fips")
  
  out <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = us_counties,
                          ggplot2::aes_(x = ~ long, y = ~ lat, group = ~ group,
                                        fill = ~ value),
                          color = "lightgray", size = 0.2) +
    ggplot2::theme_void() +
    ggplot2::coord_map()
  
  exposure_legend <- paste0("Wind speed (",
                            ifelse(wind_metric == "mps", "m / s",
                                   wind_metric),
                            ")")
  
  if(!is.null(break_point)){
    out <- out + ggplot2::scale_fill_manual(name = exposure_legend,
                                            values = c("white", "#DE2D26"),
                                            labels = levels(cut_values))
  } else{
    out <- out + ggplot2::scale_fill_manual(name = exposure_legend,
                                            values = exposure_palette)
  }
  
  return(us_counties)
}



for (nam in typhoon_names)
{  TRACK_DATA1<-TRACK_DATA[TRACK_DATA$typhoon == nam, ]
  
  
  wind_grids <- get_grid_winds(hurr_track=TRACK_DATA1,
                               grid_df=grid_points_adm3,
                               tint = 0.5,                 gust_duration_cut = 20,
                               sust_duration_cut = 20)
  
  wind_grids_ <- data.frame(typhoon_name=nam,gridid = wind_grids$gridid, 
                            glat = wind_grids$glat, glon = wind_grids$glon, 
                            vmax_gust = wind_grids$vmax_gust,
                            vmax_sust = wind_grids$vmax_sust,
                            vmax_sust = wind_grids$vmax_sust,
                            dist_track=wind_grids$dist_track)
  
  if (exists("wind_full_grid")){
    wind_full_grid<-merge(wind_full_grid,y=wind_grids_,all=TRUE)
  } else{
    wind_full_grid<-wind_grids_
  }  
}

write.table(wind_full_grid, file = "C:/SOFTWARES/stormwindmodel/wind_full_grid.csv",row.names=FALSE, na="", sep=",")


#HAIYAN from best record
TRACK_DATA1<-TRACK_DATA[TRACK_DATA$typhoon =='HAIYAN2013', ]

fulltrack<-create_full_track(hurr_track=TRACK_DATA1, tint = 0.5)

Mydata<-fulltrack %>%
  dplyr::mutate(index= 1:nrow(fulltrack))%>%
  dplyr::select(tclat,tclon,index,typhoon_name,vmax,date)

my_track <- st_as_sf(Mydata, coords = c("tclon", "tclat"), crs = "+init=epsg:4326")


wind_grids <- get_grid_winds(hurr_track=TRACK_DATA1,
                             grid_df=grid_points_adm3,
                             tint = 0.5,                 gust_duration_cut = 20,
                             sust_duration_cut = 20)

wind_grids<- wind_grids %>%
  dplyr::mutate(fips = as.numeric(substr(gridid, 3, 11))) %>%
  dplyr::select(fips,vmax_gust,vmax_sust,gust_dur,sust_dur,glon,glat)


php_admin3 <- st_read(dsn='C:\\SOFTWARES\\stormwindmodel\\data-raw',layer='phl_admin3_simpl2')



php_admin3<-php_admin3 %>%
  dplyr::mutate(long = glon,lat=glat,group=as.numeric(substr(adm2_pcode, 3, 6)),fips=as.numeric(substr(adm3_pcode, 3, 11)))%>%
  dplyr::select(adm3_en,adm3_pcode,adm2_pcode,fips,long,lat,group,geometry)



php_data2<-php_admin3 %>%
  left_join(wind_grids, by = "fips")


my_track<-my_track %>%
  dplyr::mutate(date2 =format(my_track$date,"%B %d"))


# tmap
tmap_mode(mode = "view")
tm_shape(php_data) + tm_polygons(col = "vmax_gust", border.col = "black",lwd = 0.1,lyt='dotted')+
  tm_shape(my_track) + tm_symbols(size=0.01,border.alpha = 1,col="blue") +
  tm_shape(my_track_pag) + tm_symbols(size=0.01,border.alpha = 1,col="red") +
  #tm_text("date2", col="date2",legend.col.show = F,clustering = TRUE) +
  tm_format("NLD")






##HAIYAN pagassa forecast
TRACK_DATA3 <- read.csv("C:/Users/ATeklesadik/Rode Kruis/510 - Data preparedness - Documents/[PRJ] FbF - Philippines - (PMF)/Analysis/historical_best_track_metoc_navy/best_pagasa.csv")
TRACK_DATA1<-TRACK_DATA3[TRACK_DATA3$STORMNAME =='HAIYAN', ]

fulltrack<-create_full_track(hurr_track=TRACK_DATA1, tint = 0.5)

Mydata<-fulltrack %>%
  dplyr::mutate(index= 1:nrow(fulltrack))%>%
  dplyr::select(tclat,tclon,index,typhoon_name,vmax,date)

my_track_pag <- st_as_sf(Mydata, coords = c("tclon", "tclat"), crs = "+init=epsg:4326")
wind_grids <- get_grid_winds(hurr_track=TRACK_DATA1,grid_df=grid_points_adm3,tint = 0.5, gust_duration_cut = 20,sust_duration_cut = 20)
wind_grids<- wind_grids %>%
  dplyr::mutate(fips = as.numeric(substr(gridid, 3, 11))) %>%
  dplyr::select(fips,vmax_gust,vmax_sust,gust_dur,sust_dur,glon,glat)


php_admin3 <- st_read(dsn='C:\\SOFTWARES\\stormwindmodel\\data-raw',layer='phl_admin3_simpl2')

php_admin3<-php_admin3 %>%
  dplyr::mutate(long = glon,lat=glat,group=as.numeric(substr(adm2_pcode, 3, 6)),fips=as.numeric(substr(adm3_pcode, 3, 11)))%>%
  dplyr::select(adm3_en,adm3_pcode,adm2_pcode,fips,long,lat,group,geometry)



php_data_pag<-php_admin3 %>%
  left_join(wind_grids, by = "fips")


my_track_pag<-my_track %>%
  dplyr::mutate(date2 =format(my_track$date,"%B %d"))


# tmap
tmap_mode(mode = "view")
tm_shape(php_data) + tm_polygons(col = "vmax_gust", border.col = "black",lwd = 0.1,lyt='dotted')+
  tm_shape(my_track) + tm_symbols(size=0.01,border.alpha = 1,col=) +
  #tm_text("date2", col="date2",legend.col.show = F,clustering = TRUE) +
  tm_format("NLD")




library(rgdal)
library(sp)
library(tmap)




for (nam in typhoon_names){  
  
  TRACK_DATA1<-TRACK_DATA[TRACK_DATA$typhoon == nam, ]
  fulltrack<-create_full_track(hurr_track=TRACK_DATA1, tint = 0.5)
  
  Mydata<-fulltrack %>%
    dplyr::mutate(index= 1:nrow(fulltrack))%>%
    dplyr::select(tclat,tclon,index,typhoon_name,vmax,date)
   
  #my_track <- st_as_sf(Mydata, coords = c("tclon", "tclat"), crs = "+init=epsg:4326")
  
  if (exists("track_full_grid")){
    track_full_grid<-merge(track_full_grid,y=Mydata,all=TRUE)
  } else{
    track_full_grid<-Mydata
}  


}

my_track <- st_as_sf(track_full_grid, coords = c("tclon", "tclat"), crs = "+init=epsg:4326")

#write.table(track_full_grid, file = "C:/SOFTWARES/stormwindmodel/track_full_grid.csv",row.names=FALSE, na="", sep=",")

st_write(my_track, dsn = "C:/SOFTWARES/stormwindmodel/full_track.shp", layer='full_track')

#READ OLD DATAMATRIX 
Datamatrix_typhoon<-read.csv('C:/Typhoons/model/Datamatrix_typhoon.csv',sep=';')

#READ NEW VARIABLES FOR STORME SURGE AND LANDSLID
data_matrix_stormsurge_landslide<-read.csv('C:/Typhoons/model/Datamatrix_storm_ll.csv',sep=',')%>%
  dplyr::select(Mun_Code,landslide,stormsurge)

data_matrix_new_variables<-merge(x=Datamatrix_typhoon, y=data_matrix_stormsurge_landslide, by = "Mun_Code",all.x=TRUE)

wind_full_grid<-wind_full_grid %>%
  dplyr::mutate(typoon_name2= tolower(stringr::str_sub(typhoon_name, 1, -5)))



data_matrix_new_variables<- data_matrix_new_variables %>%
  dplyr::mutate(typhoon2 = paste0(tolower(data_matrix_new_variables$typhoon_name),substr(data_matrix_new_variables$Mun_Code, 1, 14)))


wind_full_grid<- wind_full_grid %>%
  dplyr::mutate(typhoon2 = paste0(tolower(wind_full_grid$typoon_name2),substr(wind_full_grid$gridid, 1, 14)))

wind_full_grid<-wind_full_grid[wind_full_grid$vmax_gust>10,]

xx<-merge(x = data_matrix_new_variables, y = wind_full_grid, by = "typhoon2",all.x=TRUE)


datamatrix_new<- xx %>%
  dplyr::mutate(averge_speed_mhp2 = xx$vmax_gust*2.23694)
  
  
## WRITE THE VARIABLES TO A NEW MATRIX 
write.csv(datamatrix_new,'C:/Typhoons/model/datamatrix_new2.csv',row.names = FALSE)



