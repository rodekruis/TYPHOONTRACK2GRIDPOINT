R script that uses WILLOUGHBY et al. methedology to calculate windspeed at any gridpoint location in the Phillipiness area of Responsiblity(PAR), based on a typhoon track data. The script include some codes from [stormwindmodel R package](https://CRAN.R-project.org/package=stormwindmodel). This package applies the same methedology to estimate windspeed at any grid location in the US. 

The script is available under the [GPL license](LICENSE)

![alt tag](http://510.global/wp-content/uploads/2015/06/510-opengraph.png)


# Localized windspeed 
This script is developed to calculate localized wind speed and exosure duration froma a typhoon track data with information on maximum wind speed and typhoon center. This tool can be used to calculte wind speed at any grid point location in the Phillipines Responisbility Area (PAR).  In Addtion to maximum wind speed the tool will also calculate duration of exposure to a wind speed above a certain threshold. 

# Typhoon Track data
Typhoon forecast data from National Met agencies come in different formate, the following example shows a typhoon warning bulleten issued by PAGASA for the phillipiness responisiblity area. 
![A forecast information of typhoon HAIYAN issued by PAGASA](figures/pagasa.png)

The relevant information on typhoon is extracted form the text document 

![A forecast data extracted from warning issued by PAGASA(the above text document)](figures/pagasa2.png)

The typhoon track data is not suffceint inormation to estimate the possible impact of typhoon. A wind severity map, which will show the magnitude of wind at different locations is needed as an input for impact based forecasting models. For 510 typhoon model we need data at each grid location (manucipality centers) and this tool can be used to calculate wind speed at manucipality centers.

![A forecast data at a grid point calculated by typhoontrack2grid](figures/pagasa3.png)

In the figure above  
Vmax_gust maximum wind speed
Vmax_sust maximum sustaind speed
dist_track minimum distance from typhoon track
gust_dur total duration of exposure for gust wind
sust_dur total duration of exposure for sustained wind 

# Wind severty mapping
Then we can create a severity map from  the calculated gridpoint windspeed data

![Wind severity map of Haiyan](figures/haiyan.JPG)

# Referances

This tool is based on the the work of "Willoughby, HE, RWR Darling, and ME Rahn. 2006. Parametric Representation of the Primary Hurricane Vortex. Part II: A New Family of Sectionally Continuous Profiles.” Monthly Weather Review 134 (4): 1102–20.



