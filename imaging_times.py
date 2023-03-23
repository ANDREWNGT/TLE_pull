# Author: Andrew Ng 
# Date: 08/02/2023
# Script to compute imaging time using the most recent TLE for some time in the future. 
# Precision of imaging time is 1 second.
# And compute imaging time for SAR imagery on a provided spot
# Inputs: 1. Coordinates of spot as tuple 
#         2. Satellite norad cat id

import TLE_pull 
from datetime import date
from datetime import datetime, timedelta
import os
import math
import numpy as np
from astropy.time import Time
import pandas as pd
import matplotlib.pyplot as plt
from skyfield.api import load, wgs84, EarthSatellite
global earth_rotation_rate, mu
earth_rotation_rate = 7.2921150E-5 # rad/s
mu = 398600.44189 #km^3/s^2
plt.close()
def propagate_SGP4(cat_id, coords, sat_name, ts, t0, t1):

    date_format_today = "%Y%m%d"
    today = date.today()
    data_pulled_day = today.strftime(date_format_today)
    parent_output_folder = os.path.join(os.getcwd(), "pulled_data")
    if not os.path.exists(parent_output_folder):
        os.mkdir(parent_output_folder)
    output_folder = os.path.join(parent_output_folder, f"data_{data_pulled_day}")
    
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    output_file_name = TLE_pull.check_tle(cat_id, output_folder= output_folder)
    TLE_file = open(output_file_name, 'r')
    Lines = TLE_file.readlines()

    TLE_year = 2000 + int(Lines[1][18:20])
    TLE_year = datetime(int(TLE_year), 1, 1, 0, 0)
    TLE_fractional_day = float(Lines[1][20:32])
    TLE_hour = math.floor(TLE_fractional_day%1*24)
    TLE_minute = math.floor((TLE_fractional_day%1*24)%1*60)
    datetime_TLE_epoch = TLE_year + timedelta(days=int(TLE_fractional_day)-1) \
    + timedelta(hours= TLE_hour) \
    + timedelta(minutes= TLE_minute)

    

    satellite = EarthSatellite(Lines[1][:-1], Lines[2], sat_name, ts)
    spot = wgs84.latlon(coords[0], coords[1])
    number_of_seconds = (t1.utc_datetime()-t0.utc_datetime()).total_seconds()
    t_list = ts.linspace(t0, t1, num = int(number_of_seconds))
    satellite_vector = satellite.at(t_list)
    lat, lon = wgs84.latlon_of(satellite_vector)
    difference = satellite - spot 
    topocentric = difference.at(t_list)
    #1. compute euclidean distance
    #distance_between_satelite_to_spot = np.linalg.norm(topocentric.position.km , axis =0)
    alt, az, distance = topocentric.altaz()
    index_of_min_dist = np.argmin(distance.km)
    time_at_shortest_dist = t_list[index_of_min_dist].utc_datetime_and_leap_second()[0]
    elevation_angle = alt.degrees[index_of_min_dist]
    azimuth_angle = az.degrees[index_of_min_dist]
    distance_min = distance.km[index_of_min_dist]

    print(f"TLE Epoch: {datetime_TLE_epoch} UTC")
    print(f"Imaging time: {time_at_shortest_dist} UTC")
    print(f"Minimum distance: {distance_min:.3f} km ")
    print(f"Azimuth angle: {azimuth_angle:.3f} degrees")
    print(f"Elevation angle: {elevation_angle:.3f} degrees")

    #2. Match t_list indices to euclidean distance. Convert gps seconds to UTC time. 
    
    figure, ax = plt.subplots(1,1)
    plt.plot(lon.degrees[index_of_min_dist-500:index_of_min_dist+500], lat.degrees[index_of_min_dist-500:index_of_min_dist+500])
    plt.plot(coords[1], coords[0], marker='o')
    ax.set_ylim(np.min(lat.degrees[index_of_min_dist-500:index_of_min_dist+500]), np.max(lat.degrees[index_of_min_dist-500:index_of_min_dist+500]))
    ax.set_xlim(np.min(lon.degrees[index_of_min_dist-500:index_of_min_dist+500]), np.max(lon.degrees[index_of_min_dist-500:index_of_min_dist+500]))
    ax.set_ylabel("Latitude (degrees)")
    ax.set_xlabel("Longitude (degrees)")

    plt.show()

    

    #return output_df
if __name__=="__main__":

    cat_id = 43215
    coords = (1.339798, 103.90346)
    sat_name = "PAZ"
    #t0 must be prior to the expected window. Rule of thumb for leo is to put t0 as up to 45 mins before expected imaging time. 
    ts = load.timescale()
    t0 = ts.utc(2023, 2, 16, 22, 00)
    t1 = ts.utc(2023, 2, 16, 23 ,30)
    propagate_SGP4(cat_id, coords, sat_name, ts, t0, t1)

    