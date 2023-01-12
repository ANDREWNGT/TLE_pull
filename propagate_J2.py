# Author: Andrew Ng 
# Date: 11/7/2022
# Script to pull data from Celestrak's server https://celestrak.org/NORAD/elements
import TLE_pull 
from datetime import date
from datetime import datetime, timedelta
import os
import math
import numpy as np
from astropy.time import Time
import pandas as pd
import matplotlib.pyplot as plt
from compute_RAAN_drift import calc_J2_rate

global earth_rotation_rate, mu
earth_rotation_rate = 7.2921150E-5 # rad/s
mu = 398600.44189 #km^3/s^2

def propagate_J2(data_propagate_J2_list, cat_id, use_today_date = True, days_selected = 5, time_step = 5):
    # %% User inputs
    '''
    data_propagate_J2_list:  a list of strings, specifying the requested day to display the RAAN for. 
                            Given in the format "%YYYY%MM%DD"

    cat_id: int,  value of the NORAD CAT ID of selected satellite. NEUSAR --> 52937
    time_step: int, time step in minutes
    '''
    if not use_today_date: 
        past_TLE_filename = "pulled_data\\data_20221007\\221007_NEUSAR_52937.tle"
    date_format = "%Y%m%d%H%M"
    # %% Section to pull TLE data from online
    # Example to use NORAD Cat ID: 
    # TLE_pull.check_tle(cat_id = 25544)
    if use_today_date:
        date_format_today = "%Y%m%d"
        today = date.today()
        data_pulled_day = today.strftime(date_format_today)
        
        # Example to use satellite exact name. Note that name must match CelesTrak's exact format. 
        # TLE_pull.check_tle(sat_name = "STARLINK-1007")
        parent_output_folder = os.path.join(os.getcwd(), "pulled_data")
        if not os.path.exists(parent_output_folder):
            os.mkdir(parent_output_folder)
        output_folder = os.path.join(parent_output_folder, f"data_{data_pulled_day}")
        
        if not os.path.exists(output_folder):
            os.mkdir(output_folder)
        output_file_name = TLE_pull.check_tle(cat_id, output_folder= output_folder)
    else:
        output_file_name = past_TLE_filename
    # %% Section to propagate to selected day from TLE pulled
    
    # Read values from TLE. 
    TLE_file = open(output_file_name, 'r')
    Lines = TLE_file.readlines()
    current_mean_motion = float(Lines[2][52:63])
    current_a = ((mu*1E9)**(1/3) / (2*current_mean_motion*math.pi/86400)**(2/3))/1000 #km
    current_e = float(f"0.{Lines[2][26:32]}")
    current_i = float(Lines[2][8:15]) * math.pi/180 # rad
    current_RAAN = float(Lines[2][17:24]) # deg

    TLE_year = 2000 + int(Lines[1][18:20])
    TLE_year = datetime(int(TLE_year), 1, 1, 0, 0)
    TLE_fractional_day = float(Lines[1][20:32])
    TLE_hour = math.floor(TLE_fractional_day%1*24)
    TLE_minute = math.floor((TLE_fractional_day%1*24)%1*60)
    datetime_TLE_epoch = TLE_year + timedelta(days=int(TLE_fractional_day)-1) \
    + timedelta(hours= TLE_hour) \
    + timedelta(minutes= TLE_minute)

    #calcuate nodal precession using perturbed mean motion value. 
    nodal_precession_rate = calc_J2_rate(current_a, current_i, current_e)
    overall_df = pd.DataFrame()
   
    for data_propagate_J2 in data_propagate_J2_list:
        # %% Process date time objects
        #append hour and minute to data_propagate_J2
        if TLE_hour<10:
            TLE_hour_string = f"{0}{TLE_hour}"
        else:
            TLE_hour_string = str(TLE_hour)
        if TLE_minute<10:
            TLE_minute_string = f"{0}{TLE_minute}"
        else: 
            TLE_minute_string = str(TLE_minute)
        data_propagate_J2 = data_propagate_J2 + '0000' #TLE_hour_string + TLE_minute_string 
        datetime_propagate_J2_start = datetime.strptime(data_propagate_J2, date_format)
        datetime_propagate_J2_end = datetime_propagate_J2_start + timedelta(days = days_selected)

       
        #Propagate 10 days at x minute intervals, after requested epoch
        T_prop = t = np.arange(datetime_propagate_J2_start, datetime_propagate_J2_end, timedelta(minutes=time_step)).astype(datetime)
        propagation_durations = T_prop - datetime_TLE_epoch
        #nodal_precession_rate = 1.6229E-6
        final_RAAN_NS = np.zeros(propagation_durations.shape)
        angle_between_vernal_and_pm = np.zeros(propagation_durations.shape)
        hpop_correction = np.zeros(propagation_durations.shape)
        for propagation_duration in enumerate(propagation_durations):
            #end_epoch = propagation_duration[1] + datetime_propagate_J2_start
            hpop_correction[propagation_duration[0]] = correct_J2_to_HPOP(propagation_duration[1].total_seconds()/86400)
            
            final_RAAN_NS[propagation_duration[0]] = (current_RAAN + nodal_precession_rate*(propagation_duration[1].total_seconds()))%360 + hpop_correction[propagation_duration[0]]
            angle_between_vernal_and_pm[propagation_duration[0]] = output_angle_between_vernal_and_pm(T_prop[propagation_duration[0]].strftime('%Y-%m-%d %H:%M:%S'))
            #print(f"RAAN of cat id {cat_id} is {final_RAAN_NS[propagation_duration[0]]:.5f} deg on {T_prop[propagation_duration[0]].strftime('%d/%m/%Y, %H:%M')} UTC, {propagation_duration[1].days} days from today.")


        final_RAAN_T2, T_prop_T2 = RAAN_T2(final_RAAN_NS, T_prop)
        inc_degrees = current_i*180/np.pi
        print(f"Inclination of orbit is: {inc_degrees}")

    overall_df['T_prop_T2 (UTC)'] = T_prop_T2
    overall_df['final_RAAN_T2'] = final_RAAN_T2
    overall_df['inc_degrees'] = inc_degrees
    overall_df['angle_between_vernal_and_pm'] = angle_between_vernal_and_pm

    #plt.plot(angle_between_vernal_and_pm)
    #plt.show()
    return overall_df

def correct_J2_to_HPOP(days):
    '''
    Function to correct J2 perturbation with HPOP results from STK via interpolation
    Deviation = RAAN_HPOP- RAAN_J2
    '''
    df_correction = pd.read_csv("deviation_between_j2_hpop.csv")
    return np.interp(days, df_correction["days_from_TLE_epoch"], df_correction["deviation_between_J2_and_HPOP"])



def output_launch_times(overall_df, launch_site_coords, RAAN_tol):
    '''
    angle_between_vernal_and_pm tabulates the difference between vernal and prime meridian and each propagated time.
    '''
    long_ECI = output_long_ECI(overall_df["inc_degrees"], launch_site_coords, overall_df["angle_between_vernal_and_pm"])

    overall_df["long_ECI"] = long_ECI
    overall_df["del_RAAN"] = overall_df["final_RAAN_T2"] - long_ECI

    valid_launch_conditions = overall_df[(overall_df["del_RAAN"]< RAAN_tol) & (overall_df["del_RAAN"]>=0)]
    valid_launch_conditions.drop("inc_degrees", axis = 1)
    valid_launch_conditions.drop("angle_between_vernal_and_pm", axis = 1)
    return valid_launch_conditions

def RAAN_T2(RAAN_NS, T_prop):
    '''
    Function that outputs 
    RAAN_T2--> target RAAN of T2. 
    alter_timesteps --> list of times to launch accounting for orbit insertion duration.  
    '''
    RAAN_T2 = (RAAN_NS - 120)%360
    altered_timesteps = T_prop - timedelta(minutes = 19)
    return RAAN_T2, altered_timesteps

def output_long_ECI(i, launch_site_coords, angle_between_vernal_and_pm):
    # https://space.stackexchange.com/questions/21796/relation-between-launch-window-and-raan-and-argument-of-perigee
    # Outputs an array. Cannot assume that 0 degree meridian aligns with vernal point
    # Longitude of launch site is in ECEF while RAAN is ECI. 
    # This function transforms the launch site longitude into its ECI counterpart which varies with time of day. 
    # Longitude of launch site in ECI frame is given as: 
    long_ECI = np.zeros(angle_between_vernal_and_pm.shape)
    long_ECI = np.mod(-angle_between_vernal_and_pm + launch_site_coords[1] - (np.arcsin(np.tan((90-i)*np.pi/180) * np.tan(launch_site_coords[0])*np.pi/180)*180/np.pi), 360)
    
    return long_ECI 

def output_angle_between_vernal_and_pm(time):
    '''
    This angle is formally refered to as Greenwich sidereal angle (GST), the angle between the prime meridian and the vernal equinox.
    We proceed by:
    1. Obtaining the Greenwich Apparent Sidereal Time (GAST), which is the greenwich mean sidereal time 
        corrected for shift in position of vernal equinox due to nutation. Do this using astropy.
    2. Convert GAST to GST with the conversion: GST = GAST * 15 deg/hour
    Source: https://lweb.cfa.harvard.edu/~jzhao/times.html#:~:text=Greenwich%20Mean%20Sidereal%20Time%20(GMST,the%20equinox%20due%20to%20nutation.

    time: example--> '2018-03-14 23:48:00'

    '''
    from astropy.time import Time

    t = Time(time, scale='utc')
    GAST = t.sidereal_time('apparent', 'greenwich')  
    rotation_angle = GAST.value*15

    return rotation_angle
    
if __name__=="__main__":

   

    RAAN_tol = 5
    data_propagate_J2_list = ["20230319"]
    cat_id = 53297
    use_today_date = False
    overall_df = propagate_J2(data_propagate_J2_list, cat_id, use_today_date)

    launch_site_coords = (13.73204, 80.23621)
    launch_azimuth = 102 #degree
    launch_data = output_launch_times(overall_df, launch_site_coords, RAAN_tol)
    overall_df.to_csv("output_data.csv")

    launch_data.to_csv("valid_launch_times.csv", index = False)