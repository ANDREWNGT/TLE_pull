# Author: Andrew Ng 
# Date: 11/7/2022
# Script to pull data from Celestrak's server https://celestrak.org/NORAD/elements
import TLE_pull 
from datetime import date
import os
if __name__ == '__main__':
    # Example to use NORAD Cat ID: 
    # TLE_pull.check_tle(cat_id = 25544)
    
    today = date.today()
    data_pulled_day = today.strftime("%y%m%d")
    # Example to use satellite exact name. Note that name must match CelesTrak's exact format. 
    # TLE_pull.check_tle(sat_name = "STARLINK-1007")
    parent_output_folder = r"C:\Users\Ngengton\Documents"
    output_folder = os.path.join(parent_output_folder, f"data_{data_pulled_day}")
    
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
   
    TLE_pull.check_tle(cat_id = 41169, output_folder= output_folder) # Pulls data of Fengyun 1c DEB. Saves to folder that is specified
    TLE_pull.check_tle(cat_id = 25544) # Pulls data of the ISS. Saves to output folder in git repo.