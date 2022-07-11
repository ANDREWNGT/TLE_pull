# Author: Andrew Ng 
# Date: 11/7/2022
# Script to pull data from Celestrak's server https://celestrak.org/NORAD/elements
import TLE_pull 

if __name__ == '__main__':
    # Example to use NORAD Cat ID: 
    # TLE_pull.check_tle(cat_id = 25544)

    # Example to use satellite exact name. Note that name must match CelesTrak's exact format. 
    # TLE_pull.check_tle(sat_name = "STARLINK-1007")

    TLE_pull.check_tle(cat_id = 25544) # Pulls data of the ISS. User is recommended to edit the input parameter of this line. 
