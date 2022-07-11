# Author: Andrew Ng 
# Date: 11/7/2022
# Script to pull data from Celestrak's server https://celestrak.org/NORAD/elements
import TLE_pull 

if __name__ == '__main__':
    TLE_pull.check_tle(cat_id = 25544)
    TLE_pull.check_tle(sat_name = "TELEOS 1")
