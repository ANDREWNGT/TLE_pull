import requests
import os
import codecs
from datetime import date

def tle_file_naming(sat_name, cat_id):
    '''
    Purpose: Names tle file with date that data is pulled from server in format "YYMMDD_{cat_no}_{satellite name}.tle"

    Params

    sat_name: str
    cat_id: int
    '''

    today = date.today()
    data_pulled_day = today.strftime("%y%m%d")
    output_file_name = os.path.join(os.getcwd(), "output", f"{data_pulled_day}_{sat_name}_{cat_id}.tle")
    return output_file_name

def check_tle(cat_id= None, sat_name = None):
    '''
    Purpose: Pull latest tle from server at time of execution and save it to local machine using NORAD category number

    Params
    cat_no: int, NORAD category number
    
    '''
    # %% Request data from server
    if cat_id: 
        request_string = f"http://celestrak.org/NORAD/elements/gp.php?CATNR={cat_id}&FORMAT=TLE"
    elif sat_name: 
        request_string = f"http://celestrak.org/NORAD/elements/gp.php?NAME={sat_name}&FORMAT=TLE"
        

    f = requests.get(request_string)
    
    string_format = codecs.decode(f.content, 'UTF-8')
    
     # %% Determine sat name for output file
    if sat_name is None:
        sat_name = string_format.split("   ")[0]
        
    # %% Determine NORAD cat id for output file
    if cat_id is None: 
        cat_id = string_format.split("\n")[1].split(" ")[1][0:5]

    # %% Determine output file name
    output_file_name = tle_file_naming(sat_name, cat_id)
    print(output_file_name)

    # %% Save data to file
    with open(output_file_name, 'wb') as file: 
        file.write(f._content)
