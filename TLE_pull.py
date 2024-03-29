import requests
import os
import codecs
from datetime import date

def tle_file_naming(sat_name, cat_id, output_folder):
    '''
    Purpose: Names tle file with date that data is pulled from server in format "YYMMDD_{cat_no}_{satellite name}.tle"

    Params

    sat_name: str
    cat_id: int
    '''

    today = date.today()
    data_pulled_day = today.strftime("%y%m%d")
    if output_folder is None:
        path = os.path.join(os.getcwd(), "output")
        if not os.path.exists(path):
            os.mkdir(path)
    else: 
        path = output_folder

    output_file_name = os.path.join(path, f"{data_pulled_day}_{sat_name}_{cat_id}.tle")
    return output_file_name

def check_tle(cat_id= None, sat_name = None, output_folder = None):
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
    output_file_name = tle_file_naming(sat_name, cat_id, output_folder)
    print(output_file_name)

    # %% Save data to file
    with open(output_file_name, 'wb') as file: 
        file.write(f._content)
    return output_file_name

if __name__=="__main__":
    NORAD_cat_id=[52935,
                        41169,
                        32060,
                        35946,
                        40115,
                        38012,
                        39019,
                        48268,
                        49070,
                        39418,
                        40072,
                        41601,
                        41771,
                        41772,
                        41773,
                        41774,
                        42987,
                        42988,
                        42989,
                        42990,
                        42991,
                        42992,
                        43797,
                        43802,
                        45788,
                        45789,
                        45790,
                        46179,
                        46180,
                        46235,
                        ]
    date_format_today = "%Y%m%d"
    today = date.today()
    data_pulled_day = today.strftime(date_format_today)
    parent_output_folder = os.path.join(os.getcwd(), "pulled_data")
    if not os.path.exists(parent_output_folder):
        os.mkdir(parent_output_folder)
    output_folder = os.path.join(parent_output_folder, f"data_{data_pulled_day}")
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    for cat_id in NORAD_cat_id:
        check_tle(cat_id = cat_id, output_folder = output_folder)