import math


def calc_J2_rate(a, i, e):
    # J2 effect constants
    mu = 398600.44189 #km^3/s^2
    Re = 6378 #km 
    J2_const = 1.08262668* 10**(-3)
    
    mean_motion = math.sqrt(mu/(a**3))
    perturbed_mean_motion = mean_motion * (1 + 3/4*J2_const * (Re/a)**2 * (1- e**2)**(-3/2)* (3*math.cos(i)**2-1))
    nodal_precession_rate = -3/2 * J2_const * (Re/(a*(1-e**2)))**2 * perturbed_mean_motion * math.cos(i) * 180 / math.pi# deg/s
    return nodal_precession_rate
if __name__=="__main__":
    Re = 6378
    nominal_alt = 585
    for alt_offset in range(1, 50, 5):
        #alt_offset = 20 #km below
        list_a=[nominal_alt + Re, nominal_alt + Re - alt_offset]
        i=10
        e=0.01
        J2_rates=[]

        for a in list_a:
            #Time to complete RAAN maneuver via different altitudes
            J2_rates.append(calc_J2_rate(a, i, e))

        J2_diff = J2_rates[1]-J2_rates[0]
        time_to_complete_raan = 2/J2_diff/(24*60*60) #days
        print(f"{alt_offset}, {time_to_complete_raan} days")