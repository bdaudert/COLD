import sys
import numpy as np
import json

def compute_mon_day(doy):
    '''
    Reverse of compute_doy but counting every feb as having 29 days
    '''
    ndoy = int(doy)
    mon = 0
    day = 0
    if ndoy >366 or ndoy < 1:
        return None,None
    mon_day_sum = [31,60,91,121,152,182,213,244,274,305,335,366]
    for i in range(12):
        if i == 0:
            if ndoy <=31:
                mon = 1
                day = ndoy
                break
        else:
            if mon_day_sum[i-1] < ndoy and ndoy <= mon_day_sum[i]:
                mon = i+1
                day = ndoy - mon_day_sum[i-1]
                break
            else:
                continue
    return mon,day

def get_cold_days(var_name, data_dir):
    with open(data_dir + var_name + '_pca_component_1_ts.json') as f:
        ts_data = json.load(f)
    dates = []
    vals = []
    for d_idx, date_val in enumerate(ts_data):
        '''
        d_int = date_val[0]
        year = 1950 + d_int / 90
        doy = d_int % 90
        mon, day = compute_mon_day(doy)
        mon = str(mon); day = str(day)
        if len(mon) == 1: mon = '0' + mon
        if len(day) == 1: day = '0' + day
        dates.append(str(year) + '-' + mon + '-' + day)
        '''
        dates.append(date_val[0])
        vals.append(date_val[1])
    np_vals = np.array(vals)
    #find 10 coldest days
    max_idx = list(np.argsort(np_vals)[-10:])
    max_dates = []
    for m_idx in max_idx:
        max_dates.append(dates[m_idx])
    return max_dates

########
#M A I N
########
if __name__ == '__main__' :
    years = range(1951,2012)
    data_dir = 'RESULTS/livneh/'
    for var_name in ['tmin', 'tmax']:
        cold_dates = get_cold_days(var_name, data_dir)
        print cold_dates
