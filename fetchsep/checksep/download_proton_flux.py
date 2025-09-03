import datetime
import io
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests

import config

def rerequest(url, tries=0):
    """
    Runs requests.get() until a response is received to avoid chokepoints

    Parameters
    ----------
    url : string
        URL where desired data resides
    tries [untouched by user] : int
        Number of attempts
    """

    # end the requesting process if tries > 5
    if tries > 5:
        raise Exception(url + " refuses to respond after 5 attempts. Exiting.")
    try:
        output = requests.get(url, timeout=10.0)
        return output
    except requests.exceptions.Timeout as e:
        return rerequest(url, tries + 1)

def backfill(df):
    # CURE DATA
    non_time_columns = df.select_dtypes(exclude=['datetime64[ns, UTC]', 'timedelta']).columns
    df[non_time_columns] = df[non_time_columns].map(lambda x: x if x >= 0 else np.nan)
    for column in non_time_columns:
        df[column] = df[column].interpolate(method='linear')
    return df

def no_negatives(df): 
    # CURE DATA
    non_time_columns = df.select_dtypes(exclude=['datetime64[ns, UTC]', 'timedelta']).columns
    df[non_time_columns] = df[non_time_columns].map(lambda x: x if x >= 0 else np.nan)
    return df

def download_hapi_flux(flux_type_url, start_datetime, end_datetime, backfill_flag=True, no_negatives_flag=True): 
    if end_datetime < config.time.iswa_minimum_time:
        print('WARNING: ISWA does not carry data prior to 2010-04-14T00:00:00Z. Expect empty plots.')
    base_url = 'https://iswa.gsfc.nasa.gov/IswaSystemWebApp/hapi/data?'
    url = base_url + flux_type_url + '&time.min=' + start_datetime.strftime('%Y-%m-%dT%H:%M:%S.%fZ') + '&time.max=' + end_datetime.strftime('%Y-%m-%dT%H:%M:%S.%fZ') + '&format=csv'
    response = rerequest(url) 
    if response.status_code == 200:
        data = response.text
        df = pd.read_csv(io.StringIO(data))
        df.columns = ['time_tag', '>=1 MeV', '>=5 MeV', '>=10 MeV', '>=30 MeV', '>=50 MeV', '>=100 MeV', 'E>=0.8 MeV', 'E>=2 MeV', 'E>=4 MeV', '>=60 MeV', '>=500 MeV']
        df['time_tag'] = pd.to_datetime(df['time_tag'], format='%Y-%m-%dT%H:%M:%SZ', errors='coerce', utc=True)
    else:
        print(f'Failed to retrieve data. HTTP Status code: {response.status_code}')
        print('Exiting...')
        exit()
    if backfill_flag:
        df = backfill(df)
    if no_negatives_flag:
        df = no_negatives(df)
    return df




def download_goes_flux(_, start_datetime, end_datetime, backfill_flag=True, no_negatives_flag=True): 
    flux_type_url = 'id=goesp_part_flux_P5M'
    df = download_hapi_flux(flux_type_url, start_datetime, end_datetime, backfill_flag, no_negatives_flag)
    return df

def download_soho_flux(flux_type, start_datetime, end_datetime, backfill_flag=True, no_negatives_flag=True):
    pass

def download_ace_sis_flux(_, start_datetime, end_datetime, backfill_flag=True, no_negatives_flag=True):
    flux_type_url = 'id=ace_sis_P5M'
    df = download_hapi_flux(flux_type_url, start_datetime, end_datetime, backfill_flag, no_negatives_flag)
    return df

def download_flux(source, _, start_datetime, end_datetime, backfill_flag=True, no_negatives_flag=True):
    source_function_dict = {'GOES' : download_goes_flux,
                            'SOHO' : download_soho_flux,
                            'ACE SIS' : download_ace_sis_flux}
    df = source_function_dict[source](_, start_datetime, end_datetime, backfill_flag, no_negatives_flag)
    return df

    

if __name__ == '__main__':
    
    start_datetime = datetime.datetime.strptime('2024-02-01', '%Y-%m-%d').replace(tzinfo=datetime.timezone.utc)
    end_datetime = datetime.datetime.strptime('2024-03-01', '%Y-%m-%d').replace(tzinfo=datetime.timezone.utc)
    df = download_goes_flux('proton', start_datetime, end_datetime)
    plt.figure()
    plt.plot(df['time_tag'], df['>=10 MeV'], color='red')
    plt.yscale('log')
    plt.show()




