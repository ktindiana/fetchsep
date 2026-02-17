import datetime
import pandas as pd
import zulu

""" Various subroutines to convert strings to dates
    and dates to different formats.
"""

def str_to_datetime(date):
    """ String date to datetime
        
        INPUTS:
        
        :date: (string) - date as "YYYY-MM-DD" or "YYYY-MM-DD HH:MM:SS" or YYYY-MM-DDTHH:MM:SSZ or YYYY-MM-DDTHH:MM:SS
        
        OUTPUTS:
        
        :dt: (datetime) - datetime conversion of date
    """
    if isinstance(date, datetime.date):
        return date
    
    if date == '':
        return pd.NaT
    
    #If user entered zulu time or similar
    #YYYY-MM-DDTHH:MM:SSZ or YYYY-MM-DDTHH:MM:SS
    if 'T' in date:
        date = date.replace('T', ' ')
    if 'Z' in date:
        date = date.replace('Z','')
    
    if len(date) == 10: #only YYYY-MM-DD
        date = date  + ' 00:00:00'
        
    if len(date) == 16: #only YYYY-MM-DD HH:MM
        date = date  + ':00'
        
    dt = datetime.datetime.strptime(date, "%Y-%m-%d %H:%M:%S")
    return dt


def time_to_zulu(dt): #make_ccmc_zulu_time(dt):
    """ Make a datetime string in the format YYYY-MM-DDTHH:MM:SSZ
        
        INPUTS:
        
        :dt: (datetime)
        
        OUTPUTS:
        
        :zuludate: (string) in the format YYYY-MM-DDTHH:MM:SSZ
    
    """
    if dt == '':
        return ''
    if dt == None:
        return None
    if dt is pd.NaT:
        return pd.NaT
    if dt == 0:
        return 0

    if isinstance(dt,str):
        dt = str_to_datetime(dt)

    zdt = zulu.create(dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second)
    stzdt = str(zdt)
    stzdt = stzdt.split('+00:00')
    zuludate = stzdt[0] + "Z"
    return zuludate
 
 
def zulu_to_time(zt):
    """ Convert Zulu time to datetime
    
        INPUTS:
        
        :zt: (string) - date in the format "YYYY-MM-DDTHH:MM:SSZ"
        
        OUTPUTS:
        
        :dt: (datetime)
        
    """
    #Time e.g. 2014-01-08T05:05:00Z or 2012-07-12T22:25Z
    if zt == '':
        return ''
    if zt == None:
        return None
    if zt is pd.NaT:
        return pd.NaT
    if zt == 0:
        return 0
  
    strzt = zt.split('T')
    strzt[1] = strzt[1].strip('Z')
    n = strzt[1].split(':')
    stdt = strzt[0] + ' ' + strzt[1]

    if len(n) == 2:
        dt = datetime.datetime.strptime(stdt, '%Y-%m-%d %H:%M')
    if len(n) == 3:
        dt = datetime.datetime.strptime(stdt, '%Y-%m-%d %H:%M:%S')
    return dt


def fractional_year_to_datetime(fractional_year):
    """ Convert fractional year to datetime """
    # Extract the integer part as the year
    year = int(fractional_year)
    # Extract the fractional part
    remainder = fractional_year - year

    # Define the start of the year
    base_date = datetime.datetime(year, 1, 1)

    # Calculate the end of the year to determine total days
    # (base_date.replace(year=year + 1) creates the next year's Jan 1st)
    next_year_date = base_date.replace(year=year + 1)
    total_year_seconds = (next_year_date - base_date).total_seconds()

    # Calculate the number of seconds corresponding to the fractional part
    fractional_seconds = total_year_seconds * remainder

    # Add the calculated time delta to the base date
    result_date = base_date + datetime.timedelta(seconds=fractional_seconds)

    return result_date
    


def get_next_month(date):
    """ Get the date one month from the input date.
        The day will remain the same, but the month and
        possibly year will advance.
    
        INPUTS:
        
        :date: (datetime) - reference date
        
        OUTPUTS:
        
        :nextdate: (datetime) - 1 month from reference date
        
    """
    year = date.year
    month = date.month
    day = date.day
    
    if month == 12:
        year = year + 1
        month = 1
    else:
        month = month+1
        
    nextdate = datetime.datetime(year=year,month=month,day=day)
    return nextdate
