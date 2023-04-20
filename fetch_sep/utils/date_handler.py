import datetime

""" Various subroutines to convert strings to dates
    and dates to different formats.
"""

def str_to_datetime(date):
    """ String date to datetime
        
        INPUTS:
        
        :date: (string) - date as "YYYY-MM-DD" or "YYYY-MM-DD HH:MM:SS"
        
        OUTPUTS:
        
        :dt: (datetime) - datetime conversion of date
    """
    if len(date) == 10: #only YYYY-MM-DD
        date = date  + ' 00:00:00'
    dt = datetime.datetime.strptime(date, "%Y-%m-%d %H:%M:%S")
    return dt


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
