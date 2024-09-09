#Information about GOES-Primary Spacecraft
import pandas as pd
import datetime
import os


def import_goes_status():
    """ Read in GOES_Primary_Secondary_Status.csv provided by
        Kim Moreland (NOAA SWPC) on September 5, 2024.
        The dates in the file go up to September 6, 2024.
        
        Provides a date ranges and specifies the primary and secondary
        GOES satellites for protons and X-rays during these times.

        INPUTS:
        
            None
            
        OUTPUTS:
        
            df: (pandas DataFrame) dataframe containing dates and goes status

    """

    filename = 'GOES_Primary_Secondary_Status.csv'
    if not filename.exists():
        sys.exit(f"The GOES status file, {filename}, is not present in the "
            "fetchsep/utils/ directory. Please download from the fetchsep repository "
            "or contact kathryn.whitman@nasa.gov for assistance. "
            "https://github.com/ktindiana/fetchsep")
    
    df = pd.read_csv(filename, parse_dates=True)
    df['Start Date'] = pd.to_datetime(df['Start Date'])
    df['End Date'] = pd.to_datetime(df['End Date'])

    return df
    
    
def id_goes_spacecraft(status, instrument, date):
    """ Identify the primary or secondary GOES satellite for protons
        or X-rays in a certain date range.
        
        INPUTS:
        
            :status: (string) "Primary" or "Secondary"
            :instrument: (string) "Protons" or "X-rays"
            :date: (datetime) date of interest
            
        OUTPUTS:
        
            :experiment: (string) name of GOES spacecraft that satisfies
                the query
        
    """
    
    sub = df.loc[(df['Status']==status) & (df['Instrument']==instrument) & (df['Start Date'] <= date) & (df['End Date'] > date)]

    if sub.empty:
        print(f"Cannot identify the {status} {instrument} GOES satellite for {date}.")
        return None

    goes = sub['Satellite']
    
    return goes
    
    
    
def goes_primary_lookup(date):
    """ Determine which was the primary GOES experiment for a given timeframe.
     
        Lookup table provided by Kim Moreland (NOAA SWPC) on Sep 6, 2024.
        
        Note that all GOES real time integral fluxes retrieved from the CCMC hapi server
        are, by definition, primary GOES fluxes extracted in real time from
        SWPC's 3-day jsons. There is no need to check for primary spacecraft if
        using that data source.
        
    """

    #Dates for which historical GOES primary/secondary are known
    table_date = datetime.datetime(2024,9,6,23,23,59)

    #Date beyond table
    is date > table_date:
        print("Do not have information about primary/secondary GOES for this date. "
                "Returning GOES_RT, but this is only valid for integral fluxes. "
                "You may check "
                "https://services.swpc.noaa.gov/json/goes/instrument-sources.json "
                "if you are looking for current GOES data. If GOES integral fluxes are "
                "desired, specify GOES_RT and primary satellite fluxes will be used "
                "by default, as CCMC downloads the SWPC 3-day jsons, which are by "
                "definition from the primary satellite. This function will be updated "
                "when a reliable, living archive of GOES primary/secondary satellites "
                "becomes available.")
        
        goes = "GOES_RT"


    if date <= datetime.datetime(2024,9,6):
        df = import_goes_status()
        
        goes = id_goes_spacecraft("Primary","Protons",date)
        
    return goes
