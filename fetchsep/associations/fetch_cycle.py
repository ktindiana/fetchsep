from ..utils import date_handler as dh
import os
import sys
import datetime

#Solar cycle definitions used by Badhwar-O'Neill model in fractional year
solar_cycle = {
    "SC18": {   "start_time": datetime.datetime(1944, 2, 28),
                "max_time": datetime.datetime(1947, 5, 5),
                "end_time": datetime.datetime(1954, 4, 16),
                "id": 18},
    "SC19": {   "start_time": datetime.datetime(1954, 4, 16),
                "max_time": datetime.datetime(1958, 2, 28),
                "end_time": datetime.datetime(1964, 10, 23),
                "id": 19},
    "SC20": {   "start_time": datetime.datetime(1964, 10, 23),
                "max_time": datetime.datetime(1968, 10, 30),
                "end_time": datetime.datetime(1976, 5, 30),
                "id": 20},
    "SC21": {   "start_time": datetime.datetime(1976, 5, 30),
                "max_time": datetime.datetime(1979, 12, 28),
                "end_time": datetime.datetime(1986, 9, 6),
                "id": 21},
    "SC22": {   "start_time": datetime.datetime(1986, 9, 6),
                "max_time": datetime.datetime(1989, 10, 30),
                "end_time": datetime.datetime(1996, 4, 30),
                "id": 22},
    "SC23": {   "start_time": datetime.datetime(1996, 4, 30),
                "max_time": datetime.datetime(2001, 12, 2),
                "end_time": datetime.datetime(2008, 11, 18),
                "id": 23},
    "SC24": {   "start_time": datetime.datetime(2008, 11, 18),
                "max_time": datetime.datetime(2014, 4, 2),
                "end_time": datetime.datetime(2019, 12, 13),
                "id": 24},
    "SC25": {   "start_time": datetime.datetime(2019, 12, 13),
                "max_time": datetime.datetime(2025, 7, 17),
                "end_time": datetime.datetime(2031,1,1),
                "id": 25},
}


def get_solar_cycle(request_date):
    """ Identify which solar cycle contains the requested date
        and return the integer identifier, e.g. 25 for Solar Cycle 25.
        
    """
    keys = solar_cycle.keys()
    for key in keys:
        if request_date >= solar_cycle[key]['start_time'] and \
            request_date < solar_cycle[key]['end_time']:
            return solar_cycle[key]['id']

    print(f"get_solar_cycle: Could not identify solar cycle for {request_date}.")
    return None
