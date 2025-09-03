import datetime
import matplotlib.pyplot as plt
import os
import pytz
import tkinter

class Path:
    def __init__(self):
        self.base = os.path.abspath(os.path.dirname(__file__))
        self.output = os.path.join(self.base, 'output')
        self.event_list_directory = os.path.abspath(os.path.join(self.base, 'event_list_common_format'))
        self.data = os.path.abspath(os.path.join(self.base, 'data'))
path = Path()
if not os.path.exists(path.output):
    os.mkdir(path.output)

class Time:
    def __init__(self):
        self.epoch = datetime.datetime(year=1970, month=1, day=1, tzinfo=pytz.UTC)
        self.iswa_format_swap = datetime.datetime(year=2023, month=5, day=23, tzinfo=pytz.UTC)
        self.iswa_minimum_time = datetime.datetime(year=2010, month=4, day=14, hour=0, minute=0, second=0, tzinfo=pytz.UTC)
time = Time()

class Color:
    def __init__(self):
        self.associations = {'Hits'                    : '#2e7d32',
                             'Correct Negatives'       : '#0b3d91',
                             'Misses'                  : '#c62828',
                             'False Alarms'            : '#f57f17',
                             'Neutral'                 : '#0e1111',
                             'Not Evaluated'           : '#000000',
                             'Long X-Ray Flux'         : '#ffa500',
                             'Short X-Ray Flux'        : '#5400a8',
                             '>=1 MeV Proton Flux'     : '#b3b3b3',
                             '>=5 MeV Proton Flux'     : '#ffd480',
                             '>=10 MeV Proton Flux'    : '#ff0000',
                             '>=30 MeV Proton Flux'    : '#6b3d9a',
                             '>=50 MeV Proton Flux'    : '#0000ff',
                             '>=60 MeV Proton Flux'    : '#000000',
                             '>=100 MeV Proton Flux'   : '#00ff00',
                             '>=500 MeV Proton Flux'   : '#f39c12',
                             '38-53 keV Electron Flux' : '#808080',
                             '175-315 keV Electron Flux' : '#00faff',
                             'Probability Color Main'  : '#7814e3',
                             'Probability Color Sub'   : '#c90eb7',
                             'Divider'                 : '#d92906',
                             'Eruption Out of Range'   : '#000000',
                             'Trigger/Input after Observed Phenomenon' : '#000000',
                             'No Matching Threshold' : '#000000',
                             'Ongoing SEP Event' : '#000000',
                             'Unmatched' : '#000000',
                             'No Prediction Provided' : '#000000',
                             None: 'none'
                             }
        self.associations['&ge; 1'] = self.associations['> 1 MeV'] = self.associations['>=1 MeV Proton Flux']
        self.associations['&ge; 5'] = self.associations['> 5 MeV'] = self.associations['>=5 MeV Proton Flux']
        self.associations['&ge; 10'] = self.associations['> 10 MeV'] = self.associations['>=10 MeV Proton Flux']
        self.associations['&ge; 30'] = self.associations['> 30 MeV'] = self.associations['>=30 MeV Proton Flux']
        self.associations['&ge; 50'] = self.associations['> 50 MeV'] = self.associations['>=50 MeV Proton Flux']
        self.associations['&ge; 60'] = self.associations['> 60 MeV'] = self.associations['>=60 MeV Proton Flux']
        self.associations['&ge; 100'] = self.associations['> 100 MeV'] = self.associations['>=100 MeV Proton Flux']
        self.associations['&ge; 500'] = self.associations['> 500 MeV'] = self.associations['>=500 MeV Proton Flux']
        
        self.matplotlib_color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color'] * 10

color = Color()

class Window:
    def __init__(self):
        self.width = 1200
        self.height = 1200
        self.horizontal_offset = 0
        self.vertical_offset = 0
#        self.horizontal_offset = 4500
#        self.vertical_offset = 0
        self.size = str(self.width) + 'x' + str(self.height)
        self.offset = '+' + str(self.horizontal_offset) + '+' + str(self.vertical_offset)
        self.geometry = self.size + self.offset
window = Window()

class Font:
    def __init__(self):
        self.size = 16
        self.typeface = 'Arial'
        #self.button_font = tkinter.font.Font
        self.button_font_settings = {'size' : self.size, 
                                     'family' : self.typeface}
font = Font()

class Widget:
    def __init__(self):
        self.width = 20
        self.height = 1
        self.font = (font.size, font.typeface)
        self.format = {'width' : self.width,
                       'height' : self.height,
                       'font' : self.font}
widget = Widget()

class Textbox:
    def __init__(self):
        self.width = 20
        self.font = (font.size, font.typeface)
        self.format = {'width' : self.width,
                       'font' : self.font}
textbox = Textbox()

class Table:
    def __init__(self):
        self.column_width = 100
table = Table()

