import argparse
import datetime
import matplotlib
import matplotlib.backends.backend_tkagg
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pickle as pkl
import tkinter
import tkinter.ttk
import itertools

import config
import download_proton_flux


# from ..utils import read_datasets as datasets


"""
checkSEP GUI
Made by Clayton Allison and Luke Stegeman
"""


def from_fetchsep(fetchsep_filename, flux_files = None, energy_channel = '10.0--1 MeV'):
    # read idsep_file and get start/end time lists
    headers = ['Start Time', 'End Time', 'Experiment', 'Flux Type', 'Flags', 'User Experiment Name',  'User Filename', 'Options', 'BGStart', 'BGEnd', 'JSON Type', 'Spacecraft', 'IDSEP Path', 'Location', 'Species']
    
    df = pd.read_csv(fetchsep_filename) #


    # run checksep
    initial_ends, confirmed_ends = checksep(df, flux_files, energy_channel, fetchsep_filename)

    return output_filename



def checksep(df, flux_files, energy_channel, fetchsep_file):


    # INITIALIZE TKINTER ROOT WINDOW
    root = tkinter.Tk()
    geometry_string = get_window_geometry(root)
    root.geometry(geometry_string)
    
    fetchsep_df = df
    app = CheckSEPApp(root, df, flux_files, energy_channel) #, list(event_labels.values()), time_buffer=datetime.timedelta(days=time_buffer), separate_energies=separate_energies)
    root.mainloop()
    initial_ends, confirmed_ends = app()


    # get new end time lists
    print(confirmed_ends)
    print(pd.isnull(confirmed_ends))
    filter_ = np.where(pd.isnull(confirmed_ends) == False)
    print(filter_)
    filtered_ends = [confirmed_ends[x] for x in filter_[0]]
    print(filtered_ends)
    filtered_ends = list(itertools.chain.from_iterable(filtered_ends))
    print(filtered_ends)
    new_ends = [pd.to_datetime(x) for x in filtered_ends]
    print(new_ends)
    concat_ends = sorted(initial_ends + new_ends)
    print(concat_ends)

    # create output file

    output_file = create_checksep_list(fetchsep_file, concat_ends)


    # feed back to fetchsep
    return initial_ends, confirmed_ends, output_file




def create_checksep_list(fetchsep_file, checksep_ends):


    headers = ['Start Time', 'End Time', 'Experiment', 'Flux Type', 'Flags', 'User Experiment Name',  'User Filename', 'Options', 'BGStart', 'BGEnd', 'JSON Type', 'Spacecraft', 'IDSEP Path', 'Location', 'Species']
    output_filename = './output/checksep/checksep_batch_event_list.txt'
    fetchsep_dataframe = pd.read_csv(fetchsep_file)
    checksep_dataframe = pd.DataFrame()
    print(fetchsep_dataframe)
    initial_start = fetchsep_dataframe['#Start Time'][0]
    for e in range(len(checksep_ends)):
        if e == 0:
            start_entry = initial_start
            end_entry = datetime.datetime.strptime(str(checksep_ends[e]), '%Y-%m-%d %H:%M:%S')
            
        else:
            start_entry = checksep_ends[e-1]
            end_entry = datetime.datetime.strptime(str(checksep_ends[e]), '%Y-%m-%d %H:%M:%S')
        for d in range(len(fetchsep_dataframe)):
            if end_entry == fetchsep_dataframe['End Time'][d]:
                exp_entry = fetchsep_dataframe['Experiment'][d]
                flux_entry = fetchsep_dataframe['Flux Type'][d]
                flag_entry = fetchsep_dataframe['Flags'][d]
                user_exp_entry = fetchsep_dataframe['User Experiment Name'][d]
                user_file_entry = fetchsep_dataframe['User Filename'][d]
                opt_entry = fetchsep_dataframe['Options'][d]
                bgstart_entry = fetchsep_dataframe['BGStart'][d]
                bgend_entry = fetchsep_dataframe['BGEnd'][d]
                json_entry = fetchsep_dataframe['JSON Type'][d]
                space_entry = fetchsep_dataframe['Spacecraft'][d]
                path_entry = fetchsep_dataframe['IDSEP Path'][d]
                loc_entry = fetchsep_dataframe['Location'][d]
                spec_entry = fetchsep_dataframe['Species'][d]
            elif checksep_ends[e] > pd.to_datetime(fetchsep_dataframe['#Start Time'][d]) and checksep_ends[e] < pd.to_datetime(fetchsep_dataframe['End Time'][d]):
                exp_entry = fetchsep_dataframe['Experiment'][d]
                flux_entry = fetchsep_dataframe['Flux Type'][d]
                flag_entry = fetchsep_dataframe['Flags'][d]
                user_exp_entry = fetchsep_dataframe['User Experiment Name'][d]
                user_file_entry = fetchsep_dataframe['User Filename'][d]
                opt_entry = fetchsep_dataframe['Options'][d]
                bgstart_entry = fetchsep_dataframe['BGStart'][d]
                bgend_entry = fetchsep_dataframe['BGEnd'][d]
                json_entry = fetchsep_dataframe['JSON Type'][d]
                space_entry = fetchsep_dataframe['Spacecraft'][d]
                path_entry = fetchsep_dataframe['IDSEP Path'][d]
                loc_entry = fetchsep_dataframe['Location'][d]
                spec_entry = fetchsep_dataframe['Species'][d]
            else:
                pass
        row = pd.DataFrame.from_dict({'#Start Time': [start_entry],
                'End Time': [end_entry],
                'Experiment': [exp_entry],
                'Flux Type': [flux_entry],
                'Flags': [flag_entry],
                'User Experiment Name': [user_exp_entry],
                'User Filename': [user_file_entry],
                'Options': [opt_entry],
                'BGStart': [bgstart_entry],
                'BGEnd': [bgend_entry],
                'JSON Type': [json_entry],
                'Spacecraft': [space_entry],
                'IDSEP Path': [path_entry],
                'Location': [loc_entry],
                'Species': [spec_entry]
                })
        checksep_dataframe = pd.concat([checksep_dataframe, row], ignore_index = True)
    
    checksep_dataframe.to_csv(output_filename, index=False)
        



    return output_filename



def days_since_epoch(t):
    # t = pytz.timezone('UTC').localize(t)
    difference = t - config.time.epoch
    days = difference.days + difference.seconds / 60 / 60 / 24
    return days

def datetime_given_days_since_epoch(days):
    t = config.time.epoch + datetime.timedelta(days=days)
    return t


def get_proton_data(start_datetime, end_datetime, observation=None, instrument=None, download = False):
    
    
    
    if download:
        if observation is not None:
            df = pd.read_csv(observation)
            df = df.rename(columns={'dates' : 'time_tag'})
            df['time_tag'] = pd.to_datetime(df['time_tag'])
            for column in df.columns:
                if column != 'time_tag':
                    df = df.rename(columns={column : column + ' MeV'})
        else:
            if instrument in ['GOES', 'SOHO', 'ACE SIS']:
                df = download_proton_data(start_datetime, end_datetime, instrument)
            else:
                print('Sorry. This instrument is not currently supported: ', instrument)
                exit()
        return df
    else:
        if observation is not None:
            dates_df = pd.DataFrame()
            for i in range(len(observation)):
                print(observation[i])
                foo = pd.read_csv(observation[i])
                dates_df = pd.concat([foo, dates_df], ignore_index = True)
            dates_df = dates_df.rename(columns={'dates' : 'time_tag'})
            # print(dates_df.columns)
            dates_df['time_tag'] = pd.to_datetime(dates_df['time_tag'])
            for column in dates_df.columns:
                if column != 'time_tag':
                    dates_df = dates_df.rename(columns={column : column + ' MeV'})
        return dates_df

def download_proton_data(start_datetime, end_datetime, instrument='GOES', energy=[]):
    df = download_proton_flux.download_flux(instrument, 'proton', start_datetime, end_datetime, backfill_flag=False)
    return df

def format_fetchsep_list(event_list, observation=None, instrument=None, index=0, energy_channel = '10.0 - -1'):
    a = open(event_list, 'r')
    lines = a.readlines()
    a.close()
    # FIND ENERGY
    # for i in range(0, len(lines)):
    #     line = lines[i].split(':')
    #     if line[0] == '#Energy channel':
    #         energy_channel = line[1].lstrip().rstrip()
    #         break
    # print(energy_channel.split('-'))
    # if float(energy_channel.split(' - ')[1]) == -1:
    #     integral = True
    # else:
    #     integral = False

    energy = energy_channel
    print(event_list)
    headers = ['Start Time', 'End Time', 'Experiment', 'Flux Type', 'User Experiment Name',  'User Filename', 'Options', 'BGStart', 'BGEnd', 'JSON Type', 'Spacecraft', 'IDSEP Path', 'Location', 'Species']
    event_list_df = pd.read_csv(event_list)
    # event_list_df['start'] = pd.to_datetime(event_list_df['Start Time'])
    # event_list_df['end'] = pd.to_datetime(event_list_df['End Time'])
    # # event_list_df = event_list_df.drop(columns=['start_date', 'start_time', 'end_date', 'end_time'])
    event_list_df['list_index'] = index
    # event_list_df['energy'] = energy
    # event_list_df['observation'] = observation
    # event_list_df['instrument'] = 
    # event_list_df['integral'] = integral
    return event_list_df

def find_overlapping_events(dfs : list, labels : dict, separate_energies : bool, energy_channel = '10.0 - -1'):
    main_df = pd.concat(dfs, ignore_index=True)
    main_df = main_df.sort_values(by=['#Start Time', 'End Time'], ignore_index=True).reset_index(drop=True)
    # DETERMINE ALL CONTIGUOUS TIME PERIODS
    current_start = main_df.iloc[0]['#Start Time']
    current_end = main_df.iloc[0]['End Time']
    current_index = main_df.iloc[0]['list_index']
    current_energy = energy_channel
    # current_observation = main_df.iloc[0]['observation']
    current_instrument = main_df.iloc[0]['Experiment']
    current_integral_status = main_df.iloc[0]['Flux Type']
    current_group = [[current_start, current_end, current_index, labels[current_index], current_energy, current_instrument, current_integral_status]]
    groups = []
    for i in range(1, len(main_df)):
        start = main_df.iloc[i]['#Start Time']
        end = main_df.iloc[i]['End Time']
        index = main_df.iloc[i]['list_index']
        label = labels[index]
        print(label)
        # energy = main_df.iloc[i]['energy']
        # observation = main_df.iloc[i]['observation']
        instrument = main_df.iloc[i]['Experiment']
        integral_status = main_df.iloc[i]['Flux Type']
        # if separate_energies:
        #     if (start <= current_end) and (energy == current_energy):
        #         current_group.append([start, end, index, label, energy, observation, instrument, integral_status])
        #         current_end = max(current_end, end)
        #     else:
            # groups.append(current_group)
            # current_group = [[start, end, index, label, energy, observation, instrument, integral_status]]
            # current_start = start
            # current_end = end
            # current_energy = energy
        # print('where does this come from', observation)
        print(start, current_end)
        if (start <= current_end):
            current_group.append([start, end, index, label, current_energy, observation, instrument, integral_status])
            current_end = max(current_end, end)
        else:
            groups.append(current_group)
            current_group = [[start, end, index, label, current_energy, observation, instrument, integral_status]]
            current_start = start
            current_end = end
    if current_group:
        groups.append(current_group)
    dfs = pd.DataFrame()
    for group in groups:
        foo = pd.DataFrame(group, columns=['start', 'end', 'list_index', 'list', 'energy', 'observation', 'instrument', 'integral'])
        dfs = pd.concat([foo, dfs], ignore_index = True)
    return dfs

def get_window_geometry(root):
    screen_width = root.winfo_screenwidth()
    screen_height = root.winfo_screenheight()
    window_width = int(screen_width * 0.9)
    window_height = int(screen_height * 0.7)
    x_window_position = int((screen_width - window_width) / 2)
    y_window_position = int((screen_height - window_height) / 2)
    geometry_string = str(window_width) + 'x' + str(window_height) + '+' + str(x_window_position) + '+' + str(y_window_position)
    return geometry_string

class Toolbar(matplotlib.backends.backend_tkagg.NavigationToolbar2Tk):
    def __init__(self, canvas, frame, plot_click, on_click):
        super().__init__(canvas, frame)
        self.plot_click = plot_click
        self.on_click = on_click
    
    def pan(self):
        super().pan()
        if self.mode == 'pan/zoom':
            self.canvas.mpl_disconnect(self.plot_click)
        else:
            self.plot_click = self.canvas.mpl_connect('button_press_event', self.on_click)
     
    def zoom(self):
        super().zoom()
        if self.mode == 'zoom rect':
            self.canvas.mpl_disconnect(self.plot_click)
        else:
            self.plot_click = self.canvas.mpl_connect('button_press_event', self.on_click)
        


class Line:
    def __init__(self):
        self.glow = None
        self.line = None
        self.active = False

    def set_line(self, line):
        self.line = line
    
    def remove_line(self):
        self.line.remove()
        self.glow.remove()


class errorwindow:
    def __init__(self):
        self.root = root
        self.new_window = tkinter.Toplevel() 
        self.new_window.title('Warning Window')
        self.return_status = tkinter.StringVar()
        self.return_status.set("Foo")

        # SET UP CLOSING PROTOCOL
        self.root.protocol('WM_DELETE_WINDOW', self.on_close)

        self.plot_frame = tkinter.Frame(self.new_window)
        self.plot_frame.pack(side=tkinter.LEFT, fill=tkinter.BOTH, expand=True)
        
        # ADD BUTTON FRAME
        self.button_frame = tkinter.Frame(self.plot_frame)
        self.button_frame.pack(side=tkinter.BOTTOM, fill=tkinter.X)
        
        # ADD CONFIRM
        self.prev_button = tkinter.Button(self.button_frame, text="Confirm Reset", command=self.go_back, font=(config.font.typeface, config.font.size))
        self.prev_button.pack(side=tkinter.LEFT, padx=(10, 0), fill=tkinter.X, expand=True)

        # ADD REJECT
        self.rej_button = tkinter.Button(self.button_frame, text="Don't Reset", command=self.dont_go, font=(config.font.typeface, config.font.size))
        self.rej_button.pack(side=tkinter.LEFT, fill=tkinter.X, expand=True)

        # self.wait_window(self) # Wait until this window is closed
        self.new_window.wait_variable(self.return_status)
        # self.status

    def go_back(self):
        self.return_status.set("True")
        self.status
        self.new_window.destroy()                   # PROPERLY DESTROY THE TKINTER WINDOW
        
    def dont_go(self):
        self.return_status.set("False")
        self.status
        self.new_window.destroy()                   # PROPERLY DESTROY THE TKINTER WINDOW
        
    def status(self):
        return self.return_status

    def on_close(self):
        """CLEANUP AND CLOSE THE APP SAFELY."""
        self.new_window.destroy()                   # PROPERLY DESTROY THE TKINTER WINDOW



    

class PlaceholderEntry(tkinter.Entry):
    def __init__(self, master=None, placeholder='', **kwargs):
        super().__init__(master, **kwargs)
        self.placeholder = placeholder
        self.insert(0, placeholder)
        self.bind('<FocusIn>', self.on_focus_in)
        self.bind('<FocusOut>', self.on_focus_out)

    def on_focus_in(self, event):
        if self.get() == self.placeholder:
            self.delete(0, tkinter.END)
            self.config(fg='black')
    
    def on_focus_out(self, event):
        if not self.get():
            self.insert(0, self.placeholder)
            self.config(fg='gray')

class CheckSEPApp:
    def __init__(self, root, events, flux_files, energy_channel):
        self.end_times = []
        self.confirmed_times = []
        self.root = root
        self.multiple_end = False
        self.fig, self.ax = plt.subplots(figsize=(10,6))
        self.root.title('SEP Event')

        # ESTABLISH IF ENERGIES SHOULD BE SEPARATED OR NOT
        # self.separate_energies = separate_energies

        # SET UP CLOSING PROTOCOL
        self.root.protocol('WM_DELETE_WINDOW', self.on_close)

        # SAVE THE EVENT SET
        self.events = events
        self.stored_event_plots = []
        self.energies = []       
        self.confirmed_times = []

        # DEFINE THE TIME BUFFER
        self.time_buffer = datetime.timedelta(days=1)

        # SAMPLE DATA FOR MULTIPLE PLOTS (X VALUES, Y VALUES, TITLE)
        self.format_event_data(flux_files, energy_channel) 
        self.current_plot_index = 0
        self.yaxis_log_scale = True

        # CREATE A FRAME FOR THE PLOT
        self.plot_frame = tkinter.Frame(self.root)
        self.plot_frame.pack(side=tkinter.LEFT, fill=tkinter.BOTH, expand=True)

        # CREATE A FRAME FOR THE TABLE
        self.table_frame = tkinter.LabelFrame(self.root, text='CheckSEP', font=(config.font.typeface, config.font.size))
        self.table_frame.pack(side=tkinter.RIGHT, fill=tkinter.BOTH, padx=10, pady=10, expand=True)

        # ESTABLISH TABLE STYLE
        tkinter.ttk.Style().configure('Treeview', font=(config.font.typeface, config.font.size - 5))
        tkinter.ttk.Style().configure('Treeview.Heading', font=(config.font.typeface, config.font.size, 'bold')) 

        # CREATE THE TREEVIEW WIDGET FOR THE TABLE
        self.table = tkinter.ttk.Treeview(self.table_frame, columns=('Col1', 'Col2', 'Col3', 'Col4', 'Col5'), show='headings')
        self.table.heading('Col1', text='Start')
        self.table.heading('Col2', text='End')
        self.table.heading('Col3', text='Duration [hr]')
        self.table.heading('Col4', text='List')
        self.table.heading('Col5', text='Energy Channel')
        
        # SET THE COLUMNS WIDTH
        self.table.column('Col1', width=config.table.column_width)
        self.table.column('Col2', width=config.table.column_width)
        self.table.column('Col3', width=config.table.column_width)
        self.table.column('Col4', width=config.table.column_width)
        self.table.column('Col5', width=config.table.column_width)

        # ADD BUTTON FRAME
        self.button_frame = tkinter.Frame(self.plot_frame)
        self.button_frame.pack(side=tkinter.BOTTOM, fill=tkinter.X)
        
        # ADD PREVIOUS BUTTON
        self.prev_button = tkinter.Button(self.button_frame, text="<<", command=self.previous_plot, font=(config.font.typeface, config.font.size))
        self.prev_button.pack(side=tkinter.LEFT, padx=(10, 0), fill=tkinter.X, expand=True)

        # ADD A BUTTON TO TOGGLE Y-AXIS SCALE
        self.yaxis_button = tkinter.Button(self.button_frame, text="Toggle Vertical Axis Scale", command=self.toggle_yaxis_scale, font=(config.font.typeface, config.font.size))
        self.yaxis_button.pack(side=tkinter.LEFT, fill=tkinter.X, expand=True)
        
        # ADD A BUTTON TO CONFIRM YOUR MANUAL END TIME
        self.confirm_button = tkinter.Button(self.button_frame, text="Confirm New End Time", command=self.confirm_end, font=(config.font.typeface, config.font.size))
        self.confirm_button.pack(side=tkinter.LEFT, fill=tkinter.X, expand=True)

        # ADD A BUTTON TO SWITCH TO MULTIPLE ENDS MODE
        # self.mult_end_button = tkinter.Button(self.button_frame, text="Add Another End Time", command=self.add_another_end_time, font=(config.font.typeface, config.font.size))
        # self.mult_end_button.pack(side=tkinter.LEFT, padx=(0, 10), fill=tkinter.X, expand=True)

        # ADD RESET LINES BUTTON
        self.reset_lines_button = tkinter.Button(self.button_frame, text="RESET ALL MANUAL LINES FOR THIS EVENT", command=self.reset_manual_lines, font=(config.font.typeface, config.font.size))
        self.reset_lines_button.pack(side=tkinter.LEFT, padx=(0, 10), fill=tkinter.X, expand=True)

        # ADD NEXT BUTTON
        self.next_button = tkinter.Button(self.button_frame, text=">>", command=self.next_plot, font=(config.font.typeface, config.font.size))
        self.next_button.pack(side=tkinter.LEFT, padx=(0, 10), fill=tkinter.X, expand=True)

        # ADD FORCE CLOSE
        self.quit_button = tkinter.Button(self.button_frame, text="FORCE QUIT", command=self.force_quit, font=(config.font.typeface, config.font.size))
        self.quit_button.pack(side=tkinter.LEFT, padx=(0, 10), fill=tkinter.X, expand=True)

        # CREATE INITIAL PLOT
        self.create_plot()

        # PACK THE TABLE INTO THE FRAME
        self.table.pack(fill=tkinter.BOTH, expand=True)
        self.update_table()

        # SET UP GRID
        self.frame = tkinter.Frame(self.root)
        self.frame.pack(fill='both', expand=1)
        self.frame.grid_columnconfigure(0, weight=1)
        self.frame.grid_columnconfigure(1, weight=1)
        self.frame.grid_columnconfigure(2, weight=1)
        self.frame.grid_rowconfigure(1, weight=1)

        # EDIT BOX
        self.edit_frame = tkinter.Frame(self.frame)
        self.edit_frame.grid(row=0, column=2, padx=10, pady=10)
        self.edit_label = tkinter.Label(self.edit_frame, text='Marker Datetime')#, font=config.font.button_font(**config.font.button_font_settings))
        self.edit = PlaceholderEntry(self.edit_frame, placeholder='(YYYY-MM-DD HH:MM:SS)')#, font=config.font.button_font(**config.font.button_font_settings), width=config.textbox.width, state='disabled')
        self.edit.bind('<Return>', self.update_vertical_line_text)
        self.edit_label.grid(row=0, column=0, sticky='se', padx=10, pady=10) 
        self.edit.grid(row=1, column=0, sticky='se', padx=10, pady=10) 

        # BIND CLICK EVENT TO PLOT
        self.lines = []
        self.selected_line = Line()
        self.plot_click = self.fig.canvas.mpl_connect('button_press_event', self.on_click)

        # ADD THE MATPLOTLIB TOOLBAR
        self.toolbar_frame = tkinter.Frame(master=self.frame)
        self.toolbar_frame.grid(row=2, columnspan=3)
        self.toolbar = Toolbar(self.canvas, self.toolbar_frame, self.plot_click, self.on_click)
        # self.toolbar = matplotlib.backends.backend_tkagg.NavigationToolbar2Tk(self.canvas, self.plot_frame)
        self.toolbar.update()
        self.canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=True)

        # UPDATE BUTTON STATES INITIALLY
        self.update_navigation_buttons()
        

    def format_event_data(self, flux_files, energy_channel):
        counter = 0 
        for  index, event in self.events.iterrows():
            
            # print(type(self.events))
            # print('type', type(event), type(self.events))
            # print(self.events['energy'][0])
            # input()
            event_interpretations = []

            self.events['duration'] = ''# self.events.assign(['duration'][:] = '')
            # print('this is event', event)
            # for row_index, row in event.iterrows():
                # event_data_df = get_proton_data(row['Start Time'], row['End Time'], observation=row['observation'], instrument=row['instrument'])
                
                
            event_data_df = get_proton_data(event['start'], event['end'], observation=flux_files, instrument=event['instrument'])
            
            duration_timedelta = pd.to_datetime(event['end']) - pd.to_datetime(event['start'])
            duration_seconds = duration_timedelta.total_seconds()
            # self.events[counter].loc[index, 'duration'] = duration_timedelta
            event.loc['duration'] = duration_timedelta
            print(event.loc['duration'])
            # input()
            # print(type(duration_timedelta))
            # input()
            condition = (event_data_df['time_tag'] >= event['start']) & (event_data_df['time_tag'] <= event['end'])
            time = event_data_df[condition]['time_tag']
            
            flux = event_data_df[condition][energy_channel]

            event_interpretations.append((duration_seconds, duration_timedelta, time, flux)) #, row['list'], row['list_index'], row['energy'], row['observation'], row['instrument'], row['integral'])) #, buffer_time, buffer_flux))
            self.end_times.append(time.to_list()[-1])
            
            # print((duration_seconds, duration_timedelta, time, flux, row['list'], row['list_index'], row['energy'], row['observation'], row['instrument'], row['integral']))
            # input()
            try:
                event_interpretations.sort(reverse=True)
            except ValueError:
                try:
                    event_interpretations.sort(key=lambda x : x[-1])
                except ValueError:
                    pass
            self.stored_event_plots.append((event_interpretations))
                
            counter += 1

    def create_plot(self):
        
        event_interpretations = self.stored_event_plots[self.current_plot_index]
        # CREATE A NEW MATPLOTLIB FIGURE AND AXIS
        self.fig, self.ax = plt.subplots(figsize=(10,6))
        self.confirmed = False
        self.current_ends = []
        # CLEAR PREVIOUS LINES IF ANY
        self.ax.clear()
        self.lines = []
        
        self.vertical_lines = []
        counter = 0
        for event_interpretation in event_interpretations:
            
            # duration_seconds, duration_timedelta, x, y, list_name, index, energy, observation, instrument, integral = event_interpretation
            duration_seconds, duration_timedelta, x, y = event_interpretation
            print('event info')
            print(duration_seconds, duration_timedelta) #, x, y)
            x = x.tolist()
            y = y.tolist()

            # PLOT MULTIPLE LINES WITH DIFFERENT STYLES
            self.lines.append(self.ax.plot(x, y, label='line', linewidth=1, zorder=10))
          

            # PLOT VERTICAL LINES AT THE START AND END OF EACH EVENT INTERPRETATION
            self.vertical_lines.append((self.ax.axvline(x=x[0], color=config.color.matplotlib_color_cycle[counter], linewidth=2, linestyle='--'), self.ax.axvline(x=x[-1], color=config.color.matplotlib_color_cycle[counter], linewidth=2, linestyle='--')))
            
            print('counter loop', len(self.confirmed_times), counter, self.current_plot_index)
            if len(self.confirmed_times) == self.current_plot_index:
                self.confirmed_times.append(pd.NaT)
                print(self.confirmed_times)
            
            counter += 1

        
        # ENABLE GRIDLINES
        self.ax.grid(True)
       
        # FORMAT THE x-AXIS LABELS
        self.ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%Y-%m-%d %H:%M'))

        # ROTATE THE x-AXIS LABELS
        self.ax.tick_params(axis='x', rotation=45, labelright=True)
 
        # ADD A LEGEND
        self.legend = self.ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1.0))

        # ADJUST THE LAYOUT TO PREVENT CLIPPING
        self.fig.tight_layout(rect=[0.05, 0, 0.95, 1.0])

        
        print('create plot', self.current_plot_index, len(self.confirmed_times), self.confirmed_times[self.current_plot_index])            
        try:
            if pd.isnull(all(self.confirmed_times[self.current_plot_index])):
                pass
            else:
                ct = [pd.to_datetime(x) for x in self.confirmed_times[self.current_plot_index]]
                for j in range(len(ct)):
                     
                    self.confirmed_line(ct[j])
        except:
            if pd.isnull(self.confirmed_times[self.current_plot_index]):
                pass
            else:
                self.confirmed_line(pd.to_datetime(self.confirmed_times[self.current_plot_index]))
         
        
        # ADD x LABEL
        self.ax.set_xlabel('UTC')
        
        # ADD y LABEL
        # if integral:
        self.ax.set_ylabel('Integral Proton Flux [cm$^\\mathregular{-2}$ sr$^\\mathregular{-1}$ s$^\\mathregular{-1}$]')
        # else:
        #     self.ax.set_ylabel('Differential Proton Flux [cm$^\\mathregular{-2}$ sr$^\\mathregular{-1}$ s$^\\mathregular{-1}$ MeV$^\\mathregular{-1}$]')

        # EMBED THE PLOT IN THE TKINTER WINDOW
        if hasattr(self, 'canvas'):
            self.canvas.get_tk_widget().destroy()  # REMOVE THE OLD CANVAS IF IT EXISTS
        self.canvas = matplotlib.backends.backend_tkagg.FigureCanvasTkAgg(self.fig, master=self.plot_frame)
        self.canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=True)
        # BIND CLICK EVENT TO PLOT
        self.lines = []
        self.selected_line = Line()
        self.plot_click = self.fig.canvas.mpl_connect('button_press_event', self.on_click)

        if hasattr(self, 'toolbar'):
            self.toolbar.destroy()
            self.toolbar = matplotlib.backends.backend_tkagg.NavigationToolbar2Tk(self.canvas, self.plot_frame)
            self.toolbar.update()
            self.canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=True)

        # SET Y-AXIS SCALE BASED ON CURRENT STATE
        if self.yaxis_log_scale:
            self.ax.set_yscale('log')
        else:
            self.ax.set_yscale('linear')

        # MAKE LEGEND ITEMS CLICKABLE
        for i, (legend_line, legend_text) in enumerate(zip(self.legend.get_lines(), self.legend.get_texts())):
            legend_line.set_picker(True)
            legend_text.set_picker(True)
            legend_line.set_gid(i)  # ASSIGN AN ID TO IDENTIFY WHICH LINE CORRESPONDS TO WHICH LEGEND ITEM

        # CONNECT MATPLOTLIB EVENTS FOR THE NEW PLOT
        self.fig.canvas.mpl_connect('pick_event', self.on_legend_click)
        self.fig.canvas.mpl_connect('motion_notify_event', self.on_hover)
        self.canvas.draw()
        
        # UPDATE THE NAVIGATION BUTTON STATES AFTER CREATING THE PLOT
        self.update_navigation_buttons()

    def update_table(self):
        """UPDATE THE TABLE BASED ON THE CURRENT PLOT DATA."""
        # CLEAR EXISTING DATA
        for row in self.table.get_children():
            self.table.delete(row)

        # ADD NEW DATA FOR THE CURRENT PLOT
        print('update table')
        # print(self.events)
        print(self.current_plot_index)
        event = self.events.iloc[self.current_plot_index]
        print(event['duration'], type(event['duration']))
        print(event)
        print(type(event['start']))
        try:
            start_table = event['start'].strftime('%Y-%m-%d %H:%M')
            end_table = event['end'].strftime('%Y-%m-%d %H:%M')
            duration_table = event['end'].strftime('%Y-%m-%d %H:%M') - event['start'].strftime('%Y-%m-%d %H:%M')
            # output = [event['start'].strftime('%Y-%m-%d %H:%M'), event['end'].strftime('%Y-%m-%d %H:%M'), '{:.1f}'.format(event['duration'].total_seconds() / 3600), event['list'], event['energy']]
        except:
            start_table = event['start']
            end_table = event['end']
        duration_table = pd.to_datetime(end_table) - pd.to_datetime(start_table)
        output = [start_table, end_table, '{:.1f}'.format(duration_table.total_seconds() / 3600), event['list'], event['energy']]
        self.table.insert("", "end", values=output, tag=(event['list'],))

        self.highlight_rows()
 
    def highlight_rows(self):
        counter = 0
        for line in self.lines:
            label = line[0].get_label()
            self.table.tag_configure(label, foreground=config.color.matplotlib_color_cycle[counter])
            counter += 1 

    def previous_plot(self):
        """NAVIGATE TO THE PREVIOUS PLOT."""    
        if self.current_plot_index > 0:
            self.current_plot_index -= 1
            self.create_plot()
            self.update_table()

    def next_plot(self):
        """NAVIGATE TO THE NEXT PLOT."""
        if self.current_plot_index < len(self.stored_event_plots) - 1:
            self.current_plot_index += 1
            self.multiple_end = False
            self.create_plot()
            self.update_table()

    def reset_manual_lines(self):
        """NAVIGATE TO THE NEXT PLOT."""
        new_win = errorwindow()
        
        if new_win.return_status.get() == "True":
            self.current_ends = []
            self.confirmed_times[self.current_plot_index] = pd.NaT
            self.confirmed = False
            self.create_plot()
            self.update_table()

    def force_quit(self):
        """FOR WHEN SHIT DOESNT CLOSE (SAVES YOUR PROGRESS SO YOU SHOULDNT FORCE QUIT BEFORE YOURE READY)"""
        self.root.quit()                      # STOP THE MAIN TKINTER LOOP
        self.root.destroy()                   # PROPERLY DESTROY THE TKINTER WINDOW


    def update_navigation_buttons(self):
        """ENABLE OR DISABLE NAVIGATION BUTTONS BASED ON THE CURRENT PLOT INDEX."""
        if self.current_plot_index > 0:
            self.prev_button.config(state=tkinter.NORMAL)
        else:
            self.prev_button.config(state=tkinter.DISABLED)
        if self.current_plot_index < len(self.stored_event_plots) - 1:
            self.next_button.config(state=tkinter.NORMAL)
        else:
            self.next_button.config(state=tkinter.DISABLED)            

    def on_legend_click(self, event):
        """HANDLES CLICKS ON THE LEGEND TO TOGGLE LINE VISIBILITY."""
        if isinstance(event.artist, matplotlib.lines.Line2D):
            # TOGGLE THE VISIBILITY OF THE CORRESPONDING LINE
            for line in self.lines:
                if line == event.artist:
                    visible = not line.get_visible()
                    line.set_visible(visible)
                    event.artist.set_alpha(1.0 if visible else 0.2)
                    self.fig.canvas.draw()
                    break

        # CHECK IF A LEGEND TEXT WAS CLICKED
        for i, text in enumerate(self.legend.get_texts()):
            if text == event.artist:
                line = self.lines[i][0]
                visible = not line.get_visible()
                line.set_visible(visible)
                self.legend.get_lines()[i].set_alpha(1.0 if visible else 0.2)
                text.set_alpha(1.0 if visible else 0.2)
                self.fig.canvas.draw()
                break

    def on_hover(self, event):
        """HIGHLIGHTS THE LEGEND TEXT WHEN HOVERING OVER IT."""
        for text in self.legend.get_texts():
            is_hovering = text.contains(event)[0]
            if is_hovering:
                text.set_fontweight("bold")
                text.set_color("blue")
            else:
                text.set_fontweight("normal")
                text.set_color("black")
        self.fig.canvas.draw_idle()

    def add_another_end_time(self):
        """
        For When One End Time is Not Enough
        
        """
        self.multiple_end = True
    
    def toggle_yaxis_scale(self):
        """TOGGLES THE Y-AXIS SCALE BETWEEN LINEAR AND LOGARITHMIC."""
        self.yaxis_log_scale = not self.yaxis_log_scale
        self.ax.set_yscale('log' if self.yaxis_log_scale else 'linear')
        self.canvas.draw()

    def confirm_end(self):
        """Confirms your chosen end time"""
        print(self.current_plot_index)
        self.current_ends.append(self.temp_end)
        self.confirmed_times[self.current_plot_index] = self.current_ends
        self.confirmed = True
        self.confirmed_line(pd.to_datetime(self.temp_end))
        self.multiple_end = False

    def on_close(self):
        """CLEANUP AND CLOSE THE APP SAFELY."""
        
        self.root.quit()                      # STOP THE MAIN TKINTER LOOP
        self.root.destroy()                   # PROPERLY DESTROY THE TKINTER WINDOW
        print("APPLICATION CLOSED CLEANLY.")  # CONFIRMATION IN THE CONSOLE
        

    def on_click(self, event):
        if event.inaxes is not None:
            x = event.xdata
            dt = datetime_given_days_since_epoch(x)
            self.edit.delete(0, tkinter.END)
            self.edit.insert(0, dt.strftime('%Y-%m-%d %H:%M:%S'))
            
            self.edit.focus_set()
            self.edit.select_range(0, tkinter.END)
            if self.selected_line.active:
                # MOVE LINE, REMAINS ACTIVE
                self.selected_line.remove_line()
                self.temp_end = 0
              
            self.update_vertical_line(x)
            self.temp_end = dt.strftime('%Y-%m-%d %H:%M:%S')
            if self.selected_line.active: 
                self.edit.config(state='normal')
        
    def update_vertical_line(self, x):
        try:
            for i in range(len(x)):
                self.selected_line.line = self.ax.axvline(x=x[i], color='black', linestyle='--', linewidth=2)
                self.selected_line.glow = self.ax.axvline(x=x[i], color='yellow', linewidth=4, alpha=0.5)
                self.selected_line.active = True
                self.canvas.draw()
        except:
            self.selected_line.line = self.ax.axvline(x=x, color='black', linestyle='--', linewidth=2)
            self.selected_line.glow = self.ax.axvline(x=x, color='yellow', linewidth=4, alpha=0.5)
            self.selected_line.active = True
            self.canvas.draw()

    def confirmed_line(self, x):
        # self.selected_line.remove_line()
        self.edit.insert(0, x.strftime('%Y-%m-%d %H:%M:%S'))
        self.selected_line.line = self.ax.axvline(x=x, color='red', linestyle='--', linewidth=2)
        self.selected_line.glow = self.ax.axvline(x=x, color='black', linewidth=4, alpha=0.5)
        self.selected_line.active = False
        self.canvas.draw()

    def update_vertical_line_text(self, event=None):
        x = days_since_epoch(datetime.datetime.strptime(self.edit.get(), '%Y-%m-%d %H:%M:%S', ))
        self.selected_line.remove_line()
        self.selected_line.active = True
        self.update_vertical_line(x)
        
    def update_xticks(self):
        if self.ax:
            xmin, xmax = self.ax.get_xlim()
            gridlines = self.ax.get_xgridlines()
            grid_positions = [line.get_xdata()[0] for line in gridlines if line.get_xdata()]
            if not grid_positions:
                grid_positions = np.linspace(xmin, xmax, num=10)
            self.ax.set_xticks(grid_positions)
            self.ax.set_xticklabels([datetime_given_days_since_epoch(tick).strftime('%Y-%m-%d %H:%M:%S') for tick in grid_positions])
            self.canvas.draw()
        
    def __call__(self):
        return self.end_times, self.confirmed_times
    



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='GUI Program for Comparing SEP Event Definitions')
    parser.add_argument('-el', '--event-list',  type=str,   help='Full Filepath to the Event List Output from IDSEP')
    parser.add_argument('-ff',  '--flux-files', type=str,   help='Full Filepath to the Flux Files, can be an array of multiple files')
    parser.add_argument('-tb',  '--time-buffer',           type=float, help='Fractional number of days to use as a buffer; extends time range of proton flux displayed on plot.', default=1.0)
    parser.add_argument('-di',  '--instrument',            type=str,   help='Select proton measurement instrument.',                                                              default=None)
    parser.add_argument('-se',  '--separate-energies',                 help='Separates fluxes of different energy channels; they are not shown on the same plot.', action='store_true')
    parser.add_argument('-t',   '--test',                              help='Testing mode.',                                                                       action='store_true')
    args = parser.parse_args()


    # if args.test:
    #     # args_label_file = 'test/' #os.path.join('test', 'label.csv')
    #     # args_event_list_directory = 'test/'
    #      # os.path.join('test', 'event_lists')
    #     # args_observation_directory = 'test/' #os.path.join('test', 'observations')
    #     args_time_buffer = 2.0
    #     args_instrument = None
    #     args_separate_energies = False
    # else:
    
    args_event_list_directory = args.event_list
    args_observation_directory = args.flux_files
    args_time_buffer = args.time_buffer
    args_instrument = args.instrument
    args_separate_energies = args.separate_energies

    # label_df = pd.read_csv(args_label_file)
    event_lists = []
    observations = []
    instruments = []
    event_labels = {}
    # event_list_files = label_df['event_list_filename'].to_list()
    # observation_files = label_df['observation_filename'].to_list()
    # instruments = label_df['instrument'].to_list()
    # event_list_labels = label_df['label'].to_list()     
    event_list_files = ['./test/batch_event_list_GOES-07_integral.txt']
    observation_files = ['./test/fluxes_GOES-07_integral_19870301_19960831.csv']
    instruments = ['GOES']
    event_list_labels = ['GOES']
    counter = 0
    for f, o, instrument, l in zip(event_list_files, observation_files, instruments, event_list_labels):
        if f.strip()[0] == '#':
            continue
        event_list_filename = f #os.path.join(args_event_list_directory, f)
        observation_filename = o #os.path.join(args_observation_directory, o)
        
        event_lists.append(event_list_filename)
        observations.append(observation_filename)
        instruments.append(instrument)
        event_labels[counter] = l
        counter += 1

    counter = 0
    event_list_dfs = []
    energy_channel = '10.0--1 MeV'
    
    for (event_list, observation, instrument) in zip(event_lists, observations, instruments):
        event_list_df = format_fetchsep_list(event_list, observation=observation, instrument=instrument, index=counter, energy_channel = energy_channel)
        event_list_dfs.append(event_list_df)
        counter += 1

    events = find_overlapping_events(event_list_dfs, event_labels, args_separate_energies)
    print(events)
    
    # events = pd.read_csv('batch_event_list_GOES-07_integral.txt')
    # INITIALIZE TKINTER ROOT WINDOW
    initial_ends, confirmed_ends, output_file = checksep(events, observations, energy_channel, event_list_files[0])
    #  checksep(events, flux_files)
    # root = tkinter.Tk()
    # geometry_string = get_window_geometry(root)
    # root.geometry(geometry_string)
    # app = CheckSEPApp(root, events)
    # root.mainloop()
    # initial_ends, confirmed_ends = app()
    # print(initial_ends)
    # elements_to_remove = {pd.NaT}
   
    
    # Now that the end times have be selected - we now build the event
    # list to feed into opsep
    



