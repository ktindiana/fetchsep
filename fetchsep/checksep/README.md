# eventsphinx

This package is a compilation of tools used to view and compare interpretations of solar particle events.

## Dependencies
Install dependencies: `pip install -r requirements.txt`

# compare_event.py

## How to Run:

1) Make a storage directory (your choice on the name): `<data_directory>`
2) Place your SEP event lists in `<data_directory>/event_lists`.
3) Place your proton flux data in `<data_directory>/observations`.
4) Make a comma-separated file called `<data_directory>/label.csv` in the same format as the example included in the `test` directory. Specify the SEP event list files you care to display, alongside observation flux files and/or instrument specifications, and finally a label that will be used to identify the data in the plotting GUI.
5) Execute: `python compare_event.py --event-list-directory <data_directory>/event_lists --observation-directory <data_directory>/observations --label-file <data_directory>/label.csv --time_buffer <number_of_fractional_days>`

