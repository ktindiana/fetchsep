History
=======

0.2.0 (2023-09-22)
------------------

 * Fix packaging bug that failed to include the default configuration

0.2.0 (2023-09-22)
------------------

 * Add GOES-18 differential and real time integral fluxes; changes
   GOES-R integral flux source to CCMC's HAPI API
 * Add configuration by user config file in the configparser INI
   format, replacing placement of config values in config.py.  This
   fixes the pip deployment issue whereby a user could not easily
   update the code configuration.
 * Add --ExportConfig option to opsep & idsep.  Users may use this to
   prepare a custom configuration.
 * Fix: Create necessary output directories instead of crashing if
   they don't exist

0.1.4 (2023-08-22)
------------------

* Fix bug preventing command-line scripts from executing in Python interpreter

0.1.3 (2023-08-22)
------------------

* Update utils/read_datasets.py to add new version, v3-0-1, for GOES-R differential files
* Fix bug in ccmc_json_handler.py that would crash for model profiles with 0 values
* Fix deployment of command-line scripts
* Update README documentation

0.1.2 (2023-08-22)
------------------

* Fix README for PyPI distribution
* Update developer package requirements

0.1.1 (2023-08-08)
------------------

* Test bump2version procedure

0.1.0 (2023-06-01)
------------------

* First release on PyPI.
