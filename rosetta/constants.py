# A list of constants for the rosetta package.
import os, util

usrDir = os.path.expanduser('~')

# # # # #  Default files and paths:

options_filename = 'Options.ini'
output_log_filename = 'Output_Log.txt'
error_log_filename = 'Error_Log.txt'
saved_run_options_extension = '.rosetta_savefile'
docking_score_file = 'scores.fasc'
# Check that these work with unix.
searchPaths = ['.', '/', '..', '/Applications', '/Applications/Rosetta Apps', usrDir, os.path.join(usrDir,'Desktop')]

# # # # #  Default options:
default_general = [('rosetta_bundle',''), ('rosetta_database',''),
                   ('executables_dir',''), ('executables_suffix',''),
                   ('scoreView_path','')]
default_paths = [('docking_protocol',''), ('docking_prepack_protocol',''),
                 ('relax','')] # ('FloppyTail','')

default_outputDir = ''
default_keepTopDecoys = 500
default_numCPUs = max(util.determineNumCPUs() - 1, 1)

# # # # #  Minor settings:
count_completed_interval = 5 # Seconds between checking number of files finished.

