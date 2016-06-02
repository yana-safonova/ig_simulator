#!/usr/bin/env python

############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import sys
import getopt
import os
import logging
import shutil
import datetime
from time import gmtime, strftime
from matplotlib.font_manager import home

home_directory = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + '/'
ig_bin_directory = os.path.join(home_directory, "bin/")
python_src_directory = os.path.join(home_directory, "src/python_utils/")
repertoire_simulation_src_directory = os.path.join(home_directory, "src/ig_simulation/")

sys.path.append(python_src_directory)

class PathToBins:
    create_ideal_repertoire_tool = os.path.join(ig_bin_directory, "ideal_repertoire_constructor")
    paired_read_merger_tool = os.path.join(ig_bin_directory, "paired_read_merger")
    simulate_repertoire_tool = os.path.join(ig_bin_directory, "ig_simulator")

    art_illumina = os.path.join(home_directory, "bin/art_illumina")

    run_create_ideal_repertoire_tool = ig_bin_directory + "./ideal_repertoire_constructor"
    run_paired_read_merger_tool = ig_bin_directory + "./paired_read_merger"
    run_simulate_repertoire_tool = ig_bin_directory + "./ig_simulator"

    run_art_illumina = os.path.join(home_directory, "bin/art_illumina")

def PrintCommandLine(argv, log):
    command_line = " ".join([str(x) for x in argv] )
    log.info("Command line: "+ command_line)

def ReadConfig():
    if not os.path.exists(path_to_config_template):
        print("ERROR: config file " + path_to_config_template + " was not found")
    f = open(path_to_config_template, "r")
    config_params = dict()
    for line in f.readlines():
        splits = line.split()
        config_params[splits[0]] = splits[1]
    return config_params

def RunGrinder():
    config_params = ReadConfig()
    return config_params['path_to_grinder'] + "/script/./grinder"

def IgblastDirectory():
    config_params = ReadConfig()
    return config_params['path_to_igblast'] + "/"

def RunIgblast():
    config_params = ReadConfig()
    return config_params['path_to_igblast'] + "/bin/igblastn"

class Unbuffered:
   def __init__(self, stream):
           self.stream = stream
   def write(self, data):
           self.stream.write(data)
           self.stream.flush()
   def __getattr__(self, attr):
           return getattr(self.stream, attr)

def ErrorMsg(log):
    log.error("Something goes wrong. Please contact us and send .log file")
    sys.exit(1)

def AbnormalFinishMsg(log, program_name):
    log.info("Script " + program_name + " finished abnormally. Please contact us and send .log file")
