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
ig_bin_directory = os.path.join(home_directory, "bin/ig_tools/")
python_src_directory = os.path.join(home_directory, "src/python_utils/")
repertoire_simulation_src_directory = os.path.join(home_directory, "src/ig_simulation/")

sys.path.append(python_src_directory)

class PathToBins:
    create_ideal_repertoire_tool = os.path.join(ig_bin_directory, "create_ideal_repertoire")
    paired_read_merger_tool = os.path.join(ig_bin_directory, "paired_read_merger")
    simulate_repertoire_tool = os.path.join(ig_bin_directory, "ig_simulator")
    fastq_to_fasta_tool = os.path.join(ig_bin_directory, "fastq_to_fasta")
    merged_reads_stats_calc_tool = os.path.join(ig_bin_directory, "compute_merged_reads_stats")
    create_clusters_tool = os.path.join(ig_bin_directory, "create_clusters")
    singleton_gluer = os.path.join(ig_bin_directory, "singleton_gluer")
    filtered_reads_processer = os.path.join(ig_bin_directory, "add_filtered_reads")
    aa_verificator = os.path.join(ig_bin_directory, "verify_by_aa")

    hammer = os.path.join(spades_release_bins, "hammer")
    iterative_BH_tool = home_directory

    art_illumina = os.path.join(home_directory, "src/ig_tools/art_bin_VanillaIceCream/art_illumina")
    art_454 = os.path.join(home_directory, "src/ig_tools/art_bin_VanillaIceCream/art_454")

    run_create_ideal_repertoire_tool = ig_bin_directory + "./create_ideal_repertoire"
    run_paired_read_merger_tool = ig_bin_directory + "./paired_read_merger"
    run_repertoire_evaluator_tool = ig_bin_directory + "./repertoire_evaluator"
    run_simulate_repertoire_tool = ig_bin_directory + "./simulate_repertoire"
    run_split_paired_fastq_reads_tool = ig_bin_directory + "./split_paired_fastq_reads"
    run_fastq_to_fasta_tool = ig_bin_directory + "./fastq_to_fasta"
    run_merged_reads_stats_calc_tool = ig_bin_directory + "./compute_merged_reads_stats"
    run_create_clusters = ig_bin_directory + "./create_clusters"
    run_singleton_gluer = ig_bin_directory + './singleton_gluer'
    run_filtered_reads_processer = ig_bin_directory + "./add_filtered_reads"
    run_aa_verificator = ig_bin_directory + "./verify_by_aa"

    run_hammer = os.path.join(spades_release_bins, "hammer")
    run_iterative_BH = home_directory + './iterativeBH.py'
    run_changer_of_K = home_directory + './changer_of_K_bh.py'

    run_art_illumina = os.path.join(home_directory, "src/ig_tools/art_bin_VanillaIceCream/./art_illumina")
    run_art_454 = os.path.join(home_directory, "src/ig_tools/art_bin_VanillaIceCream/./art_454")

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
