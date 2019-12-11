#!/usr/bin/env python3

## \file fsi_computation.py
#  \brief Orchestrator for Dynresp CFD coupled approach. It launches SU2 for a steady and 1 unsteady simulations in series. Each physical timestep SU2
#  \ reads modal displacements and delivers generalized forces
#  \author  Rocco Bombardieri (UC3M)
#  \version 7.0.0
#
# SU2 Original Developers: Dr. Francisco D. Palacios.
#                          Dr. Thomas D. Economon.
#
# SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
#                 Prof. Piero Colonna's group at Delft University of Technology.
#                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
#                 Prof. Rafael Palacios' group at Imperial College London.
#                 Prof. Edwin van der Weide's group at the University of Twente.
#                 Prof. Vincent Terrapon's group at the University of Liege.
#
# Copyright (C) 2012-2017 SU2, the open-source CFD code.
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os, sys, shutil, copy
import time as timer
from math import *	# use mathematical expressions
from optparse import OptionParser	# use a parser for configuration
import subprocess
# -------------------------------------------------------------------
#  Main
# -------------------------------------------------------------------

def main():

    # --- Get the FSI conig file name form the command line options --- #
    parser = OptionParser()
    parser.add_option("-c", "--conf_file", dest="conf_file",
		      help="read config from conf_file", metavar="CONF_FILE")
    parser.add_option("-e", "--echo", dest="echo",
		      help="write ECHO file", metavar="ECHO")
    parser.add_option("-d", "--dynresp", dest="dyn_filename",
		      help="read config from FILE", metavar="FILE")
    parser.add_option("-n", "--partitions", dest="partitions", default=1,
		      help="number of PARTITIONS", metavar="PARTITIONS")

    # imputs handling
    (options, args)=parser.parse_args()
    options.partitions  = int( options.partitions )
    # THE NAME OF THE FILE IS GIVEN: 'FSICoupler_config.cfg'
    pos = options.conf_file.find('.',1,-1)
    # Identifying dynresp file
    filename_dynresp = options.dyn_filename
    # Identifying steady and unsteady file
    filename_unsteady = options.conf_file


    #=== Assembly of the commands to be run in sequence
    # Dynresp
    # Format: dynresp_exe  dynresp_file
    # (thhe locatino needs to be changed according to the machine)
    dynresp_location = '/media/rocco/290CF0EB732D1122/Technion_project/Dynresp/Linux_Exe/bin/dynrespCFD'
    command_dynresp =  dynresp_locatison + '   ' + filename_dynresp

    # SU2
    # Format: fsi_steady_unsteady_dynresp.py -d GWTS_CFD_1_.inp -c FSICoupler_config.cfg -n 22 -e echo
    command_SU2 = 'fsi_steady_unsteady_dynresp.py ' + ' -d ' + filename_dynresp + ' -c ' + filename_unsteady + ' -n ' + str(options.partitions) + ' -e ' + options.echo

    try:
       # Now launching the programs
       print('Initializing Dynresp...')
       print(command_dynresp)
       proc1 = subprocess.Popen([command_dynresp],shell=True)
       #os.system(command_dynresp)
       timer.sleep(1)
       print('Initializing SU2...')
       print(command_SU2)
       proc2 = subprocess.Popen([command_SU2], shell=True)
       #os.system(command_SU2)
    except KeyboardInterrupt:
       proc1.terminate()  # <-- terminate the process 1
       proc2.terminate()  # <-- terminate the process 2


# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# --- This is only accessed if running from command prompt --- #
if __name__ == '__main__':
    main()