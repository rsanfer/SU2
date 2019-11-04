#!/usr/bin/env python3

## \file fsi_computation.py
#  \brief Python wrapper code for NITRO_FSI. It launches SU2 for a steady and 2 unsteady simulations in series.
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
from libFSI import FSI_config as io  # imports FSI python tools
from libFSI import Interface as FSI # imports FSI python tools

# -------------------------------------------------------------------
#  Main
# -------------------------------------------------------------------

def main():

   # --- Get the FSI conig file name form the command line options --- #
   parser=OptionParser()
   parser.add_option("-c", "--conf_file",       dest="conf_file",
                      help="read config from conf_file", metavar="CONF_FILE")
   parser.add_option("-e", "--echo",       dest="echo",
                      help="write ECHO file", metavar="ECHO") 
   parser.add_option("-n", "--partitions", dest="partitions", default=1,
                      help="number of PARTITIONS", metavar="PARTITIONS")       

   # NUmber of blended step times simulated in simulation 1
   tq_mult = 3

   # imputs handling
   (options, args)=parser.parse_args()
   options.partitions  = int( options.partitions )
   # THE NAME OF THE FILE IS GIVEN: 'FSICoupler_config.cfg'
   pos = options.conf_file.find('.',1,-1)
   # Identifying steady and unsteady file
   filename_unsteady = options.conf_file
   filename_steady = options.conf_file[0:pos] + '_steady' + options.conf_file[pos:]
  
   # I get something from the unsteady conf file (ONLY!!!!) like the time scheme and the output directory where to operate
   FSI_config = io.FSIConfig(filename_unsteady)  # FSI configuration file
  
   # I create the output folder in which all the results are going to be saved
   try:
      os.system('mkdir ' + FSI_config['OUTPUT_DIRECTORY'])
   except:
      pass
  
  
   # temporary file which contains the bash commands (to be executed from shell)
   # it is positioned in the input file folder
   pos = filename_unsteady.find('FSICoupler_config.cfg')
   Command_file = filename_unsteady[0:pos] + 'Command_steady_unsteady.txt'
  
   # Open command file
   temp = open(Command_file , "w")
  
   # Command From steady simulation
   string1 = 'mpirun -np ' +  str(options.partitions) + ' run_fsi.py -f ' + filename_steady + ' --parallel ... | tee ' + FSI_config['OUTPUT_DIRECTORY']  + '/' + options.echo   # ' ../fsi_computation_imposed.py -f '
   temp.write(string1 + '\n')
  

   # Rename the restart file in case of following unsteady
   # I need to know if the dual time method is single step or dual
   # FSI configuration file
   if FSI_config['UNSTEADY_SCHEME'] == 'DUAL_TIME_STEPPING-1ST_ORDER':
     #string2 = 'mv ./' + options.output_folder + '/restart_flow.dat   '  + './' + options.output_folder + '/restart_flow_00000.dat' # this is for my local
     string2 = 'mv ' +  'zrestart_flow.dat   '   + 'zrestart_flow_00000.dat'
   else:
     string2 = '';
     print('ERROR: Chosen time Advancing scheme not implemented yet')
     sys.exit("Good bye!")
    
   temp.write(string2 + '\n')
      
   # Command to execute the following unsteady simulation
  
   for i in range(1, FSI_config['UNST_TOTAL_SIMUL_NUMBER'] +1  ):
     if i==1:
       index = ''
     else:
         index = str(i)
     string3 = 'mpirun -np ' +  str(options.partitions) + ' run_fsi.py -f ' + filename_unsteady[0:-4] + index + '.cfg' + ' --parallel ... | tee ' + FSI_config['OUTPUT_DIRECTORY'] + '/' + options.echo + '_unsteady' + str(i)
     string3_bis= 'mv ' + 'zrestart_flow_' +str(tq_mult*int(FSI_config['BS_TIMESTEP_1'])-1).zfill(5) + '.dat   '  + 'zrestart_flow_' + str(tq_mult*int(FSI_config['BS_TIMESTEP_2'])-1).zfill(5) + '.dat'
     temp.write(string3 + '\n')
     if i != FSI_config['UNST_TOTAL_SIMUL_NUMBER']:
        temp.write(string3_bis + '\n')
  
   temp.close()

   command = 'bash   ' + Command_file
  
   #command_del = 'rm ' + Command_file
  
   os.system (command)
  
   #os.system (command_del)









# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# --- This is only accessed if running from command prompt --- #
if __name__ == '__main__':
    main()
