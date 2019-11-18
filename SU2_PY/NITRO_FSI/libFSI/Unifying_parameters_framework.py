#!/usr/bin/env python3

## \file UnifyingParameters_framework.py
#  \brief Function that unifies parameters for the input files FSI and SU2
#  \author Rocco Bombardieri
#  \version 7.0.0
#
# The current SU2 release has been coordinated by the
# SU2 International Developers Society <www.su2devsociety.org>
# with selected contributions from the open-source community.
#
# The main research teams contributing to the current release are:
#  - Prof. Juan J. Alonso's group at Stanford University.
#  - Prof. Piero Colonna's group at Delft University of Technology.
#  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
#  - Prof. Rafael Palacios' group at Imperial College London.
#  - Prof. Vincent Terrapon's group at the University of Liege.
#  - Prof. Edwin van der Weide's group at the University of Twente.
#  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
#
# Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
#                      Tim Albring, and the SU2 contributors.
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

from math import *	# use mathematical expressions
import os, sys, shutil, copy
import fileinput


def UnifyFSIConfig(DYN_config,confFile):
    '''This routine reads some info from the dynresp input file and unifies cards in SU2 and FSICoupler files'''

    configfile2 = open(confFile + '_temp',"w")
    with open(confFile, 'r') as configfile:
      while 1:
        line = configfile.readline()
        string = line
        if not line:
          break

        # remove line returns
        line = line.strip('\r\n')
        # make sure it has useful data
        if (not "=" in line) or (line[0] == '%'):
         configfile2.write(string)
        else:
         # split across equal sign
         line = line.split("=",1)
         this_param = line[0].strip()
         this_value = line[1].strip()

         #float values
         if this_param == "V_INF":
                    stringalt = 'V_INF = '+ str(DYN_config['V']) + '   \r\n'
                    configfile2.write(stringalt)
         if this_param == "FREESTREAM_DENSITY":
                    stringalt = 'FREESTREAM_DENSITY = '+ str(DYN_config['RHO']) + '   \r\n'
                    configfile2.write(stringalt)

         #float values
         if this_param == "START_TIME":
                    stringalt = 'START_TIME = '+ str(DYN_config['TW0']) + '   \r\n'
                    configfile2.write(stringalt)
         if this_param == "UNST_TIME":
                    stringalt = 'UNST_TIME = '+ str(DYN_config['TWF']) + '   \r\n'
                    configfile2.write(stringalt)
         if this_param == "UNST_TIMESTEP":
                    stringalt = 'UNST_TIMESTEP = '+ str(DYN_config['DT']) + '   \r\n'
                    configfile2.write(stringalt)

         else:
             configfile2.write(string)

    configfile.close()
    configfile2.close()
    # the file is now replaced
    os.remove(confFile)
    os.rename(confFile + '_temp', confFile)

def UnifyingParameters_framework(FSI_config,confFile,myid ):
          
    rootProcess = 0
    
    #FREESTREAM_PRESSURE_default_SU2= 101325.0 # [ N/m^2]
    FREESTREAM_TEMPERATURE_default_SU2= 288.15 # [K]
    GAMMA_VALUE_default_SU2= 1.4
    GAS_CONSTANT_default_SU2= 287.058 # [J/kg*K]
    if FSI_config['CSD_SOLVER'] == 'NITRO_FRAMEWORK':
        sound_speed = sqrt(GAMMA_VALUE_default_SU2*FSI_config['FREESTREAM_PRESSURE']/FSI_config['FREESTREAM_DENSITY'])
    else:    
       sound_speed = sqrt(GAMMA_VALUE_default_SU2*GAS_CONSTANT_default_SU2*FREESTREAM_TEMPERATURE_default_SU2) #[m/s]
    #    V_12inf = sound_speed*Mach_inf_SU2;
    MACH_NUMBER = FSI_config['V_INF'] / sound_speed
    
    tq_mult = 3
    
    # Restart options
    
    if FSI_config['UNSTEADY_SIMULATION'] == 'NO':
       FSI_config['RESTART_SOL'] = 'NO'
       if myid == rootProcess:
          print("Simulation {}, FSI conf file {}, SU2 conf file {}, MLS conf file {}\n".format('STEADY', confFile,FSI_config['CFD_CONFIG_FILE_NAME'],FSI_config['MLS_CONFIG_FILE_NAME'] ))

          
    if FSI_config['MOTION_TYPE'] == 'BLENDED_STEP' and FSI_config['UNSTEADY_SIMULATION'] == 'YES':                                
       # Unsteady timestep suggestion (BLENDED_STEP type)    
       t_q = 3.14*FSI_config['L_REF']/ FSI_config['V_INF']/(FSI_config['K_MAX']/2)
       FSI_config['BLE_STEP_LENGTH'] = t_q
       #if myid == rootProcess: print("FSI_config['BLE_STEP_LENGTH'] = {}".format(FSI_config['BLE_STEP_LENGTH']))
       if FSI_config['UNST_NR'] == 1:
          FSI_config['UNST_TIMESTEP'] = t_q/(FSI_config['BS_TIMESTEP_1'])
          if myid == rootProcess: print("HC: FSI_config['UNST_TIMESTEP'] = {}".format(FSI_config['UNST_TIMESTEP']))
       elif FSI_config['UNST_NR'] == 2:   
          FSI_config['UNST_TIMESTEP'] = t_q/(FSI_config['BS_TIMESTEP_2'])
       else:   
          print('ERROR: More than two simulations are not implemented yet: change UNST_NR and/or UNST_TOTAL_SIMUL_NUMBER ')
          sys.exit("Good bye!")       
       
       # If we have unsteady and blended step and restart and second order and restart iter =2 then we start the movement at the first timestep 
       if FSI_config['RESTART_SOL'] == 'YES' and FSI_config['UNSTEADY_SCHEME'] == 'DUAL_TIME_STEPPING-2ND_ORDER':
          # start motion time for blended step depends on the dual time method chosen (1st-2nd)  
          FSI_config['START_MOTION_TIME'] = FSI_config['UNST_TIMESTEP']*1    
            
    # we have to play in two different time advancing lines cause the timestep is different  
       if FSI_config['RESTART_SOL'] == 'YES' and FSI_config['UNST_NR'] == 1: 
          FSI_config['START_TIME'] = FSI_config['UNST_TIMESTEP']*(FSI_config['RESTART_ITER']); # restart iter is chosen manually 
          FSI_config['UNST_TIME'] = tq_mult*FSI_config['BLE_STEP_LENGTH'] - FSI_config['START_TIME']
          if myid == rootProcess:
             print("Simulation number {}, FSI conf file {}, SU2 conf file {}, MLS conf file {}\n".format(FSI_config['UNST_NR'], confFile,FSI_config['CFD_CONFIG_FILE_NAME'],FSI_config['MLS_CONFIG_FILE_NAME'] ))
             print("HC: FSI_config['START_TIME'] = {}".format(FSI_config['START_TIME']))
             print("HC: FSI_config['UNST_TIME'] = {}".format(FSI_config['UNST_TIME']))
          
       elif FSI_config['RESTART_SOL'] == 'YES' and FSI_config['UNST_NR'] == 2:
          FSI_config['START_TIME'] = tq_mult*FSI_config['BLE_STEP_LENGTH'] + FSI_config['UNST_TIMESTEP'] ; 
          FSI_config['RESTART_ITER'] = int(tq_mult*FSI_config['BS_TIMESTEP_2'] )  #(the total number of the old simulation is tq_mult*FSI_config['BS_TIMESTEP_2'] (the first timetep is always 0))
          if myid == rootProcess:
             print("Simulation number {}, FSI conf file {}, SU2 conf file {}, MLS conf file {}".format(FSI_config['UNST_NR'], confFile,FSI_config['CFD_CONFIG_FILE_NAME'],FSI_config['MLS_CONFIG_FILE_NAME'] ))
             print("HC: FSI_config['START_TIME'] = {}".format(FSI_config['START_TIME']))
             print("HC: FSI_config['RESTART_ITER'] = {}".format(FSI_config['RESTART_ITER']))

          # unsteady time is given manually and corresponds to the end of the simulation
       else:   
          print('ERROR: More than two simulations are not implemented yet: change UNST_NR# and/or UNST_TOTAL_SIMUL_NUMBER ')
          sys.exit("Good bye!")          
           
           
    if myid == rootProcess and FSI_config['UNSTEADY_SIMULATION'] == 'YES' and FSI_config['MOTION_TYPE'] == 'BLENDED_STEP':
       print("\nBlended step length (plus Start Motion Time, in case): {} [sec]".format(t_q + FSI_config['START_MOTION_TIME']))
       print("\nBlended step timesteps nr. (starting from 0): {} [-]".format(( float(FSI_config['BLE_STEP_LENGTH']) )/(float(FSI_config['UNST_TIMESTEP']))))

    if myid == rootProcess: 
       UnifyFluid(FSI_config, FREESTREAM_TEMPERATURE_default_SU2, GAMMA_VALUE_default_SU2, GAS_CONSTANT_default_SU2, MACH_NUMBER)
       #UnifyStructure(FSI_config)
    
def UnifyFluid(FSI_config, FREESTREAM_TEMPERATURE_default_SU2, GAMMA_VALUE_default_SU2, GAS_CONSTANT_default_SU2, MACH_NUMBER):   
    
    
    configfile2 = open(FSI_config['CFD_CONFIG_FILE_NAME'] + '_temp',"w")
    with open(FSI_config['CFD_CONFIG_FILE_NAME'], 'r') as configfile:
      while 1:          
        line = configfile.readline()
        string = line
        if not line:
          break

        # remove line returns
        line = line.strip('\r\n')
        # make sure it has useful data
        if (not "=" in line) or (line[0] == '%'):  
         configfile2.write(string)
        else: 
         # split across equal sign
         line = line.split("=",1)
         this_param = line[0].strip()
         this_value = line[1].strip()

         #float values
         if this_param == "MACH_NUMBER":
                    stringalt = 'MACH_NUMBER = '+ str(MACH_NUMBER) + '   \r\n'
                    configfile2.write(stringalt)                            
         elif this_param == "MACH_MOTION":
             if str(FSI_config['UNSTEADY_SIMULATION']) == 'YES':
                 stringalt = 'MACH_MOTION = ' + str(MACH_NUMBER) + '   \r\n'
                 configfile2.write(stringalt)
         elif this_param == "FREESTREAM_TEMPERATURE":
                    stringalt = 'FREESTREAM_TEMPERATURE = '+ str(FREESTREAM_TEMPERATURE_default_SU2) + '   \r\n'
                    configfile2.write(stringalt)
         elif this_param == "GAMMA_VALUE":
                    stringalt = 'GAMMA_VALUE = '+ str(GAMMA_VALUE_default_SU2) + '   \n'
                    configfile2.write(stringalt)
         elif this_param == "GAS_CONSTANT":
                    stringalt = 'GAS_CONSTANT = '+ str(GAS_CONSTANT_default_SU2) + '   \n'
                    configfile2.write(stringalt)
         elif this_param == "TIME_STEP":
                    stringalt = 'TIME_STEP = '+ str(FSI_config['UNST_TIMESTEP']) + '   \n'
                    configfile2.write(stringalt)
         elif this_param == "MAX_TIME":
                    stringalt = 'MAX_TIME = '+ str(FSI_config['UNST_TIME']) + '   \n'
                    configfile2.write(stringalt)
         elif this_param == "UNST_RESTART_ITER":
                    stringalt = 'UNST_RESTART_ITER = '+ str(FSI_config['RESTART_ITER']) + '   \n'
                    configfile2.write(stringalt)
         elif this_param == "TIME_ITER" and FSI_config['UNST_NR'] == 1:
                    stringalt = 'TIME_ITER = ' + str(int(FSI_config['BS_TIMESTEP_1'] -1)) + '   \n'
                    configfile2.write(stringalt)
         elif this_param == "TIME_ITER" and FSI_config['UNST_NR'] == 2:
                    stringalt = 'TIME_ITER = ' + str(int(FSI_config['BS_TIMESTEP_2'] -1)) + '   \n'
                    configfile2.write(stringalt)
         elif this_param == "WRT_SOL_FREQ_DUALTIME": 
             if FSI_config['UNSTEADY_SIMULATION'] == 'YES': 
                #if   FSI_config['UNST_NR']  == 1:
                if  int(FSI_config['BS_TIMESTEP_1']/10) == 0: 
                   stringalt = 'WRT_SOL_FREQ_DUALTIME = '+ str(1) + '   \n'  
                else:
                   stringalt = 'WRT_SOL_FREQ_DUALTIME = '+ str(int(FSI_config['BS_TIMESTEP_1']/10)) + '   \n'  
                configfile2.write(stringalt)
                    
         #string values
         elif this_param == "TIME_MARCHING":
                    if str(FSI_config['UNSTEADY_SIMULATION']) == 'YES':
                       stringalt = 'TIME_MARCHING = '+ FSI_config['UNSTEADY_SCHEME'] + '   \n'
                       configfile2.write(stringalt)  
                    else:   
                       stringalt = 'TIME_MARCHING = '+ 'NO' + '   \n'
                       configfile2.write(stringalt)                     
         elif this_param == "GRID_MOVEMENT":
                    if str(FSI_config['UNSTEADY_SIMULATION']) == 'YES':
                       stringalt = 'GRID_MOVEMENT = '+ 'YES' + '   \n'                    
                       configfile2.write(stringalt)  
                    else:   
                       stringalt = 'GRID_MOVEMENT = '+ 'NO' + '   \n'                    
                       configfile2.write(stringalt)    
         elif this_param == "RESTART_SOL":
                    stringalt = 'RESTART_SOL = ' + str(FSI_config['RESTART_SOL']) + '   \n' 
                    configfile2.write(stringalt)
         #For the NITRO_FRAMEWORK approach it is important to memorize the output file in precise locations                     
         elif this_param == "SURFACE_FILENAME":
                    if FSI_config['CSD_SOLVER'] == 'NITRO_FRAMEWORK':
                       #stringalt = 'SURFACE_FILENAME = ' + FSI_config['OUTPUT_DIRECTORY'] + '/surface_flow'     + '   \n'
                       stringalt = 'SURFACE_FILENAME = '  + 'surface_flow'     + '   \n'
                       # I just want to memorize the surface-flow location since it is not possible to avoide the writing
                       #FSI_config._ConfigContent[this_param] = str(FSI_config['OUTPUT_DIRECTORY'] + '/surface_flow')
                       FSI_config._ConfigContent[this_param] =  'surface_flow'
                       configfile2.write(stringalt)
                    else: # 'NITRO'
                       # surf_flow_filename = this_value
                       FSI_config._ConfigContent[this_param] = this_value
                       configfile2.write(string)
         elif this_param == "VOLUME_FILENAME":
                    #stringalt = 'VOLUME_FILENAME = ' + FSI_config['OUTPUT_DIRECTORY'] + '/flow'     + '   \n'
                    stringalt = 'VOLUME_FILENAME = '  + 'flow'     + '   \n'
                    configfile2.write(stringalt)
         elif this_param == "RESTART_FLOW_FILENAME":                    
                    stringalt = 'RESTART_FLOW_FILENAME = '  + 'zrestart_flow.dat'  + '   \n'
                    configfile2.write(stringalt)           
         elif this_param == "BREAKDOWN_FILENAME":                    
                    stringalt = 'BREAKDOWN_FILENAME = ' + FSI_config['OUTPUT_DIRECTORY'] + '/forces_breakdown.dat'   + '   \n'                   
                    configfile2.write(stringalt)      
         elif this_param == "CONV_FILENAME":                    
                    stringalt = 'CONV_FILENAME = ' + FSI_config['OUTPUT_DIRECTORY'] + '/history'           + '   \n'           
                    configfile2.write(stringalt)        
         elif this_param == "SOLUTION_FILENAME":
                    stringalt = 'SOLUTION_FILENAME = '  + 'zrestart_flow.dat'    + '   \n'
                    configfile2.write(stringalt)         
         elif this_param == "MESH_OUT_FILENAME":                    
                    stringalt = 'MESH_OUT_FILENAME = ' + FSI_config['OUTPUT_DIRECTORY'] + '/mesh_out.su2'        + '   \n'              
                    configfile2.write(stringalt)     
         #For the NITRO_FRAMEWORK approach it is important to set correctly the pressure and density for the standard air                    
         elif this_param == "FREESTREAM_OPTION": 
                    if FSI_config['CSD_SOLVER'] == 'NITRO_FRAMEWORK': 
                       stringalt = 'FREESTREAM_OPTION = ' + FSI_config['FREESTREAM_OPTION']        + '   \n'              
                       configfile2.write(stringalt)  
         elif this_param == "FREESTREAM_PRESSURE":  
                    if FSI_config['CSD_SOLVER'] == 'NITRO_FRAMEWORK': 
                       stringalt = 'FREESTREAM_PRESSURE = ' + str(FSI_config['FREESTREAM_PRESSURE'])        + '   \n'              
                       configfile2.write(stringalt)                     
         elif this_param == "FREESTREAM_DENSITY":
                    if FSI_config['CSD_SOLVER'] == 'NITRO_FRAMEWORK': 
                       stringalt = 'FREESTREAM_DENSITY = ' + str(FSI_config['FREESTREAM_DENSITY'])        + '   \n'   
                       configfile2.write(stringalt)
         else:
                    configfile2.write(string) 

         
         
    configfile.close()    
    configfile2.close()
    # the file is now replaced
    os.remove(FSI_config['CFD_CONFIG_FILE_NAME'] )
    os.rename(FSI_config['CFD_CONFIG_FILE_NAME'] + '_temp', FSI_config['CFD_CONFIG_FILE_NAME'] )
        
        
def UnifyStructure(FSI_config):   
    
    configfile2 = open(FSI_config['CSD_CONFIG_FILE_NAME'] + '_temp',"w")
    with open(FSI_config['CSD_CONFIG_FILE_NAME'], 'r') as configfile:
      while 1:          
        line = configfile.readline()
        string = line
        if not line:
          break

        # remove line returns
        line = line.strip('\r\n')
        # make sure it has useful data
        if (not "=" in line) or (line[0] == '%'):  
         configfile2.write(string)
        else: 
         # split across equal sign
         line = line.split("=",1)
         this_param = line[0].strip()
         this_value = line[1].strip()

         #float values
         if str(FSI_config['UNSTEADY_SIMULATION']) == 'YES':
                if this_param == "UNSTEADY_SIMULATION":
                    stringalt = 'UNSTEADY_SIMULATION = '+ str(FSI_config['UNSTEADY_SIMULATION']) + '   \n'
                    configfile2.write(stringalt)                    
                #elif this_param == "DELTA_T":
                #    stringalt = 'DELTA_T = '+ str(FSI_config['UNST_TIMESTEP']) + '   \n'
                #    configfile2.write(stringalt)
                #elif this_param == "START_TIME":
                #    stringalt = 'START_TIME = '+ str(FSI_config['START_TIME']) + '   \n'                    
                #    configfile2.write(stringalt)        
                #elif this_param == "STOP_TIME":
                #    stringalt = 'STOP_TIME = '+ str(FSI_config['STOP_TIME']) + '   \n'                    
                #    configfile2.write(stringalt)    
                else:
                    configfile2.write(string)                    
         #string values
         elif str(FSI_config['UNSTEADY_SIMULATION']) == 'NO':
                if this_param == "UNSTEADY_SIMULATION":
                    stringalt = 'UNSTEADY_SIMULATION = '+ str(FSI_config['UNSTEADY_SIMULATION']) + '   \n'
                    configfile2.write(stringalt)  
                else:
                    configfile2.write(string) 

         
         
    configfile.close()    
    configfile2.close()
    # the file is now replaced
    os.remove(FSI_config['CSD_CONFIG_FILE_NAME'] )
    os.rename(FSI_config['CSD_CONFIG_FILE_NAME'] + '_temp', FSI_config['CSD_CONFIG_FILE_NAME'] )        