#!/usr/bin/env python3

## \file NITRO_Tester.py
#  \brief NITRO Tester solver (for the NITRO approach involving forced moving boundary condition) used for testing the Py wrapper for external FSI coupling.
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

import os, sys, shutil, copy
from libFSI.util.switch import switch
import numpy as np
#import scipy as sp
#import scipy.linalg as linalg
from math import *
#from cmath import *
#from Read_decomp_coeff import *

# ----------------------------------------------------------------------
#  Config class
# ----------------------------------------------------------------------

class Point:
  """ Description. """

  def __init__(self):
    self.Coord0 = np.zeros((3,1)) # marks the initial position of the nodes: to this we superimpose the mode shape      
    self.Coord = np.zeros((3,1))
    self.Coord_prec = np.zeros((3,1))
    self.Vel = np.zeros((3,1))
    self.Force = np.zeros((3,1))
    
  def GetCoord0(self):
    return self.Coord0    

  def GetCoord(self):
    return self.Coord

  def GetCoord_prec(self):
    return self.Coord_prec

  def GetVel(self):
    return self.Vel

  def GetForce(self):
    return self.Force

  def SetCoord0(self, val_Coord):
    x, y, z = val_Coord
    self.Coord0[0] = x
    self.Coord0[1] = y
    self.Coord0[2] = z

  def SetCoord(self, val_Coord):
    x, y, z = val_Coord
    self.Coord[0] = x
    self.Coord[1] = y
    self.Coord[2] = z

  def SetCoord_prec(self, val_Coord):
    x, y, z = val_Coord
    self.Coord_prec[0] = x
    self.Coord_prec[1] = y
    self.Coord_prec[2] = z

  def SetVel(self, val_Vel):
    vx, vy, vz = val_Vel
    self.Vel[0] = vx
    self.Vel[1] = vy
    self.Vel[2] = vz

  def SetForce(self, val_Force):
    fx, fy, fz = val_Force
    self.Force[0] = fx
    self.Force[1] = fy
    self.Force[2] = fz

class NITRO:
  """Description"""

  def __init__(self, config_fileName = None):
    """ Description. """

    self.Config_file = config_fileName
    self.Config = {}

    print("\n------------------------------ Configuring the structural tester solver for FSI simulation ------------------------------")
    if config_fileName !=None:
       self.__readConfig()

    #self.Mesh_file = self.Config['MESH_FILE']
    #self.FSI_marker = self.Config['MOVING_MARKER']
    self.nDof = 3
    print("Structural model : wing for NITRO framework.")


    self.nDim= int()
    self.nElem = int()
    self.nVolumeFluidPoint = int()
    self.nPoint = 0
    self.nMarker = int()
    self.node = []
    self.markers = {}
    self.tracking_pos = []
    self.tracking_vel = []
    self.tracking_force_x =[]
    self.tracking_force_y =[]
    self.tracking_force_z =[]
    self.time = [];
    self.Flutter_mode_fluid_x = []
    self.Flutter_mode_fluid_y = []
    self.Flutter_mode_fluid_z = []
    self.Aq = float()  # Scaled Amplitude of the blended step 
    self.nModes = 0  # in case of CSD= NITRO/NITRO_FRAMEWORK helps understanding how many modes we are using for the computation
    #print("\n------------------------------ Reading the SU2 mesh ------------------------------")
    #self.__readSU2Mesh()
    self.GlobalCoordinates0 =[]

  def SetNodesProperties(self):
        for iPoint in range(0, self.nPoint):
          self.node.append(Point())
          self.node[iPoint].SetCoord((self.GlobalCoordinates0[iPoint][0],self.GlobalCoordinates0[iPoint][1], self.GlobalCoordinates0[iPoint][2]))
          self.node[iPoint].SetCoord0((self.GlobalCoordinates0[iPoint][0],self.GlobalCoordinates0[iPoint][1],self.GlobalCoordinates0[iPoint][2]))
          self.node[iPoint].SetCoord_prec((self.GlobalCoordinates0[iPoint][0],self.GlobalCoordinates0[iPoint][1],self.GlobalCoordinates0[iPoint][2]))

          #self.node[iPoint].SetCoord((self.GlobalCoordinates0[iPoint][0],self.GlobalCoordinates0[iPoint][1], self.GlobalCoordinates0[iPoint][2]))
          #self.node[iPoint].SetCoord0((self.GlobalCoordinates0[iPoint][0],self.GlobalCoordinates0[iPoint][1], self.GlobalCoordinates0[iPoint][2]))
          #self.node[iPoint].SetCoord_prec((self.GlobalCoordinates0[iPoint][0],self.GlobalCoordinates0[iPoint][1], self.GlobalCoordinates0[iPoint][2]))

  def ExtractDisplacements(self, iVertex):
      Coord0 = self.node[iVertex].GetCoord0()
      Coord_prec = self.node[iVertex].GetCoord_prec()
      Coord  = self.node[iVertex].GetCoord()
      diff = Coord - Coord0
      dispX = diff[0]; dispY = diff[1]; dispZ = diff[2]
      #print("{} {} {}".format(dispX, dispY, dispZ))
      return dispX, dispY, dispZ

  def __readConfig(self):
    """ Read structural tester config file. """

    with open(self.Config_file) as configfile:
      while 1:
        line = configfile.readline()
        if not line:
          break

        # remove line returns
        line = line.strip('\r\n')
        # make sure it has useful data
        if (not "=" in line) or (line[0] == '%'):
          continue
        # split across equal sign
        line = line.split("=",1)
        this_param = line[0].strip()
        this_value = line[1].strip()

        for case in switch(this_param):
          #float values
          if case("PITCHING_AMPL")          : pass
          if case("INITIAL_DISP")           : pass
          if case("INITIAL_ANGLE")          : pass
          if case("PLUNGING_AMPL")          : pass
          if case("PITCHING_FREQ")          : pass
          if case("PLUNGING_FREQ")          : pass
          if case("PITCHING_0")             : pass
          if case("PLUNGING_0")             :
            self.Config[this_param] = float(this_value)
            break

          #string values
          if case("MESH_FILE")		    : pass
          if case("SOLID_FORCES_OUTPUT")    : pass
          if case("CSD_SOLVER")	       	    : pass
          if case("MOVING_MARKER")	    : pass
          if case("STRUCT_TYPE")            :
            self.Config[this_param] = this_value
            break

          if case():
            print(this_param + " is an invalid option !")
            break

  def __computeInterfacePosVel(self,time,FSI_config, MLS_Spline):
    """ Description. """

    newCoord = np.zeros((3,1))
    newVel = np.zeros((3,1))
    rotCoorddot = np.zeros((3,1))

    if not FSI_config['START_MOTION_TIME']:
        FSI_config['START_MOTION_TIME'] = float(0.0)
     
    if (time < FSI_config['START_MOTION_TIME'] or MLS_Spline == None):  #self.startTime:
      if FSI_config['RESTART_SOL'] == 'NO': 
        #print(self.markers.keys())
        for iPoint in range(0, self.nPoint):
                Coord0 = self.node[iPoint].GetCoord0()
                newCoord = Coord0
                
                if iPoint == int(FSI_config["TRACKING_NODE"]):
                       
                   #print("Marker id {} phase and module: {} {} ".format(iPoint, phase(Flutter_mode_fluid_z[i]), scaling_coeff*abs(Flutter_mode_fluid_z[i]) ))
                   #print("Marker id {} z position: {}".format(iPoint, newCoord[2] ))
                       self.tracking_pos.append(float(newCoord[2]));
                       self.tracking_vel.append(float(newVel[2]));
                       self.time.append(float(time));

                self.node[iPoint].SetCoord_prec((newCoord[0], newCoord[1], newCoord[2]))
                self.node[iPoint].SetCoord((newCoord[0], newCoord[1], newCoord[2]))
                self.node[iPoint].SetVel((newVel[0], newVel[1], newVel[2]))


    
    else:
        
     if FSI_config['MOTION_TYPE'] == 'HARMONIC': 
        time = time - FSI_config['START_MOTION_TIME'] #self.startTime # this trick should allow the current formulation for  which at t=0 the movement starts from the initial conditions

        # the amplitude of the modal shape displacements is scaled in torder to reduced the appearent velocity given by the mode (Romanelli)
        conc = np.concatenate((np.absolute(self.Flutter_mode_fluid_x),np.absolute(self.Flutter_mode_fluid_y),np.absolute(self.Flutter_mode_fluid_z)) , axis=0)
        maxl = np.amax(conc, axis=0)
        if maxl < 1.e-16:
            scaling_coeff = 0
        else:    
            scaling_coeff = tan(radians(1))*(2*FSI_config["L_REF"]/FSI_config["K_MAX"])/maxl      #FSI_config["V_INF"]/FSI_config["FREQ"]
        '''
        Now, aeropoints to generate the interpolation matrix (query points) are ordered into a vector follwing the list of markers[iMarker]
        so I just need to extract them in order starting from 1 till the end during the cycle for iPoint in vertexList:
        '''
        i = 0
        for iPoint in range(0, self.nPoint):
                Coord0 = self.node[iPoint].GetCoord0()
                Coord = self.node[iPoint].GetCoord()
                newCoord[0] = Coord0[0]  + scaling_coeff*abs(self.Flutter_mode_fluid_x[i])*sin(FSI_config["OMEGA"]*time + phase(self.Flutter_mode_fluid_x[i]))
                newCoord[1] = Coord0[1]  + scaling_coeff*abs(self.Flutter_mode_fluid_y[i])*sin(FSI_config["OMEGA"]*time + phase(self.Flutter_mode_fluid_y[i]))
                newCoord[2] = Coord0[2]  + scaling_coeff*abs(self.Flutter_mode_fluid_z[i])*sin(FSI_config["OMEGA"]*time + phase(self.Flutter_mode_fluid_z[i]))
  
                newVel[0] = scaling_coeff*FSI_config["OMEGA"]*abs(self.Flutter_mode_fluid_x[i])*cos(FSI_config["OMEGA"]*time + phase(self.Flutter_mode_fluid_x[i]))
                newVel[1] = scaling_coeff*FSI_config["OMEGA"]*abs(self.Flutter_mode_fluid_y[i])*cos(FSI_config["OMEGA"]*time + phase(self.Flutter_mode_fluid_y[i]))
                newVel[2] = scaling_coeff*FSI_config["OMEGA"]*abs(self.Flutter_mode_fluid_z[i])*cos(FSI_config["OMEGA"]*time + phase(self.Flutter_mode_fluid_z[i]))                
                
                if iPoint == int(FSI_config["TRACKING_NODE"]):
                       
                   #print("Marker id {} phase and module: {} {} ".format(iPoint, phase(Flutter_mode_fluid_z[i]), scaling_coeff*abs(Flutter_mode_fluid_z[i]) ))
                   #print("Marker id {} z position: {}".format(iPoint, newCoord[2] ))
                       self.tracking_pos.append(float(newCoord[2]));
                       self.tracking_vel.append(float(newVel[2]));
                       self.time.append(float(time + FSI_config['START_MOTION_TIME']));

                self.node[iPoint].SetCoord_prec((Coord[0], Coord[1], Coord[2]))
                self.node[iPoint].SetCoord((newCoord[0], newCoord[1], newCoord[2]))
                self.node[iPoint].SetVel((newVel[0], newVel[1], newVel[2]))
                i = i+1
                
     elif FSI_config['MOTION_TYPE'] == 'BLENDED_STEP':
         
        time = time - FSI_config['START_MOTION_TIME'] #self.startTime # this trick should allow the current formulation for  which at t=0 the movement starts from the initial conditions

        ## the amplitude of the modal shape displacements is scaled in torder to reduced the appearent velocity given by the mode (Romanelli)
        #conc = np.concatenate((np.absolute(self.Flutter_mode_fluid_x),np.absolute(self.Flutter_mode_fluid_y),np.absolute(self.Flutter_mode_fluid_z)) , axis=0)
        #maxl = np.amax(conc, axis=0)
        ## Amplitude of the blended step
        #if maxl < 1.e-16:
        #    Aq = 0
        #else:            
        #    Aq = tan(radians(1))*(4*FSI_config["L_REF"]/FSI_config["K_MAX"])/maxl
        
        tau = time*FSI_config["V_INF"]/FSI_config["L_REF"]
        tau_q = pi*2/FSI_config["K_MAX"]
        
        if tau < tau_q:

            q = self.Aq/2*(1-cos(FSI_config["K_MAX"]/2*tau))
            q_dot = self.Aq*FSI_config["K_MAX"]/4 * sin(FSI_config["K_MAX"]/2*tau)  
        else:    
            q = self.Aq
            q_dot = 0

        ''' 
        print("DEBUG: q = {0:1.20f}".format(q))
        print("DEBUG: FSI_config[K_MAX] = {0:1.20f}".format(FSI_config["K_MAX"]))
        print("DEBUG: VINF = {0:1.20f}".format(FSI_config["V_INF"]))
        print("DEBUG: L_REF = {0:1.20f}".format(FSI_config["L_REF"]))
        print("DEBUG: tau = {0:1.20f}".format(tau))
        print("DEBUG: self.Flutter_mode_fluid_x[i] = {0:1.20f}".format(self.Flutter_mode_fluid_x[15114]))
        print("DEBUG: self.Flutter_mode_fluid_y[i] = {0:1.20f}".format(self.Flutter_mode_fluid_y[15114]))
        print("DEBUG: self.Flutter_mode_fluid_z[i] = {0:1.20f}".format(self.Flutter_mode_fluid_z[15114]))
        '''
        '''
        Now, aeropoints to generate the interpolation matrix (query points) are ordered into a vector follwing the list of markers[iMarker]
        so I just need to extract them in order starting from 1 till the end during the cycle for iPoint in vertexList:
        '''


        for iPoint in range(0, self.nPoint):
                Coord0 = self.node[iPoint].GetCoord0()
                Coord = self.node[iPoint].GetCoord()
                newCoord[0] = Coord0[0]  + q * self.Flutter_mode_fluid_x[iPoint];
                newCoord[1] = Coord0[1]  + q * self.Flutter_mode_fluid_y[iPoint];
                newCoord[2] = Coord0[2]  + q * self.Flutter_mode_fluid_z[iPoint];

                newVel[0] = q_dot * self.Flutter_mode_fluid_x[iPoint];
                newVel[1] = q_dot * self.Flutter_mode_fluid_y[iPoint];
                newVel[2] = q_dot * self.Flutter_mode_fluid_z[iPoint];
                #print("HD - Marker id {} Coord0[x] Coord0[y] Coord0[z]: {} {} {}".format(iPoint, Coord0[0], Coord0[1], Coord0[2] ))
                if iPoint == int(FSI_config["TRACKING_NODE"]):
                       
                   #print("Marker id {} Coord0[2]: {}".format(iPoint, Coord0[2] ))
                   #print("Marker id {} Flutter_mode_fluid_z: {}".format(iPoint, Flutter_mode_fluid_z[i] ))
                   #print("Marker id {} tau: {}".format(iPoint, tau ))
                   #print("Marker id {} q: {}".format(iPoint, q ))
                   #print("Marker id {} q_dot: {}".format(iPoint, q_dot ))
                   #print("Marker id {} Aq: {}".format(iPoint, Aq ))
                   #print("Marker id {} z position: {}".format(iPoint, newCoord[2] ))
                   #print("Marker id {} z vel: {}".format(iPoint, newVel[2] ))
                   self.tracking_pos.append(float(newCoord[2]));
                   self.tracking_vel.append(float(newVel[2]));
                   self.time.append(float(time + FSI_config['START_MOTION_TIME']));

                self.node[iPoint].SetCoord_prec((Coord[0], Coord[1], Coord[2]))
                self.node[iPoint].SetCoord((newCoord[0], newCoord[1], newCoord[2]))
                self.node[iPoint].SetVel((newVel[0], newVel[1], newVel[2]))


  def initialize_OutputForces(self, NbTimeIter,FSI_config):
    """ Description. """

    TimeIter_Simul = NbTimeIter - FSI_config['RESTART_ITER']
    
    #print("HD: NbTimeIter = {}".format(NbTimeIter))
    #print("HD: FSI_config['RESTART_ITER'] = {}".format(FSI_config['RESTART_ITER']))
    MovingVertex = 0
    if FSI_config['MEMO_FORCE_OUTPUT'] == 'YES':
       self.NodalForces = np.zeros((self.nPoint,5,NbTimeIter))
    #if FSI_config['MEMO_GEN_FORCE_OUTPUT'] == 'YES': 
    #   self.GenForces = np.zeros((NbTimeIter,2))

  def exit(self):
    """ Description. """

    print("\n**************** Exiting the structural tester solver ****************")

  def run(self,time,FSI_config, MLS_Spline):
    """ Description. """

    #self.__temporalIteration()  NOT NEEDED ANYMORE AS THE DISPLACEMENT IS IMPOSED

    print("Time")#\tDisp 1\tDisp2\tVel 1\tVel2\tAcc 1\tAcc 2")
    print(str(time))# + '\t' + str(float(self.q[0])) + '\t' + str(float(self.q[1])) + '\t' + str(float(self.qdot[0])) + '\t' + str(float(self.qdot[1])) + '\t' + str(float(self.qddot[0])) + '\t' + str(float(self.qddot[1])))

    self.__computeInterfacePosVel(time,FSI_config, MLS_Spline)

  def run_restart(self,restart_time,FSI_config, MLS_Spline):
    """ Description. """

    print("Solid Node updated positions from restart point")#\tDisp 1\tDisp2\tVel 1\tVel2\tAcc 1\tAcc 2")
    #print('Restart time =' + str(restart_time))# + '\t' + str(float(self.q[0])) + '\t' + str(float(self.q[1])) + '\t' + str(float(self.qdot[0])) + '\t' + str(float(self.qdot[1])) + '\t' + str(float(self.qddot[0])) + '\t' + str(float(self.qddot[1])))

    self.__computeInterfacePosVel(restart_time, FSI_config, MLS_Spline)  
    # A spurious tracking value is added now that has to be removed
    self.time = []; self.tracking_pos = []; self.tracking_vel = [];
    
    
  def setInitialDisplacements(self,FSI_config, MLS_Spline):
    """ Description. """
    time = 0;
    print("Set Initial Displacements")
    self.__computeInterfacePosVel(time,FSI_config, MLS_Spline)

  def writeSolution(self, time_iter, time, FSI_config):  #TimeIter, NbTimeIter):
    """ Description. """
    #wait = input("PRESS ENTER TO CONTINUE.")
    
    cont_time_iter = time_iter - FSI_config['RESTART_ITER']
    if FSI_config['CSD_SOLVER']	== 'IMPOSED': # at the moment can't be used (towards unification: IMPOSED works for airfoil pitching and  plunging)
     a=0
     #if time == 0:
     #  histFile = open(FSI_config['NODAL_FORCE_FILE'], "w")
     #  histFile.write("Time\tDisp 1\tDisp2\tVel 1\tVel2\tAcc 1\tAcc 2\tAccVar 1\tAccVar 2\n")
     #else:
     #  histFile = open(FSI_config['NODAL_FORCE_FILE'], "a")
     #  histFile.write(str(time) + '\t' + str(float(self.q[0])) + '\t' + str(float(self.q[1])) + '\t' + str(float(self.qdot[0])) + '\t' + str(float(self.qdot[1])) + '\t' + str(float(self.qddot[0])) + '\t' + str(float(self.qddot[1])) + '\t' + str(float(self.a[0])) + '\t' + str(float(self.a[1])) + '\n')
     #histFile.close()
    elif ((FSI_config['CSD_SOLVER'] == 'NITRO' or FSI_config['CSD_SOLVER'] == 'NITRO_FRAMEWORK') and FSI_config['UNSTEADY_SIMULATION'] == 'YES'  ):
     if FSI_config['WRITE_FORCE_OUTPUT'] == 'YES':   
        # Opens the force file and writes the header 
        histFile = open(FSI_config['NODAL_FORCE_FILE'] + '_' + str(time_iter) + '.dat' , "w")   
        histFile.write("Node\tTime\tForce x\tForce y\tForce z\n")
     for iPoint in range(0, self.nPoint):
            Force = self.node[iPoint].GetForce()  #FX += float(Force[0])  #FY += float(Force[1])  #FZ += float(Force[2])
            # Writes in the force file the current line for the node
            if FSI_config['WRITE_FORCE_OUTPUT'] == 'YES':
               histFile.write(str(iPoint) + '\t' + str(float(time)) +'\t' + str(float(Force[0])) + '\t' + str(float(Force[1])) + '\t' + str(float(Force[2])) + '\n')               
            # Memorizes in the force matrix the current line for the node
            if FSI_config['MEMO_FORCE_OUTPUT'] == 'YES':
               self.NodalForces[a,0,cont_time_iter] =  iPoint; self.NodalForces[a,1,0] =  time; self.NodalForces[a,2,cont_time_iter] =  Force[0]; self.NodalForces[a,3,cont_time_iter] =  Force[1]; self.NodalForces[a,4,cont_time_iter] =  Force[2]; 
            # Memorizes for the tracking node the quantities to track
            #print("HD - Marker id {} Force[x] Force[y] Force[z]: {} {} {}".format(iPoint, Force[0], Force[1], Force[2] ))
            if iPoint == int(FSI_config["TRACKING_NODE"]):
               self.tracking_force_x.append(float(Force[0])); 
               self.tracking_force_y.append(float(Force[1])); 
               self.tracking_force_z.append(float(Force[2])); 
     # Close the file force if open       
     if FSI_config['WRITE_FORCE_OUTPUT'] == 'YES':       
        histFile.close()
    elif ((FSI_config['CSD_SOLVER'] == 'NITRO' or FSI_config['CSD_SOLVER'] == 'NITRO_FRAMEWORK') and FSI_config['UNSTEADY_SIMULATION'] == 'NO'  ):
     if FSI_config['WRITE_FORCE_OUTPUT'] == 'YES':
        # Opens the force file and writes the header 
        histFile = open(FSI_config['NODAL_FORCE_FILE'] + '_Steady' + '.dat', "w")   
        histFile.write("Node\tTime\tForce x\tForce y\tForce z\n")


     for iPoint in range(0, self.nPoint):
            Force = self.node[iPoint].GetForce()  #FX += float(Force[0])  #FY += float(Force[1])  #FZ += float(Force[2])
            # Writes in the force file the current line for the node
            if FSI_config['WRITE_FORCE_OUTPUT'] == 'YES':
               histFile.write(str(iPoint) + '\t' + str(0) + '\t' + str(float(Force[0])) + '\t' + str(float(Force[1])) + '\t' + str(float(Force[2])) + '\n')               
            # Memorizes in the force matrix the current line for the node
            if FSI_config['MEMO_FORCE_OUTPUT'] == 'YES':            
               self.NodalForces[a,0,0] =  iPoint; self.NodalForces[a,1,0] =  0; self.NodalForces[a,2,0] =  Force[0]; self.NodalForces[a,3,0] = Force[1]; self.NodalForces[a,4,0] = Force[2]; 
        
     if FSI_config['WRITE_FORCE_OUTPUT'] == 'YES':       
        histFile.close()

    if FSI_config['WRITE_GEN_FORCE_OUTPUT'] == 'YES':
      for mode in range(self.nModes):
         self.writeGenForceFile( time_iter, time, FSI_config, mode)    
         

  def writeGenForceFile(self, time_iter, time, FSI_config, mode):
    if FSI_config['CSD_SOLVER']	== 'IMPOSED': # at the moment can't be used (towards unification: IMPOSED works for airfoil pitching and  plunging)
     a=0
    elif ((FSI_config['CSD_SOLVER'] == 'NITRO' or FSI_config['CSD_SOLVER'] == 'NITRO_FRAMEWORK') and FSI_config['UNSTEADY_SIMULATION'] == 'YES'  ):
     # Opens the generalized force file and writes the header 
     if time_iter == 0:
           genForceHistFile = open(FSI_config['GENERALIZED_FORCE_FILE'] + str(mode) + '.dat' , "w") 
           genForceHistFile.write("Time\tGeneralized Force\n")
     elif time_iter == 2 and FSI_config['MOTION_TYPE'] == 'BLENDED_STEP' and FSI_config['RESTART_SOL'] == 'YES' and FSI_config['RESTART_ITER'] == 2:
           genForceHistFile = open(FSI_config['GENERALIZED_FORCE_FILE'] + str(mode) + '.dat' , "r+")
           genForceHistFile.readline(); line_to_get = genForceHistFile.readline().split(); # the second line is the one to be repeated
           Generalized_force_i = line_to_get[1]         
           genForceHistFile.write(str(float(FSI_config['UNST_TIMESTEP'])) +'\t' + str(  float(  Generalized_force_i) ) + '\n' ) # we extend the force on the first timestep and we start the movement here  
     else:
           genForceHistFile = open(FSI_config['GENERALIZED_FORCE_FILE'] + str(mode) + '.dat' , "a")     

     Generalized_force_i=0

     for iPoint in range(0, self.nPoint):
            Force = self.node[iPoint].GetForce()  #FX += float(Force[0])  #FY += float(Force[1])  #FZ += float(Force[2])            
            # Generalized forces
            Generalized_force_i =  Generalized_force_i + (1/self.Aq)*Force[0]*self.mode_fluid_x[iPoint,mode] + (1/self.Aq)*Force[1]*self.mode_fluid_y[iPoint,mode] + (1/self.Aq)*Force[2]*self.mode_fluid_z[iPoint,mode]

     # appends to the generalized force file the value relative to the current timestep       
    
     genForceHistFile.write( str(float(time) ) +'\t' + str( float( Generalized_force_i)  ) + '\n' )  
     genForceHistFile.close()
        
    elif ((FSI_config['CSD_SOLVER'] == 'NITRO' or FSI_config['CSD_SOLVER'] == 'NITRO_FRAMEWORK') and FSI_config['UNSTEADY_SIMULATION'] == 'NO'  ): 
     # Opens the generalized force file and writes the header 
     if time_iter == 0: # forced to be zero in input in case of steady simulation
           genForceHistFile = open(FSI_config['GENERALIZED_FORCE_FILE'] + str(mode) + '.dat' , "w") 
           genForceHistFile.write("Time\tGeneralized Force\n")
     else:    # not really used
           genForceHistFile = open(FSI_config['GENERALIZED_FORCE_FILE'] + str(mode) + '.dat' , "a")        # not really used

     Generalized_force_i=0

     for iPoint in range(0, self.nPoint):
            Force = self.node[iPoint].GetForce()  #FX += float(Force[0])  #FY += float(Force[1])  #FZ += float(Force[2]) 
            # Generalized forces   
            Generalized_force_i =  Generalized_force_i + (1/self.Aq)*Force[0]*self.mode_fluid_x[iPoint,mode] + (1/self.Aq)*Force[1]*self.mode_fluid_y[iPoint,mode] + (1/self.Aq)*Force[2]*self.mode_fluid_z[iPoint,mode]

     # appends to the generalized force file the value relative to the current timestep       
     #if FSI_config['WRITE_GEN_FORCE_OUTPUT'] == 'YES':        
     genForceHistFile.write(str(float(time)) +'\t' + str(  float(Generalized_force_i))  + '\n' )  
     genForceHistFile.close()         
         
  def EvaluateIntefaceFluidDisplacements(self, FSI_config, MLS_Spline):  #TimeIter, NbTimeIter):
    """ Description. """
    Interf_matrix = MLS_Spline.interpolation_matrix;

    # format modes
    if FSI_config['FORMAT_MODES'] == 'CSHELL':
       idx = 0; idy = 1; idz = 2;
    elif FSI_config['FORMAT_MODES'] == 'NASTRANF06':
       idx = 1;  idy = 2; idz = 3;
    else:
       print("Format {} not known !!".format(FORMAT_MODES))
       sys.exit("Goodbye!")

    if FSI_config['CSD_SOLVER'] == 'NITRO':
       self.gen_coord_arr = FSI_config["MODAL_COEFF"]; 
       self.nModes = len(self.gen_coord_arr)
       #Modes[i].GetMode()[:,0]
       Flutter_mode_str = np.zeros((MLS_Spline.nStructNodes, 6),dtype=np.complex_)  # I need to define numpy matrix as complex
       for i in range(self.nModes):

          Flutter_mode_str = Flutter_mode_str + self.gen_coord_arr[i] * MLS_Spline.Modes[i].GetMode();

         
       self.Flutter_mode_fluid_x =  Interf_matrix.dot(Flutter_mode_str[:,0])   # let's try to evaluate the displacemeny only to sum to the aero boundary nodes
       self.Flutter_mode_fluid_y =  Interf_matrix.dot(Flutter_mode_str[:,1])
       self.Flutter_mode_fluid_z =  Interf_matrix.dot(Flutter_mode_str[:,2])    
     
       # HARD CODED we take only the real part of the mode as the immaginary HAS to be zero!!!!
       self.Flutter_mode_fluid_x = self.Flutter_mode_fluid_x.real
       self.Flutter_mode_fluid_y = self.Flutter_mode_fluid_y.real
       self.Flutter_mode_fluid_z = self.Flutter_mode_fluid_z.real
    
       # Need to memorize also all the other modes in order to calculate the relativ egeneralized forces
       self.mode_fluid_x = np.zeros((self.nPoint,self.nModes))
       self.mode_fluid_y = np.zeros((self.nPoint,self.nModes))
       self.mode_fluid_z = np.zeros((self.nPoint,self.nModes))
    
       for i in range(self.nModes):
          self.mode_fluid_x[:,i] = Interf_matrix.dot( MLS_Spline.Modes[i].GetMode()[:,idx] )
          self.mode_fluid_y[:,i] = Interf_matrix.dot( MLS_Spline.Modes[i].GetMode()[:,idy] )
          self.mode_fluid_z[:,i] = Interf_matrix.dot( MLS_Spline.Modes[i].GetMode()[:,idz] )

    else: # int this case I mean NITRO_FRAMEWORK   
       
       self.nModes = FSI_config["NMODES"]
       #Modes[i].GetMode()[:,0]
       Flutter_mode_str = np.zeros((MLS_Spline.nStructNodes, 6))  # I need to define numpy matrix as complex

       Flutter_mode_str = MLS_Spline.Modes[FSI_config["MODE_TO_SIMULATE"]].GetMode();

         
       self.Flutter_mode_fluid_x =  Interf_matrix.dot(Flutter_mode_str[:,0])   # let's try to evaluate the displacemeny only to sum to the aero boundary nodes
       self.Flutter_mode_fluid_y =  Interf_matrix.dot(Flutter_mode_str[:,1])
       self.Flutter_mode_fluid_z =  Interf_matrix.dot(Flutter_mode_str[:,2])    
      
       # Need to memorize also all the other modes in order to calculate the relativ egeneralized forces
       self.mode_fluid_x = np.zeros((self.nPoint,self.nModes))
       self.mode_fluid_y = np.zeros((self.nPoint,self.nModes))
       self.mode_fluid_z = np.zeros((self.nPoint,self.nModes))

       for i in range(self.nModes):
          self.mode_fluid_x[:,i] = Interf_matrix.dot( MLS_Spline.Modes[i].GetMode()[:,idx] )
          self.mode_fluid_y[:,i] = Interf_matrix.dot( MLS_Spline.Modes[i].GetMode()[:,idy] )
          self.mode_fluid_z[:,i] = Interf_matrix.dot( MLS_Spline.Modes[i].GetMode()[:,idz] )
          
          
    if FSI_config['MOTION_TYPE'] == 'BLENDED_STEP':    
    
     # the amplitude of the modal shape displacements is scaled in order to reduced the appearent velocity given by the mode (Romanelli)
     conc = np.concatenate((np.absolute(self.Flutter_mode_fluid_x),np.absolute(self.Flutter_mode_fluid_y),np.absolute(self.Flutter_mode_fluid_z)) , axis=0)
     maxl = np.amax(conc, axis=0)
     # Amplitude of the blended step
     if maxl < 1.e-16:
         self.Aq = 0
     else:
         self.Aq = tan(radians(1))*(4*FSI_config["L_REF"]/FSI_config["K_MAX"])/maxl
         #self.Aq = round(self.Aq,5)
     print("Aq = {}".format(self.Aq))
    
  def updateSolution(self):
    """ Description. """

    makerID = self.markers.keys()[0]
    nodeList = self.markers[makerID]

    for iPoint in nodeList:
      self.node[iPoint].updateCoordVel()

  def applyload(self, iVertex, fx, fy, fz):
    """ Description """

    self.node[iVertex].SetForce((fx,fy,fz))

  def getFSIMarkerID(self):
    """ Description. """
    
    list = self.markers.keys()
    return list[0]

  def getNumberOfSolidInterfaceNodes(self, markerID):
    """ Description. """

    return len(self.markers[markerID])

  def getInterfaceNodeGlobalIndex(self, markerID, iVertex):
    """ Description. """

    return self.markers[markerID][iVertex]

  def getInterfaceNodePosX(self, markerID, iVertex):
    """ Desciption. """

    iPoint = self.markers[markerID][iVertex]
    Coord = self.node[iPoint].GetCoord()
    return float(Coord[0])

  def getInterfaceNodePosY(self, markerID, iVertex):
    """ Desciption. """

    iPoint = self.markers[markerID][iVertex]
    Coord = self.node[iPoint].GetCoord()
    return float(Coord[1])

  def getInterfaceNodePosZ(self, markerID, iVertex):
    """ Desciption. """

    iPoint = self.markers[markerID][iVertex]
    Coord = self.node[iPoint].GetCoord()
    return float(Coord[2])

  def getInterfaceNodeDispX(self, markerID, iVertex):
    """ Desciption. """

    iPoint = self.markers[markerID][iVertex]
    Coord = self.node[iPoint].GetCoord()
    Coord0 = self.node[iPoint].GetCoord0()
    return float(Coord[0]-Coord0[0])

  def getInterfaceNodeDispY(self, markerID, iVertex):
    """ Desciption. """

    iPoint = self.markers[markerID][iVertex]
    Coord = self.node[iPoint].GetCoord()
    Coord0 = self.node[iPoint].GetCoord0()
    return float(Coord[1]-Coord0[1])

  def getInterfaceNodeDispZ(self, markerID, iVertex):
    """ Desciption. """

    iPoint = self.markers[markerID][iVertex]
    Coord = self.node[iPoint].GetCoord()
    Coord0 = self.node[iPoint].GetCoord0()
    return float(Coord[2]-Coord0[2])



  def printForceDispl(self,TimeIter):

       outC = open("./Output/Coord.txt" + str(TimeIter), "w")
       outF = open("./Output/Forces.txt" + str(TimeIter), "w")
       outD = open("./Output/Disp.txt"+ str(TimeIter), "w")
       outG = open("./Output/Gen_forces.txt"+ str(TimeIter), "w")
       for iPoint in range(0, self.nPoint):
           Force = self.node[iPoint].GetForce()
           outF.write(str(iPoint))
           outF.write("\t")
           outF.write(str(float(Force[0])))
           outF.write("\t")
           outF.write(str(float(Force[1])))
           outF.write("\t")
           outF.write(str(float(Force[2])))
           outF.write("\n")

           POS = self.node[iPoint].GetCoord()
           outC.write(str(iPoint))
           outC.write("\t")
           outC.write(str(float(POS[0])))
           outC.write("\t")
           outC.write(str(float(POS[1])))
           outC.write("\t")
           outC.write(str(float(POS[2] )))
           outC.write("\n")

           outD.write(str(iPoint))
           outD.write("\t")
           outD.write(str(float(self.mode_fluid_x[iPoint,0] )))
           outD.write("\t")
           outD.write(str(float(self.mode_fluid_y[iPoint,0] )))
           outD.write("\t")
           outD.write(str(float(self.mode_fluid_z[iPoint,0]) ))
           outD.write("\n")

           outG.write(str(float((1/self.Aq)*Force[0]*self.mode_fluid_x[iPoint,0] + (1/self.Aq)*Force[1]*self.mode_fluid_y[iPoint,0] + (1/self.Aq)*Force[2]*self.mode_fluid_z[iPoint,0])))
           outG.write("\n")

       outC.close()
       outF.close()
       outD.close()

  def printNode(self, TimeIter,steady_unsteady):

      iPoint = 15114
      if steady_unsteady == 0:
         outC = open("./Output/Resume"  + ".txt", "w")

         outC.write(str(self.Aq))
         outC.write("\t")
         outC.write(str(iPoint))
         outC.write("\n")

      else:
         outC = open("./Output/Resume"  + ".txt", "a")

      POS = self.node[iPoint].GetCoord()
      POS0 = self.node[iPoint].GetCoord0()
      POS_prec = self.node[iPoint].GetCoord_prec()

      outC.write(str(float(TimeIter)))
      outC.write("\n")

      outC.write(str(float(POS0[0])))
      outC.write("\t")
      outC.write(str(float(POS0[1])))
      outC.write("\t")
      outC.write(str(float(POS0[2])))
      outC.write("\n")

      outC.write(str(float(POS_prec[0])))
      outC.write("\t")
      outC.write(str(float(POS_prec[1])))
      outC.write("\t")
      outC.write(str(float(POS_prec[2])))
      outC.write("\n")

      outC.write(str(float(POS[0])))
      outC.write("\t")
      outC.write(str(float(POS[1])))
      outC.write("\t")
      outC.write(str(float(POS[2])))
      outC.write("\n")


      outC.write(str(float(self.mode_fluid_x[iPoint, 0])))
      outC.write("\t")
      outC.write(str(float(self.mode_fluid_y[iPoint, 0])))
      outC.write("\t")
      outC.write(str(float(self.mode_fluid_z[iPoint, 0])))
      outC.write("\n")
      
      outC.close()
