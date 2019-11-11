#!/usr/bin/env python3

## \file ReadModes.py
#  \brief Function that reads structural modes according to different formats.
#  \author Rocco Bombardieri
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


import pdb
import os, sys, shutil, copy
import numpy as np
import scipy as sp
import scipy.linalg as linalg
from math import *

class CModes:  # for modes

  def __init__(self):    
    self.Mode = np.empty((0,6),float)# one row only... to add more lines with vstack
    
  def SetModeLine(self,node_l):
    self.Mode = np.append(self.Mode, np.array([node_l]), axis=0)   # es. A = numpy.vstack([A, newrow]) 
    
  def GetMode(self):
    return self.Mode  

  def GetL(self):
    a = self.Mode.shape()[0]
    return a


class CModes_f06:  # for modes

    def __init__(self):
        self.Mode = np.empty((0, 7), float)  # one row only... to add more lines with vstack

    def SetModeLine(self, node_l):
        self.Mode = np.append(self.Mode, np.array([node_l]), axis=0)  # es. A = numpy.vstack([A, newrow])

    def GetMode(self):
        return self.Mode

    def GetL(self):
        a = self.Mode.shape()[0]
        return a

def ReadModes(Modes, Mode_file, FORMAT_MODES, nModes):
    index = -1
    with open(Mode_file, 'r') as modefile:
        if FORMAT_MODES == "CSHELL":
            print('Opened mode file ' + Mode_file + '.' + " Format = " + FORMAT_MODES)
            line = modefile.readline()
            while 1:
               if not line:
                  break
               pos = line.find('T1')  
               if pos != -1:
                  index = index+1
                  Modes.append(CModes())
                  line = modefile.readline()
                  if not line:
                   break
                  while line.find('T1') ==-1:
                    line = line.split()
                    #print("index ={}".format(index))
                    node_disp = [float(line[0]), float(line[1]), float(line[2]), float(line[3]), float(line[4]), float(line[5])]
                    Modes[index].SetModeLine(node_disp)
                    #print(node_disp)
                    line = modefile.readline()
               
        elif FORMAT_MODES == "NASTRANF06":
            print('Opened mode file ' + Mode_file + '.' + " Format = " + FORMAT_MODES)
            line = modefile.readline()
            while 1:
               if not line:
                  break
               # find the string
               pos = line.find('R E A L   E I G E N V E C T O R   N O .')
               if pos != -1:
                  # Memorize the mode
                  print(line)
                  ModeNr = int( line.split('R E A L   E I G E N V E C T O R   N O .',1)[1] )
                  # now go ahead till you find
                  if index < ModeNr:
                     index = ModeNr
                     Modes.append(CModes_f06())
                  # look for the line starting with the modes
                  while line.find('T1') == -1:
                    if not line:
                       break
                    line = modefile.readline()
                  line = modefile.readline()
                  # starting with the modes
                  while line.find('  G  ') != -1:
                    print(line)
                    line = line.split()
                    node_disp = [float(line[0]), float(line[2]), float(line[3]), float(line[4]), float(line[5]), float(line[6]), float(line[7])]
                    print(node_disp)
                    Modes[index-1].SetModeLine(node_disp)
                    if not line:
                       break
                    line = modefile.readline()
               line = modefile.readline()
        else:
           print("Format not known !!".format(FORMAT_MODES))
           sys.exit("Goodbye!")
    nModes.append(int(index))
    print("Total number of available structural modes from input file is: {}".format(int(nModes[0])))


'''
if __name__ == "__main__":

      Modes = [];
      Mode_file = "./GTA_modes.f06"
      FORMAT_MODES = "NASTRANF06"
      nModes = []
      ReadModes(Modes, Mode_file, FORMAT_MODES, nModes)
      print("Mode 1")
      print(Modes[0].GetMode())
'''






















