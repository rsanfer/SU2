#!/usr/bin/env python3

## \file ReadModes.py
#  \brief Function that reads structural modes according to different formats.
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


import pdb
import os, sys, shutil, copy
import numpy as np
import scipy as sp
import scipy.linalg as linalg
from math import *

class CModes:  # for modes

  def __init__(self):    
    self.Mode = np.empty((0,6),float)# one row only... to add more lines
    
  def SetModeLine(self,node_l):
    self.Mode = np.append(self.Mode, np.array([node_l]), axis=0)
    
  def GetMode(self):
    return self.Mode  

  def GetL(self):
    a = self.Mode.shape()[0]
    return a


class CModes_f06:  # for modes

    def __init__(self):
        self.Mode = np.empty((0, 7), float)  # one row only... to add more lines

    def SetModeLine(self, node_l):
        self.Mode = np.append(self.Mode, np.array([node_l]), axis=0)

    def GetMode(self):
        return self.Mode

    def GetL(self):
        a = self.Mode.shape()[0]
        return a

    def OrderMode(self):
        # order the modes in ascending order with respect of the node grid
      self.Mode = self.Mode[self.Mode[:,0].argsort()]

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
           print("Format {} not known !!".format(FORMAT_MODES))
           sys.exit("Goodbye!")
    nModes.append(int(index))
    if FORMAT_MODES == 'NASTRANF06':
        # order modes in asxcending node order (to ensure compatibility with structural grid independently of the reading order)
        for i in range(0, nModes[0]):
            Modes[i].OrderMode()

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






















