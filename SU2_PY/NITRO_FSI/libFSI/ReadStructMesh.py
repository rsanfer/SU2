#!/usr/bin/env python3

## \file ReadModes.py
#  \brief Function that reads structural mesh according to different formats.
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

# ----------------------------------------------------------------------
#  Input file loading
# ----------------------------------------------------------------------
class Point:  # for nodes
  """ Description. """
  
  def __init__(self):
    self.Id = int(0)
    self.Coord = np.zeros((3,1))

  def SetId(self,Id):
      self.Id = Id

  def SetCoord(self,val_Coord):
    x, y, z = val_Coord
    self.Coord[0] = x
    self.Coord[1] = y
    self.Coord[2] = z
    
  def GetCoord(self):
    return self.Coord

  def GetId(self):
    return self.Id


def ReadStructMesh(Mesh_file, Mesh_format, node, nPoint):
    """ Read the mesh and saves the structural mode   
    """
    with open(Mesh_file, 'r') as Meshfile:
        iPoint = -1
        if Mesh_format == "CSHELL":
            print('Opened mesh file ' + Mesh_file + '.' + " Format = " + Mesh_format)
            line = Meshfile.readline()
            line1 = line.split()
            nPoint.append(int(line1[0]) -1)
            while 1:
               line = Meshfile.readline()
               if line:
                  while iPoint < nPoint[0]:
                     iPoint = iPoint + 1                      
                     line = Meshfile.readline() 
                     line = line.split()
                     node.append(Point())
                     node[iPoint].SetId(iPoint)
                     node[iPoint].SetCoord([float(line[0]), float(line[1]), float(line[2])])
                     if not line:
                        break
               if not line:
                  break
            nPoint[0] = nPoint[0] +1                
        elif Mesh_format == "NASTRAN":
            print('Opened mesh file ' + Mesh_file + '.' + " Format = " + Mesh_format)
            while 1:
               line = Meshfile.readline()
               if not line:
                  break
               pos = line.find('GRID')  
               if pos != -1:                  
                  while line.find('GRID') !=-1:
                       #line1 = line.split()
                       iPoint = iPoint +1
                       node.append(Point())
                       iGRID =line[9:18].strip();  x = line[24:32].strip(); y = line[32:40].strip(); z = line[40:48].strip();
                       posx = x.find('-'); posy = y.find('-'); posz = z.find('-');  
                       # === In case the nastran file has exponential (-10 instead of e-10) 
                       # X coord
                       if posx != -1:
                          if posx == 0:
                             xx = x[1:]; posxx = xx.find('-')
                             if posxx != -1:
                                x_bis = '-' + xx[:posxx] + 'e' + xx[posxx:]
                             else:
                                x_bis = '-' + xx
                          else:   
                             x_bis = x[:posx] + 'e' + x[posx:]
                       else:
                          x_bis = x
                       # Y coord      
                       if posy != -1:
                          if posy == 0:
                             yy = y[1:]; posyy = yy.find('-')
                             if posyy != -1:
                                y_bis = '-' + yy[:posyy] + 'e' + yy[posyy:]
                             else:
                                y_bis = '-' + yy
                          else:   
                             y_bis = y[:posy] + 'e' + y[posy:]
                       else:
                          y_bis = y                         
                       # Z coord      
                       if posz != -1:
                          if posz == 0:
                             zz = z[1:]; poszz = zz.find('-')
                             if poszz != -1:
                                z_bis = '-' + zz[:poszz] + 'e' + zz[poszz:]
                             else:
                                z_bis = '-' + zz   
                          else:   
                             z_bis = z[:posz] + 'e' + z[posz:]
                       else:
                          z_bis = z                       
                       
                       
                       #x_bis = x[:posx] + 'e' + x[posx:] if (posx != -1 and posx != 0) else x
                       #y_bis = y[:posy] + 'e' + y[posy:] if (posy != -1 and posy != 0) else y
                       #z_bis = z[:posz] + 'e' + z[posz:] if (posz != -1 and posz != 0) else z
                       node[iPoint].SetId(int(iGRID))
                       node[iPoint].SetCoord([float(x_bis), float(y_bis), float(z_bis)])
                       #print("Node nr. {}".format(line[8:16]))
                       line = Meshfile.readline()
                       if not line:
                          break
            nPoint.append(int(iPoint +1))
            #nPoint = int(iPoint +1)
        else:
           print("Reading of format {} not yet implemented. Check routine ReadStructMesh".format(Mesh_format))
    print("Total number of structural nodes from input file is: {} \n".format(nPoint[0]))
     # readStructMesh(Mesh_file, Mesh_format, node, nPoint)
    
'''
if __name__ == "__main__":
    
    print("Storing Structural mesh information from the input file")
    Mesh_file = "./Modal_20_span.bdf"#MLS_conf['STRUCTURAL_NODES_FILE_NAME']
    Mesh_format = "NASTRAN" # MLS_conf['FORMAT_SRUCT_NODES']
    nStrPoint = []
    StructNodes = []   
    node = []
    nPoint = []
    ReadStructMesh(Mesh_file, Mesh_format, node, nPoint)
    nPoint = nPoint[0]

    if Mesh_format == 'NASTRAN' :

        idList = [node[i].GetId() for i in range( nPoint)]
        print(idList)
        Order = np.argsort(np.argsort(idList))

        for i in range(nPoint):
            print(node[i].GetCoord())
        print(Order)
        reorder(node, Order, nPoint)

        for i in range(nPoint):
            print(node[i].GetCoord())
'''

