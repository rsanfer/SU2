#!/usr/bin/env python3

## \file ReadModes.py
#  \brief Function that reads Modal displacements according to different formats.
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
import math
import numpy as np
from numpy import float32, float64, int32, int64, product
import matplotlib.pyplot as plt

# decrease output precision
np.set_printoptions(precision=3, threshold=20)


class CModalDisp:  # for modes

    def __init__(self):
        self.Mode = np.empty((0, 7), float)  # one row only... to add more lines with vstack

    def SetModeLine(self, node_l):
        self.Mode = np.append(self.Mode, np.array([node_l]), axis=0)  # es. A = numpy.vstack([A, newrow])

    def GetMode(self):
        return self.Mode

    def GetL(self):
        a = self.Mode.shape()[0]
        return a

def ReadModalDisplacements(FSIConfig):
    T0 = FSIConfig['START_TIME']
    TF = FSIConfig['UNST_TIME']
    TW = TF - T0
    DT = FSIConfig['UNST_TIMESTEP']
    NT = TW / DT + 1

    op4_filename = FSIConfig['MODAL_DISPLACEMENT']
    # Reading op4 file with time domain response
    matrices = readOutput4File(op4_filename)
    print(matrices[0][5])
    print(matrices[0][5].shape)
    # extract the modal values
    nModes = matrices[0][1]  # Total number of modes provided by dynresp
    mod_displ = np.zeros((nModes,int(NT)))
    mod_displ = matrices[0][5][:,::3]  # matrix nModes x NT

    return mod_displ

def readOutput4File(fileName):
    Out = []
    nextHeader = 1;
    length = 16
    matrix = []
    with open(fileName, 'r') as fid:
       print('Opened modal coordinates file ' + fileName + '.' + " Format = OUTPUT4")
       while 1:
          crtLine = fid.readline()
          #print(crtLine)
          if not crtLine:
              break
          elif not crtLine.strip():
              continue

          if nextHeader:
             ncols, nrows, dformat, dtype, name = getHeaderData(crtLine);
          del matrix
          matrix = []
          matrix.append(name);matrix.append(ncols);matrix.append(nrows);matrix.append(dformat);matrix.append(dtype)
          Vals = np.zeros((nrows,ncols))
          for iCol in range(1,ncols+2):
              crtLine = fid.readline()
              colIdx, idxFirstNZRow, NoNZElems = getColHeader(crtLine)
              if dtype == 'RealSimple':
                 NRowVals = nrows
              elif dtype == 'RealDouble':
                 NRowVals = nrows
              #elif dtype == 'ComplexSimple':
              #   NRowVals = 2*nrows
              #elif dtype == 'ComplexDouble':
              #   NRowVals = 2*nrows
              else:
                  continue
              ColLines = math.ceil(NoNZElems / 5)
              FirstNZ = idxFirstNZRow
              NoNZ = NoNZElems
              vals = np.zeros((NRowVals))
              for iCL in range(1,ColLines+1):
                  crtLine = fid.readline()
                  if iCL < ColLines:
                      #range_d = np.linspace(FirstNZ - 1 + ((iCL - 1) * 5 + 1),FirstNZ - 1 + iCL * 5, iCL * 5 - ((iCL - 1) * 5 + 1) +1, endpoint = True)
                      range_d = range(FirstNZ-1+((iCL-1)*5+1),FirstNZ-1+ iCL*5+1)
                  else:
                      #range_d = np.linspace(FirstNZ - 1 + ((iCL - 1) * 5) + 1, FirstNZ - 1 + NoNZ, NoNZ - ((iCL - 1) * 5) + 1 +1, endpoint = True)
                      range_d = range(FirstNZ-1+((iCL-1)*5+1),FirstNZ-1+ NoNZ+1  )
                  range_d = [x - 1 for x in range_d]
                  chunks = [crtLine[i:i + length] for i in range(0, len(crtLine), length)]
                  if chunks[-1] == '\n':
                      chunks.remove(chunks[-1])
                  #print(chunks)
                  chunks_map = [float(x) for x in chunks]
                  #print(chunks_map)
                  #print(range_d)
                  for i in range(0,len(range_d)):
                     vals[range_d[i]] = chunks_map[i]
              if colIdx <= ncols:
                  if dtype == 'RealSimple':
                     Vals[:, colIdx - 1] = vals
                  if dtype == 'RealDouble':
                     Vals[:, colIdx - 1] = vals
                  #if dtype == 'ComplexSimple':
                  #   Vals[:, colIdx - 1] = vals
                  #if dtype == 'ComplexDouble':
                  #   Vals[:, colIdx - 1] = vals
              else:
                  #print('Break')
                  break
          matrix.append(Vals)
          Out.append(matrix)
    fid.close()
    return Out


def getHeaderData(line):

    fields = [None]*6
    for i in range(1,6):
        fields[i-1] = line[ ((i - 1) * 8 + 1)-1:i * 8 ]
    #print(fields)
    fields[5] = line[41:-1]
    ncols = int(fields[0].strip())
    nrows = int(fields[1].strip())
    format = int(fields[2].strip())
    if format ==1:
        dformat = 'square'
    elif format ==2:
        dformat = 'rectangular'
    elif format ==3:
        dformat = 'diagonal'
    elif format ==4:
        dformat = 'lowerTriangular'
    elif format ==5:
        dformat = 'upperTriangular'
    elif format ==6:
        dformat = 'symmetric'
    elif format ==7:
        dformat = 'symmetric'
    elif format ==8:
        dformat = 'rowVector'
    else:
        print("Invalid format {}".format(format))
        sys.exit("Goodbye!")

    type = int(fields[3].strip())
    if type ==1:
        dtype = 'RealSimple'
    elif type ==2:
        dtype = 'RealDouble'
    elif type ==3:
        dtype = 'ComplexSimple'; print("Warning: cannot read complex OUTPUT4 matrices yet"); sys.exit("Goodbye!")
    elif type ==4:
        dtype = 'ComplexDouble'; print("Warning: cannot read complex OUTPUT4 matrices yet"); sys.exit("Goodbye!")
    else:
        print("Invalid format {}".format(type))
        sys.exit("Goodbye!")

    name = fields[4].strip()
    return ncols, nrows, dformat, dtype, name


def getColHeader(line):

    data = []
    fields = [None] * 3
    for ix in range(0,3+1):
        fields[ix-1] = line[ ((ix-1)*8+1):ix*8 ];
    colIdx = int(fields[0].strip())
    idxFirstNZRow = int(fields[1].strip())
    NoNZElems = int(fields[2].strip())

    return colIdx, idxFirstNZRow, NoNZElems



if __name__ == "__main__":


      ReadModalDisplacements()
      print("Mode 1")
























