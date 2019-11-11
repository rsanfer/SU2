#!/usr/bin/env python

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
import pyNastran
pkg_path = pyNastran.__path__[0]

from pyNastran.utils import print_bad_path
from pyNastran.op4.op4 import read_op4
import numpy as np
from numpy import float32, float64, int32, int64, product

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

def ReadModalDisplacements():  # FSIConfig
    T0 = 0 #FSIConfig['START_TIME']
    TF = 20 #FSIConfig['UNST_TIME']
    TW = TF - T0
    DT = 0.005 #FSIConfig['UNST_TIMESTEP']
    NT = TW / DT + 1

    op4_filename = './gta_d_d_dff__5_old.f48'
    name ='MODALXVG'
    # read martrices
    #help(read_op4)
    matrices =read_op4(op4_filename)#,matrix_names='MODALXVG',  precision='double', debug=True, log=None)



if __name__ == "__main__":


      ReadModalDisplacements()
      print("Mode 1")
























