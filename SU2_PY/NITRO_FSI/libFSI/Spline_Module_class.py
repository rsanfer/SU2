#!/usr/bin/env python

## \file FSI_config.py
#  \brief Python class for handling configuration file for FSI computation.
#  \author Rocco Bombardieri, Ruben Sanchez
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

import pyMLS as Spline
import libFSI.pyMLSConfig as io
from libFSI.ReadModes import ReadModes
from libFSI.Plot_modes import Plot_modes
from libFSI.ReadStructMesh import ReadStructMesh
#import scipy.io
import numpy as np
import sys

np.set_printoptions(threshold=sys.maxsize)


# ----------------------------------------------------------------------
#  MLS_Spline Interface Class
# ----------------------------------------------------------------------

class MLS_Spline:

    """
    MLS_Spline class that handles fluid/solid solver synchronisation and communication
    """

    def __init__(self, MLS_Config_File, nDim, AeroNodes, FSI_config, DYN_config = None): # NrAeroElem, BoundElem,
        """
         Class constructor. Declare some variables and do some screen outputs.
        """
        # AeroPoint is numpy matrix nAeropoint*3 (AeroPoint[markers[FSI_marker])
        # BoundElem is an object. here we need the function .GetNodes() or something similar
        # The MLS configurations parameters are stored from the MLS input file
        print("Storing MLS parameters from input file ")
        MLS_conf = io.pyMLSConfig(MLS_Config_File)

        # Storing structural modes from relative input file
        print("Storing structural modes from the input file ")
        #print("NB: Remember we want structural modes to be mass normalized!")
        self.Modes = []; # It's an object and further elements will be "appended"
        Mode_file = FSI_config['STRUCTURAL_MODES_FILE_NAME']
        print("Mode_file = {}".format(Mode_file))
        FORMAT_MODES = FSI_config['FORMAT_MODES']
        self.nModes = []
        ReadModes(self.Modes, Mode_file, FORMAT_MODES, self.nModes)
        self.nModes = int(self.nModes[0])

        # If dynresp
        if DYN_config != None:
           # First select the max number of modes to be used
           if DYN_config['NMODE'] != int(0):
              nModes =  DYN_config['NMODE']
              list = [i for i in range(nModes ,self.nModes)]
              for index in sorted(list, reverse=True):
                  del self.Modes[index]
              self.nModes = nModes
              print(" Max number of modes before elimination by MLIST1 is: {}".format(int(self.nModes)))
           # Then, elimination of selected modes in SMODES card
           if DYN_config.Meliminated != None:
              nModes = self.nModes - len(DYN_config.Meliminated)
              for index in sorted(DYN_config.Meliminated, reverse=True):
                  del self.Modes[index - 1]
              self.nModes = nModes
              print(" Further elimination of modes: {}".format(DYN_config.Meliminated))

        # aerodynamic boundary nodes
        self.nAeroNodes = np.shape(AeroNodes)[0]
        lenAeroNodes = self.nAeroNodes * 3


        print("Storing Structural mesh information from the input file")
        Mesh_file = MLS_conf['STRUCTURAL_NODES_FILE_NAME']
        Mesh_format = MLS_conf['FORMAT_SRUCT_NODES']
        self.nStructNodes = []
        StructNodes = []
        ReadStructMesh(Mesh_file, Mesh_format, StructNodes, self.nStructNodes)
        self.nStructNodes = int(self.nStructNodes[0])
        lenStructNodes = self.nStructNodes * 3

        # if structural mesh comes from Nastran, nodes need to be reordered in ascending order given their Id. (to match with modes displacements)
        if (MLS_conf['FORMAT_SRUCT_NODES'] == 'NASTRAN' and FSI_config['FORMAT_MODES'] == 'NASTRANF06' ):
           idList = [StructNodes[i].GetId() for i in range(0,self.nStructNodes)]
           Order = np.argsort(np.argsort(idList))
           reorder(StructNodes, Order, self.nStructNodes)

        # Performing the meshless method
        print("Performing the Meshless Method")
        # Arrange structural nodes in the wrapped standard vector

        str_data_std = Spline.DoubleVector(lenStructNodes)
        l = 0
        for i in range(0, 3):
            for j in range(0, self.nStructNodes):
                str_data_std[l] = float(StructNodes[j].GetCoord()[i])  # str_data[j,i]
                l = l + 1

        # Arrange aerodynamic nodes in the wrapped standard vector
        aero_data_std = Spline.DoubleVector(lenAeroNodes)
        l = 0
        for i in range(0, 3):
            for j in range(0, self.nAeroNodes):
                aero_data_std[l] = float(AeroNodes[j][i])
                l = l + 1

        interpolation_matrix_std = Spline.DoubleVector(self.nAeroNodes * self.nStructNodes)
        norm_err_std = Spline.DoubleVector(self.nAeroNodes)

        # check
        if MLS_conf['POINTS'] > self.nStructNodes:
           print("Warning! Nr. of query points requested for MLS is {}. Nr. of availavle points is {}. Set POINTS to {}".format(MLS_conf['POINTS'],self.nStructNodes,self.nStructNodes ))
           MLS_conf['POINTS'] = self.nStructNodes
        '''
        /print("str_Matrix = {}".format(StructNodes))
        print("Aero_Matrix = {}".format(AeroNodes))
        print("NrAeroPoint = {}".format(self.nAeroNodes))
        print("self.nStrPoint = {}".format(self.nStructNodes))
        print("MLS_conf['POLY'] = {}".format(MLS_conf['POLY']))
        print("MLS_conf['WEIGHT'] = {}".format(MLS_conf['WEIGHT']))
        print("MLS_conf['RMAX'] = {}".format(MLS_conf['RMAX']))
        print("MLS_conf['DELTA'] = {}".format(MLS_conf['DELTA']))     
        print("MLS_conf['TOLL_SVD'] = {}".format(MLS_conf['TOLL_SVD'])) 
        '''
        Spline.mls_interface(interpolation_matrix_std, norm_err_std, self.nStructNodes, self.nAeroNodes, str_data_std,
                             aero_data_std, MLS_conf['POLY'], MLS_conf['WEIGHT'], MLS_conf['POINTS'],
                             MLS_conf['RMAX'], MLS_conf['DELTA'], MLS_conf['TOLL_SVD'])

        # --- OUTPUT ----------------------------------------------------------------
        self.interpolation_matrix = np.zeros((self.nAeroNodes, self.nStructNodes))
        l = 0
        for i in range(0, self.nStructNodes):
            for j in range(0, self.nAeroNodes):
                self.interpolation_matrix[j][i] = interpolation_matrix_std[l]
                l = l + 1

        print("Splining: norm of interpolation error over nodes position = {}".format(np.linalg.norm(norm_err_std)))

        if MLS_conf['DEBUG'] == "YES":

            # --- PLOTTING ---DEBUGGING options------------------------------------------------------------
            Aero_Matrix = np.zeros((self.nAeroNodes, 3))
            l = 0
            for i in range(0, nDim):
                for j in range(0, self.nAeroNodes):
                    Aero_Matrix[j][i] = aero_data_std[l]
                    l = l + 1

            str_Matrix = np.zeros((self.nStructNodes, 3))
            l = 0
            for i in range(0, nDim):
                for j in range(0, self.nStructNodes):
                    str_Matrix[j][i] = str_data_std[l]
                    l = l + 1

            # ------ Plotting options   -------------------

            # error on the interpolation plot

            fig = plt.figure(0)
            X = np.linspace(0, self.nAeroNodes - 1, self.nAeroNodes)

            plt.plot(X, np.asarray(norm_err_std), 'ro')

            plt.xlabel('Query nodes')
            plt.ylabel('[%] Error')
            plt.title('Interpolation error');
            plt.grid()
            plt.draw()

            # plotting base configuration
            # plotmodes(Aero_Matrix[:,0], Aero_Matrix[:,1], Aero_Matrix[:,2])

            # Plotting modes

            for i in range(0, FSI_config['NMODES']):
                X_mode = self.interpolation_matrix.dot(
                    str_Matrix[:, 0] + self.Modes[i].GetMode()[:, 0] * MLS_conf['MAGNIF_FACTOR'])  #
                Y_mode = self.interpolation_matrix.dot(
                    str_Matrix[:, 1] + self.Modes[i].GetMode()[:, 1] * MLS_conf['MAGNIF_FACTOR'])  #
                Z_mode = self.interpolation_matrix.dot(
                    str_Matrix[:, 2] + self.Modes[i].GetMode()[:, 2] * MLS_conf['MAGNIF_FACTOR'])  #
                Plot_modes(X_mode, Y_mode, Z_mode, i)

            print("PRESS ENTER TO END PROGRAM.")
            wait = input("PROGRAM TERMINATED CORRECTLY.")


def reorder(arr, index, n):
    temp = [0] * n;

    # arr[i] should be
    # present at index[i] index
    for i in range(0, n):
        temp[index[i]] = arr[i]

        # Copy temp[] to arr[]
    for i in range(0, n):
        arr[i] = temp[i]
        index[i] = i