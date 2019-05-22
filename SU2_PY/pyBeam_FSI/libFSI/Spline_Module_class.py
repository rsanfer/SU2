#!/usr/bin/env python

## \file FSI_config.py
#  \brief Python class for handling configuration file for FSI computation.
#  \author Rocco Bombardieri, Ruben Sanchez
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


from libFSI import ReadModes as ReadModes
from libFSI import MLS_config as MLS_config
from libFSI import Plot_modes as Plot_modes
import scipy.io
import numpy as np
import sys
import pyMLS_Cpp as Spline

np.set_printoptions(threshold=sys.maxsize)


# ----------------------------------------------------------------------
#  MLS_Spline Interface Class
# ----------------------------------------------------------------------

class MLS_Spline:

    """
    MLS_Spline class that handles fluid/solid solver synchronisation and communication
    """

    def __init__(self, MLS_Config_File, nDim, AeroNodes, StructNodes, FSI_config): # NrAeroElem, BoundElem,
        """
         Class constructor. Declare some variables and do some screen outputs.
        """
        # AeroPoint is numpy matrix nAeropoint*3 (AeroPoint[markers[FSI_marker])
        # BoundElem is an object. here we need the function .GetNodes() or something similar
        # The MLS configurations parameters are stored from the MLS input file
        print("Storing MLS parameters from input file ")
        MLS_conf = MLS_config.MLSConfig(MLS_Config_File)

        self.nAeroNodes = np.shape(AeroNodes)[0]
        self.nStructNodes = np.shape(StructNodes)[0]

        lenAeroNodes = self.nAeroNodes * 3
        lenStructNodes = self.nStructNodes * 3

        # Performing the meshless method
        print("Performing the Meshless Method")
        # Arrange structural nodes in the wrapped standard vector

        str_data_std = Spline.DoubleVector(lenStructNodes)
        l = 0
        for i in range(0, 3):
            for j in range(0, self.nStructNodes):
                str_data_std[l] = float(StructNodes[j][i])  # str_data[j,i]
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

        # Print norm error
        print("Splining: norm of interpolation error over nodes position = {}".format(np.linalg.norm(norm_err_std)))

        if MLS_conf['DEBUG'] == "YES":

            # Storing structural modes from relative input file (this can be used for validation purposes)
            print("Storing structural modes from the input file ")
            print("NB: Remember we want structural modes to be mass normalized!")
            self.Modes = []  # It's an object and further elements will be "appended"
            print("Mode_file = {}".format(MLS_conf['STRUCTURAL_MODES_FILE_NAME']))
            self.nModes = []
            readModes(self.Modes, MLS_conf['STRUCTURAL_MODES_FILE_NAME'], MLS_conf['FORMAT_MODES'], self.nModes)
            self.nModes = int(self.nModes[0])

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
            # print(str_Matrix)

            # Connectivity = np.zeros((NrAeroElem, nDim))
            #
            # for i in range(0, NrAeroElem):
            #     for j in range(0, nDim):
            #         Connectivity[i][j] = int(BoundElem[i].GetNodes()[j])  # this has to be reviewed in

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
                plotmodes(X_mode, Y_mode, Z_mode, i)

            print("PRESS ENTER TO END PROGRAM.")
            wait = input("PROGRAM TERMINATED CORRECTLY.")
