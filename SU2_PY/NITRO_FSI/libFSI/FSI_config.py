#!/usr/bin/env python

## \file FSI_config.py
#  \brief Python class for handling configuration file for FSI computation.
#  \author David Thomas
#  \version 6.2.0 "Falcon"
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

from util.switch import switch


def decomp_coeff_read(File_conf, array):
    f = open(File_conf["MODAL_COEFF_FILE_NAME"], "r+")

    for line in f:
        line = line.replace('i', 'j')  # .split()
        if line:
            # array.append(Point())
            # array[i].SetCoeff([ float(line[0]) float(line[1])  ])
            array.append(complex(line))


# ----------------------------------------------------------------------
#  FSI Configuration Class
# ----------------------------------------------------------------------

class FSIConfig:
    """
    Class that contains all the parameters coming from the FSI configuration file.
    Read the file and store all the options into a dictionary.
    """

    def __init__(self, FileName):
        self.ConfigFileName = FileName
        self._ConfigContent = {}
        self.readConfig()

    def __str__(self):
        tempString = str()
        for key, value in self._ConfigContent.items():
            tempString += "{} = {}\n".format(key, value)
        return tempString

    def __getitem__(self, key):
        return self._ConfigContent[key]

    def __setitem__(self, key, value):
        self._ConfigContent[key] = value

    def readConfig(self):
        input_file = open(self.ConfigFileName)
        while 1:
            line = input_file.readline()
            if not line:
                break
            # remove line returns
            line = line.strip('\r\n')
            # make sure it has useful data
            if (not "=" in line) or (line[0] == '%'):
                continue
            # split across equal sign
            line = line.split("=", 1)
            this_param = line[0].strip()
            this_value = line[1].strip()

            for case in switch(this_param):
                # integer values
                if case("DYN_COUPLED_TIMESTEP_NR")	: pass  # in case: CSD_SOLVER = NITRO_FRAMEWORK
                if case("BS_NR")		: pass  # in case: CSD_SOLVER = NITRO_FRAMEWORK
                if case("NMODES")		      : pass  # in case: CSD_SOLVER = NITRO_FRAMEWORK
                if case("MODE_TO_SIMULATE")	      : pass  # in case: CSD_SOLVER = NITRO_FRAMEWORK
                if case("ITER")	                      : pass  # in case: CSD_SOLVER = NITRO_FRAMEWORK
                if case("NDIM")			      : pass
                if case("TRACKING_NODE")              : pass
                if case("RESTART_ITER")		      : pass
                if case("UNST_NR")		      : pass
                if case("UNST_TOTAL_SIMUL_NUMBER")    : pass
                if case("NB_EXT_ITER")		      :
                    self._ConfigContent[this_param] = int(this_value)
                    break

                # float values
                if case("FREESTREAM_DENSITY")         : pass
                if case("FREESTREAM_PRESSURE")        : pass
                if case("K_MAX")                      : pass
                if case("L_REF")                      : pass
                if case("OMEGA")                      : pass
                if case("V_INF")                      : pass
                if case("START_MOTION_TIME")          : pass
                if case("START_TIME")		      : pass
                if case("UNST_TIMESTEP")	      : pass
                if case("UNST_TIME")		      : pass
                if case("BS_TIMESTEP_1")	      : pass
                if case("BS_TIMESTEP_2")	      :
                    self._ConfigContent[this_param] = float(this_value)
                    break

                # string values  MEMO_GEN_FORCE_OUTPUT
                if case("MODAL_DISPLACEMENT")         : pass
                if case("MODAL_DISPLACEMENT_FORMAT")  : pass
                if case("FREESTREAM_OPTION")          : pass
                if case("OUTPUT_DIRECTORY")           : pass
                if case("UNSTEADY_SCHEME")            : pass
                if case("MEMO_GEN_FORCE_OUTPUT")      : pass
                if case("MEMO_FORCE_OUTPUT")          : pass
                if case("WRITE_GEN_FORCE_OUTPUT")     : pass
                if case("GENERALIZED_FORCE_FILE")     : pass
                if case("MOTION_TYPE")                : pass
                if case("REAL_TIME_TRACKING")         : pass
                if case("WRITE_FORCE_OUTPUT")         : pass
                if case("NODAL_FORCE_FILE")           : pass
                if case("MODAL_COEFF_FILE_NAME")      : pass
                if case("MLS_CONFIG_FILE_NAME")       : pass
                if case("STRUCTURAL_MODES_FILE_NAME") : pass
                if case("FORMAT_MODES")	              : pass
                if case("CFD_CONFIG_FILE_NAME")	      : pass
                if case("CSD_SOLVER")		          : pass
                if case("CSD_CONFIG_FILE_NAME")	      : pass
                if case("RESTART_SOL")		          : pass
                if case("MATCHING_MESH")	          : pass
                if case("MESH_INTERP_METHOD")         : pass
                if case("DISP_PRED")		          : pass
                if case("UNSTEADY_SIMULATION")	      : pass
                if case("MESH_DEF_METHOD")	      : pass
                if case("INTERNAL_FLOW")	      :
                    self._ConfigContent[this_param] = this_value
                    break

                if case("MEMO_GEN_FORCE_OUTPUT")      :
                    print("\nWarning! MEMO_GEN_FORCE_OUTPUT is obsolete and not used any more" )

                if case():
                    print(this_param + " is an invalid option !")
                    break
                # end for

        if self._ConfigContent["CSD_SOLVER"] == 'NITRO':
             array = []
             decomp_coeff_read(self,array)
             if array:
                self._ConfigContent["MODAL_COEFF"] = array
