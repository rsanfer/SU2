#!/usr/bin/env python3

## \file Dynresp_config.py
#  \brief Python class for handling Dynresp configuration file for FSI computation.
#  \author Rocco Bombardieri
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


# ----------------------------------------------------------------------
#  FSI Configuration Class
# ----------------------------------------------------------------------

class DYNData:
    """
    Class that contains all the parameters coming from the FSI configuration file.
    Read the file and store all the options into a dictionary.
    """

    def __init__(self, FileName):
        self.ConfigFileName = FileName
        self._ConfigContent = {}
        self.Meliminated = None
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
        length = 8
        input_file = open(self.ConfigFileName)
        print(input_file)
        while 1:
            line = input_file.readline()
            #print(line)
            # commented line
            if not line:
                break
            if (line[0] == '$'):
                continue
            # remove line returns
            line = line.strip('\r\n')
            if line[0:7].strip() == 'TIMEF':
               chunks = [line[i:i + length] for i in range(0, len(line), length)]
               # === In case the nastran file has exponential (-10 instead of e-10)
               for i in range(2,4+1):
                   posx = chunks[i].find('-')
                   if posx != -1:
                       if posx == 0:
                           xx = chunks[i][1:]
                           posxx = xx.find('-')
                           if posxx != -1:
                               chunks[i] = '-' + xx[:posxx] + 'e' + xx[posxx:]
                           else:
                               chunks[i] = '-' + xx
                       else:
                           chunks[i] = chunks[i][:posx] + 'e' + chunks[i][posx:]
                   else:
                       continue

               if chunks[2].strip():
                  self._ConfigContent['TW0'] = float(chunks[2].strip())
               else:
                  self._ConfigContent['TW0'] = float(0)

               if chunks[3].strip():
                  self._ConfigContent['TWF'] = float(chunks[3].strip())
               else:
                  self._ConfigContent['TWF'] = float(0)

               self._ConfigContent['DT'] = float(chunks[4].strip())

            if line[0:7].strip() == 'AERO':
               chunks = [line[i:i + length] for i in range(0, len(line), length)]
               # === In case the nastran file has exponential (-10 instead of e-10)
               for i in range(2,3+1):
                   posx = chunks[i].find('-')
                   if posx != -1:
                       if posx == 0:
                           xx = chunks[i][1:]
                           posxx = xx.find('-')
                           if posxx != -1:
                               chunks[i] = '-' + xx[:posxx] + 'e' + xx[posxx:]
                           else:
                               chunks[i] = '-' + xx
                       else:
                           chunks[i] = chunks[i][:posx] + 'e' + chunks[i][posx:]
                   else:
                       continue
               self._ConfigContent['V'] = float(chunks[2].strip())
               self._ConfigContent['RHO'] = float(chunks[3].strip())

            if line[0:7].strip() == 'SMODES':
               chunks = [line[i:i + length] for i in range(0, len(line), length)]
               if not chunks[1].strip():
                   self._ConfigContent['NMODE'] = int(0)
               else:
                   self._ConfigContent['NMODE'] = int(chunks[1].strip())

               if len(chunks) >= 4 and chunks[3].strip():
                  self._ConfigContent['MLIST1'] = int(chunks[3].strip())
               else:
                  self._ConfigContent['MLIST1'] = None


        # If there is a list of modes to eliminate (MLIST1)
        if self._ConfigContent['MLIST1'] is not None:
           next_line = True
           modes = []
           input_file.seek(0)
           input_file.seek(0)
           while 1:
            line = input_file.readline()
            if not line:
                break
            if (line[0] == '$'):
                continue
            # remove line returns
            line = line.strip('\r\n')
            # Look for SET1 card with MLIST1 index
            if line[0:7].strip() == 'SET1':
               chunks = [line[i:i + length] for i in range(0, len(line), length)]
               print(chunks)
               # if this is the SET1 card I'm looking for I read the modes
               if int(chunks[1].strip()) == self._ConfigContent['MLIST1']:
                  # reading modes
                  # Needs to check if there are more lines and read the modes of all the lines
                  pos = len(chunks)
                  for i in range (2,min(pos,9)): # in case there are more lines
                      modes.append(int(chunks[i]))
                  while 1:
                      # If there is s asecond line, read that one also
                      if pos != 10:
                         break

                      line = input_file.readline()
                      chunks = [line[i:i + length] for i in range(0, len(line), length)]
                      pos = len(chunks)
                      for i in range(1, min(pos ,9)):
                          modes.append(int(chunks[i]))

               else:
                   continue
           self.Meliminated = modes
        else:
           self.Meliminated = []

'''
if __name__ == "__main__":

      DYNConfig = None;
      file = "./GTA_gust_maneu_D_.inp"

      Dyn_config = DYNData(file)

      print("Dyn_config[TW0] = {}".format(Dyn_config['TW0']))
      print("Dyn_config[TWF] = {}".format(Dyn_config['TWF']))
      print("Dyn_config[DT] = {}".format(Dyn_config['DT']))
      print("Dyn_config[V] = {}".format(Dyn_config['V']))
      print("Dyn_config[RHO] = {}".format(Dyn_config['RHO']))
      print("Dyn_config[NMODE] = {}".format(Dyn_config['NMODE']))
      print("Dyn_config[MLIST1] = {}".format(Dyn_config['MLIST1']))
      print("Dyn_config[.Neliminated = {}".format(Dyn_config.Meliminated))
'''