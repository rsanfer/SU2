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
                  index = int( test_string.plit('R E A L   E I G E N V E C T O R   N O .',1)[1] )
                  # now go ahead till you find
                  Modes.append(CModes_f06())
                  while line.find('G') == -1:
                    line = modefile.readline()
                    if not line:
                       break
                  while line.find('  G  ') == 1:
                    line = line.split()
                    node_disp = [float(line[0]), float(line[2]), float(line[3]), float(line[4]), float(line[5]), float(line[6]), float(line[7])]
                    Modes[index].SetModeLine(node_disp)
        else:
           print("Format not known !!".format(FORMAT_MODES))
           sys.exit("Goodbye!")
    nModes.append(int(index))
    print("Total number of available structural modes from input file is: {}".format(int(nModes[0])))
    
























