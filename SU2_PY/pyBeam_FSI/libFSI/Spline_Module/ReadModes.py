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
    return void

def readModes(Modes, Mode_file, FORMAT_MODES, nModes):
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
               
        elif FORMAT_MODES == "NASTRAN":
            print('Opened mode file ' + Mode_file + '.' + " Format = " + FORMAT_MODES)
        else:
           print "CHE CAZZO DI FORMATO SAREBBE {}? Babbeo !!".format(FORMAT_MODES)    
    nModes.append(int(index))
    print("Total number of available structural modes from input file is: {}".format(int(nModes[0])))
    
'''   
  if __name__ == "__main__":
   
      Modes = [];
      Mode_file = "./Input/MODES_ICASE111.str"
      FORMAT_MODES = "CSHELL"
      nModes = []
      readModes(Modes, Mode_file, FORMAT_MODES, nModes)      
      '''