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
    self.Coord = np.zeros((3,1))
    
  def SetCoord(self,val_Coord):
    x, y, z = val_Coord
    self.Coord[0] = x
    self.Coord[1] = y
    self.Coord[2] = z
    
  def GetCoord(self):
    return self.Coord  


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
                       x = line[24:32].strip(); y = line[32:40].strip(); z = line[40:48].strip();
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
    Mesh_file = "./Input/Modal_20_span.bdf"#MLS_conf['STRUCTURAL_NODES_FILE_NAME']
    Mesh_format = "NASTRAN" # MLS_conf['FORMAT_SRUCT_NODES']
    nStrPoint = []
    StructNodes = []   
    node = []
    nPoint = []
    readStructMesh(Mesh_file, Mesh_format, node, nPoint)    
    '''
