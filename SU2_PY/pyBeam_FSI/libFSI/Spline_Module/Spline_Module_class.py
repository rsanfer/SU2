from petsc4py import PETSc
from ReadModes import *
from ReadSU2Mesh import *
from MLS_config import *
from ReadStructMesh import *
from Plot_modes import *
import scipy.io 
import pyMLS_Cpp as Spline
np.set_printoptions(threshold=sys.maxsize)

# ----------------------------------------------------------------------
#  MLS_Spline Interface Class
# ----------------------------------------------------------------------

class MLS_Spline:
    """ 
    MLS_Spline class that handles fluid/solid solver synchronisation and communication
    """
   
    def __init__(self, MLS_Config_File, FSI_config, SolidSolver):
	""" 
	Class constructor. Declare some variables and do some screen outputs.
	"""
        # DEBUG
        #MLS_Config_File = "./Input/MLS_config.cfg"
        #Mesh_file = "./Input/mesh_ONERAM6.su2"   # this has to be given from the main FSI cfg file  TO BE REMOVED
        #self.omega = MLS_conf['FREQ']
        #Mode_file = MLS_conf['STRUCTURAL_MODES_FILE_NAME']
        #FORMAT_MODES = MLS_conf['FORMAT_MODES']
        
        # The MLS configurations parameters are stored from the MLS input file 
        print("Storing MLS parameters from input file ")
        MLS_conf = MLSConfig(MLS_Config_File)
        
    
        # Storing structural modes from relative input file
        Mesh_file = SolidSolver.Config['MESH_FILE']
        print("Storing structural modes from the input file ")
        print("NB: Remember we want structural modes to be mass normalized!")
        self.Modes = [] # It's an object and further elements will be "appended"
        Mode_file = FSI_config['STRUCTURAL_MODES_FILE_NAME']
        print("Mode_file = {}".format(Mode_file))
        FORMAT_MODES = FSI_config['FORMAT_MODES']
        self.nModes = []
        readModes(self.Modes, Mode_file, FORMAT_MODES, self.nModes)
        self.nModes = int(self.nModes[0])
            
        # Storing aerodynamic mesh information from the SU2 mesh file    
        print("Storing aerodynamic mesh information from the SU2 mesh file ")
        # SolidSolver.Config["MOVING_MARKER"]
        FSI_marker = SolidSolver.Config["MOVING_MARKER"]#"WING"                    # this has to be given from the main structural cfg file
        AeroPoint = []                        # It's an object and further elements will be "appended"
        BoundElem = []                        # It's an object and further elements will be "appended"
        markers = {}                          # It's an object and further elements will be "appended"
        SU2mesh_info = []
         # readSU2Mesh(Mesh_file, AeroPoint, BoundElem, markers, FSI_marker, SU2mesh_info):
        readSU2Mesh(Mesh_file, AeroPoint, BoundElem, markers, FSI_marker, SU2mesh_info)
        NrAeroPoint = SU2mesh_info[0].nPointBound
        NrAeroElem = SU2mesh_info[0].nElem
        nDim = SU2mesh_info[0].nDim
        # Storing structural nodes
        print("Storing Structural mesh information from the input file")
        Mesh_file = MLS_conf['STRUCTURAL_NODES_FILE_NAME']
        Mesh_format = MLS_conf['FORMAT_SRUCT_NODES']
        self.nStrPoint = []
        StructNodes = []
        readStructMesh(Mesh_file, Mesh_format, StructNodes, self.nStrPoint)
        self.nStrPoint = int(self.nStrPoint[0])
        # Performing the meshless method
        print("Performing the Meshless Method")
        # Arrange structural nodes in the wrapped standard vector
        str_data_std = Spline.DoubleVector(self.nStrPoint*3) 
        l=0
        for i in range(0,3): 
          for j in range(0,self.nStrPoint):   
           str_data_std[l] = round(float(StructNodes[j].GetCoord()[i])*pow(10,5))/pow(10,5) #  str_data[j,i]
           l = l+1
        # Arrange aerodynamic nodes in the wrapped standard vector
        aero_data_std = Spline.DoubleVector(NrAeroPoint*3)
        l=0
        for i in range(0,3):
           for j in range(0,NrAeroPoint):
              aero_data_std[l] = float(AeroPoint[markers[FSI_marker][j]].GetCoord()[i]) # aero_data[j,i]     The j-th node of the boundary is given by AeroPoint[markers[FSI_marker][j]      
              l = l+1
          
         
        interpolation_matrix_std = Spline.DoubleVector(NrAeroPoint*self.nStrPoint)
        norm_err_std = Spline.DoubleVector(NrAeroPoint)

        #print(markers[FSI_marker]) 
        
        Spline.MLS(interpolation_matrix_std, norm_err_std, self.nStrPoint, NrAeroPoint, str_data_std, aero_data_std,
                   MLS_conf['POLY'], MLS_conf['WEIGHT'], MLS_conf['POINTS'],
                   MLS_conf['RMAX'], MLS_conf['DELTA'], MLS_conf['TOLL_SVD'])

        # --- OUTPUT ----------------------------------------------------------------
        self.interpolation_matrix = np.zeros((NrAeroPoint,self.nStrPoint))
        l=0
        for i in range(0,self.nStrPoint):
           for j in range(0,NrAeroPoint):
               self.interpolation_matrix[j][i] = interpolation_matrix_std[l]
               l=l+1
      
        if MLS_conf['DEBUG'] == "YES":
         
         # Print norm error
         print("Splining: max of interpolation error over nodes position = {}".format(np.asarray(norm_err_std)))

           
         # --- PLOTTING ---------------------------------------------------------------
         Aero_Matrix = np.zeros((NrAeroPoint,3))
         l=0
         for i in range(0,3):
           for j in range(0,NrAeroPoint):
              Aero_Matrix[j][i] = aero_data_std[l]
              l=l+1
         #print(Aero_Matrix)
         #scipy.io.savemat('GRID_python.mat', {'GRID' : Aero_Matrix})

         str_Matrix = np.zeros((self.nStrPoint,3))
         l=0
         for i in range(0,3):
           for j in range(0,self.nStrPoint):
              str_Matrix[j][i] = str_data_std[l]
              l=l+1
         #print(str_Matrix)
    
         Connectivity = np.zeros((NrAeroElem,nDim))

         for i in range(0,NrAeroElem):
            for j in range(0,nDim):
                Connectivity[i][j] = int(BoundElem[i].GetNodes()[j]) #
         #print(Connectivity)
         '''
         scipy.io.savemat('GRID_python.mat', {'GRID' : Aero_Matrix})
         scipy.io.savemat('R_python.mat', {'R' : str_Matrix})
         scipy.io.savemat('Interpolation_matrix.mat', {'interpolation_matrix' : interpolation_matrix})
         '''
         #------ Plotting options   -------------------

         # error on the interpolation plot
     
         fig = plt.figure(0)
         X = np.linspace(0,NrAeroPoint-1,NrAeroPoint)
     
         plt.plot( X, np.asarray(norm_err_std), 'ro')

         plt.xlabel('Query nodes')
         plt.ylabel('[%] Error')
         plt.title('Interpolation error');
         plt.grid()
         plt.draw()

         # plotting base configuration
         #plotmodes(Aero_Matrix[:,0], Aero_Matrix[:,1], Aero_Matrix[:,2])
    
         # Plotting modes
    
         for i in range(0,FSI_config['NMODES']):
           X_mode = self.interpolation_matrix.dot( str_Matrix[:,0] + self.Modes[i].GetMode()[:,0]*MLS_conf['MAGNIF_FACTOR'] )   # 
           Y_mode = self.interpolation_matrix.dot( str_Matrix[:,1] + self.Modes[i].GetMode()[:,1]*MLS_conf['MAGNIF_FACTOR'] )   #  
           Z_mode = self.interpolation_matrix.dot( str_Matrix[:,2] + self.Modes[i].GetMode()[:,2]*MLS_conf['MAGNIF_FACTOR'] )   # 
           plotmodes(X_mode, Y_mode, Z_mode,i)

         print("PRESS ENTER TO END PROGRAM.")
         wait = input("PROGRAM TERMINATED CORRECTLY.")
         '''  DEBUGGING
         scipy.io.savemat('X_modepy.mat', {'X_modepy' : X_mode})
         scipy.io.savemat('Y_modepy.mat', {'Y_modepy' : Y_mode})
         scipy.io.savemat('Z_modepy.mat', {'Z_modepy' : Z_mode})
         '''
    

    

