#!/usr/bin/env python

## \file Interface.py
#  \brief Interface class that handles fluid/solid solvers synchronisation and communication.
#  \author Rocco Bombardieri, Ruben Sanchez based on previous work by David Thomas.
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


# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import numpy as np
import shutil
import os
# ----------------------------------------------------------------------
#  FSI Interface Class
# ----------------------------------------------------------------------


class Interface:

    """
    FSI interface class that handles fluid/solid solvers synchronisation and communication
    """

    def __init__(self, FSI_config, have_MPI):
        """
        Class constructor. Declare some variables and do some screen outputs.
        """

        if have_MPI:
            from mpi4py import MPI
            self.MPI = MPI
            self.comm = MPI.COMM_WORLD  # MPI World communicator
            self.have_MPI = True
            myid = self.comm.Get_rank()
        else:
            self.comm = 0
            self.have_MPI = False
            myid = 0

        self.rootProcess = 0  # the root process is chosen to be MPI rank = 0 (NOTE: it's always the same which holds the structural solver)

        self.nDim = FSI_config['NDIM']  # problem dimension

        self.haveFluidSolver = False  # True if the fluid solver is initialized on the current rank
        self.haveSolidSolver = False  # True if the solid solver is initialized on the current rank
        self.haveFluidInterface = False  # True if the current rank owns at least one fluid interface node
        self.haveSolidInterface = False  # True if the current rank owns at least one solid interface node

        self.fluidSolverProcessors = list()  # list of partitions where the fluid solver is initialized
        #self.solidSolverProcessors = list()  # list of partitions where the solid solver is initialized
        self.fluidInterfaceProcessors = list()  # list of partitions where there are fluid interface nodes
        self.solidInterfaceProcessors = list()  # list of partitions where there are solid interface nodes

        self.fluidInterfaceIdentifier = None  # object that can identify the f/s interface within the fluid solver
        self.solidInterfaceIdentifier = None  # object that can identify the f/s interface within the solid solver

        self.fluidGlobalIndexRange = {}  # contains the global FSI indexing of each fluid interface node for all partitions
        self.solidGlobalIndexRange = {}  # contains the global FSI indexing of each solid interface node for all partitions

        self.FluidHaloNodeList = {}  # contains the the indices (fluid solver indexing) of the halo nodes for each partition
        self.fluidIndexing = {}  # links between the fluid solver indexing and the FSI indexing for the interface nodes
        self.SolidHaloNodeList = {}  # contains the the indices (solid solver indexing) of the halo nodes for each partition
        self.solidIndexing = {}  # links between the solid solver indexing and the FSI indexing for the interface nodes

        self.nLocalFluidInterfaceNodes = 0  # number of nodes (halo nodes included) on the fluid interface, on each partition
        self.nLocalFluidInterfaceHaloNode = 0  # number of halo nodes on the fluid intrface, on each partition
        self.nLocalFluidInterfacePhysicalNodes = 0  # number of physical (= non halo) nodes on the fluid interface, on each partition
        self.nFluidInterfaceNodes = 0  # number of nodes on the fluid interface, sum over all the partitions
        self.nFluidInterfacePhysicalNodes = 0  # number of physical nodes on the fluid interface, sum over all partitions

        self.nLocalSolidInterfaceNodes = 0  # number of physical nodes on the solid interface, on each partition
        self.nLocalSolidInterfaceHaloNode = 0  # number of halo nodes on the solid intrface, on each partition
        self.nLocalSolidInterfacePhysicalNodes = 0  # number of physical (= non halo) nodes on the solid interface, on each partition
        self.nSolidInterfaceNodes = 0  # number of nodes on the solid interface, sum over all partitions
        self.nSolidInterfacePhysicalNodes = 0  # number of physical nodes on the solid interface, sum over all partitions

        self.MappingMatrix = None  # interpolation/mapping matrix for meshes interpolation/mapping - interpolation happening with MLS module
        self.MappingMatrix_T = None  # transposed interpolation/mapping matrix for meshes interpolation/mapping - interpolation happening with MLS module
        self.d_RBF = 0

        self.localFluidInterface_vertex_indices = None

        self.globalFluidInterfaceXcoor = None
        self.globalFluidInterfaceYcoor = None
        self.globalFluidInterfaceZcoor = None

        self.globalFluidCoordinates0 = None
        self.globalFluidCoordinates0X = None
        self.globalFluidCoordinates0Y = None
        self.globalFluidCoordinates0Z = None


        self.sendCounts = None
        self.globalFluidDispX = None
        self.globalFluidDispY = None
        self.globalFluidDispZ = None

        self.globalSolidDispX = None
        self.globalSolidDispY = None
        self.globalSolidDispZ = None

        self.globalSolidDispXOld = None
        self.globalSolidDispYOld = None
        self.globalSolidDispZOld = None

        self.globalFluidLoadX = None
        self.globalFluidLoadY = None
        self.globalFluidLoadZ = None

        self.globalSolidLoadX = None
        self.globalSolidLoadY = None
        self.globalSolidLoadZ = None

        self.haloNodesPositionsInit = {}  # initial position of the halo nodes (fluid side only)

        self.solidInterface_array_DispX = None  # solid interface displacement
        self.solidInterface_array_DispY = None
        self.solidInterface_array_DispZ = None

        self.solidInterfaceResidual_array_X = None  # solid interface position residual
        self.solidInterfaceResidual_array_Y = None
        self.solidInterfaceResidual_array_Z = None

        self.solidInterfaceResidualnM1_array_X = None  # solid interface position residual at the previous BGS iteration
        self.solidInterfaceResidualnM1_array_Y = None
        self.solidInterfaceResidualnM1_array_Z = None

        self.fluidInterface_array_DispX = None  # fluid interface displacement
        self.fluidInterface_array_DispY = None
        self.fluidInterface_array_DispZ = None

        self.fluidLoads_array_X = None  # loads on the fluid side of the f/s interface
        self.fluidLoads_array_Y = None
        self.fluidLoads_array_Z = None

        self.solidLoads_array_X = None  # loads on the solid side of the f/s interface
        self.solidLoads_array_Y = None
        self.solidLoads_array_Z = None

        self.FSIIter = 0  # current FSI iteration
        self.unsteady = False  # flag for steady or unsteady simulation (default is steady)

        # ---Some screen output ---
        self.MPIPrint('Fluid solver : SU2_CFD')
        self.MPIPrint('Matching fluid-solid interface using Moving Least Squares method')
        self.MPIPrint('FSI interface is set')

    def MPIPrint(self, message):
        """
        Print a message on screen only from the master process.
        """

        if self.have_MPI:
            myid = self.comm.Get_rank()
        else:
            myid = 0

        if myid == self.rootProcess:
            print(message)

    def MPIBarrier(self):
        """
        Perform a synchronization barrier in case of parallel run with MPI.
        """

        if self.have_MPI:
            self.comm.barrier()

    def checkMPI(self):
        """
        Return the MPI characteristics of the problem
        """

        if self.have_MPI:
            myid = self.comm.Get_rank()
            MPIsize = self.comm.Get_size()
        else:
            myid = 0
            MPIsize = 1

        return myid, MPIsize

    def connect(self, FSI_config, FluidSolver, SolidSolver):
        """
        Connection between solvers.
        Creates the communication support between the two solvers.
        Gets information about f/s interfaces from the two solvers.
        """

        # Recover the process and the size of the parallelization
        myid, MPIsize = self.checkMPI()

        # --- Identify the fluid interface and store the number of nodes for each partition ---#
        self.fluidInterfaceIdentifier = None
        self.nLocalFluidInterfaceNodes = 0
        if FluidSolver != None:
            print('Fluid solver is initialized on process {}'.format(myid))
            self.haveFluidSolver = True
            allInterfaceMarkersTags = FluidSolver.GetAllDeformMeshMarkersTag()
            allMarkersID = FluidSolver.GetAllBoundaryMarkers()
            if not allInterfaceMarkersTags:
                raise Exception('No moving marker was defined in SU2.')
            else:
                if allInterfaceMarkersTags[0] in allMarkersID.keys():
                    self.fluidInterfaceIdentifier = allMarkersID[allInterfaceMarkersTags[0]]
            if self.fluidInterfaceIdentifier != None:
                self.nLocalFluidInterfaceNodes = FluidSolver.GetNumberVertices(self.fluidInterfaceIdentifier)
            if self.nLocalFluidInterfaceNodes != 0:
                self.haveFluidInterface = True
                print('Number of interface fluid nodes (halo nodes included) on proccess {} and marker {}: {}'\
                      .format(myid,allInterfaceMarkersTags[0],self.nLocalFluidInterfaceNodes))
        else:
            pass

        # --- Exchange information about processors --- #
        if self.have_MPI:
            if self.haveFluidSolver:
                sendBufFluid = np.array(int(1))
            else:
                sendBufFluid = np.array(int(0))
            if self.haveFluidInterface:
                sendBufFluidInterface = np.array(int(1))
            else:
                sendBufFluidInterface = np.array(int(0))
            rcvBufFluid = np.zeros(MPIsize, dtype=int)
            rcvBufFluidInterface = np.zeros(MPIsize, dtype=int)
            self.comm.Allgather(sendBufFluid, rcvBufFluid)
            self.comm.Allgather(sendBufFluidInterface, rcvBufFluidInterface)

            for iProc in range(MPIsize):
                if rcvBufFluid[iProc] == 1:
                    self.fluidSolverProcessors.append(iProc)
                if rcvBufFluidInterface[iProc] == 1:
                    self.fluidInterfaceProcessors.append(iProc)

            del sendBufFluid, rcvBufFluid, sendBufFluidInterface, rcvBufFluidInterface
        else:
            self.fluidSolverProcessors.append(0)
            self.fluidInterfaceProcessors.append(0)

        self.MPIBarrier()

        # --- Calculate the total number of nodes at the fluid interface (sum over all the partitions)
        # --- Calculate the number of halo nodes on each partition
        self.nLocalFluidInterfaceHaloNode = 0
        for iVertex in range(self.nLocalFluidInterfaceNodes):
            if FluidSolver.IsAHaloNode(self.fluidInterfaceIdentifier, iVertex) == True:
                GlobalIndex = FluidSolver.GetVertexGlobalIndex(self.fluidInterfaceIdentifier, iVertex)
                self.FluidHaloNodeList[GlobalIndex] = iVertex
                self.nLocalFluidInterfaceHaloNode += 1

        # Calculate the number of physical (= not halo) nodes on each partition
        self.nLocalFluidInterfacePhysicalNodes = self.nLocalFluidInterfaceNodes - self.nLocalFluidInterfaceHaloNode
        if self.have_MPI == True:
            self.FluidHaloNodeList = self.comm.allgather(self.FluidHaloNodeList)
        else:
            self.FluidHaloNodeList = [{}]

        # --- Calculate the total number of nodes (with and without halo) at the fluid interface (sum over all the partitions) and broadcast the number accross all processors ---
        sendBuffHalo = np.array(int(self.nLocalFluidInterfaceNodes))
        sendBuffPhysical = np.array(int(self.nLocalFluidInterfacePhysicalNodes))
        rcvBuffHalo = np.zeros(1, dtype=int)
        rcvBuffPhysical = np.zeros(1, dtype=int)
        if self.have_MPI:
            self.comm.barrier()
            self.comm.Allreduce(sendBuffHalo, rcvBuffHalo, op=self.MPI.SUM)
            self.comm.Allreduce(sendBuffPhysical, rcvBuffPhysical, op=self.MPI.SUM)
            self.nFluidInterfaceNodes = rcvBuffHalo[0]
            self.nFluidInterfacePhysicalNodes = rcvBuffPhysical[0]
        else:
            self.nFluidInterfaceNodes = np.copy(sendBuffHalo)
            self.nFluidInterfacePhysicalNodes = np.copy(sendBuffPhysical)
        del sendBuffHalo, rcvBuffHalo, sendBuffPhysical, rcvBuffPhysical

        # --- Store the number of physical interface nodes on each processor and allgather the information ---
        self.fluidPhysicalInterfaceNodesDistribution = np.zeros(MPIsize, dtype=int)
        if self.have_MPI:
            sendBuffPhysical = np.array(int(self.nLocalFluidInterfacePhysicalNodes))
            self.comm.Allgather(sendBuffPhysical, self.fluidPhysicalInterfaceNodesDistribution)
            del sendBuffPhysical
        else:
            self.fluidPhysicalInterfaceNodesDistribution[0] = self.nFluidInterfacePhysicalNodes

        # --- Calculate and store the global indexing of interface physical nodes on each processor and allgather the information ---
        if self.have_MPI:
            if myid in self.fluidInterfaceProcessors:
                globalIndexStart = 0
                for iProc in range(myid):
                    globalIndexStart += self.fluidPhysicalInterfaceNodesDistribution[iProc]
                globalIndexStop = globalIndexStart + self.nLocalFluidInterfacePhysicalNodes - 1
            else:
                globalIndexStart = 0
                globalIndexStop = 0
            self.fluidGlobalIndexRange[myid] = [globalIndexStart, globalIndexStop]
            self.fluidGlobalIndexRange = self.comm.allgather(self.fluidGlobalIndexRange)
        else:
            temp = {}
            temp[0] = [0, self.nLocalFluidInterfacePhysicalNodes - 1]
            self.fluidGlobalIndexRange = list()
            self.fluidGlobalIndexRange.append(temp)

        self.MPIBarrier()

        # --- Get, for the fluid interface on each partition:
        # --- The vertex indices, which are stored locally on each processor
        # --- The coordinates X, Y and Z, which are stored in a global list in processor 0
        GlobalIndex = int()
        localIndex = 0
        fluidIndexing_temp = {}
        self.localFluidInterface_vertex_indices = np.zeros(self.nLocalFluidInterfacePhysicalNodes, dtype=int)
        localFluidInterface_array_X_init = np.zeros(self.nLocalFluidInterfacePhysicalNodes)
        localFluidInterface_array_Y_init = np.zeros(self.nLocalFluidInterfacePhysicalNodes)
        localFluidInterface_array_Z_init = np.zeros(self.nLocalFluidInterfacePhysicalNodes)

        for iVertex in range(self.nLocalFluidInterfaceNodes):
            GlobalIndex = FluidSolver.GetVertexGlobalIndex(self.fluidInterfaceIdentifier, iVertex)
            posx, posy, posz = FluidSolver.GetVertex_UndeformedCoord(self.fluidInterfaceIdentifier, iVertex)

            if GlobalIndex in self.FluidHaloNodeList[myid].keys():
                self.haloNodesPositionsInit[GlobalIndex] = (posx, posy, posz)
            else:
                fluidIndexing_temp[GlobalIndex] = self.__getGlobalIndex('fluid', myid, localIndex)
                localFluidInterface_array_X_init[localIndex] = posx
                localFluidInterface_array_Y_init[localIndex] = posy
                localFluidInterface_array_Z_init[localIndex] = posz
                self.localFluidInterface_vertex_indices[localIndex] = int(iVertex)
                localIndex += 1

        #print("rank: {}, local_vertex_indices: {}".format(myid, self.localFluidInterface_vertex_indices))

        if self.have_MPI:
            fluidIndexing_temp = self.comm.allgather(fluidIndexing_temp)
            for ii in range(len(fluidIndexing_temp)):
                for key, value in fluidIndexing_temp[ii].items():
                    self.fluidIndexing[key] = value

             # --- Collect local array sizes using gather in root process
            bufXCoor = np.array(localFluidInterface_array_X_init)
            bufYCoor = np.array(localFluidInterface_array_Y_init)
            bufZCoor = np.array(localFluidInterface_array_Z_init)
            self.sendCounts = np.array(self.comm.gather(self.nLocalFluidInterfacePhysicalNodes, 0))

            if myid == self.rootProcess:
                print("sendCounts: {}, total: {}".format(self.sendCounts, sum(self.sendCounts)))
                self.globalFluidInterfaceXcoor = np.empty(sum(self.sendCounts))
                self.globalFluidInterfaceYcoor = np.empty(sum(self.sendCounts))
                self.globalFluidInterfaceZcoor = np.empty(sum(self.sendCounts))

            self.comm.Gatherv(sendbuf=bufXCoor, recvbuf=(self.globalFluidInterfaceXcoor, self.sendCounts), root=0)
            self.comm.Gatherv(sendbuf=bufYCoor, recvbuf=(self.globalFluidInterfaceYcoor, self.sendCounts), root=0)
            self.comm.Gatherv(sendbuf=bufZCoor, recvbuf=(self.globalFluidInterfaceZcoor, self.sendCounts), root=0)
            #if myid == 0:
                #print("Gathered array X: {}".format(self.globalFluidInterfaceXcoor))
                #print("Gathered array Y: {}".format(self.globalFluidInterfaceYcoor))
                #print("Gathered array Z: {}".format(self.globalFluidInterfaceZcoor))

        else:
            self.fluidIndexing = fluidIndexing_temp.copy()
            self.globalFluidInterfaceXcoor = localFluidInterface_array_X_init.copy()
            self.globalFluidInterfaceYcoor = localFluidInterface_array_Y_init.copy()
            self.globalFluidInterfaceZcoor = localFluidInterface_array_Z_init.copy()

        self.MPIBarrier()

        # Store the global fluid coordinates
        if myid == self.rootProcess:
            self.globalFluidCoordinates0 = np.zeros((self.nFluidInterfacePhysicalNodes, 3))
            self.globalFluidCoordinates0X = np.zeros(self.nFluidInterfacePhysicalNodes)
            self.globalFluidCoordinates0Y = np.zeros(self.nFluidInterfacePhysicalNodes)
            self.globalFluidCoordinates0Z = np.zeros(self.nFluidInterfacePhysicalNodes)

            for i in range(0, self.nFluidInterfacePhysicalNodes):
                self.globalFluidCoordinates0[i][0] = self.globalFluidInterfaceXcoor[i]
                self.globalFluidCoordinates0X[i]   = self.globalFluidInterfaceXcoor[i]

                self.globalFluidCoordinates0[i][1] = self.globalFluidInterfaceYcoor[i]
                self.globalFluidCoordinates0Y[i]   = self.globalFluidInterfaceYcoor[i]

                self.globalFluidCoordinates0[i][2] = self.globalFluidInterfaceZcoor[i]
                self.globalFluidCoordinates0Z[i]   = self.globalFluidInterfaceZcoor[i]

        del fluidIndexing_temp, localFluidInterface_array_X_init, \
            localFluidInterface_array_Y_init, localFluidInterface_array_Z_init

        ###################################################################################################
        # Initialize the local load array
        # The initial displacements are set to 0 (this might be changed when we have restart capabilities)
        ###################################################################################################

        self.globalFluidLoadX = np.zeros(self.nFluidInterfacePhysicalNodes)
        self.globalFluidLoadY = np.zeros(self.nFluidInterfacePhysicalNodes)
        self.globalFluidLoadZ = np.zeros(self.nFluidInterfacePhysicalNodes)

        # --- Identify the solid interface and store the number of nodes (single core) ---#
        if SolidSolver != None:
           print('Solid solver is initialized on process {}'.format(myid))
           SolidSolver.nPoint = self.nFluidInterfacePhysicalNodes
           self.haveSolidSolver = True
           self.nSolidInterfaceNodes = SolidSolver.nPoint
           self.nSolidInterfacePhysicalNodes = SolidSolver.nPoint
           self.nLocalSolidInterfaceNodes = SolidSolver.nPoint

           SolidSolver.GlobalCoordinates0 = np.zeros((self.nSolidInterfacePhysicalNodes, 3))
           SolidSolver.GlobalCoordinates0 = self.globalFluidCoordinates0

           # Further addition to the list inside SolidSolver for future usage
           SolidSolver.SetNodesProperties()

        else:
            pass

        if self.haveSolidSolver:

            self.globalSolidLoadX = np.zeros(self.nSolidInterfaceNodes)
            self.globalSolidLoadY = np.zeros(self.nSolidInterfaceNodes)
            self.globalSolidLoadZ = np.zeros(self.nSolidInterfaceNodes)

            self.globalSolidDispX = np.zeros(self.nSolidInterfaceNodes)
            self.globalSolidDispY = np.zeros(self.nSolidInterfaceNodes)
            self.globalSolidDispZ = np.zeros(self.nSolidInterfaceNodes)

        self.MPIBarrier()


        self.solidInterface_array_DispX = np.zeros(self.nSolidInterfacePhysicalNodes)
        self.solidInterface_array_DispY = np.zeros(self.nSolidInterfacePhysicalNodes)
        self.solidInterface_array_DispZ = np.zeros(self.nSolidInterfacePhysicalNodes)

        self.MPIPrint(
            'Total number of fluid interface nodes (halo nodes included) : {}'.format(self.nFluidInterfaceNodes))
        self.MPIPrint(
            'Total number of physical fluid interface nodes : {}'.format(self.nFluidInterfacePhysicalNodes))
        self.MPIPrint(
            'Total number of beam interface nodes : {}'.format(self.nSolidInterfaceNodes))
        self.MPIPrint('Total number of fluid interface nodes : {}'.format(self.nFluidInterfacePhysicalNodes))
        self.MPIPrint('Total number of solid interface nodes : {}'.format(self.nSolidInterfacePhysicalNodes))

        self.MPIBarrier()



    def transferFluidTractions(self, FluidSolver, SolidSolver):
        """
        Transfer fluid tractions.
        Gathers the fluid tractions from the interface into the root process.
        Interpolates the tractions using the transposed matrix.
        Applies the tractions into the solid solver.
        """

        # Recover the process and the size of the parallelization
        myid, MPIsize = self.checkMPI()

        ################################################################################################################
        # --- STEP 1: Retrieve the fluid loads
        # --- Get, for the fluid interface on each partition:
        # --- The vertex indices, which are stored locally on each processor
        # --- The coordinates X, Y and Z, which are stored in a global list in processor 0
        ################################################################################################################

        # Initialize the local load array
        localFluidLoadX = np.zeros(self.nLocalFluidInterfacePhysicalNodes)
        localFluidLoadY = np.zeros(self.nLocalFluidInterfacePhysicalNodes)
        localFluidLoadZ = np.zeros(self.nLocalFluidInterfacePhysicalNodes)

        localIndex = 0
        # For the vertices that belong to the interface
        for iVertex in self.localFluidInterface_vertex_indices:
            # Compute the vertex forces on the fluid solver
            loadX, loadY, loadZ = FluidSolver.GetFlowLoad(self.fluidInterfaceIdentifier, int(iVertex))
            #print("Boundary node {}. loadX, loadY, loadZ = {} {} {}".format(iVertex, loadX,loadY,loadZ))
            # Store them in the local load array
            localFluidLoadX[localIndex] = loadX
            localFluidLoadY[localIndex] = loadY
            localFluidLoadZ[localIndex] = loadZ
            localIndex += 1

        if self.have_MPI:

            # Store the local loads in buffers in the form of numpy arrays
            bufXLoad = np.array(localFluidLoadX)
            bufYLoad = np.array(localFluidLoadY)
            bufZLoad = np.array(localFluidLoadZ)

            # Initialize the global load array
            if myid == self.rootProcess:
                print("sendCounts: {}, total: {}".format(self.sendCounts, sum(self.sendCounts)))
                self.globalFluidLoadX = np.empty(sum(self.sendCounts))
                self.globalFluidLoadY = np.empty(sum(self.sendCounts))
                self.globalFluidLoadZ = np.empty(sum(self.sendCounts))

            # Gatherv using self.sendCounts maintains the ordering of the coordinates
            self.comm.Gatherv(sendbuf=bufXLoad, recvbuf=(self.globalFluidLoadX, self.sendCounts), root=0)
            self.comm.Gatherv(sendbuf=bufYLoad, recvbuf=(self.globalFluidLoadY, self.sendCounts), root=0)
            self.comm.Gatherv(sendbuf=bufZLoad, recvbuf=(self.globalFluidLoadZ, self.sendCounts), root=0)

            # if myid == 0:
            #     print("Gathered array X: {}".format(self.globalFluidLoadX))
            #     print("Gathered array Y: {}".format(self.globalFluidLoadY))
            #     print("Gathered array Z: {}".format(self.globalFluidLoadZ))

        else:
            self.globalFluidLoadX = localFluidLoadX.copy()
            self.globalFluidLoadY = localFluidLoadY.copy()
            self.globalFluidLoadZ = localFluidLoadZ.copy()

        # Delete local variables
        del localFluidLoadX, localFluidLoadY, localFluidLoadZ

        self.MPIBarrier()

        ################################################################################################################
        # --- STEP 2: Send it to the corresponding "ghost" nodes in the structural solver
        ################################################################################################################

        if myid == self.rootProcess:
            self.globalSolidLoadX = self.globalFluidLoadX
            self.globalSolidLoadY = self.globalFluidLoadY
            self.globalSolidLoadZ = self.globalFluidLoadZ

            for iVertex in range (0, self.nSolidInterfaceNodes):
                SolidSolver.applyload(iVertex, self.globalSolidLoadX[iVertex], self.globalSolidLoadY[iVertex], self.globalSolidLoadZ[iVertex])

        #     print("Drag force: ", FluidSolver.Get_Drag())
        #     print("Lift force: ", FluidSolver.Get_Lift())
        #     print("Drag coefficient: ", FluidSolver.Get_DragCoeff())
        #     print("Lift coefficient: ", FluidSolver.Get_LiftCoeff())
        #
            outF = open("loadsFlow.txt", "w")
            index = 0
            for i in self.globalFluidLoadX:
                outF.write(str(index))
                outF.write("\t")
                outF.write(str(self.globalFluidLoadX[index]))
                outF.write("\t")
                outF.write(str(self.globalFluidLoadY[index]))
                outF.write("\t")
                outF.write(str(self.globalFluidLoadZ[index]))
                outF.write("\n")
                index += 1
            outF.close()

            outF = open("loadsFEA.txt", "w")
            index = 0
            for i in self.globalSolidLoadZ:
                outF.write(str(index))
                outF.write("\t")
                outF.write(str(self.globalSolidLoadX[index]))
                outF.write("\t")
                outF.write(str(self.globalSolidLoadY[index]))
                outF.write("\t")
                outF.write(str(self.globalSolidLoadZ[index]))
                outF.write("\n")
                index += 1
            outF.close()
        #
        # exit()

        # ---> Output: self.globalSolidLoadX, self.globalSolidLoadY, self.globalSolidLoadZ

        ################################################################################################################
        # --- STEP 3: Check conservation
        ################################################################################################################

        # --- Check for total force conservation after interpolation
        if myid == self.rootProcess:

            # Total loads before interpolation, fluid side
            FFX = self.globalFluidLoadX.sum()
            FFY = self.globalFluidLoadY.sum()
            FFZ = self.globalFluidLoadZ.sum()

            # Total loads after interpolation, solid side
            FX = self.globalSolidLoadX.sum()
            FY = self.globalSolidLoadY.sum()
            FZ = self.globalSolidLoadZ.sum()

            self.MPIPrint("Checking f/s interface total force... --- For NITRO it's pointless (to be removed) ")
            self.MPIPrint('Solid side (Fx, Fy, Fz) = ({}, {}, {})'.format(FX, FY, FZ))
            self.MPIPrint('Fluid side (Fx, Fy, Fz) = ({}, {}, {})'.format(FFX, FFY, FFZ))

            force_file = open("history_forces.dat", "a")
            force_file.write(str(FFX) + "\t" + str(FFY) + "\t" + str(FFZ) + "\n" + str(FX) + "\t" + str(FY) + "\t" + str(FZ) + "\n")
            force_file.close()

            #f = open('pyBeam_Loads_Iter' + str(self.FSIIter) + '.dat', "w+")
            #for iVertex in range(0, self.nSolidInterfaceNodes):
            #    f.write('beam.SetLoads(' + str(iVertex) +',' + str(self.globalSolidLoadX[iVertex]) +',' + str(self.globalSolidLoadY[iVertex]) +',' + str(self.globalSolidLoadZ[iVertex]) + ')\n' )
            #    print(str(iVertex) +',' + str(self.globalSolidLoadX[iVertex]) +',' + str(self.globalSolidLoadY[iVertex]) +',' + str(self.globalSolidLoadZ[iVertex]) )
            #f.close()
            #os.rename('surface_flow.vtk', 'surface_flow_' + str(self.FSIIter) + '.vtk')
    def transferStructuralDisplacements(self, FluidSolver, SolidSolver ):
        """
        Transfer structural displacements.
        Gathers the structural displacements.
        Applies the fluid displacements by scattering them into the correct position for the fluid solver.
        """

        # Recover the process and the size of the parallelization
        myid, MPIsize = self.checkMPI()

        ################################################################################################################
        # --- STEP 1: Retrieve the structural displacements from pyBeam and apply the relaxation
        # --- pyBeam runs in single core, so there is no need to deal with the parallelization
        ################################################################################################################

        if self.haveSolidSolver:

            # Recover the relaxation parameter

            # Store the old displacements
            self.globalSolidDispXOld = self.globalSolidDispX.copy()
            self.globalSolidDispYOld = self.globalSolidDispY.copy()
            self.globalSolidDispZOld = self.globalSolidDispZ.copy()


            # For the vertices that belong to the interface
            for iVertex in range(0, self.nSolidInterfaceNodes):

                # Store the new displacements in the global load array directly
                dispX, dispY, dispZ = SolidSolver.ExtractDisplacements(iVertex)
                self.globalSolidDispX[iVertex] = dispX
                self.globalSolidDispY[iVertex] = dispY
                self.globalSolidDispZ[iVertex] = dispZ

            outF = open("DispStr.txt", "w")
            for iVertex in range(0, self.nSolidInterfaceNodes):
                outF.write(str(iVertex))
                outF.write("\t")
                outF.write(str(self.globalSolidDispX[iVertex]))
                outF.write("\t")
                outF.write(str(self.globalSolidDispY[iVertex]))
                outF.write("\t")
                outF.write(str(self.globalSolidDispZ[iVertex]))
                outF.write("\n")
            outF.close()



        ################################################################################################################
        # --- STEP 2: Interpolate
        ################################################################################################################

        # ---> Input: relaxedSolidDispX, relaxedSolidDispY, relaxedSolidDispZ

        if myid == self.rootProcess:

            self.globalFluidDispX = self.globalSolidDispX
            self.globalFluidDispY = self.globalSolidDispY
            self.globalFluidDispZ = self.globalSolidDispZ

        # ---> Output: self.globalFluidDispX, self.globalFluidDispY, self.globalFluidDispZ

        ################################################################################################################
        # --- STEP 3: Check conservation
        ################################################################################################################

        # --- Checking conservation ---
        '''
        if myid == self.rootProcess:

            WSX = self.globalSolidLoadX.dot(self.globalSolidDispX)
            WSY = self.globalSolidLoadY.dot(self.globalSolidDispY)
            WSZ = self.globalSolidLoadZ.dot(self.globalSolidDispZ)

            WFX = self.globalFluidLoadX.dot(self.globalFluidDispX)
            WFY = self.globalFluidLoadY.dot(self.globalFluidDispY)
            WFZ = self.globalFluidLoadZ.dot(self.globalFluidDispZ)

            self.MPIPrint("Checking f/s interface conservation...")
            self.MPIPrint('Solid side (Wx, Wy, Wz) = ({}, {}, {})'.format(WSX, WSY, WSZ))
            self.MPIPrint('Fluid side (Wx, Wy, Wz) = ({}, {}, {})'.format(WFX, WFY, WFZ))
        '''
        ################################################################################################################
        # --- STEP 4: Transfer to the fluid solver
        ################################################################################################################

        # --- Recover them from the interpolated vectors
        localFluidDispX = np.zeros(self.nLocalFluidInterfacePhysicalNodes)
        localFluidDispY = np.zeros(self.nLocalFluidInterfacePhysicalNodes)
        localFluidDispZ = np.zeros(self.nLocalFluidInterfacePhysicalNodes)

        localFluidCoord0X = np.zeros(self.nLocalFluidInterfacePhysicalNodes)
        localFluidCoord0Y = np.zeros(self.nLocalFluidInterfacePhysicalNodes)
        localFluidCoord0Z = np.zeros(self.nLocalFluidInterfacePhysicalNodes)

        if self.have_MPI:
            # here it communicates with the root (apparently)
            self.comm.Scatterv(sendbuf=(self.globalFluidDispX, self.sendCounts), recvbuf=localFluidDispX, root=0)
            self.comm.Scatterv(sendbuf=(self.globalFluidDispY, self.sendCounts), recvbuf=localFluidDispY, root=0)
            self.comm.Scatterv(sendbuf=(self.globalFluidDispZ, self.sendCounts), recvbuf=localFluidDispZ, root=0)

            self.comm.Scatterv(sendbuf=(self.globalFluidCoordinates0X, self.sendCounts), recvbuf=localFluidCoord0X, root=0)
            self.comm.Scatterv(sendbuf=(self.globalFluidCoordinates0Y, self.sendCounts), recvbuf=localFluidCoord0Y, root=0)
            self.comm.Scatterv(sendbuf=(self.globalFluidCoordinates0Z, self.sendCounts), recvbuf=localFluidCoord0Z, root=0)

            # print("rank: {}, local_array X: {}".format(myid, localFluidDispX))
            # print("rank: {}, local_array Y: {}".format(myid, localFluidDispY))
            # print("rank: {}, local_array Z: {}".format(myid, localFluidDispZ))

        else:
            localFluidDispX = self.globalFluidDispX.copy()
            localFluidDispY = self.globalFluidDispY.copy()
            localFluidDispZ = self.globalFluidDispZ.copy()

            localFluidCoord0X = self.globalFluidCoordinates0X.copy()
            localFluidCoord0Y = self.globalFluidCoordinates0Y.copy()
            localFluidCoord0Z = self.globalFluidCoordinates0Z.copy()

        # For the vertices that belong to the interface
        localIndex = 0
        for iVertex in self.localFluidInterface_vertex_indices:
            # Store them in the mesh displacement routine
            FluidSolver.SetMeshDisplacement(self.fluidInterfaceIdentifier, int(iVertex), localFluidDispX[localIndex],
                                            localFluidDispY[localIndex], localFluidDispZ[localIndex])

            #FluidSolver.SetVertexCoordX(self.fluidInterfaceIdentifier, int(iVertex), localFluidCoord0X[localIndex] + localFluidDispX[localIndex])
            #FluidSolver.SetVertexCoordY(self.fluidInterfaceIdentifier, int(iVertex), localFluidCoord0Y[localIndex] + localFluidDispY[localIndex])
            #FluidSolver.SetVertexCoordZ(self.fluidInterfaceIdentifier, int(iVertex), localFluidCoord0Z[localIndex] + localFluidDispZ[localIndex])
            #print("Processor {}. Local index {} for node {} and displacements: {}, {}, {}.".format(myid,localIndex,iVertex, localFluidDispX[localIndex], localFluidDispY[localIndex], localFluidDispZ[localIndex] ))
            # Increment the local index
            localIndex += 1
            #nodalVarCoordNorm = FluidSolver.SetVertexVarCoord(self.fluidInterfaceIdentifier, int(iVertex))
            #print("nodalVarCoordNorm = {}".format(nodalVarCoordNorm) )

        # Delete local variables
        del localFluidDispX, localFluidDispY, localFluidDispZ, localFluidCoord0X, localFluidCoord0Y, localFluidCoord0Z

    def getSolidInterfaceDisplacement(self, SolidSolver):
        """
        Gets the current solid interface position from the solid solver.
        """
        print("getSolidInterfaceDisplacement is replaced by transferStructuralDisplacements")
        print("self.solidInterface_array_DispX/Y/Z are replaced by self.globalFluidDispX/Y/Z")
        '''
        if self.have_MPI == True:
          myid = self.comm.Get_rank()
        else:
          myid = 0

        # --- Get the solid interface position from the solid solver and directly fill the corresponding PETSc vector ---
        GlobalIndex = int()
        localIndex = 0
        #for iVertex in self.localFluidInterface_vertex_indices:
        for iVertex in range(self.nLocalSolidInterfaceNodes):
          GlobalIndex = SolidSolver.getInterfaceNodeGlobalIndex(self.solidInterfaceIdentifier, iVertex)
          if GlobalIndex in self.SolidHaloNodeList[myid].keys():
            pass
          else:
            newDispx = SolidSolver.getInterfaceNodeDispX(self.solidInterfaceIdentifier, iVertex)
            newDispy = SolidSolver.getInterfaceNodeDispY(self.solidInterfaceIdentifier, iVertex)
            newDispz = SolidSolver.getInterfaceNodeDispZ(self.solidInterfaceIdentifier, iVertex)
            iGlobalVertex = self.__getGlobalIndex('solid', myid, localIndex)
            self.solidInterface_array_DispX[iGlobalVertex] = newDispx
            self.solidInterface_array_DispY[iGlobalVertex] = newDispy
            self.solidInterface_array_DispZ[iGlobalVertex] = newDispz
            localIndex += 1
        '''
    def setFluidInterfaceVarCoord(self, FluidSolver):
        """
        Communicate the change of coordinates of the fluid interface to the fluid solver.
        Prepare the fluid solver for mesh deformation.
        """
        if self.have_MPI == True:
          myid = self.comm.Get_rank()
        else:
          myid = 0

        # --- Send the new fluid interface position to the fluid solver (on each partition, halo nodes included) ---
        localIndex = 0
        for iVertex in range(self.nLocalFluidInterfaceNodes):
            GlobalIndex = FluidSolver.GetVertexGlobalIndex(self.fluidInterfaceIdentifier, iVertex)
            if GlobalIndex in self.FluidHaloNodeList[myid].keys():
              posX0, posY0, posZ0 = self.haloNodesPositionsInit[GlobalIndex]
              DispX, DispY, DispZ = self.haloNodesDisplacements[GlobalIndex]
              #if posY0 == 0.0:
              #  posX = posX0
              #  posY = posY0
              #  posZ = posZ0
              #else:
              posX = posX0 + DispX
              posY = posY0 + DispY
              posZ = posZ0 + DispZ
              FluidSolver.SetVertexCoordX(self.fluidInterfaceIdentifier, iVertex, posX)
              FluidSolver.SetVertexCoordY(self.fluidInterfaceIdentifier, iVertex, posY)
              FluidSolver.SetVertexCoordZ(self.fluidInterfaceIdentifier, iVertex, posZ)
            else:
              #if self.localFluidInterface_array_Y_init[localIndex] == 0.0:				# !!! This is temporary and case dependent, it should be removed ASAP !!!
                #posX = self.localFluidInterface_array_X_init[localIndex]
                #posY = self.localFluidInterface_array_Y_init[localIndex]
                #posZ = self.localFluidInterface_array_Z_init[localIndex]
              #else:
              posX = self.localFluidInterface_array_DispX[localIndex] + self.localFluidInterface_array_X_init[localIndex]
              posY = self.localFluidInterface_array_DispY[localIndex] + self.localFluidInterface_array_Y_init[localIndex]
              posZ = self.localFluidInterface_array_DispZ[localIndex] + self.localFluidInterface_array_Z_init[localIndex]
              FluidSolver.SetVertexCoordX(self.fluidInterfaceIdentifier, iVertex, posX)
              FluidSolver.SetVertexCoordY(self.fluidInterfaceIdentifier, iVertex, posY)
              FluidSolver.SetVertexCoordZ(self.fluidInterfaceIdentifier, iVertex, posZ)
              localIndex += 1
            # Prepares the mesh deformation in the fluid solver
            nodalVarCoordNorm = FluidSolver.SetVertexVarCoord(self.fluidInterfaceIdentifier, iVertex)
            #print nodalVarCoordNorm


    def SteadyFSI(self, FSI_config, FluidSolver, SolidSolver, MLSSolver):
        """
        Runs the steady FSI computation
        Synchronizes the fluid and solid solver with data exchange at the f/s interface.
        """

        # Recover the process and the size of the parallelization
        myid, MPIsize = self.checkMPI()

        # --- Set some general variables for the steady computation --- #
        NbIter = FSI_config['NB_EXT_ITER']  # number of fluid iteration at each FSI step


        if myid == self.rootProcess:
            cd_file = open("history_CD.dat", "w")
            cd_file.write("Drag Coefficient\n")
            cd_file.close()
            cl_file = open("history_CL.dat", "w")
            cl_file.write("Lift Coefficient\n")
            cl_file.close()
            force_file = open("history_forces.dat", "w")
            force_file.write("Forces Flow (X, Y, Z) \t Forces FEA (X, Y, Z)\n")
            force_file.close()

        # --- Initialize matrix of boundary nodal forces  --- #
        #if myid == self.rootProcess:
        #    SolidSolver.EvaluateIntefaceFluidDisplacements(FSI_config, MLSSolver) # Flutter_mode_fluid_x/y/z are stored (root) once and for all
        #    SolidSolver.initialize_OutputForces(1,FSI_config)


        self.MPIPrint('\n********************************')
        self.MPIPrint('* Begin steady FSI computation *')
        self.MPIPrint('********************************\n')

        self.MPIBarrier()
        self.MPIPrint('\nLaunching fluid solver for a steady computation...')
        # --- Fluid solver call for FSI subiteration ---#




        FluidSolver.ResetConvergence()  # Make sure the solver starts convergence from 0
        FluidSolver.Preprocess(0)  # Time iteration pre-processing
        FluidSolver.Run()  # Run one time-step (static: one simulation)
        FluidSolver.Postprocess()
        FluidSolver.Update()  # Update the solver for the next time iteration
        FluidSolver.Monitor(0)  # Monitor the solver and output solution to file if required
        FluidSolver.Output(0)  # Output the solution to file

        # uncomment
        # --- Surface fluid loads interpolation and communication ---#
        self.MPIPrint('\nProcessing interface fluid loads...\n')
        self.MPIBarrier()
        self.transferFluidTractions(FluidSolver, SolidSolver)

        # --- Solid solver call for FSI subiteration --- #
        self.MPIPrint('\nLaunching solid solver for a static computation and Generalized force evaluation...\n')
        if myid == self.rootProcess:
                #SolidSolver.printForceDispl(0)
                SolidSolver.printNode(0, 0)
                SolidSolver.writeSolution(0, 0, FSI_config)

        # Move the restart file to a solution file
        if myid == self.rootProcess:
                new_name_flow = "./Output/flow_00000"  + ".vtk"
                new_name_surf = "./Output/surface_flow_00000" + ".vtk"
                shutil.move("flow.vtk", new_name_flow)
                shutil.move("surface_flow.vtk", new_name_surf)

                cd_file = open("history_CD.dat", "a")
                cd_file.write(str(FluidSolver.Get_DragCoeff()) + "\n")
                cd_file.close()
                cl_file = open("history_CL.dat", "a")
                cl_file.write(str(FluidSolver.Get_LiftCoeff()) + "\n")
                cl_file.close()

        self.printMeshCoord_bis(FluidSolver, 0)
        self.MPIBarrier()
        self.MPIPrint(' ')
        self.MPIPrint('*************************')
        self.MPIPrint('*  End FSI computation  *')
        self.MPIPrint('*************************')
        self.MPIPrint(' ')

    def UnsteadyFSI(self, FSI_config, FluidSolver, SolidSolver, MLS_Spline):
        """
	Run the unsteady FSI computation by synchronizing the fluid and solid solvers.
	F/s interface data are exchanged through interface mapping and interpolation (if non mathcing meshes).
	"""

        if self.have_MPI == True:
            myid = self.comm.Get_rank()
            numberPart = self.comm.Get_size()
        else:
            myid = 0
            numberPart = 1

        # write some data
        if myid == self.rootProcess:
            cd_file = open("history_CD.dat", "w")
            cd_file.write("Drag Coefficient\n")
            cd_file.close()
            cl_file = open("history_CL.dat", "w")
            cl_file.write("Lift Coefficient\n")
            cl_file.close()
            force_file = open("history_forces.dat", "w")
            force_file.write("Forces Flow (X, Y, Z) \t Forces FEA (X, Y, Z)\n")
            force_file.close()


        # --- Set some general variables for the unsteady computation --- #
        deltaT = FSI_config['UNST_TIMESTEP']  # physical time step
        totTime = FSI_config['UNST_TIME']  # physical simulation time
        TimeIterTreshold = 0  # time iteration from which we allow the solid to deform
        if FSI_config['MOTION_TYPE'] == 'BLENDED_STEP':
            blended_step_lenght = FSI_config['BLE_STEP_LENGTH']  # time required to perform the blended step

        if FSI_config['RESTART_SOL'] == 'YES':
            startTime = FSI_config['START_TIME']
            NbTimeIter = ((totTime) / deltaT) + 1
            time = startTime
            TimeIter = FSI_config['RESTART_ITER']
        else:
            NbTimeIter = (totTime / deltaT) + 1  # number of time iterations
            time = 0.0  # initial time
            TimeIter = 0  # initial time iteration

        #NbTimeIter = 15 ; print("NbTimeIter = 15 forced")

        self.MPIPrint('\n**********************************')
        self.MPIPrint('* Begin unsteady FSI computation *')
        self.MPIPrint('**********************************\n')

        # --- Initialize matrix of boundary nodal forces  --- #
        #if myid == self.rootProcess:
        #    SolidSolver.initialize_OutputForces(NbTimeIter, FSI_config)  ## NEEDS TO BE REBUILD IN CASE OF RESTART!!!!!
        # --- Initialize the coupled solution --- #
        if FSI_config['RESTART_SOL'] == 'YES':
            TimeIterTreshold = -1
            if myid == self.rootProcess:  # HARD CODE
                SolidSolver.run(FSI_config['UNST_TIMESTEP'] * (FSI_config['RESTART_ITER']), FSI_config, MLS_Spline)
            if self.have_MPI == True:
                self.comm.barrier()
            #self.transferStructuralDisplacements( FluidSolver, SolidSolver)
            #self.getSolidInterfaceDisplacement(SolidSolver)
        # If no restart
        else:
            self.MPIPrint('Setting FSI initial conditions')
            if myid in self.solidSolverProcessors:
                #SolidSolver.EvaluateIntefaceFluidDisplacements(FSI_config,MLS_Spline)  # Flutter_mode_fluid_x/y/z are stored (root) once and for all
                SolidSolver.setInitialDisplacements(FSI_config, MLS_Spline)
            if self.have_MPI == True:
                self.comm.barrier()
            self.transferStructuralDisplacements(FluidSolver, SolidSolver)
            #self.interpolateSolidPositionOnFluidMesh(FSI_config) #OLD VERSION
            #self.setFluidInterfaceVarCoord(FluidSolver)
            FluidSolver.Preprocess(0)  # if there is an initial deformation in the solid, it has to be communicated to the fluid solver
            self.MPIPrint('\nFSI initial conditions are set')
            self.MPIPrint('Beginning time integration\n')


        # --- External temporal loop --- #
        while TimeIter <= NbTimeIter:

            self.FSIIter = 0
            FSIConv = False
            self.MPIPrint("Timeiter = {}".format(TimeIter))

            #self.MPIPrint("\n>>>> Time iteration {} / FSI iteration {} <<<<".format(TimeIter, self.FSIIter))

            if TimeIter != 0:
                if FSI_config['MOTION_TYPE'] == 'BLENDED_STEP':
                    if time <= FSI_config['START_MOTION_TIME'] + blended_step_lenght and time >= FSI_config[
                        'START_MOTION_TIME']:
                        # --- Mesh morphing step (displacements interpolation, displacements communication, and mesh morpher call) --- #
                        self.transferStructuralDisplacements( FluidSolver, SolidSolver)
                        #self.interpolateSolidPositionOnFluidMesh(FSI_config) #OLD VERSION
                        self.MPIPrint('\nPerforming dynamic mesh deformation (ALE)...\n')
                        #self.setFluidInterfaceVarCoord(FluidSolver)
                        #FluidSolver.DynamicMeshUpdate(TimeIter)
                        FluidSolver.Preprocess(TimeIter)
                else:
                    # --- Mesh morphing step (displacements interpolation, displacements communication, and mesh morpher call) --- #
                    self.transferStructuralDisplacements(FluidSolver, SolidSolver)
                    #self.interpolateSolidPositionOnFluidMesh(FSI_config) #OLD VERSION
                    self.MPIPrint('\nPerforming dynamic mesh deformation (ALE)...\n')
                    #self.setFluidInterfaceVarCoord(FluidSolver)
                    #FluidSolver.DynamicMeshUpdate(TimeIter)
                    FluidSolver.Preprocess(TimeIter)

            self.printMeshCoord_bis(FluidSolver, TimeIter)
            # --- Fluid solver call for FSI subiteration --- #
            self.MPIPrint('\nLaunching fluid solver for one single dual-time iteration...')
            self.MPIPrint("Time Iter = {}".format(FluidSolver.GetTime_Iter()))
            self.MPIBarrier()

            FluidSolver.ResetConvergence()
            FluidSolver.Run()
            FluidSolver.Postprocess()
            # --- Update, monitor and output the fluid solution before the next time step  ---#
            FluidSolver.Update()
            FluidSolver.Monitor(TimeIter)
            FluidSolver.Output(TimeIter)
            self.MPIBarrier()

            # --- Surface fluid loads interpolation and communication --- #
            self.MPIPrint('\nProcessing interface fluid loads...\n')
            self.MPIBarrier()
            self.transferFluidTractions(FluidSolver, SolidSolver)
            self.MPIBarrier()


            if myid == self.rootProcess:
                    # --- Output the solid solution before thr next time step --- #
                    #SolidSolver.printForceDispl(TimeIter)
                    SolidSolver.printNode(TimeIter, 1)
                    SolidSolver.writeSolution(TimeIter, time, FSI_config)

            if myid == self.rootProcess and FSI_config['REAL_TIME_TRACKING'] == 'YES':
                monitor(FSI_config, SolidSolver)

            # Move the restart file to a solution file
            if myid == self.rootProcess:


                new_name_flow = "./Output/flow_"  + str(TimeIter).zfill(5) + ".vtk"
                new_name_surf = "./Output/surface_flow_" + str(TimeIter).zfill(5) + ".vtk"
                #shutil.move("flow_" + str(0).zfill(5) + ".vtk", new_name_flow)
                #shutil.move("surface_flow_" + str(0).zfill(5) + ".vtk", new_name_surf)
                os.remove("surface_flow_" + str(TimeIter).zfill(5) + ".vtk")
                os.remove("flow_" + str(TimeIter).zfill(5) + ".vtk")
                cd_file = open("history_CD.dat", "a")
                cd_file.write(str(FluidSolver.Get_DragCoeff()) + "\n")
                cd_file.close()
                cl_file = open("history_CL.dat", "a")
                cl_file.write(str(FluidSolver.Get_LiftCoeff()) + "\n")
                cl_file.close()

            TimeIter += 1
            time += deltaT

            ## --- Solid solver call for FSI subiteration --- #
            self.MPIPrint('\nEvaluating the solid displacement for the considerd time-step')
            if myid == self.rootProcess:
                SolidSolver.run(time, FSI_config, MLS_Spline)
            #self.transferStructuralDisplacements(FluidSolver, SolidSolver)
            #self.getSolidInterfaceDisplacement(
            #    SolidSolver)  # this has to be done for every processor (not only the one of the SolidSolver!!)

            # --- End of the temporal loop --- #

            self.MPIBarrier()


        self.MPIPrint('\n*************************')
        self.MPIPrint('*  End FSI computation  *')
        self.MPIPrint('*************************\n')

    def __getGlobalIndex(self, physics, iProc, iLocalVertex):
        """
        Calculate the global indexing of interface nodes accross all the partitions. This does not include halo nodes.
        """

        if physics == 'fluid':
            globalStartIndex = self.fluidGlobalIndexRange[iProc][iProc][0]
        elif physics == 'solid':
            globalStartIndex = self.solidGlobalIndexRange[iProc][iProc][0]

        globalIndex = globalStartIndex + iLocalVertex

        return globalIndex


    def printMeshCoord_bis(self, FluidSolver, TimeIter):

        if self.have_MPI == True:
            myid = self.comm.Get_rank()
            numberPart = self.comm.Get_size()
        else:
            myid = 0
            numberPart = 1

        outC = open("./Output/FluidBoundCoord_" + str(TimeIter) +"Proc_" +str(myid) + ".dat", "w")

        self.MPIBarrier()

        print("DEBUG")

        # For the vertices that belong to the interface
        for iVertex in self.localFluidInterface_vertex_indices:
                outC.write(str(int(iVertex)))
                outC.write("\t")
                outC.write(str(float(FluidSolver.GetVertexCoordX(self.fluidInterfaceIdentifier, int(iVertex)))))
                outC.write("\t")
                outC.write(str(float(FluidSolver.GetVertexCoordY(self.fluidInterfaceIdentifier, int(iVertex)))))
                outC.write("\t")
                outC.write(str(float(FluidSolver.GetVertexCoordZ(self.fluidInterfaceIdentifier, int(iVertex)))))
                outC.write("\n")



        self.MPIBarrier()
        outC.close()


        if myid == self.rootProcess:
          filenames = ['Hi' for i in range(numberPart)]
          for i in range(0,numberPart):
               filenames[i] = "./Output/FluidBoundCoord_" + str(TimeIter) +"Proc_" +str(i) + ".dat"


          with open( "./Output/FluidBoundCoord_" + str(TimeIter) + ".dat", 'w') as outfile:
            for fname in filenames:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)

                os.remove(fname)
        self.MPIBarrier()


    def SteadyFSI_dyn(self, FSI_config, FluidSolver, SolidSolver, MLSSolver):
        """
        Runs the steady FSI computation from dynresp
        Synchronizes the fluid and solid solver with data exchange at the f/s interface.
        """

        # first it is necessary to read the modal displacement file



        # Recover the process and the size of the parallelization
        myid, MPIsize = self.checkMPI()

        # --- Set some general variables for the steady computation --- #
        NbIter = FSI_config['NB_EXT_ITER']  # number of fluid iteration at each FSI step


        self.MPIPrint('\n********************************')
        self.MPIPrint('* Begin steady FSI computation *')
        self.MPIPrint('********************************\n')

        self.MPIBarrier()
        self.MPIPrint('\nLaunching fluid solver for a steady computation...')
        # --- Fluid solver call for FSI subiteration ---#


        if myid == self.rootProcess:
            # first it is necessary to read the modal displacement file
            SolidSolver.mod_displ = ReadModalDisplacements(FSIConfig)
            # all modes are evaluated on the fluid boundary: PHI_CFD
            #SolidSolver.EvaluateIntefaceFluidDisplacements(FSI_config, MLSSolver) # Flutter_mode_fluid_x/y/z are stored (root) once and for all
            SolidSolver.run(0, FSI_config, MLS_Spline,0)


        FluidSolver.ResetConvergence()  # Make sure the solver starts convergence from 0
        FluidSolver.Preprocess(0)  # Time iteration pre-processing
        FluidSolver.Run()  # Run one time-step (static: one simulation)
        FluidSolver.Postprocess()
        FluidSolver.Update()  # Update the solver for the next time iteration
        FluidSolver.Monitor(0)  # Monitor the solver and output solution to file if required
        FluidSolver.Output(0)  # Output the solution to file

        # uncomment
        # --- Surface fluid loads interpolation and communication ---#
        self.MPIPrint('\nProcessing interface fluid loads...\n')
        self.MPIBarrier()
        self.transferFluidTractions(FluidSolver, SolidSolver)

        # --- Solid solver call for FSI subiteration --- #
        self.MPIPrint('\nLaunching solid solver for a static computation and Generalized force evaluation...\n')
        if myid == self.rootProcess:
                #SolidSolver.printForceDispl(0)
                SolidSolver.printNode(0, 0)
                SolidSolver.writeSolution(0, 0, FSI_config)

        # Move the restart file to a solution file
        if myid == self.rootProcess:
                new_name_flow = "./Output/flow_00000"  + ".vtk"
                new_name_surf = "./Output/surface_flow_00000" + ".vtk"
                shutil.move("flow.vtk", new_name_flow)
                shutil.move("surface_flow.vtk", new_name_surf)

                cd_file = open("history_CD.dat", "a")
                cd_file.write(str(FluidSolver.Get_DragCoeff()) + "\n")
                cd_file.close()
                cl_file = open("history_CL.dat", "a")
                cl_file.write(str(FluidSolver.Get_LiftCoeff()) + "\n")
                cl_file.close()

        self.printMeshCoord_bis(FluidSolver, 0)
        self.MPIBarrier()
        self.MPIPrint(' ')
        self.MPIPrint('*************************')
        self.MPIPrint('*  End FSI computation  *')
        self.MPIPrint('*************************')
        self.MPIPrint(' ')

    def UnsteadyFSI_dyn_sequential(self, FSI_config, FluidSolver, SolidSolver, MLS_Spline):
        """
	Run the unsteady FSI computation by synchronizing the fluid and solid solvers.
	F/s interface data are exchanged through interface mapping and interpolation (if non mathcing meshes).
	"""
        print('Not Ready yet!')
        if self.have_MPI == True:
            myid = self.comm.Get_rank()
            numberPart = self.comm.Get_size()
        else:
            myid = 0
            numberPart = 1

        # write some data
        if myid == self.rootProcess:
            cd_file = open("history_CD.dat", "w")
            cd_file.write("Drag Coefficient\n")
            cd_file.close()
            cl_file = open("history_CL.dat", "w")
            cl_file.write("Lift Coefficient\n")
            cl_file.close()
            force_file = open("history_forces.dat", "w")
            force_file.write("Forces Flow (X, Y, Z) \t Forces FEA (X, Y, Z)\n")
            force_file.close()


        # --- Set some general variables for the unsteady computation --- #
        deltaT = FSI_config['UNST_TIMESTEP']  # physical time step
        totTime = FSI_config['UNST_TIME']  # physical simulation time
        TimeIterTreshold = 0  # time iteration from which we allow the solid to deform
        if FSI_config['MOTION_TYPE'] == 'BLENDED_STEP':
            blended_step_lenght = FSI_config['BLE_STEP_LENGTH']  # time required to perform the blended step

        if FSI_config['RESTART_SOL'] == 'YES':
            startTime = FSI_config['START_TIME']
            NbTimeIter = ((totTime) / deltaT) + 1
            time = startTime
            TimeIter = FSI_config['RESTART_ITER']
        else:
            NbTimeIter = (totTime / deltaT) + 1  # number of time iterations
            time = 0.0  # initial time
            TimeIter = 0  # initial time iteration

        #NbTimeIter = 15 ; print("NbTimeIter = 15 forced")

        self.MPIPrint('\n**********************************')
        self.MPIPrint('* Begin unsteady FSI computation *')
        self.MPIPrint('**********************************\n')

        # Read modal displacement file (for sequential approach has to be done only once (beginning) as it is not changing every timestep)
        if myid == self.rootProcess:
           SolidSolver.mod_displ = ReadModalDisplacements(FSIConfig)

        if self.have_MPI == True:
            self.comm.barrier()

        # --- Initialize the coupled solution --- #
        if FSI_config['RESTART_SOL'] == 'YES':
            TimeIterTreshold = -1
            if myid == self.rootProcess:  # HARD CODE
                SolidSolver.run(FSI_config['UNST_TIMESTEP'] * (FSI_config['RESTART_ITER']), FSI_config, MLS_Spline, FSI_config['RESTART_ITER'])
            if self.have_MPI == True:
                self.comm.barrier()
            #self.transferStructuralDisplacements( FluidSolver, SolidSolver)
            #self.getSolidInterfaceDisplacement(SolidSolver)
        # If no restart
        else:
            self.MPIPrint('Setting FSI initial conditions')
            if myid in self.solidSolverProcessors:
                SolidSolver.setInitialDisplacements(FSI_config, MLS_Spline)
            if self.have_MPI == True:
                self.comm.barrier()
            self.transferStructuralDisplacements(FluidSolver, SolidSolver)
            #self.interpolateSolidPositionOnFluidMesh(FSI_config) #OLD VERSION
            #self.setFluidInterfaceVarCoord(FluidSolver)
            FluidSolver.Preprocess(0)  # if there is an initial deformation in the solid, it has to be communicated to the fluid solver
            self.MPIPrint('\nFSI initial conditions are set')
            self.MPIPrint('Beginning time integration\n')


        # --- External temporal loop --- #
        while TimeIter <= NbTimeIter:

            self.FSIIter = 0
            FSIConv = False
            self.MPIPrint("Timeiter = {}".format(TimeIter))

            #self.MPIPrint("\n>>>> Time iteration {} / FSI iteration {} <<<<".format(TimeIter, self.FSIIter))

            if TimeIter != 0:
                if FSI_config['MOTION_TYPE'] == 'BLENDED_STEP':
                    if time <= FSI_config['START_MOTION_TIME'] + blended_step_lenght and time >= FSI_config[
                        'START_MOTION_TIME']:
                        # --- Mesh morphing step (displacements interpolation, displacements communication, and mesh morpher call) --- #
                        self.transferStructuralDisplacements( FluidSolver, SolidSolver)
                        #self.interpolateSolidPositionOnFluidMesh(FSI_config) #OLD VERSION
                        self.MPIPrint('\nPerforming dynamic mesh deformation (ALE)...\n')
                        #self.setFluidInterfaceVarCoord(FluidSolver)
                        #FluidSolver.DynamicMeshUpdate(TimeIter)
                        FluidSolver.Preprocess(TimeIter)
                else:
                    # --- Mesh morphing step (displacements interpolation, displacements communication, and mesh morpher call) --- #
                    self.transferStructuralDisplacements(FluidSolver, SolidSolver)
                    #self.interpolateSolidPositionOnFluidMesh(FSI_config) #OLD VERSION
                    self.MPIPrint('\nPerforming dynamic mesh deformation (ALE)...\n')
                    #self.setFluidInterfaceVarCoord(FluidSolver)
                    #FluidSolver.DynamicMeshUpdate(TimeIter)
                    FluidSolver.Preprocess(TimeIter)

            self.printMeshCoord_bis(FluidSolver, TimeIter)
            # --- Fluid solver call for FSI subiteration --- #
            self.MPIPrint('\nLaunching fluid solver for one single dual-time iteration...')
            self.MPIPrint("Time Iter = {}".format(FluidSolver.GetTime_Iter()))
            self.MPIBarrier()

            FluidSolver.ResetConvergence()
            FluidSolver.Run()
            FluidSolver.Postprocess()
            # --- Update, monitor and output the fluid solution before the next time step  ---#
            FluidSolver.Update()
            FluidSolver.Monitor(TimeIter)
            FluidSolver.Output(TimeIter)
            self.MPIBarrier()

            # --- Surface fluid loads interpolation and communication --- #
            self.MPIPrint('\nProcessing interface fluid loads...\n')
            self.MPIBarrier()
            self.transferFluidTractions(FluidSolver, SolidSolver)
            self.MPIBarrier()


            if myid == self.rootProcess:
                    # --- Output the solid solution before thr next time step --- #
                    #SolidSolver.printForceDispl(TimeIter)
                    SolidSolver.printNode(TimeIter, 1)
                    SolidSolver.writeSolution(TimeIter, time, FSI_config)

            if myid == self.rootProcess and FSI_config['REAL_TIME_TRACKING'] == 'YES':
                monitor(FSI_config, SolidSolver)

            # Move the restart file to a solution file
            if myid == self.rootProcess:


                new_name_flow = "./Output/flow_"  + str(TimeIter).zfill(5) + ".vtk"
                new_name_surf = "./Output/surface_flow_" + str(TimeIter).zfill(5) + ".vtk"
                #shutil.move("flow_" + str(0).zfill(5) + ".vtk", new_name_flow)
                #shutil.move("surface_flow_" + str(0).zfill(5) + ".vtk", new_name_surf)
                os.remove("surface_flow_" + str(TimeIter).zfill(5) + ".vtk")
                os.remove("flow_" + str(TimeIter).zfill(5) + ".vtk")
                cd_file = open("history_CD.dat", "a")
                cd_file.write(str(FluidSolver.Get_DragCoeff()) + "\n")
                cd_file.close()
                cl_file = open("history_CL.dat", "a")
                cl_file.write(str(FluidSolver.Get_LiftCoeff()) + "\n")
                cl_file.close()

            TimeIter += 1
            time += deltaT

            ## --- Solid solver call for FSI subiteration --- #
            self.MPIPrint('\nEvaluating the solid displacement for the considerd time-step')
            if myid == self.rootProcess:
                SolidSolver.run(time, FSI_config, MLS_Spline,TimeIter)
            #self.transferStructuralDisplacements(FluidSolver, SolidSolver)
            #self.getSolidInterfaceDisplacement(
            #    SolidSolver)  # this has to be done for every processor (not only the one of the SolidSolver!!)

            # --- End of the temporal loop --- #

            self.MPIBarrier()


        self.MPIPrint('\n*************************')
        self.MPIPrint('*  End FSI computation  *')
        self.MPIPrint('*************************\n')
