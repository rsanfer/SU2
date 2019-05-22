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

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os, sys, shutil, copy
import numpy as np
import scipy as sp
import scipy.spatial.distance as spdist
from math import *
import shelve
import random

# ----------------------------------------------------------------------
#  FSI Interface Class
# ----------------------------------------------------------------------

class Interface:
    """
    FSI interface class that handles fluid/solid solvers synchronisation and communication
    """

    def __init__(self, FSI_config, FluidSolver, SolidSolver, MLS_Solver, have_MPI):
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

        self.rootProcess = 0  # the root process is chosen to be MPI rank = 0

        self.nDim = FSI_config['NDIM']  # problem dimension

        self.MPI_trhld = 10000  # Appearently I have to divide the size of the data for Send/Rcv

        self.haveFluidSolver = False  # True if the fluid solver is initialized on the current rank
        self.haveSolidSolver = False  # True if the solid solver is initialized on the current rank
        self.haveFluidInterface = False  # True if the current rank owns at least one fluid interface node
        self.haveSolidInterface = False  # True if the current rank owns at least one solid interface node

        self.fluidSolverProcessors = list()  # list of partitions where the fluid solver is initialized
        self.solidSolverProcessors = list()  # list of partitions where the solid solver is initialized
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

        self.globalFluidCoordinates = None
        self.globalSolidCoordinates = None

        self.sendCounts = None
        self.globalFluidDispX = None
        self.globalFluidDispY = None
        self.globalFluidDispZ = None

        self.globalSolidDispX = None
        self.globalSolidDispY = None
        self.globalSolidDispZ = None

        self.globalFluidLoadX = None
        self.globalFluidLoadY = None
        self.globalFluidLoadZ = None

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
        self.MPIPrint('Solid solver : pyBeam')
        self.MPIPrint('Steady coupled simulation')
        self.MPIPrint('Matching fluid-solid interface using Moving Least Squares method')
        self.MPIPrint('Maximum number of FSI iterations : {}'.format(FSI_config['NB_FSI_ITER']))
        self.MPIPrint('FSI tolerance : {}'.format(FSI_config['FSI_TOLERANCE']))
        self.MPIPrint('Static under-relaxation with constant parameter {}'.format(FSI_config['RELAX_PARAM']))
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

    def connect(self, FSI_config, FluidSolver, SolidSolver):
        """
        Connection between solvers.
        Creates the communication support between the two solvers.
        Gets information about f/s interfaces from the two solvers.
        """
        if self.have_MPI:
            myid = self.comm.Get_rank()
            MPIsize = self.comm.Get_size()
        else:
            myid = 0
            MPIsize = 1

        # --- Identify the fluid interface and store the number of nodes for each partition ---#
        self.fluidInterfaceIdentifier = None
        self.nLocalFluidInterfaceNodes = 0
        if FluidSolver != None:
            print('Fluid solver is initialized on process {}'.format(myid))
            self.haveFluidSolver = True
            allMovingMarkersTags = FluidSolver.GetAllMovingMarkersTag()
            allMarkersID = FluidSolver.GetAllBoundaryMarkers()
            if not allMovingMarkersTags:
                raise Exception('No moving marker was defined in SU2.')
            else:
                if allMovingMarkersTags[0] in allMarkersID.keys():
                    self.fluidInterfaceIdentifier = allMarkersID[allMovingMarkersTags[0]]
            if self.fluidInterfaceIdentifier != None:
                self.nLocalFluidInterfaceNodes = FluidSolver.GetNumberVertices(self.fluidInterfaceIdentifier)
            if self.nLocalFluidInterfaceNodes != 0:
                self.haveFluidInterface = True
                print('Number of interface fluid nodes (halo nodes included) on proccess {} and marker {}: {}'\
                      .format(myid,allMovingMarkersTags[0],self.nLocalFluidInterfaceNodes))
        else:
            pass

        # --- Identify the solid interface and store the number of nodes (single core) ---#
        if SolidSolver != None:
            print('Solid solver is initialized on process {}'.format(myid))
            self.haveSolidSolver = True
            self.nSolidInterfaceNodes = SolidSolver.nPoint
            self.nSolidInterfacePhysicalNodes = SolidSolver.nPoint
            self.nLocalSolidInterfaceNodes = SolidSolver.nPoint
            self.globalSolidCoordinates = np.zeros((SolidSolver.nPoint, 3))
            for iPoint in range(0, SolidSolver.nPoint):
                coordX, coordY, coordZ = SolidSolver.getInitialCoordinates(iPoint)
                self.globalSolidCoordinates[iPoint, 0] = coordX
                self.globalSolidCoordinates[iPoint, 1] = coordY
                self.globalSolidCoordinates[iPoint, 2] = coordZ
            # self.haveSolidSolver = True
            # self.solidInterfaceIdentifier = SolidSolver.getFSIMarkerID()
            # self.nLocalSolidInterfaceNodes = SolidSolver.getNumberOfSolidInterfaceNodes(self.solidInterfaceIdentifier)
            # if self.nLocalSolidInterfaceNodes != 0:
            #     self.haveSolidInterface = True
            #     print('Number of interface solid nodes (halo nodes included) on proccess {} : {}'.format(myid,
            #                                                                             self.nLocalSolidInterfaceNodes))
            #self.nSolidInterfaceNodes = np.copy(sendBuffHalo)
            #self.nSolidInterfacePhysicalNodes = np.copy(sendBuffPhysical)
        else:
            pass

        # --- Exchange information about processors on which the solvers are defined and where the interface nodes are lying --- #
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


        self.MPIPrint(
            'Total number of fluid interface nodes (halo nodes included) : {}'.format(self.nFluidInterfaceNodes))
        self.MPIPrint(
            'Total number of physical fluid interface nodes : {}'.format(self.nFluidInterfacePhysicalNodes))
        self.MPIPrint(
            'Total number of beam interface nodes : {}'.format(self.nSolidInterfaceNodes))
        self.MPIPrint('Total number of fluid interface nodes : {}'.format(self.nFluidInterfacePhysicalNodes))
        self.MPIPrint('Total number of solid interface nodes : {}'.format(self.nSolidInterfacePhysicalNodes))

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
            posx = FluidSolver.GetVertexCoordX(self.fluidInterfaceIdentifier, iVertex)
            posy = FluidSolver.GetVertexCoordY(self.fluidInterfaceIdentifier, iVertex)
            posz = FluidSolver.GetVertexCoordZ(self.fluidInterfaceIdentifier, iVertex)

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
                self.globalFluidInterface_vertex_indices = np.empty(sum(self.sendCounts))
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
            self.globalFluidInterface_vertex_indices = fluidIndexing_temp.copy()
            self.globalFluidInterfaceXcoor = localFluidInterface_array_X_init.copy()
            self.globalFluidInterfaceYcoor = localFluidInterface_array_Y_init.copy()
            self.globalFluidInterfaceZcoor = localFluidInterface_array_Z_init.copy()

        # Store the global fluid coordinates
        if myid == self.rootProcess:
            self.globalFluidCoordinates = np.zeros((self.nFluidInterfacePhysicalNodes, 3))
            for i in range(0, self.nFluidInterfacePhysicalNodes):
                self.globalFluidCoordinates[i][0] = self.globalFluidInterfaceXcoor[i]
                self.globalFluidCoordinates[i][1] = self.globalFluidInterfaceYcoor[i]
                self.globalFluidCoordinates[i][2] = self.globalFluidInterfaceZcoor[i]

            print(self.globalFluidCoordinates.shape)
            print(self.globalSolidCoordinates.shape)

        del fluidIndexing_temp, localFluidInterface_array_X_init, \
            localFluidInterface_array_Y_init, localFluidInterface_array_Z_init


        # # --- Scatter the values of the coordinates (to test)
        # # --- The vertex indices
        # # --- The coordinates X, Y and Z
        # GlobalIndex = int()
        # localIndex = 0
        # fluidIndexing_temp = {}
        # localCoorX = np.zeros(self.nLocalFluidInterfacePhysicalNodes)
        # localCoorY = np.zeros(self.nLocalFluidInterfacePhysicalNodes)
        # localCoorZ = np.zeros(self.nLocalFluidInterfacePhysicalNodes)
        #
        # if self.have_MPI:
        #
        #     self.comm.Scatterv(sendbuf=(self.globalFluidInterfaceXcoor, self.sendCounts), recvbuf=localCoorX, root=0)
        #     self.comm.Scatterv(sendbuf=(self.globalFluidInterfaceYcoor, self.sendCounts), recvbuf=localCoorY, root=0)
        #     self.comm.Scatterv(sendbuf=(self.globalFluidInterfaceZcoor, self.sendCounts), recvbuf=localCoorZ, root=0)
        #
        #     print("rank: {}, local_array X: {}".format(myid, localCoorX))
        #     print("rank: {}, local_array Y: {}".format(myid, localCoorY))
        #     print("rank: {}, local_array Z: {}".format(myid, localCoorZ))
        #
        # else:
        #     self.globalFluidInterface_vertex_indices = fluidIndexing_temp.copy()
        #     self.globalFluidInterfaceXcoor = self.localFluidInterface_array_X_init.copy()
        #     self.globalFluidInterfaceYcoor = self.localFluidInterface_array_Y_init.copy()
        #     self.globalFluidInterfaceZcoor = self.localFluidInterface_array_Z_init.copy()



    def transferFluidTractions(self, FluidSolver, SolidSolver):
        """
        Transfer fluid tractions.
        Gathers the fluid tractions from the interface into the root process.
        Interpolates the tractions using the transposed matrix.
        Applies the tractions into the solid solver.
        """
        if self.have_MPI:
            myid = self.comm.Get_rank()
            MPIsize = self.comm.Get_size()
        else:
            myid = 0
            MPIsize = 1

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
             # Store them in the local load array
            localFluidLoadX[localIndex] = FluidSolver.GetVertexForceX(self.fluidInterfaceIdentifier, iVertex)
            localFluidLoadY[localIndex] = FluidSolver.GetVertexForceY(self.fluidInterfaceIdentifier, iVertex)
            localFluidLoadZ[localIndex] = FluidSolver.GetVertexForceZ(self.fluidInterfaceIdentifier, iVertex)
            localIndex += 1

        if self.have_MPI:

            # Store the local loads in buffers in the form of numpy arrays
            bufXLoad = np.array(localFluidLoadX)
            bufYLoad = np.array(localFluidLoadY)
            bufZLoad = np.array(localFluidLoadZ)

            # Initialize the global load array
            if myid == self.rootProcess:
                print("sendCounts: {}, total: {}".format(self.sendCounts, sum(self.sendCounts)))
                globalFluidLoadX = np.empty(sum(self.sendCounts))
                globalFluidLoadY = np.empty(sum(self.sendCounts))
                globalFluidLoadZ = np.empty(sum(self.sendCounts))

            # Gatherv using self.sendCounts maintains the ordering of the coordinates
            self.comm.Gatherv(sendbuf=bufXLoad, recvbuf=(globalFluidLoadX, self.sendCounts), root=0)
            self.comm.Gatherv(sendbuf=bufYLoad, recvbuf=(globalFluidLoadY, self.sendCounts), root=0)
            self.comm.Gatherv(sendbuf=bufZLoad, recvbuf=(globalFluidLoadZ, self.sendCounts), root=0)

            if myid == 0:
                print("Gathered array X: {}".format(globalFluidLoadX))
                print("Gathered array Y: {}".format(globalFluidLoadY))
                print("Gathered array Z: {}".format(globalFluidLoadZ))

        else:
            self.globalFluidLoadX = localFluidLoadX.copy()
            self.globalFluidLoadY = localFluidLoadY.copy()
            self.globalFluidLoadZ = localFluidLoadZ.copy()

        # Delete local variables
        del localFluidLoadX, localFluidLoadY, localFluidLoadZ

        ################################################################################################################
        # --- STEP 2: Interpolate
        ################################################################################################################

        # ---> Input: self.globalSolidDispX, self.globalSolidDispX, self.globalSolidDispX

        # ---> Output: self.globalSolidLoadX, self.globalSolidLoadY, self.globalSolidLoadZ

        ################################################################################################################
        # --- STEP 3: Check conservation
        ################################################################################################################

        #     FY = 0.0  # solid-side resultant forces
        #     FX = 0.0
        #     FZ = 0.0
        #     FFX = 0.0  # fluid-side resultant forces
        #     FFY = 0.0
        #     FFZ = 0.0
        #
        #     # --- Check for total force conservation after interpolation
        #     FFX = self.fluidLoads_array_X.sum()
        #     FFY = self.fluidLoads_array_Y.sum()
        #     FFZ = self.fluidLoads_array_Z.sum()
        #
        #     for iVertex in range(self.nLocalSolidInterfaceNodes):
        #         FX += self.localSolidLoads_array_X[iVertex]
        #         FY += self.localSolidLoads_array_Y[iVertex]
        #         FZ += self.localSolidLoads_array_Z[iVertex]
        #
        #     if self.have_MPI == True:
        #         FX = self.comm.allreduce(FX)
        #         FY = self.comm.allreduce(FY)
        #         FZ = self.comm.allreduce(FZ)
        #
        #     self.MPIPrint("Checking f/s interface total force...")
        #     self.MPIPrint('Solid side (Fx, Fy, Fz) = ({}, {}, {})'.format(FX, FY, FZ))
        #     self.MPIPrint('Fluid side (Fx, Fy, Fz) = ({}, {}, {})'.format(FFX, FFY, FFZ))

        ################################################################################################################
        # --- STEP 4: Transfer to the structural solver
        # --- pyBeam runs in single core, so there is no need to deal with the parallelization
        ################################################################################################################

        # For the vertices that belong to the interface
        for iVertex in range(0, self.nSolidInterfaceNodes):
             # Store them in the solid solver directly
            SolidSolver.SetLoads(iVertex, 0, self.globalSolidDispX[iVertex])
            SolidSolver.SetLoads(iVertex, 1, self.globalSolidDispY[iVertex])
            SolidSolver.SetLoads(iVertex, 2, self.globalSolidDispZ[iVertex])


    def transferStructuralDisplacements(self, FluidSolver, SolidSolver, MLSSolver):
        """
        Transfer structural displacements tractions.
        Gathers the structural displacements from the interface.
        Interpolates the displacements using the transposed matrix.
        Applies the fluid displacements by scattering them into the correct position for the fluid solver.
        """
        if self.have_MPI:
            myid = self.comm.Get_rank()
            MPIsize = self.comm.Get_size()
        else:
            myid = 0
            MPIsize = 1

        ################################################################################################################
        # --- STEP 1: Retrieve the structural displacements from pyBeam
        # --- pyBeam runs in single core, so there is no need to deal with the parallelization
        ################################################################################################################

        # Initialize the local load array
        self.globalSolidDispX = np.zeros(self.nSolidInterfaceNodes)
        self.globalSolidDispY = np.zeros(self.nSolidInterfaceNodes)
        self.globalSolidDispZ = np.zeros(self.nSolidInterfaceNodes)

        # For the vertices that belong to the interface
        for iVertex in range(0, self.nSolidInterfaceNodes):
             # Store them in the global load array directly
            self.globalSolidDispX[iVertex] = SolidSolver.ExtractDisplacements(iVertex, 0)
            self.globalSolidDispY[iVertex] = SolidSolver.ExtractDisplacements(iVertex, 1)
            self.globalSolidDispZ[iVertex] = SolidSolver.ExtractDisplacements(iVertex, 2)

        ################################################################################################################
        # --- STEP 2: Interpolate
        ################################################################################################################
        # ---> Input: self.globalSolidDispX, self.globalSolidDispX, self.globalSolidDispX

        self.globalFluidDispX = MLSSolver.interpolation_matrix.dot(self.globalSolidDispX)
        self.globalFluidDispY = MLSSolver.interpolation_matrix.dot(self.globalSolidDispY)
        self.globalFluidDispZ = MLSSolver.interpolation_matrix.dot(self.globalSolidDispZ)

        # ---> Output: self.globalFluidDispX, self.globalFluidDispY, self.globalFluidDispZ

        ################################################################################################################
        # --- STEP 3: Check conservation
        ################################################################################################################

        #     # --- Checking conservation ---
        #     WSX = self.solidLoads_array_X.dot(self.solidInterface_array_DispX)
        #     WSY = self.solidLoads_array_Y.dot(self.solidInterface_array_DispY)
        #     WSZ = self.solidLoads_array_Z.dot(self.solidInterface_array_DispZ)
        #
        #     WFX = self.fluidLoads_array_X.dot(self.fluidInterface_array_DispX)
        #     WFY = self.fluidLoads_array_Y.dot(self.fluidInterface_array_DispY)
        #     WFZ = self.fluidLoads_array_Z.dot(self.fluidInterface_array_DispZ)
        #
        #     self.MPIPrint("Checking f/s interface conservation...")
        #     self.MPIPrint('Solid side (Wx, Wy, Wz) = ({}, {}, {})'.format(WSX, WSY, WSZ))
        #     self.MPIPrint('Fluid side (Wx, Wy, Wz) = ({}, {}, {})'.format(WFX, WFY, WFZ))

        ################################################################################################################
        # --- STEP 4: Transfer to the fluid solver
        ################################################################################################################

        # --- Recover them from the interpolated vectors
        localFluidDispX = np.zeros(self.nLocalFluidInterfacePhysicalNodes)
        localFluidDispY = np.zeros(self.nLocalFluidInterfacePhysicalNodes)
        localFluidDispZ = np.zeros(self.nLocalFluidInterfacePhysicalNodes)

        if self.have_MPI:

            self.comm.Scatterv(sendbuf=(self.globalFluidDispX, self.sendCounts), recvbuf=localFluidDispX, root=0)
            self.comm.Scatterv(sendbuf=(self.globalFluidDispY, self.sendCounts), recvbuf=localFluidDispY, root=0)
            self.comm.Scatterv(sendbuf=(self.globalFluidDispZ, self.sendCounts), recvbuf=localFluidDispZ, root=0)

            print("rank: {}, local_array X: {}".format(myid, localFluidDispX))
            print("rank: {}, local_array Y: {}".format(myid, localFluidDispY))
            print("rank: {}, local_array Z: {}".format(myid, localFluidDispZ))

        else:
            localFluidDispX = self.globalFluidDispX.copy()
            localFluidDispY = self.globalFluidDispY.copy()
            localFluidDispZ = self.globalFluidDispZ.copy()

        # For the vertices that belong to the interface
        localIndex = 0
        for iVertex in self.localFluidInterface_vertex_indices:
            # Store them in the mesh displacement routine
            FluidSolver.SetMeshDisplacement(self.fluidInterfaceIdentifier, int(iVertex), localFluidDispX[localIndex],
                                            localFluidDispY[localIndex], localFluidDispZ[localIndex])
            # Increment the local index
            localIndex += 1

        # Delete local variables
        del localFluidDispX, localFluidDispY, localFluidDispZ


    #
    # def interfaceMapping(self, FluidSolver, SolidSolver, FSI_config):
    #     """
    #     Creates the one-to-one mapping between interfaces in case of matching meshes.
    #     Creates the interpolation rules between interfaces in case of non-matching meshes.
    #     """
    #     if self.have_MPI:
    #         myid = self.comm.Get_rank()
    #         MPIsize = self.comm.Get_size()
    #     else:
    #         myid = 0
    #         MPIsize = 1
    #
    #     # --- Get the fluid interface from fluid solver on each partition ---
    #     GlobalIndex = int()
    #     localIndex = 0
    #     fluidIndexing_temp = {}
    #     # fluidnode = open( "./Output/fluid_nodes" + '.dat' , "w")  #     HARD CODED
    #     self.localFluidInterface_array_X_init = np.zeros((self.nLocalFluidInterfacePhysicalNodes))
    #     self.localFluidInterface_array_Y_init = np.zeros((self.nLocalFluidInterfacePhysicalNodes))
    #     self.localFluidInterface_array_Z_init = np.zeros((self.nLocalFluidInterfacePhysicalNodes))
    #
    #     for iVertex in range(self.nLocalFluidInterfaceNodes):
    #         GlobalIndex = FluidSolver.GetVertexGlobalIndex(self.fluidInterfaceIdentifier, iVertex)
    #         posx = FluidSolver.GetVertexCoordX(self.fluidInterfaceIdentifier, iVertex)
    #         posy = FluidSolver.GetVertexCoordY(self.fluidInterfaceIdentifier, iVertex)
    #         posz = FluidSolver.GetVertexCoordZ(self.fluidInterfaceIdentifier, iVertex)
    #
    #         if GlobalIndex in self.FluidHaloNodeList[myid].keys():
    #             self.haloNodesPositionsInit[GlobalIndex] = (posx, posy, posz)
    #         else:
    #             fluidIndexing_temp[GlobalIndex] = self.__getGlobalIndex('fluid', myid, localIndex)
    #             self.localFluidInterface_array_X_init[localIndex] = posx
    #             self.localFluidInterface_array_Y_init[localIndex] = posy
    #             self.localFluidInterface_array_Z_init[localIndex] = posz
    #             localIndex += 1
    #
    #     if self.have_MPI:
    #         fluidIndexing_temp = self.comm.allgather(fluidIndexing_temp)
    #         for ii in range(len(fluidIndexing_temp)):
    #             for key, value in fluidIndexing_temp[ii].items():
    #                 self.fluidIndexing[key] = value
    #     else:
    #         self.fluidIndexing = fluidIndexing_temp.copy()
    #     del fluidIndexing_temp
    #
    #
    #     # --- Get the solid interface from solid solver on each partition ---
    #     localIndex = 0
    #     solidIndexing_temp = {}
    #     self.localSolidInterface_array_X = np.zeros(self.nLocalSolidInterfaceNodes)
    #     self.localSolidInterface_array_Y = np.zeros(self.nLocalSolidInterfaceNodes)
    #     self.localSolidInterface_array_Z = np.zeros(self.nLocalSolidInterfaceNodes)
    #
    #     for iVertex in range(self.nLocalSolidInterfaceNodes):
    #         GlobalIndex = SolidSolver.getInterfaceNodeGlobalIndex(self.solidInterfaceIdentifier, iVertex)
    #
    #         posx = SolidSolver.getInterfaceNodePosX(self.solidInterfaceIdentifier, iVertex)
    #         posy = SolidSolver.getInterfaceNodePosY(self.solidInterfaceIdentifier, iVertex)
    #         posz = SolidSolver.getInterfaceNodePosZ(self.solidInterfaceIdentifier, iVertex)
    #
    #         if GlobalIndex in self.SolidHaloNodeList[myid].keys():
    #             pass
    #         else:
    #             solidIndexing_temp[GlobalIndex] = self.__getGlobalIndex('solid', myid, localIndex)
    #             self.localSolidInterface_array_X[localIndex] = posx
    #             self.localSolidInterface_array_Y[localIndex] = posy
    #             self.localSolidInterface_array_Z[localIndex] = posz
    #             localIndex += 1
    #     if self.have_MPI:
    #         solidIndexing_temp = self.comm.allgather(solidIndexing_temp)
    #         for ii in range(len(solidIndexing_temp)):
    #             for key, value in solidIndexing_temp[ii].items():
    #                 self.solidIndexing[key] = value
    #     else:
    #         self.solidIndexing = solidIndexing_temp.copy()
    #     del solidIndexing_temp
    #     # print("DEBUG MESSAGE from proc {}, solidIndexing = {}".format(myid, self.solidIndexing))
    #
    #     if self.have_MPI == True:
    #         self.MappingMatrix = PETSc.Mat().create(self.comm)
    #         self.MappingMatrix_T = PETSc.Mat().create(self.comm)
    #         self.MappingMatrix.setType('mpiaij')
    #         self.MappingMatrix_T.setType('mpiaij')
    #     else:
    #         self.MappingMatrix = PETSc.Mat().create()
    #         self.MappingMatrix_T = PETSc.Mat().create()
    #         self.MappingMatrix.setType('aij')
    #         self.MappingMatrix_T.setType('aij')
    #
    #     self.MappingMatrix.setSizes((self.nFluidInterfacePhysicalNodes, self.nSolidInterfacePhysicalNodes))
    #     self.MappingMatrix.setUp()
    #     self.MappingMatrix.setOption(PETSc.Mat().Option.NEW_NONZERO_ALLOCATION_ERR, False)
    #     self.MappingMatrix_T.setSizes((self.nSolidInterfacePhysicalNodes, self.nFluidInterfacePhysicalNodes))
    #     self.MappingMatrix_T.setUp()
    #     self.MappingMatrix_T.setOption(PETSc.Mat().Option.NEW_NONZERO_ALLOCATION_ERR, False)
    #
    #     # --- Fill the interpolation matrix in parallel (working in serial too) --- #
    #     self.MPIPrint("Building interpolation matrix...")
    #     self.MPIBarrier()
    #
    #     # Appearently I have to divide the size of the data for Send/Rcv: the treshold size is defined in init()
    #
    #     if self.have_MPI == True:
    #         for iProc in self.solidInterfaceProcessors:
    #             if myid == iProc:
    #                 for jProc in self.fluidInterfaceProcessors:
    #                     if self.localSolidInterface_array_X.size > self.MPI_trhld:
    #                         iMPI = 0
    #                         while iMPI < self.localSolidInterface_array_X.size:
    #                             # print("HARD CODE: Check. self.localSolidInterface_array_X.size = {}, iMPI = {}".format(self.localSolidInterface_array_X.size, iMPI))
    #                             if self.localSolidInterface_array_X.size - iMPI > self.MPI_trhld:
    #                                 self.comm.Send(self.localSolidInterface_array_X[iMPI:iMPI + self.MPI_trhld],
    #                                                dest=jProc, tag=iMPI + 1)
    #                                 self.comm.Send(self.localSolidInterface_array_Y[iMPI:iMPI + self.MPI_trhld],
    #                                                dest=jProc, tag=iMPI + 2)
    #                                 self.comm.Send(self.localSolidInterface_array_Z[iMPI:iMPI + self.MPI_trhld],
    #                                                dest=jProc, tag=iMPI + 3)
    #                             else:
    #                                 self.comm.Send(
    #                                     self.localSolidInterface_array_X[iMPI:self.localSolidInterface_array_X.size],
    #                                     dest=jProc, tag=iMPI + 1)
    #                                 self.comm.Send(
    #                                     self.localSolidInterface_array_Y[iMPI:self.localSolidInterface_array_X.size],
    #                                     dest=jProc, tag=iMPI + 2)
    #                                 self.comm.Send(
    #                                     self.localSolidInterface_array_Z[iMPI:self.localSolidInterface_array_X.size],
    #                                     dest=jProc, tag=iMPI + 3)
    #                             iMPI = iMPI + self.MPI_trhld
    #                             # print(iMPI)
    #                     else:
    #                         self.comm.Send(self.localSolidInterface_array_X, dest=jProc, tag=1)
    #                         self.comm.Send(self.localSolidInterface_array_Y, dest=jProc, tag=2)
    #                         self.comm.Send(self.localSolidInterface_array_Z, dest=jProc, tag=3)
    #                     # print("DEBUG MESSAGE From proc {}, sent to {}".format(myid, jProc))    # HARD CODED
    #             if myid in self.fluidInterfaceProcessors:
    #                 sizeOfBuff = self.solidPhysicalInterfaceNodesDistribution[iProc]
    #                 solidInterfaceBuffRcv_X = np.zeros(sizeOfBuff)
    #                 solidInterfaceBuffRcv_Y = np.zeros(sizeOfBuff)
    #                 solidInterfaceBuffRcv_Z = np.zeros(sizeOfBuff)
    #                 if sizeOfBuff > self.MPI_trhld:  # Can't compare to the thershold in case of receiving the self.localSolidInterface_array_X.size as some processors don't have such info the equivalent is sizeOfBuff
    #                     iMPI = 0
    #                     while iMPI < sizeOfBuff:
    #                         if sizeOfBuff - iMPI > self.MPI_trhld:
    #                             self.comm.Recv(solidInterfaceBuffRcv_X[iMPI:iMPI + self.MPI_trhld], iProc, tag=iMPI + 1)
    #                             self.comm.Recv(solidInterfaceBuffRcv_Y[iMPI:iMPI + self.MPI_trhld], iProc, tag=iMPI + 2)
    #                             self.comm.Recv(solidInterfaceBuffRcv_Z[iMPI:iMPI + self.MPI_trhld], iProc, tag=iMPI + 3)
    #                         else:
    #                             self.comm.Recv(solidInterfaceBuffRcv_X[iMPI:sizeOfBuff], iProc, tag=iMPI + 1)
    #                             self.comm.Recv(solidInterfaceBuffRcv_Y[iMPI:sizeOfBuff], iProc, tag=iMPI + 2)
    #                             self.comm.Recv(solidInterfaceBuffRcv_Z[iMPI:sizeOfBuff], iProc, tag=iMPI + 3)
    #                         iMPI = iMPI + self.MPI_trhld
    #                         # print(iMPI)
    #                 else:
    #                     self.comm.Recv(solidInterfaceBuffRcv_X, iProc, tag=1)
    #                     self.comm.Recv(solidInterfaceBuffRcv_Y, iProc, tag=2)
    #                     self.comm.Recv(solidInterfaceBuffRcv_Z, iProc, tag=3)
    #                 # print((self.localSolidInterface_array_X - solidInterfaceBuffRcv_X))                    # HARD CODED
    #                 # print("DEBUG MESSAGE : From processor {}, received from {}".format(myid, iProc))     # HARD CODED
    #                 if FSI_config['MATCHING_MESH'] == 'NO':
    #                     if FSI_config['MESH_INTERP_METHOD'] == 'RBF':
    #                         self.RBFMeshMapping_B(solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y,
    #                                               solidInterfaceBuffRcv_Z, iProc, self.RBF_rad)
    #                     elif FSI_config['MESH_INTERP_METHOD'] == 'TPS':
    #                         self.TPSMeshMapping_B(solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y,
    #                                               solidInterfaceBuffRcv_Z, iProc)
    #                     else:
    #                         self.NearestNeighboorMeshMapping(solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y,
    #                                                          solidInterfaceBuffRcv_Z, iProc)
    #                 else:
    #                     self.matchingMeshMapping(solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y,
    #                                              solidInterfaceBuffRcv_Z, iProc)
    #     else:
    #         self.matchingMeshMapping(self.localSolidInterface_array_X, self.localSolidInterface_array_Y,
    #                                  self.localSolidInterface_array_Z, 0)
    #
    #     self.MappingMatrix.assemblyBegin()
    #     self.MappingMatrix.assemblyEnd()
    #     self.MappingMatrix_T.assemblyBegin()
    #     self.MappingMatrix_T.assemblyEnd()
    #     self.MPIPrint("Interpolation matrix is built.")
    #
    #     self.MPIBarrier()
    #
    #     # print("DEBUG MESSAGE From proc {} : MAPPING DONE !".format(myid))
    #     del self.localSolidInterface_array_X
    #     del self.localSolidInterface_array_Y
    #     del self.localSolidInterface_array_Z
    #
    # def matchingMeshMapping(self, solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, iProc):
    #     """
    #     Fill the mapping matrix in case of matching meshes at the f/s interface.
    #     """
    #     if self.have_MPI == True:
    #         myid = self.comm.Get_rank()
    #     else:
    #         myid = 0
    #
    #     # --- Instantiate the spatial indexing ---
    #     prop_index = index.Property()
    #     prop_index.dimension = self.nDim
    #     SolidSpatialTree = index.Index(properties=prop_index)
    #
    #     nSolidNodes = solidInterfaceBuffRcv_X.shape[0]
    #
    #     for jVertex in range(nSolidNodes):
    #         posX = solidInterfaceBuffRcv_X[jVertex]
    #         posY = solidInterfaceBuffRcv_Y[jVertex]
    #         posZ = solidInterfaceBuffRcv_Z[jVertex]
    #         if self.nDim == 2:
    #             SolidSpatialTree.add(jVertex, (posX, posY))
    #         else:
    #             SolidSpatialTree.add(jVertex, (posX, posY, posZ))
    #
    #     if self.nFluidInterfacePhysicalNodes != self.nSolidInterfacePhysicalNodes:
    #         raise Exception("Fluid and solid interface must have the same number of nodes for matching meshes ! ")
    #
    #     # --- For each fluid interface node, find the nearest solid interface node and fill the boolean mapping matrix ---
    #     for iVertexFluid in range(self.nLocalFluidInterfacePhysicalNodes):
    #         posX = self.localFluidInterface_array_X_init[iVertexFluid]
    #         posY = self.localFluidInterface_array_Y_init[iVertexFluid]
    #         posZ = self.localFluidInterface_array_Z_init[iVertexFluid]
    #         if self.nDim == 2:
    #             neighboors = list(SolidSpatialTree.nearest((posX, posY), 1))
    #         elif self.nDim == 3:
    #             neighboors = list(SolidSpatialTree.nearest((posX, posY, posZ), 1))
    #         jVertexSolid = neighboors[0]
    #         # Check if the distance is small enough to ensure coincidence
    #         NodeA = np.array([posX, posY, posZ])
    #         NodeB = np.array([solidInterfaceBuffRcv_X[jVertexSolid], solidInterfaceBuffRcv_Y[jVertexSolid],
    #                           solidInterfaceBuffRcv_Z[jVertexSolid]])
    #         distance = spdist.euclidean(NodeA, NodeB)
    #         iGlobalVertexFluid = self.__getGlobalIndex('fluid', myid, iVertexFluid)
    #         jGlobalVertexSolid = self.__getGlobalIndex('solid', iProc, jVertexSolid)
    #         if distance > 1e-6:
    #             print(
    #                 "WARNING : Tolerance for matching meshes is not matched between node F{} and S{} : ({}, {}, {})<-->({}, {}, {}) , DISTANCE : {} !".format(
    #                     iGlobalVertexFluid, jGlobalVertexSolid, posX, posY, posZ, solidInterfaceBuffRcv_X[jVertexSolid],
    #                     solidInterfaceBuffRcv_Y[jVertexSolid], solidInterfaceBuffRcv_Z[jVertexSolid], distance))
    #         # print("DEBUG MESSAGE From proc {} : Local fluid : {} , Global fluid : {}".format(myid, iVertexFluid, iGlobalVertexFluid))
    #         # print("DEBUG MESSAGE From proc {} : Local solid : {} , Global solid : {}".format(myid, jVertexSolid, jGlobalVertexSolid))
    #         # print("DEBUG MESSAGE From Proc {} with globalFluid = {}-({},{},{}) and globalSolid = {}-({},{},{}) ".format(myid, iGlobalVertexFluid, posX, posY, posZ, jGlobalVertexSolid,solidInterfaceBuffRcv_X[jVertexSolid], solidInterfaceBuffRcv_Y[jVertexSolid], solidInterfaceBuffRcv_Z[jVertexSolid] ))
    #         self.MappingMatrix.setValue(iGlobalVertexFluid, jGlobalVertexSolid, 1.0)
    #         self.MappingMatrix_T.setValue(jGlobalVertexSolid, iGlobalVertexFluid, 1.0)
    #
    #     del solidInterfaceBuffRcv_X
    #     del solidInterfaceBuffRcv_Y
    #     del solidInterfaceBuffRcv_Z
    #
    # def interpolateSolidPositionOnFluidMesh(self, FSI_config):
    #     """
    #     Applies the one-to-one mapping or the interpolaiton rules from solid to fluid mesh.
    #     """
    #     if self.have_MPI == True:
    #         myid = self.comm.Get_rank()
    #         MPIsize = self.comm.Get_size()
    #     else:
    #         myid = 0
    #         MPIsize = 1
    #
    #     # --- Interpolate (or map) in parallel the solid interface displacement on the fluid interface ---
    #     # print("DEBUG MESSAGE From proc {} : Ready for interpolation !".format(myid))
    #     if FSI_config['MATCHING_MESH'] == 'NO' and (
    #             FSI_config['MESH_INTERP_METHOD'] == 'RBF' or FSI_config['MESH_INTERP_METHOD'] == 'TPS'):
    #         if self.have_MPI == True:
    #             gamma_array_DispX = PETSc.Vec().create(self.comm)
    #             gamma_array_DispY = PETSc.Vec().create(self.comm)
    #             gamma_array_DispZ = PETSc.Vec().create(self.comm)
    #             gamma_array_DispX.setType('mpi')
    #             gamma_array_DispY.setType('mpi')
    #             gamma_array_DispZ.setType('mpi')
    #             KSP_solver = PETSc.KSP().create(self.comm)
    #         else:
    #             gamma_array_DispX = PETSc.Vec().create()
    #             gamma_array_DispY = PETSc.Vec().create()
    #             gamma_array_DispZ = PETSc.Vec().create()
    #             gamma_array_DispX.setType('seq')
    #             gamma_array_DispY.setType('seq')
    #             gamma_array_DispZ.setType('seq')
    #             KSP_solver = PETSc.KSP().create()
    #         gamma_array_DispX.setSizes(self.nSolidInterfacePhysicalNodes + self.d_RBF)
    #         gamma_array_DispY.setSizes(self.nSolidInterfacePhysicalNodes + self.d_RBF)
    #         gamma_array_DispZ.setSizes(self.nSolidInterfacePhysicalNodes + self.d_RBF)
    #         gamma_array_DispX.set(0.0)
    #         gamma_array_DispY.set(0.0)
    #         gamma_array_DispZ.set(0.0)
    #         KSP_solver.setType('fgmres')
    #         KSP_solver.getPC().setType('jacobi')
    #         KSP_solver.setOperators(self.MappingMatrixA)
    #         KSP_solver.setFromOptions()
    #         # print(KSP_solver.getInitialGuessNonzero())
    #         KSP_solver.setInitialGuessNonzero(True)
    #         # print(KSP_solver.getInitialGuessNonzero())
    #         KSP_solver.solve(self.solidInterface_array_DispX, gamma_array_DispX)
    #         KSP_solver.solve(self.solidInterface_array_DispY, gamma_array_DispY)
    #         KSP_solver.solve(self.solidInterface_array_DispZ, gamma_array_DispZ)
    #         self.MappingMatrixB.mult(gamma_array_DispX, self.fluidInterface_array_DispX)
    #         self.MappingMatrixB.mult(gamma_array_DispY, self.fluidInterface_array_DispY)
    #         self.MappingMatrixB.mult(gamma_array_DispZ, self.fluidInterface_array_DispZ)
    #         gamma_array_DispX.destroy()
    #         gamma_array_DispY.destroy()
    #         gamma_array_DispZ.destroy()
    #         KSP_solver.destroy()
    #         del gamma_array_DispX
    #         del gamma_array_DispY
    #         del gamma_array_DispZ
    #         del KSP_solver
    #     else:
    #         # print(self.solidInterface_array_DispX)  # HARD CODE
    #         # print('\n')                             # HARD CODE
    #         # print(self.solidInterface_array_DispY)  # HARD CODE
    #         # print('\n')                             # HARD CODE
    #         # print(self.solidInterface_array_DispZ)  # HARD CODE
    #         self.MappingMatrix.mult(self.solidInterface_array_DispX,
    #                                 self.fluidInterface_array_DispX)  # linea problema Freeze
    #         self.MappingMatrix.mult(self.solidInterface_array_DispY, self.fluidInterface_array_DispY)
    #         self.MappingMatrix.mult(self.solidInterface_array_DispZ, self.fluidInterface_array_DispZ)
    #
    #     # --- Checking conservation ---
    #     WSX = self.solidLoads_array_X.dot(self.solidInterface_array_DispX)
    #     WSY = self.solidLoads_array_Y.dot(self.solidInterface_array_DispY)
    #     WSZ = self.solidLoads_array_Z.dot(self.solidInterface_array_DispZ)
    #
    #     WFX = self.fluidLoads_array_X.dot(self.fluidInterface_array_DispX)
    #     WFY = self.fluidLoads_array_Y.dot(self.fluidInterface_array_DispY)
    #     WFZ = self.fluidLoads_array_Z.dot(self.fluidInterface_array_DispZ)
    #
    #     self.MPIPrint("Checking f/s interface conservation...")
    #     self.MPIPrint('Solid side (Wx, Wy, Wz) = ({}, {}, {})'.format(WSX, WSY, WSZ))
    #     self.MPIPrint('Fluid side (Wx, Wy, Wz) = ({}, {}, {})'.format(WFX, WFY, WFZ))
    #
    #     # --- Redistribute the interpolated fluid interface according to the partitions that own the fluid interface ---
    #     # Gather the fluid interface on the master process
    #     if self.have_MPI:
    #         sendBuff_X = None
    #         sendBuff_Y = None
    #         sendBuff_Z = None
    #         self.fluidInterface_array_DispX_recon = None
    #         self.fluidInterface_array_DispY_recon = None
    #         self.fluidInterface_array_DispZ_recon = None
    #
    #         if myid == self.rootProcess:
    #             self.fluidInterface_array_DispX_recon = np.zeros(self.nFluidInterfacePhysicalNodes)
    #             self.fluidInterface_array_DispY_recon = np.zeros(self.nFluidInterfacePhysicalNodes)
    #             self.fluidInterface_array_DispZ_recon = np.zeros(self.nFluidInterfacePhysicalNodes)
    #
    #         myNumberOfNodes = self.fluidInterface_array_DispX.getArray().shape[0]
    #         sendBuffNumber = np.array([myNumberOfNodes], dtype=int)
    #         rcvBuffNumber = np.zeros(MPIsize, dtype=int)
    #         self.comm.Allgather(sendBuffNumber, rcvBuffNumber)
    #
    #         counts = tuple(rcvBuffNumber)
    #         displ = np.zeros(MPIsize, dtype=int)
    #         for ii in range(rcvBuffNumber.shape[0]):
    #             displ[ii] = rcvBuffNumber[0:ii].sum()
    #         displ = tuple(displ)
    #
    #         del sendBuffNumber, rcvBuffNumber
    #
    #         # print("DEBUG MESSAGE From proc {}, counts = {}".format(myid, counts))
    #         # print("DEBUG MESSAGE From proc {}, displ = {}".format(myid, displ))
    #
    #         self.comm.Gatherv(self.fluidInterface_array_DispX.getArray(),
    #                           [self.fluidInterface_array_DispX_recon, counts, displ, self.MPI.DOUBLE],
    #                           root=self.rootProcess)
    #         self.comm.Gatherv(self.fluidInterface_array_DispY.getArray(),
    #                           [self.fluidInterface_array_DispY_recon, counts, displ, self.MPI.DOUBLE],
    #                           root=self.rootProcess)
    #         self.comm.Gatherv(self.fluidInterface_array_DispZ.getArray(),
    #                           [self.fluidInterface_array_DispZ_recon, counts, displ, self.MPI.DOUBLE],
    #                           root=self.rootProcess)
    #
    #         # Send the partitioned interface to the right fluid partitions
    #         if myid == self.rootProcess:
    #             for iProc in self.fluidInterfaceProcessors:
    #                 sendBuff_X = np.zeros(self.fluidPhysicalInterfaceNodesDistribution[iProc])
    #                 sendBuff_Y = np.zeros(self.fluidPhysicalInterfaceNodesDistribution[iProc])
    #                 sendBuff_Z = np.zeros(self.fluidPhysicalInterfaceNodesDistribution[iProc])
    #                 globalIndex = self.fluidGlobalIndexRange[iProc][iProc][0]
    #                 for iVertex in range(self.fluidPhysicalInterfaceNodesDistribution[iProc]):
    #                     sendBuff_X[iVertex] = self.fluidInterface_array_DispX_recon[globalIndex]
    #                     sendBuff_Y[iVertex] = self.fluidInterface_array_DispY_recon[globalIndex]
    #                     sendBuff_Z[iVertex] = self.fluidInterface_array_DispZ_recon[globalIndex]
    #                     globalIndex += 1
    #                 self.comm.Send(sendBuff_X, dest=iProc, tag=1)
    #                 self.comm.Send(sendBuff_Y, dest=iProc, tag=2)
    #                 self.comm.Send(sendBuff_Z, dest=iProc, tag=3)
    #         if myid in self.fluidInterfaceProcessors:
    #             self.localFluidInterface_array_DispX = np.zeros(self.nLocalFluidInterfacePhysicalNodes)
    #             self.localFluidInterface_array_DispY = np.zeros(self.nLocalFluidInterfacePhysicalNodes)
    #             self.localFluidInterface_array_DispZ = np.zeros(self.nLocalFluidInterfacePhysicalNodes)
    #             self.comm.Recv(self.localFluidInterface_array_DispX, source=self.rootProcess, tag=1)
    #             self.comm.Recv(self.localFluidInterface_array_DispY, source=self.rootProcess, tag=2)
    #             self.comm.Recv(self.localFluidInterface_array_DispZ, source=self.rootProcess, tag=3)
    #         del sendBuff_X
    #         del sendBuff_Y
    #         del sendBuff_Z
    #     else:
    #         self.localFluidInterface_array_DispX = self.fluidInterface_array_DispX.getArray().copy()
    #         self.localFluidInterface_array_DispY = self.fluidInterface_array_DispY.getArray().copy()
    #         self.localFluidInterface_array_DispZ = self.fluidInterface_array_DispZ.getArray().copy()
    #
    #     # Special treatment for the halo nodes on the fluid interface
    #     self.haloNodesDisplacements = {}
    #     sendBuff = {}
    #     if self.have_MPI == True:
    #         if myid == self.rootProcess:
    #             for iProc in self.fluidInterfaceProcessors:
    #                 sendBuff = {}
    #                 for key in self.FluidHaloNodeList[iProc].keys():
    #                     globalIndex = self.fluidIndexing[key]
    #                     DispX = self.fluidInterface_array_DispX_recon[globalIndex]
    #                     DispY = self.fluidInterface_array_DispY_recon[globalIndex]
    #                     DispZ = self.fluidInterface_array_DispZ_recon[globalIndex]
    #                     sendBuff[key] = (DispX, DispY, DispZ)
    #                 self.comm.send(sendBuff, dest=iProc, tag=4)
    #         if myid in self.fluidInterfaceProcessors:
    #             self.haloNodesDisplacements = self.comm.recv(source=self.rootProcess, tag=4)
    #         del sendBuff
    #
    #     # print("DEBUG MESSAGE From proc {} , received haloNodesDisplacements = {}".format(myid, self.haloNodesDisplacements))
    #
    # def interpolateFluidLoadsOnSolidMesh(self, FSI_config):
    #     """
    #     Applies the one-to-one mapping or the interpolation rules from fluid to solid mesh.
    #     """
    #     if self.have_MPI == True:
    #         myid = self.comm.Get_rank()
    #         MPIsize = self.comm.Get_size()
    #     else:
    #         myid = 0
    #         MPIsize = 1
    #
    #     self.MappingMatrix_T.mult(self.fluidLoads_array_X, self.solidLoads_array_X)
    #     self.MappingMatrix_T.mult(self.fluidLoads_array_Y, self.solidLoads_array_Y)
    #     self.MappingMatrix_T.mult(self.fluidLoads_array_Z, self.solidLoads_array_Z)
    #
    #     # --- Checking conservation ---
    #     # FSX = self.solidLoads_array_X.sum()
    #     # FSY = self.solidLoads_array_Y.sum()
    #     # FSZ = self.solidLoads_array_Z.sum()
    #     # FFX = self.fluidLoads_array_X.sum()
    #     # FFY = self.fluidLoads_array_Y.sum()
    #     # FFZ = self.fluidLoads_array_Z.sum()
    #
    #     # if myid == MPIsize-1:
    #     # myLocalNumberOfElem = self.solidLoads_array_X.getArray().shape[0]
    #     # print(self.solidLoads_array_X.getArray()[myLocalNumberOfElem-4:])
    #     # FSX -= self.solidLoads_array_X.getArray()[myLocalNumberOfElem-4:].sum()
    #     # FSY -= self.solidLoads_array_Y.getArray()[myLocalNumberOfElem-4:].sum()
    #     # FSZ -= self.solidLoads_array_Z.getArray()[myLocalNumberOfElem-4:].sum()
    #     # print("Total fluid force = ({}, {}, {})".format(FFX, FFY, FFZ))
    #     # print("Total solid force = ({}, {}, {})".format(FSX, FSY, FSZ))
    #
    #     # --- Redistribute the interpolated solid loads according to the partitions that own the solid interface ---
    #     # Gather the solid loads on the master process
    #     if self.have_MPI:
    #         sendBuff_X = None
    #         sendBuff_Y = None
    #         sendBuff_Z = None
    #         self.solidLoads_array_X_recon = None
    #         self.solidLoads_array_Y_recon = None
    #         self.solidLoads_array_Z_recon = None
    #         if myid == self.rootProcess:
    #             self.solidLoads_array_X_recon = np.zeros(self.nSolidInterfacePhysicalNodes + self.d_RBF)
    #             self.solidLoads_array_Y_recon = np.zeros(self.nSolidInterfacePhysicalNodes + self.d_RBF)
    #             self.solidLoads_array_Z_recon = np.zeros(self.nSolidInterfacePhysicalNodes + self.d_RBF)
    #         myNumberOfNodes = self.solidLoads_array_X.getArray().shape[0]
    #         sendBuffNumber = np.array([myNumberOfNodes], dtype=int)
    #         rcvBuffNumber = np.zeros(MPIsize, dtype=int)
    #         self.comm.Allgather(sendBuffNumber, rcvBuffNumber)
    #
    #         counts = tuple(rcvBuffNumber)
    #         displ = np.zeros(MPIsize, dtype=int)
    #         for ii in range(rcvBuffNumber.shape[0]):
    #             displ[ii] = rcvBuffNumber[0:ii].sum()
    #         displ = tuple(displ)
    #
    #         del sendBuffNumber, rcvBuffNumber
    #
    #         self.comm.Gatherv(self.solidLoads_array_X.getArray(),
    #                           [self.solidLoads_array_X_recon, counts, displ, self.MPI.DOUBLE], root=self.rootProcess)
    #         self.comm.Gatherv(self.solidLoads_array_Y.getArray(),
    #                           [self.solidLoads_array_Y_recon, counts, displ, self.MPI.DOUBLE], root=self.rootProcess)
    #         self.comm.Gatherv(self.solidLoads_array_Z.getArray(),
    #                           [self.solidLoads_array_Z_recon, counts, displ, self.MPI.DOUBLE], root=self.rootProcess)
    #
    #         # Send the partitioned loads to the right solid partitions
    #         if myid == self.rootProcess:
    #             for iProc in self.solidInterfaceProcessors:
    #                 sendBuff_X = np.zeros(self.solidPhysicalInterfaceNodesDistribution[iProc])
    #                 sendBuff_Y = np.zeros(self.solidPhysicalInterfaceNodesDistribution[iProc])
    #                 sendBuff_Z = np.zeros(self.solidPhysicalInterfaceNodesDistribution[iProc])
    #                 globalIndex = self.solidGlobalIndexRange[iProc][iProc][0]
    #                 for iVertex in range(self.solidPhysicalInterfaceNodesDistribution[iProc]):
    #                     sendBuff_X[iVertex] = self.solidLoads_array_X_recon[globalIndex]
    #                     sendBuff_Y[iVertex] = self.solidLoads_array_Y_recon[globalIndex]
    #                     sendBuff_Z[iVertex] = self.solidLoads_array_Z_recon[globalIndex]
    #                     globalIndex += 1
    #                 if sendBuff_X.size > self.MPI_trhld:
    #                     iMPI = 0
    #                     while iMPI < sendBuff_X.size:
    #                         if sendBuff_X.size - iMPI > self.MPI_trhld:
    #                             self.comm.Send(sendBuff_X[iMPI:iMPI + self.MPI_trhld], dest=iProc, tag=iMPI + 1)
    #                             self.comm.Send(sendBuff_Y[iMPI:iMPI + self.MPI_trhld], dest=iProc, tag=iMPI + 2)
    #                             self.comm.Send(sendBuff_Z[iMPI:iMPI + self.MPI_trhld], dest=iProc, tag=iMPI + 3)
    #                         else:
    #                             self.comm.Send(sendBuff_X[iMPI:sendBuff_X.size], dest=iProc, tag=iMPI + 1)
    #                             self.comm.Send(sendBuff_Y[iMPI:sendBuff_X.size], dest=iProc, tag=iMPI + 2)
    #                             self.comm.Send(sendBuff_Z[iMPI:sendBuff_X.size], dest=iProc, tag=iMPI + 3)
    #                         iMPI = iMPI + self.MPI_trhld
    #                 else:
    #                     self.comm.Send(sendBuff_X, dest=iProc, tag=1)
    #                     self.comm.Send(sendBuff_Y, dest=iProc, tag=2)
    #                     self.comm.Send(sendBuff_Z, dest=iProc, tag=3)
    #         if myid in self.solidInterfaceProcessors:
    #             self.localSolidLoads_array_X = np.zeros(self.nLocalSolidInterfaceNodes)
    #             self.localSolidLoads_array_Y = np.zeros(self.nLocalSolidInterfaceNodes)
    #             self.localSolidLoads_array_Z = np.zeros(self.nLocalSolidInterfaceNodes)
    #             if self.nLocalSolidInterfaceNodes > self.MPI_trhld:  # Can't compare to the thershold in case of receiving the sendBuff_X.size as not all the structural processors may have the info of the fluid
    #                 iMPI = 0
    #                 while iMPI < self.nLocalSolidInterfaceNodes:
    #                     if self.nLocalSolidInterfaceNodes - iMPI > self.MPI_trhld:
    #                         self.comm.Recv(self.localSolidLoads_array_X[iMPI:iMPI + self.MPI_trhld],
    #                                        source=self.rootProcess, tag=iMPI + 1)
    #                         self.comm.Recv(self.localSolidLoads_array_Y[iMPI:iMPI + self.MPI_trhld],
    #                                        source=self.rootProcess, tag=iMPI + 2)
    #                         self.comm.Recv(self.localSolidLoads_array_Z[iMPI:iMPI + self.MPI_trhld],
    #                                        source=self.rootProcess, tag=iMPI + 3)
    #                     else:
    #                         self.comm.Recv(self.localSolidLoads_array_X[iMPI:self.nLocalSolidInterfaceNodes],
    #                                        source=self.rootProcess, tag=iMPI + 1)
    #                         self.comm.Recv(self.localSolidLoads_array_Y[iMPI:self.nLocalSolidInterfaceNodes],
    #                                        source=self.rootProcess, tag=iMPI + 2)
    #                         self.comm.Recv(self.localSolidLoads_array_Z[iMPI:self.nLocalSolidInterfaceNodes],
    #                                        source=self.rootProcess, tag=iMPI + 3)
    #                     iMPI = iMPI + self.MPI_trhld
    #             else:
    #                 self.comm.Recv(self.localSolidLoads_array_X, source=self.rootProcess, tag=1)
    #                 self.comm.Recv(self.localSolidLoads_array_Y, source=self.rootProcess, tag=2)
    #                 self.comm.Recv(self.localSolidLoads_array_Z, source=self.rootProcess, tag=3)
    #         del sendBuff_X
    #         del sendBuff_Y
    #         del sendBuff_Z
    #     else:
    #         self.localSolidLoads_array_X = self.solidLoads_array_X.getArray().copy()
    #         self.localSolidLoads_array_Y = self.solidLoads_array_Y.getArray().copy()
    #         self.localSolidLoads_array_Z = self.solidLoads_array_Z.getArray().copy()
    #
    # def getSolidInterfaceDisplacement(self, SolidSolver):
    #     """
    #     Gets the current solid interface position from the solid solver.
    #     """
    #     if self.have_MPI == True:
    #         myid = self.comm.Get_rank()
    #     else:
    #         myid = 0
    #
    #     # --- Get the solid interface position from the solid solver and directly fill the corresponding PETSc vector ---
    #     GlobalIndex = int()
    #     localIndex = 0
    #     for iVertex in range(self.nLocalSolidInterfaceNodes):
    #         GlobalIndex = SolidSolver.getInterfaceNodeGlobalIndex(self.solidInterfaceIdentifier, iVertex)
    #         if GlobalIndex in self.SolidHaloNodeList[myid].keys():
    #             pass
    #         else:
    #             newDispx = SolidSolver.getInterfaceNodeDispX(self.solidInterfaceIdentifier, iVertex)
    #             newDispy = SolidSolver.getInterfaceNodeDispY(self.solidInterfaceIdentifier, iVertex)
    #             newDispz = SolidSolver.getInterfaceNodeDispZ(self.solidInterfaceIdentifier, iVertex)
    #             iGlobalVertex = self.__getGlobalIndex('solid', myid, localIndex)
    #             self.solidInterface_array_DispX.setValues([iGlobalVertex], newDispx)
    #             self.solidInterface_array_DispY.setValues([iGlobalVertex], newDispy)
    #             self.solidInterface_array_DispZ.setValues([iGlobalVertex], newDispz)
    #             localIndex += 1
    #             # print("DEBUG MESSAGE From proc {} : Y Displacement : {}".format(myid, newDispy))
    #
    #     # print("DEBUG MESSAGE From proc {} : Prepare for assembly !".format(myid))
    #
    #     self.solidInterface_array_DispX.assemblyBegin()
    #     self.solidInterface_array_DispX.assemblyEnd()
    #     self.solidInterface_array_DispY.assemblyBegin()
    #     self.solidInterface_array_DispY.assemblyEnd()
    #     self.solidInterface_array_DispZ.assemblyBegin()
    #     self.solidInterface_array_DispZ.assemblyEnd()
    #
    # def getFluidInterfaceNodalForce(self, FSI_config, FluidSolver):
    #     """
    #     Gets the fluid interface loads from the fluid solver.
    #     """
    #     if self.have_MPI == True:
    #         myid = self.comm.Get_rank()
    #     else:
    #         myid = 0
    #
    #     localIndex = 0
    #     FX = 0.0
    #     FY = 0.0
    #     FZ = 0.0
    #
    #     # --- Get the fluid interface loads from the fluid solver and directly fill the corresponding PETSc vector ---
    #     for iVertex in range(self.nLocalFluidInterfaceNodes):
    #         halo = FluidSolver.ComputeVertexForces(self.fluidInterfaceIdentifier,
    #                                                iVertex)  # !!we have to ignore halo node coming from mesh partitioning because they introduice non-physical forces
    #         if halo == False:
    #             if FSI_config['CSD_SOLVER'] == 'GETDP':
    #                 newFx = FluidSolver.GetVertexForceDensityX(self.fluidInterfaceIdentifier, iVertex)
    #                 newFy = FluidSolver.GetVertexForceDensityY(self.fluidInterfaceIdentifier, iVertex)
    #                 newFz = FluidSolver.GetVertexForceDensityZ(self.fluidInterfaceIdentifier, iVertex)
    #             else:
    #                 newFx = FluidSolver.GetVertexForceX(self.fluidInterfaceIdentifier, iVertex)
    #                 newFy = FluidSolver.GetVertexForceY(self.fluidInterfaceIdentifier, iVertex)
    #                 newFz = FluidSolver.GetVertexForceZ(self.fluidInterfaceIdentifier, iVertex)
    #             iGlobalVertex = self.__getGlobalIndex('fluid', myid, localIndex)
    #             self.fluidLoads_array_X.setValues([iGlobalVertex], newFx)
    #             self.fluidLoads_array_Y.setValues([iGlobalVertex], newFy)
    #             self.fluidLoads_array_Z.setValues([iGlobalVertex], newFz)
    #             FX += newFx
    #             FY += newFy
    #             FZ += newFz
    #             localIndex += 1
    #
    #     if self.have_MPI == True:
    #         FX = self.comm.allreduce(FX)
    #         FY = self.comm.allreduce(FY)
    #         FZ = self.comm.allreduce(FZ)
    #
    #     self.fluidLoads_array_X.assemblyBegin()
    #     self.fluidLoads_array_X.assemblyEnd()
    #     self.fluidLoads_array_Y.assemblyBegin()
    #     self.fluidLoads_array_Y.assemblyEnd()
    #     self.fluidLoads_array_Z.assemblyBegin()
    #     self.fluidLoads_array_Z.assemblyEnd()
    #
    #     FX_b = self.fluidLoads_array_X.sum()
    #     FY_b = self.fluidLoads_array_Y.sum()
    #     FZ_b = self.fluidLoads_array_Z.sum()
    #
    # def setFluidInterfaceVarCoord(self, FluidSolver):
    #     """
    #     Communicate the change of coordinates of the fluid interface to the fluid solver.
    #     Prepare the fluid solver for mesh deformation.
    #     """
    #     if self.have_MPI == True:
    #         myid = self.comm.Get_rank()
    #     else:
    #         myid = 0
    #
    #     # --- Send the new fluid interface position to the fluid solver (on each partition, halo nodes included) ---
    #     localIndex = 0
    #     for iVertex in range(self.nLocalFluidInterfaceNodes):
    #         GlobalIndex = FluidSolver.GetVertexGlobalIndex(self.fluidInterfaceIdentifier, iVertex)
    #         if GlobalIndex in self.FluidHaloNodeList[myid].keys():
    #             posX0, posY0, posZ0 = self.haloNodesPositionsInit[GlobalIndex]
    #             DispX, DispY, DispZ = self.haloNodesDisplacements[GlobalIndex]
    #             # if posY0 == 0.0:
    #             #  posX = posX0
    #             #  posY = posY0
    #             #  posZ = posZ0
    #             # else:
    #             posX = posX0 + DispX
    #             posY = posY0 + DispY
    #             posZ = posZ0 + DispZ
    #             FluidSolver.SetVertexCoordX(self.fluidInterfaceIdentifier, iVertex, posX)
    #             FluidSolver.SetVertexCoordY(self.fluidInterfaceIdentifier, iVertex, posY)
    #             FluidSolver.SetVertexCoordZ(self.fluidInterfaceIdentifier, iVertex, posZ)
    #         else:
    #             # if self.localFluidInterface_array_Y_init[localIndex] == 0.0:				# !!! This is temporary and case dependent, it should be removed ASAP !!!
    #             # posX = self.localFluidInterface_array_X_init[localIndex]
    #             # posY = self.localFluidInterface_array_Y_init[localIndex]
    #             # posZ = self.localFluidInterface_array_Z_init[localIndex]
    #             # else:
    #             posX = self.localFluidInterface_array_DispX[localIndex] + self.localFluidInterface_array_X_init[localIndex]
    #             posY = self.localFluidInterface_array_DispY[localIndex] + self.localFluidInterface_array_Y_init[localIndex]
    #             posZ = self.localFluidInterface_array_DispZ[localIndex] + self.localFluidInterface_array_Z_init[localIndex]
    #             FluidSolver.SetVertexCoordX(self.fluidInterfaceIdentifier, iVertex, posX)
    #             FluidSolver.SetVertexCoordY(self.fluidInterfaceIdentifier, iVertex, posY)
    #             FluidSolver.SetVertexCoordZ(self.fluidInterfaceIdentifier, iVertex, posZ)
    #             localIndex += 1
    #         # Prepares the mesh deformation in the fluid solver
    #         nodalVarCoordNorm = FluidSolver.SetVertexVarCoord(self.fluidInterfaceIdentifier, iVertex)
    #         # print nodalVarCoordNorm
    #
    # def setSolidInterfaceLoads(self, SolidSolver, FSI_config, time):
    #     """
    #     Communicates the new solid interface loads to the solid solver.
    #     In case of rigid body motion, calculates the new resultant forces (lift, drag, ...).
    #     """
    #     if self.have_MPI == True:
    #         myid = self.comm.Get_rank()
    #     else:
    #         myid = 0
    #
    #     FY = 0.0  # solid-side resultant forces
    #     FX = 0.0
    #     FZ = 0.0
    #     FFX = 0.0  # fluid-side resultant forces
    #     FFY = 0.0
    #     FFZ = 0.0
    #
    #     # --- Check for total force conservation after interpolation
    #     FFX = self.fluidLoads_array_X.sum()
    #     FFY = self.fluidLoads_array_Y.sum()
    #     FFZ = self.fluidLoads_array_Z.sum()
    #
    #     for iVertex in range(self.nLocalSolidInterfaceNodes):
    #         FX += self.localSolidLoads_array_X[iVertex]
    #         FY += self.localSolidLoads_array_Y[iVertex]
    #         FZ += self.localSolidLoads_array_Z[iVertex]
    #
    #     if self.have_MPI == True:
    #         FX = self.comm.allreduce(FX)
    #         FY = self.comm.allreduce(FY)
    #         FZ = self.comm.allreduce(FZ)
    #
    #     self.MPIPrint("Checking f/s interface total force...")
    #     self.MPIPrint('Solid side (Fx, Fy, Fz) = ({}, {}, {})'.format(FX, FY, FZ))
    #     self.MPIPrint('Fluid side (Fx, Fy, Fz) = ({}, {}, {})'.format(FFX, FFY, FFZ))
    #
    #     # --- Send the new solid interface loads to the solid solver (on each partition, halo nodes included) ---
    #     GlobalIndex = int()
    #     localIndex = 0
    #     if myid in self.solidInterfaceProcessors:
    #         for iVertex in range(self.nLocalSolidInterfaceNodes):
    #             GlobalIndex = SolidSolver.getInterfaceNodeGlobalIndex(self.solidInterfaceIdentifier, iVertex)
    #             if GlobalIndex in self.SolidHaloNodeList[myid].keys():
    #                 pass
    #             else:
    #                 Fx = self.localSolidLoads_array_X[localIndex]
    #                 Fy = self.localSolidLoads_array_Y[localIndex]
    #                 Fz = self.localSolidLoads_array_Z[localIndex]
    #                 SolidSolver.applyload(iVertex, Fx, Fy, Fz, time)
    #                 # SolidSolver.applyload(GlobalIndex, Fx, Fy, Fz, time)
    #                 localIndex += 1
    #         if FSI_config['CSD_SOLVER'] == 'NATIVE':
    #             SolidSolver.setGeneralisedForce()
    #             SolidSolver.setGeneralisedMoment()
    #
    # def computeSolidInterfaceResidual(self, SolidSolver):
    #     """
    #     Computes the solid interface FSI displacement residual.
    #     """
    #
    #     if self.have_MPI == True:
    #         myid = self.comm.Get_rank()
    #     else:
    #         myid = 0
    #
    #     normInterfaceResidualSquare = 0.0
    #
    #     # --- Create and fill the PETSc vector for the predicted solid interface position (predicted by the solid computation) ---
    #     if self.have_MPI == True:
    #         predDisp_array_X = PETSc.Vec().create(self.comm)
    #         predDisp_array_X.setType('mpi')
    #         predDisp_array_Y = PETSc.Vec().create(self.comm)
    #         predDisp_array_Y.setType('mpi')
    #         predDisp_array_Z = PETSc.Vec().create(self.comm)
    #         predDisp_array_Z.setType('mpi')
    #     else:
    #         predDisp_array_X = PETSc.Vec().create()
    #         predDisp_array_X.setType('seq')
    #         predDisp_array_Y = PETSc.Vec().create()
    #         predDisp_array_Y.setType('seq')
    #         predDisp_array_Z = PETSc.Vec().create()
    #         predDisp_array_Z.setType('seq')
    #     predDisp_array_X.setSizes(self.nSolidInterfacePhysicalNodes + self.d_RBF)
    #     predDisp_array_Y.setSizes(self.nSolidInterfacePhysicalNodes + self.d_RBF)
    #     predDisp_array_Z.setSizes(self.nSolidInterfacePhysicalNodes + self.d_RBF)
    #
    #     if myid in self.solidSolverProcessors:
    #         for iVertex in range(self.nLocalSolidInterfaceNodes):
    #             predDispx = SolidSolver.getInterfaceNodeDispX(self.solidInterfaceIdentifier, iVertex)
    #             predDispy = SolidSolver.getInterfaceNodeDispY(self.solidInterfaceIdentifier, iVertex)
    #             predDispz = SolidSolver.getInterfaceNodeDispZ(self.solidInterfaceIdentifier, iVertex)
    #             iGlobalVertex = self.__getGlobalIndex('solid', myid, iVertex)
    #             predDisp_array_X.setValues([iGlobalVertex], predDispx)
    #             predDisp_array_Y.setValues([iGlobalVertex], predDispy)
    #             predDisp_array_Z.setValues([iGlobalVertex], predDispz)
    #
    #     predDisp_array_X.assemblyBegin()
    #     predDisp_array_X.assemblyEnd()
    #     predDisp_array_Y.assemblyBegin()
    #     predDisp_array_Y.assemblyEnd()
    #     predDisp_array_Z.assemblyBegin()
    #     predDisp_array_Z.assemblyEnd()
    #
    #     # --- Calculate the residual (vector and norm) ---
    #     self.solidInterfaceResidual_array_X = predDisp_array_X - self.solidInterface_array_DispX
    #     self.solidInterfaceResidual_array_Y = predDisp_array_Y - self.solidInterface_array_DispY
    #     self.solidInterfaceResidual_array_Z = predDisp_array_Z - self.solidInterface_array_DispZ
    #
    #     normInterfaceResidual_X = self.solidInterfaceResidual_array_X.norm()
    #     normInterfaceResidual_Y = self.solidInterfaceResidual_array_Y.norm()
    #     normInterfaceResidual_Z = self.solidInterfaceResidual_array_Z.norm()
    #
    #     # print("DEBUG MESSAGE FROM PROC {} : NormInterfaceResidual_X = {}".format(myid, normInterfaceResidual_X))
    #     # print("DEBUG MESSAGE FROM PROC {} : NormInterfaceResidual_Y = {}".format(myid, normInterfaceResidual_Y))
    #     # print("DEBUG MESSAGE FROM PROC {} : NormInterfaceResidual_Z = {}".format(myid, normInterfaceResidual_Z))
    #     # if myid == 0:
    #     #  print predPos_array_DispX.getArray()
    #     #  print "*************"
    #     #  print self.solidInterface_array_DispX.getArray()
    #
    #     normInterfaceResidualSquare = normInterfaceResidual_X ** 2 + normInterfaceResidual_Y ** 2 + normInterfaceResidual_Z ** 2
    #
    #     predDisp_array_X.destroy()
    #     predDisp_array_Y.destroy()
    #     predDisp_array_Z.destroy()
    #     del predDisp_array_X
    #     del predDisp_array_Y
    #     del predDisp_array_Z
    #
    #     return sqrt(normInterfaceResidualSquare)
    #
    # def relaxSolidPosition(self, FSI_config):
    #     """
    #     Apply solid displacement under-relaxation.
    #     """
    #     if self.have_MPI == True:
    #         myid = self.comm.Get_rank()
    #     else:
    #         myid = 0
    #
    #     # --- Set the Aitken coefficient for the relaxation ---
    #     if FSI_config['AITKEN_RELAX'] == 'STATIC':
    #         self.aitkenParam = FSI_config['AITKEN_PARAM']
    #     elif FSI_config['AITKEN_RELAX'] == 'DYNAMIC':
    #         self.setAitkenCoefficient(FSI_config)
    #     else:
    #         self.aitkenParam = 1.0
    #
    #     self.MPIPrint('Aitken under-relaxation step with parameter {}'.format(self.aitkenParam))
    #
    #     # --- Relax the solid interface position ---
    #     self.solidInterface_array_DispX += self.aitkenParam * self.solidInterfaceResidual_array_X
    #     self.solidInterface_array_DispY += self.aitkenParam * self.solidInterfaceResidual_array_Y
    #     self.solidInterface_array_DispZ += self.aitkenParam * self.solidInterfaceResidual_array_Z
    #
    #
    # def writeFSIHistory(self, TimeIter, time, varCoordNorm, FSIConv):
    #     """
    #     Write the FSI history file of the computaion.
    #     """
    #
    #     if self.have_MPI == True:
    #         myid = self.comm.Get_rank()
    #     else:
    #         myid = 0
    #
    #     if myid == self.rootProcess:
    #         if self.unsteady:
    #             if TimeIter == 0:
    #                 histFile = open('FSIhistory.dat', "w")
    #                 histFile.write("TimeIter\tTime\tFSIRes\tFSINbIter\n")
    #             else:
    #                 histFile = open('FSIhistory.dat', "a")
    #             if FSIConv:
    #                 histFile.write(str(TimeIter) + '\t' + str(time) + '\t' + str(varCoordNorm) + '\t' + str(
    #                     self.FSIIter + 1) + '\n')
    #             else:
    #                 histFile.write(
    #                     str(TimeIter) + '\t' + str(time) + '\t' + str(varCoordNorm) + '\t' + str(self.FSIIter) + '\n')
    #             histFile.close()
    #         else:
    #             if self.FSIIter == 0:
    #                 histFile = open('FSIhistory.dat', "w")
    #                 histFile.write("FSI Iter\tFSIRes\n")
    #             else:
    #                 histFile = open('FSIhistory.dat', "a")
    #             histFile.write(str(self.FSIIter) + '\t' + str(varCoordNorm) + '\n')
    #             histFile.close()
    #
    #     self.MPIBarrier()
    #
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
    #
    # def SteadyFSI(self, FSI_config, FluidSolver, SolidSolver, MLS_Spline):
    #     """
    #     Runs the steady FSI computation by synchronizing the fluid and solid solver with data exchange at the f/s interface.
    #     """
    #
    #     if self.have_MPI == True:
    #         myid = self.comm.Get_rank()
    #         numberPart = self.comm.Get_size()
    #     else:
    #         myid = 0
    #         numberPart = 1
    #
    #     # --- Set some general variables for the steady computation --- #
    #     NbIter = FSI_config['NB_EXT_ITER']  # number of fluid iteration at each FSI step
    #     NbFSIIterMax = FSI_config['NB_FSI_ITER']  # maximum number of FSI iteration (for each time step)
    #     FSITolerance = FSI_config['FSI_TOLERANCE']  # f/s interface tolerance
    #     varCoordNorm = 0.0
    #
    #     # --- Initialize matrix of boundary nodal forces  --- #
    #     if myid in self.solidSolverProcessors:
    #         SolidSolver.EvaluateIntefaceFluidDisplacements(FSI_config,
    #                                                        MLS_Spline)  # Flutter_mode_fluid_x/y/z are stored (root) once and for all
    #         SolidSolver.initialize_OutputForces(1, FSI_config)
    #
    #     self.MPIPrint('\n********************************')
    #     self.MPIPrint('* Begin steady FSI computation *')
    #     self.MPIPrint('********************************\n')
    #     self.MPIPrint('\n*************** Enter Block Gauss Seidel (BGS) method for strong coupling FSI ***************')
    #
    #     self.getSolidInterfaceDisplacement(SolidSolver)
    #
    #     # --- External FSI loop --- #
    #     self.FSIIter = 0
    #     while self.FSIIter < NbFSIIterMax:
    #         self.MPIPrint("\n>>>> FSI iteration {} <<<<".format(self.FSIIter))
    #         self.MPIPrint('\nLaunching fluid solver for a steady computation...')
    #         # --- Fluid solver call for FSI subiteration ---#
    #         Iter = 0
    #         FluidSolver.ResetConvergence()
    #         while Iter < NbIter:
    #             FluidSolver.PreprocessExtIter(Iter)
    #             FluidSolver.Run()
    #             StopIntegration = FluidSolver.Monitor(Iter)
    #             FluidSolver.Output(Iter)
    #             if StopIntegration:
    #                 break;
    #             Iter += 1
    #
    #         # --- Surface fluid loads interpolation and communication ---#
    #         self.MPIPrint('\nProcessing interface fluid loads...\n')
    #         self.MPIBarrier()
    #         self.getFluidInterfaceNodalForce(FSI_config, FluidSolver)
    #         self.MPIBarrier()
    #         self.interpolateFluidLoadsOnSolidMesh(FSI_config)
    #         self.setSolidInterfaceLoads(SolidSolver, FSI_config, 0.05)
    #
    #         # --- Solid solver call for FSI subiteration --- #
    #         self.MPIPrint('\nLaunching solid solver for a static computation...\n')
    #         if myid in self.solidSolverProcessors:
    #             if FSI_config['CSD_SOLVER'] == 'NATIVE':
    #                 SolidSolver.staticComputation()
    #             else:
    #                 # SolidSolver.run(0.0, 0.05)  OLD ONE
    #                 SolidSolver.run(0.0, FSI_config, MLS_Spline)
    #             # if FSI_config['WRITE_FORCE_OUTPUT']   == 'YES':
    #             SolidSolver.writeSolution(0, 0, FSI_config)
    #             # print(SolidSolver.NodalForces)
    #         # --- Compute and monitor the FSI residual --- #
    #         varCoordNorm = self.computeSolidInterfaceResidual(SolidSolver)
    #         self.MPIPrint('\nFSI displacement norm : {}\n'.format(varCoordNorm))
    #         self.writeFSIHistory(0, 0.0, varCoordNorm, False)
    #         ## remove spurious surface flow file
    #         # the location is only in memory of the structural CPU node!!
    #         if myid in self.solidSolverProcessors:
    #             command_remove = "rm -r " + FSI_config['SURFACE_FLOW_FILENAME'] + "*"
    #             os.system(command_remove)
    #         if varCoordNorm < FSITolerance:
    #             break
    #
    #         # --- Relaxe the solid displacement and update the solid solution --- #
    #         self.MPIPrint('\nProcessing interface displacements...\n')
    #         self.relaxSolidPosition(FSI_config)
    #         if myid in self.solidSolverProcessors:
    #             SolidSolver.updateSolution()
    #
    #         # --- Mesh morphing step (displacement interpolation, displacements communication, and mesh morpher call) --- #
    #         self.interpolateSolidPositionOnFluidMesh(FSI_config)
    #         self.MPIPrint('\nPerforming static mesh deformation...\n')
    #         self.setFluidInterfaceVarCoord(FluidSolver)
    #         FluidSolver.StaticMeshUpdate()
    #         self.FSIIter += 1
    #
    #     self.MPIBarrier()
    #
    #     self.MPIPrint('\nBGS is converged (strong coupling)')
    #     self.MPIPrint(' ')
    #     self.MPIPrint('*************************')
    #     self.MPIPrint('*  End FSI computation  *')
    #     self.MPIPrint('*************************')
    #     self.MPIPrint(' ')
