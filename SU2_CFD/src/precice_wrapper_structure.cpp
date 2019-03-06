/*!
 * \file precice_wrapper_structure.cpp
 * \brief Subroutines to run precice with SU2, based on the precice adaptor: https://github.com/precice/su2-adapter.git
 * \author R. Sanchez, based on the work by Alexander Rusch and the preCICE community
 * \version 6.1.0 "Falcon"
 */

#include "../include/precice_structure.hpp"

/*!
* \brief Constructor of the class CPrecice.
* \param[in] processRank - Index of this solver Process.
* \param[in] processSize - Overall number of solver Processes.
*/
CPrecice::CPrecice( int processRank, int processSize, CGeometry**** geometry_container, CSolver***** solver_container,
                    CConfig** config_container, CVolumetricMovement*** grid_movement)
:
processRank(processRank),
processSize(processSize),
solverInterface(config_container[ZONE_0]->GetpreCICE_Participant(), processRank, processSize),
coric(precice::constants::actionReadIterationCheckpoint()),
cowic(precice::constants::actionWriteIterationCheckpoint())
{
  /*--- Store the pointers of the problem ---*/
  geometry = geometry_container[ZONE_0][INST_0][MESH_0];
  solver = solver_container[ZONE_0][INST_0][MESH_0];
  config = config_container[ZONE_0];

  /*--- Store the dimensionality ---*/
  nDim = geometry->GetnDim();
  nPoint = geometry->GetnPoint();

  /*--- Initialize quantities ---*/
  activeProcess = true;
  verbose = config->GetpreCICE_VerbosityLevel_High();
  nWetSurfaces = config->GetpreCICE_NumberWetSurfaces();
  nWetSurfacesDomain = 0;


  /*--- There must be at least one wet surface ---*/
  if(nWetSurfaces < 1)
    SU2_MPI::Error("There must be at least one wet surface! Now exiting...", CURRENT_FUNCTION);

};

/*!
* \brief Destructor of the class CPrecice.
*/
CPrecice::~CPrecice() {

};

CPreciceFlow::CPreciceFlow( int processRank, int processSize, CGeometry**** geometry_container, CSolver***** solver_container,
                            CConfig** config_container, CVolumetricMovement*** grid_movement )
                           : CPrecice(processRank, processSize, geometry_container, solver_container, config_container, grid_movement)
{

  unsigned long iPoint;
  unsigned short iDim;

  /*--- Number of variables ---*/
  nVar = solver[FLOW_SOL]->GetnVar();

  dt_savedState = 0.0;
  StopCalc_savedState = false;

  /*--- Initialize structures ---*/
  vertexIDs = NULL;
  forceID = NULL;
  displDeltaID = NULL;
  meshID = NULL;

  Coord_Saved = new su2double*[nPoint];
  Coord_n_Saved = new su2double*[nPoint];
  Coord_n1_Saved = new su2double*[nPoint];
  Coord_p1_Saved = new su2double*[nPoint];
  GridVel_Saved = new su2double*[nPoint];
  GridVel_Grad_Saved = new su2double**[nPoint];
  solution_Saved = new su2double*[nPoint];
  solution_time_n_Saved = new su2double*[nPoint];
  solution_time_n1_Saved = new su2double*[nPoint];
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    Coord_Saved[iPoint] = new su2double[nDim];
    Coord_n_Saved[iPoint] = new su2double[nDim];
    Coord_n1_Saved[iPoint] = new su2double[nDim];
    Coord_p1_Saved[iPoint] = new su2double[nDim];
    GridVel_Saved[iPoint] = new su2double[nDim];
    GridVel_Grad_Saved[iPoint] = new su2double*[nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      GridVel_Grad_Saved[iPoint][iDim] = new su2double[nDim];
    solution_Saved[iPoint] = new su2double[nVar];
    solution_time_n_Saved[iPoint] = new su2double[nVar];
    solution_time_n1_Saved[iPoint] = new su2double[nVar];
  }

  /*--- Initialize magnitudes ---*/
  Velocity_Real = NULL;
  Velocity_ND = NULL;
  Density_Real = 0.0;
  Density_ND = 0.0;
  Velocity2_Real = 0.0;
  Velocity2_ND = 0.0;

  /*--- Configure preCICE with its configuration file ---*/
  Configure(config->GetpreCICE_ConfigFileName());

};

CPreciceFlow::~CPreciceFlow()
{

  unsigned long iPoint;
  unsigned short iDim;

  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    if (Coord_Saved[iPoint]  != NULL) delete [] Coord_Saved[iPoint];
    if (Coord_n_Saved[iPoint]  != NULL) delete [] Coord_n_Saved[iPoint];
    if (Coord_n1_Saved[iPoint] != NULL) delete [] Coord_n1_Saved[iPoint];
    if (Coord_p1_Saved[iPoint] != NULL) delete [] Coord_p1_Saved[iPoint];
    if (GridVel_Saved[iPoint]  != NULL) delete [] GridVel_Saved[iPoint];

    for (int iDim = 0; iDim < nDim; iDim++) {
      if (GridVel_Grad_Saved[iPoint][iDim] != NULL) delete [] GridVel_Grad_Saved[iPoint][iDim];
    }

    if (GridVel_Grad_Saved[iPoint]     != NULL) delete [] GridVel_Grad_Saved[iPoint];
    if (solution_Saved[iPoint]         != NULL) delete [] solution_Saved[iPoint];
    if (solution_time_n_Saved[iPoint]  != NULL) delete [] solution_time_n_Saved[iPoint];
    if (solution_time_n1_Saved[iPoint] != NULL) delete [] solution_time_n1_Saved[iPoint];
  }

  if (Coord_Saved            != NULL) delete [] Coord_Saved;
  if (Coord_n_Saved          != NULL) delete [] Coord_n_Saved;
  if (Coord_n1_Saved         != NULL) delete [] Coord_n1_Saved;
  if (Coord_p1_Saved         != NULL) delete [] Coord_p1_Saved;
  if (GridVel_Saved          != NULL) delete [] GridVel_Saved;
  if (GridVel_Grad_Saved     != NULL) delete [] GridVel_Grad_Saved;
  if (solution_Saved         != NULL) delete [] solution_Saved;
  if (solution_time_n_Saved  != NULL) delete [] solution_time_n_Saved;
  if (solution_time_n1_Saved != NULL) delete [] solution_time_n1_Saved;

  if (vertexIDs              != NULL) delete [] vertexIDs;
  if (forceID                != NULL) delete [] forceID;
  if (displDeltaID           != NULL) delete [] displDeltaID;
  if (meshID                 != NULL) delete [] meshID;

};

su2double CPreciceFlow::Initialize() {

  unsigned long iSurface, jSurface;
  unsigned long iPoint;
  unsigned short iDim;
  int iVertex;
  su2double precice_dt;

  /*--- Allocate structures to interchange data with preCICE ---*/
  meshID = new int[nWetSurfaces];
  forceID = new int[nWetSurfaces];
  displDeltaID = new int[nWetSurfaces];

  /*--- Recover preCICE mesh IDs ---*/
  for (iSurface = 0; iSurface < nWetSurfaces; iSurface++) {
    meshID[iSurface] = solverInterface.getMeshID("SU2_Mesh_Fluid" + to_string(iSurface));
  }

  /*--- Store the number of surfaces belonging to this rank ---*/
  for (iSurface = 0; iSurface < nWetSurfaces; iSurface++) {
    if (config->GetMarker_All_TagBound(config->GetpreCICE_WetSurfaceMarkerName() + to_string(iSurface)) != -1) {
      nWetSurfacesDomain++;
    }
  }

  /*--- Disable the processes in this rank if there are no wet surfaces that belong to it ---*/
  if (nWetSurfacesDomain < 1)  activeProcess = false;

  /*--- When this rank is active (that is, it holds a piece of the wet boundary) ---*/
  if (activeProcess) {

    /*--- Store the wet surface marker IDs in an array ---*/
    valueMarkerWet = new short[nWetSurfacesDomain];
    markerLocalToGlobal = new unsigned long[nWetSurfacesDomain];
    jSurface = 0;
    for (iSurface = 0; iSurface < nWetSurfaces; iSurface++) {
      if (config->GetMarker_All_TagBound(config->GetpreCICE_WetSurfaceMarkerName() + to_string(iSurface)) != -1) {
        valueMarkerWet[jSurface] = config->GetMarker_All_TagBound(config->GetpreCICE_WetSurfaceMarkerName() + to_string(iSurface));
        markerLocalToGlobal[jSurface] = iSurface;
        jSurface++;
      }
    }

    vertexIDs = new int*[nWetSurfacesDomain];
    nVertex = new int[nWetSurfacesDomain];

    for (iSurface = 0; iSurface < nWetSurfacesDomain; iSurface++) {

      nVertex[iSurface] = geometry->nVertex[valueMarkerWet[iSurface]];

      /*--- Allocate memory to store all the node coordinates in the surface ---*/
      su2double **nodeCoord;
      nodeCoord = new su2double*[nVertex[iSurface]];
      for (iVertex = 0; iVertex < nVertex[iSurface]; iVertex++)
        nodeCoord[iVertex] = new su2double[nDim];

      /*--- Loop over the nodes of each surface ---*/
      for (iVertex = 0; iVertex < nVertex[iSurface]; iVertex++) {

        /*--- Store the point index for the vertex nodes ---*/
        iPoint = geometry->vertex[valueMarkerWet[iSurface]][iVertex]->GetNode();

        /*--- Retrieve the coordinates ---*/
        for (iDim = 0; iDim < nDim; iDim++) {
          nodeCoord[iVertex][iDim] = geometry->node[iPoint]->GetCoord(iDim);
          if (verbose) cout << "Process #" << processRank << "/" << processSize-1 << ": Initial coordinates of node (local index, global index, node color): (" << iVertex << ", " << iPoint << ", " << geometry->node[iPoint]->GetColor() << "): " << nodeCoord[iVertex][iDim] << endl; /*--- for debugging purposes ---*/
        }
      }

      /*--- Passive double structure to pass on the vector of coordinates into preCICE ---*/
      passivedouble *coord;
      coord = new passivedouble[nVertex[iSurface]*nDim];

      /*--- Convert the coordinates to a passive double structure that will be communicated to preCICE ---*/
      for (iVertex = 0; iVertex < nVertex[iSurface]; iVertex++) {
        for (iDim = 0; iDim < nDim; iDim++) {
          coord[iVertex*nDim + iDim] = SU2_TYPE::GetValue(nodeCoord[iVertex][iDim]);
        }
      }

      /*--- Variable internal to preCICE ---*/
      vertexIDs[iSurface] = new int[nVertex[iSurface]];

      /*--- Interface to preCICE ---*/
      solverInterface.setMeshVertices(meshID[markerLocalToGlobal[iSurface]], nVertex[iSurface], coord, vertexIDs[iSurface]);
      forceID[markerLocalToGlobal[iSurface]] = solverInterface.getDataID("Forces" + to_string(markerLocalToGlobal[iSurface]), meshID[markerLocalToGlobal[iSurface]]);
      displDeltaID[markerLocalToGlobal[iSurface]] = solverInterface.getDataID("DisplacementDeltas" + to_string(markerLocalToGlobal[iSurface]), meshID[markerLocalToGlobal[iSurface]]);

      /*--- Deallocate the containers of the coordinates ---*/
      if (coord != NULL) delete [] coord;
      for (iVertex = 0; iVertex < nVertex[iSurface]; iVertex++){
        if (nodeCoord[iVertex]!= NULL) delete [] nodeCoord[iVertex];
      }
      if (nodeCoord!= NULL) delete [] nodeCoord;

    }
    for (iSurface = 0; iSurface < nWetSurfaces; iSurface++) {
      for (jSurface = 0; jSurface < nWetSurfacesDomain; jSurface++) {
        if (markerLocalToGlobal[jSurface] != iSurface) {
          solverInterface.setMeshVertices(meshID[iSurface], 0, NULL, NULL);
          forceID[iSurface] = solverInterface.getDataID("Forces" + to_string(iSurface), meshID[iSurface]);
          displDeltaID[iSurface] = solverInterface.getDataID("DisplacementDeltas" + to_string(iSurface), meshID[iSurface]);
        }
      }
    }
  } else {
    /*--- If the surface is not active ---*/
    for (iSurface = 0; iSurface < nWetSurfaces; iSurface++) {
      solverInterface.setMeshVertices(meshID[iSurface], 0, NULL, NULL);
      forceID[iSurface] = solverInterface.getDataID("Forces" + to_string(iSurface), meshID[iSurface]);
      displDeltaID[iSurface] = solverInterface.getDataID("DisplacementDeltas" + to_string(iSurface), meshID[iSurface]);
    }
  }

  /*--- Finally, initiallize preCICE with the information that we have sent to it ---*/
  precice_dt = solverInterface.initialize();

  /*--- readData ---*/
  if (solverInterface.isReadDataAvailable()) {
  /*--- Retrieve the displacements from preCICE for the wet surfaces that belong to the active process---*/
    for (iSurface = 0; iSurface < nWetSurfacesDomain; iSurface++) {

      /*--- Allocate memory to store all the displacements in the wet surface ---*/
      su2double **displacementDeltas_su2;
      displacementDeltas_su2 = new su2double*[nVertex[iSurface]];
      for (iVertex = 0; iVertex < nVertex[iSurface]; iVertex++)
        displacementDeltas_su2[iVertex] = new su2double[nDim];

      /*--- Allocate memory to retrieve all the displacements in the wet surface coming from preCICE ---*/
      passivedouble *displacementDeltas;
      displacementDeltas = new passivedouble[nVertex[iSurface]*nDim];

      /*--- Read the displacement data from preCICE ---*/
      solverInterface.readBlockVectorData(displDeltaID[markerLocalToGlobal[iSurface]], nVertex[iSurface], vertexIDs[iSurface], displacementDeltas);

      /*--- Store the displacement coming from preCICE in the container for SU2 ---*/
      for (iVertex = 0; iVertex < nVertex[iSurface]; iVertex++) {
        for (iDim = 0; iDim < nDim; iDim++) {
          displacementDeltas_su2[iVertex][iDim] = displacementDeltas[iVertex*nDim + iDim];
//if (geometry->node[geometry->vertex[valueMarkerWet[iSurface]][iVertex]->GetNode()]->GetGlobalIndex() == 300){
//cout << "DispInFlow 300  " <<  geometry->node[geometry->vertex[valueMarkerWet[iSurface]][iVertex]->GetNode()]->GetGlobalIndex()  << "   " << displacementDeltas[iVertex*nDim + iDim] << endl;
//}
	}
      }


      /*--- Set the incremental displacement to the mesh solver ---*/
      for (iVertex = 0; iVertex < nVertex[iSurface]; iVertex++) {
        geometry->vertex[valueMarkerWet[iSurface]][iVertex]->SetVarCoord(displacementDeltas_su2[iVertex]);
      }

      /*--- Deallocate the containers of the delta-displacement ---*/
      if (displacementDeltas != NULL)  delete [] displacementDeltas;
      for (iVertex = 0; iVertex < nVertex[iSurface]; iVertex++){
        if (displacementDeltas_su2[iVertex]!= NULL) delete [] displacementDeltas_su2[iVertex];
      }
      if (displacementDeltas_su2!= NULL) delete [] displacementDeltas_su2;

    }
  }

  return precice_dt;

};

su2double CPreciceFlow::Advance( su2double computedTimestep ) {

  unsigned long iSurface;
  unsigned short iDim, jDim;
  int iVertex;

  su2double max_precice_dt;

  if (activeProcess) {

    /*--- Determine if the flow is viscous ---*/
    bool viscous_flow = ((config->GetKind_Solver() == NAVIER_STOKES) || (config->GetKind_Solver() == RANS));

    /*--- Compute the redimensionalization parameter to compute the forces in SI ---*/
    Velocity_Real = config->GetVelocity_FreeStream();
    Density_Real = config->GetDensity_FreeStream();

    Velocity_ND = config->GetVelocity_FreeStreamND();
    Density_ND = config->GetDensity_FreeStreamND();

    Velocity2_Real = 0.0;
    Velocity2_ND = 0.0;
    for (iDim = 0; iDim < nDim; iDim++){
      Velocity2_Real += Velocity_Real[iDim]*Velocity_Real[iDim];
      Velocity2_ND += Velocity_ND[iDim]*Velocity_ND[iDim];
    }
    factorForces = Density_Real*Velocity2_Real/(Density_ND*Velocity2_ND);

    /*--- Compute the forces for the wet surfaces that belong to the active process ---*/
    for (iSurface = 0; iSurface < nWetSurfacesDomain; iSurface++) {

      /*--- Allocate memory to store all the flow forces in the surface ---*/
      su2double **forces_su2;
      forces_su2 = new su2double*[nVertex[iSurface]];
      for (iVertex = 0; iVertex < nVertex[iSurface]; iVertex++)
        forces_su2[iVertex] = new su2double[nDim];

      /*--- Allocate a passive double structure to store the forces as they are transferred to preCICE ---*/
      passivedouble *forces;
      forces = new passivedouble[nVertex[iSurface]*nDim];

      /*--- Allocate memory to store the intermediate forces for load calculation ---*/
      unsigned long Point_Flow;
      su2double Pn = 0.0, Pinf = 0.0, div_vel = 0.0, Viscosity = 0.0;
      su2double Tau[3][3] = {{0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0}};
      su2double Grad_Vel[3][3] = {{0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0}};
      su2double delta[3][3] = {{1.0, 0.0, 0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};
      su2double *Normal_Flow;

      /*--- Loop over vertices of coupled boundary ---*/
      for (iVertex = 0; iVertex < nVertex[iSurface]; iVertex++) {

        /*--- Get the vertex index and the normal at the vertex ---*/
        Point_Flow = geometry->vertex[valueMarkerWet[iSurface]][iVertex]->GetNode();
        Normal_Flow = geometry->vertex[valueMarkerWet[iSurface]][iVertex]->GetNormal();

        /*--- Retrieve the values of pressure and viscosity ---*/
        Pn = solver[FLOW_SOL]->node[Point_Flow]->GetPressure();
        Pinf = solver[FLOW_SOL]->GetPressure_Inf();

        /*--- Compute the inviscid term of the traction ---*/
        for (iDim = 0; iDim < nDim; iDim++)
          forces_su2[iVertex][iDim] = -(Pn-Pinf)*Normal_Flow[iDim];

        /*--- Compute the viscous term of the traction at the surface ---*/
        if (viscous_flow) {

          /*--- Read the viscosity ---*/
          Viscosity =solver[FLOW_SOL]->node[Point_Flow]->GetLaminarViscosity();

          /*--- Read the velocity gradient ---*/
          for (iDim = 0; iDim < nDim; iDim++) {
            for (jDim = 0 ; jDim < nDim; jDim++) {
              Grad_Vel[iDim][jDim] = solver[FLOW_SOL]->node[Point_Flow]->GetGradient_Primitive(iDim+1, jDim);
            }
          }

          /*--- Compute the divergence of the velocity ---*/
          div_vel = 0.0; for (iDim = 0; iDim < nDim; iDim++) div_vel += Grad_Vel[iDim][iDim];

          /*--- Add the viscous component to the traction ---*/
          for (iDim = 0; iDim < nDim; iDim++) {
            for (jDim = 0 ; jDim < nDim; jDim++) {
              Tau[iDim][jDim] = Viscosity*(Grad_Vel[jDim][iDim] + Grad_Vel[iDim][jDim]) - TWO3*Viscosity*div_vel*delta[iDim][jDim];
              forces_su2[iVertex][iDim] += Tau[iDim][jDim]*Normal_Flow[jDim];
            }
          }
        }

        /*--- Redimensionalize forces to SI units ---*/
        for (iDim = 0; iDim < nDim; iDim++) {
          forces_su2[iVertex][iDim] = forces_su2[iVertex][iDim]*factorForces;
        }

      }

      /*--- Store the forces in a passive double structure, for transfer into preCICE ---*/
      for (iVertex = 0; iVertex < nVertex[iSurface]; iVertex++) {
        for (iDim = 0; iDim < nDim; iDim++) {
          if (geometry->node[Point_Flow]->GetColor() == processRank)          /* TODO PointFlow defined??*/
            forces[iVertex*nDim + iDim] = SU2_TYPE::GetValue(forces_su2[iVertex][iDim]);
          else
            forces[iVertex*nDim + iDim] = 0.0;
//cout << "ForcesOutFlow   " <<  geometry->node[geometry->vertex[valueMarkerWet[iSurface]][iVertex]->GetNode()]->GetGlobalIndex()  << "   " << forces[iVertex*nDim + iDim] << endl;
//if (geometry->node[geometry->vertex[valueMarkerWet[iSurface]][iVertex]->GetNode()]->GetGlobalIndex() == 300){
//cout << "ForcesOutFlow 300  " <<  geometry->node[geometry->vertex[valueMarkerWet[iSurface]][iVertex]->GetNode()]->GetGlobalIndex()  << "   " << forces[iVertex*nDim + iDim] << endl;
//         }
        }
      }

      /*--- Apply a linear ramp to the forces ---*/
      if (config->GetpreCICE_LoadRamping() && ((config->GetExtIter() - config->GetUnst_RestartIter()) < config->GetpreCICE_LoadRampingDuration())) {
        *forces = *forces * ((config->GetExtIter() - config->GetUnst_RestartIter()) + 1) / config->GetpreCICE_LoadRampingDuration();
      }

      /*--- Send the data to preCICE ---*/
      solverInterface.writeBlockVectorData(forceID[markerLocalToGlobal[iSurface]], nVertex[iSurface], vertexIDs[iSurface], forces);

      /*--- Deallocate the containers of the forces ---*/
      if (forces != NULL) delete [] forces;
      for (iVertex = 0; iVertex < nVertex[iSurface]; iVertex++){
        if (forces_su2[iVertex]!= NULL) delete [] forces_su2[iVertex];
      }
      if (forces_su2!= NULL) delete [] forces_su2;

    }

    /*--- Advance the interface in preCICE ---*/
    max_precice_dt = solverInterface.advance( SU2_TYPE::GetValue(computedTimestep) );

    /*--- Retrieve the displacements from preCICE for the wet surfaces that belong to the active process---*/
    for (iSurface = 0; iSurface < nWetSurfacesDomain; iSurface++) {

      /*--- Allocate memory to store all the displacements in the wet surface ---*/
      su2double **displacementDeltas_su2;
      displacementDeltas_su2 = new su2double*[nVertex[iSurface]];
      for (iVertex = 0; iVertex < nVertex[iSurface]; iVertex++)
        displacementDeltas_su2[iVertex] = new su2double[nDim];

      /*--- Allocate memory to retrieve all the displacements in the wet surface coming from preCICE ---*/
      passivedouble *displacementDeltas;
      displacementDeltas = new passivedouble[nVertex[iSurface]*nDim];

      /*--- Read the displacement data from preCICE ---*/
      solverInterface.readBlockVectorData(displDeltaID[markerLocalToGlobal[iSurface]], nVertex[iSurface], vertexIDs[iSurface], displacementDeltas);

      /*--- Store the displacement coming from preCICE in the container for SU2 ---*/
      for (iVertex = 0; iVertex < nVertex[iSurface]; iVertex++) {
        for (iDim = 0; iDim < nDim; iDim++) {
          displacementDeltas_su2[iVertex][iDim] = displacementDeltas[iVertex*nDim + iDim];
//if (geometry->node[geometry->vertex[valueMarkerWet[iSurface]][iVertex]->GetNode()]->GetGlobalIndex() == 300){
//cout << "DispInFlow 300   " <<  geometry->node[geometry->vertex[valueMarkerWet[iSurface]][iVertex]->GetNode()]->GetGlobalIndex()  << "   " << displacementDeltas[iVertex*nDim + iDim] << endl;
//}
	}
      }

      /*--- Set the incremental displacement to the mesh solver ---*/
      for (iVertex = 0; iVertex < nVertex[iSurface]; iVertex++) {
        geometry->vertex[valueMarkerWet[iSurface]][iVertex]->SetVarCoord(displacementDeltas_su2[iVertex]);
      }
      /*--- Deallocate the containers of the delta-displacement ---*/
      if (displacementDeltas != NULL)  delete [] displacementDeltas;
      for (iVertex = 0; iVertex < nVertex[iSurface]; iVertex++){
        if (displacementDeltas_su2[iVertex]!= NULL) delete [] displacementDeltas_su2[iVertex];
      }
      if (displacementDeltas_su2!= NULL) delete [] displacementDeltas_su2;

    }

    return max_precice_dt;
  }
  else {
    /*--- Advance the interface (as nothing is transferred from this process) ---*/
    max_precice_dt = solverInterface.advance( SU2_TYPE::GetValue(computedTimestep) );
    return max_precice_dt;
  }

};

void CPreciceFlow::Set_OldState( bool *StopCalc, double *dt ) {

  unsigned long iPoint;
  unsigned short iVar, iDim, jDim;

  for (iPoint = 0; iPoint < nPoint; iPoint++) {
if (geometry->node[iPoint]->GetGlobalIndex()==11){
/*cout << "sol node 300 x before set  " << solver[FLOW_SOL]->node[iPoint]->GetSolution(0)<< "   " << solver[FLOW_SOL]->node[iPoint]->GetSolution(1) <<endl;
cout << "soltn node 300 x before set  " << solver[FLOW_SOL]->node[iPoint]->GetSolution_time_n(0)<< "   " << solver[FLOW_SOL]->node[iPoint]->GetSolution_time_n(1) <<endl;
cout << "soltn1 node 300 x before set  " << (solver[FLOW_SOL]->node[iPoint]->GetSolution_time_n1())[0] << "   " <<(solver[FLOW_SOL]->node[iPoint]->GetSolution_time_n1())[1] <<endl;*/
//cout << "coord node 11 before set  " << geometry->node[iPoint]->GetCoord(0)<< "   " << geometry->node[iPoint]->GetCoord(1) <<endl;
/*cout << "coordn node 300 x before set  " << geometry->node[iPoint]->GetCoord_n(0)<< "   " << geometry->node[iPoint]->GetCoord_n(1) <<endl;
cout << "coordn1 node 300 x before set  " << geometry->node[iPoint]->GetCoord_n1(0)<< "   " << geometry->node[iPoint]->GetCoord_n1(1) <<endl;
cout << "coordp1 node 300 x before set  " << geometry->node[iPoint]->GetCoord_p1()<< "   " << geometry->node[iPoint]->GetCoord_p1() <<endl;
cout << "grid_vel node 300 x before set  " << geometry->node[iPoint]->GetGridVel()<< "   " << geometry->node[iPoint]->GetGridVel() <<endl;*/
}
    for (iVar = 0; iVar < nVar; iVar++) {
      /*--- Save solutions at last and current time step ---*/
      solution_Saved[iPoint][iVar] = (solver[FLOW_SOL]->node[iPoint]->GetSolution())[iVar];
      if (config->GetUnsteady_Simulation() != NO){
        solution_time_n_Saved[iPoint][iVar] = (solver[FLOW_SOL]->node[iPoint]->GetSolution_time_n())[iVar];
        solution_time_n1_Saved[iPoint][iVar] = (solver[FLOW_SOL]->node[iPoint]->GetSolution_time_n1())[iVar];
      } 
    }
      /*-- Save coordinates at last, current and next time step ---*/
    for (iDim = 0; iDim < nDim; iDim++) {
      Coord_Saved[iPoint][iDim] = (geometry->node[iPoint]->GetCoord())[iDim];
      if (config->GetUnsteady_Simulation() != NO){
        Coord_n_Saved[iPoint][iDim] = (geometry->node[iPoint]->GetCoord_n())[iDim];
        Coord_n1_Saved[iPoint][iDim] = (geometry->node[iPoint]->GetCoord_n1())[iDim];
        Coord_p1_Saved[iPoint][iDim] = (geometry->node[iPoint]->GetCoord_p1())[iDim];

        /*--- Save grid velocity and gradient ---*/
        GridVel_Saved[iPoint][iDim] = (geometry->node[iPoint]->GetGridVel())[iDim];
        for (jDim = 0; jDim < nDim; jDim++) {
          GridVel_Grad_Saved[iPoint][iDim][jDim] = (geometry->node[iPoint]->GetGridVel_Grad())[iDim][jDim];
        }
      }
    }
if (geometry->node[iPoint]->GetGlobalIndex()==11){
/*cout << "sol node 300 x after set  " << solver[FLOW_SOL]->node[iPoint]->GetSolution(0)<< "   " << solver[FLOW_SOL]->node[iPoint]->GetSolution(1) <<endl;
cout << "soltn node 300 x after set  " << solver[FLOW_SOL]->node[iPoint]->GetSolution_time_n(0)<< "   " << solver[FLOW_SOL]->node[iPoint]->GetSolution_time_n(1) <<endl;
cout << "soltn1 node 300 x after set  " << (solver[FLOW_SOL]->node[iPoint]->GetSolution_time_n1())[0] << "   " <<(solver[FLOW_SOL]->node[iPoint]->GetSolution_time_n1())[1] <<endl;*/
//cout << "coord node 11 after set  " << geometry->node[iPoint]->GetCoord(0)<< "   " << geometry->node[iPoint]->GetCoord(1) <<endl;
/*cout << "coordn node 300 x after set  " << geometry->node[iPoint]->GetCoord_n(0)<< "   " << geometry->node[iPoint]->GetCoord_n(1) <<endl;
cout << "coordn1 node 300 x after set  " << geometry->node[iPoint]->GetCoord_n1(0)<< "   " << geometry->node[iPoint]->GetCoord_n1(1) <<endl;
cout << "coordp1 node 300 x after set  " << geometry->node[iPoint]->GetCoord_p1()<< "   " << geometry->node[iPoint]->GetCoord_p1() <<endl;
cout << "grid_vel node 300 x after set  " << geometry->node[iPoint]->GetGridVel()<< "   " << geometry->node[iPoint]->GetGridVel() <<endl;*/
}
  }

  /*--- Store the convergence state and the time step of the current simulation ---*/
  StopCalc_savedState = *StopCalc;
  dt_savedState = *dt;

  /*--- Inform preCICE that the writing task has been fulfilled ---*/
  solverInterface.fulfilledAction(cowic);

};

void CPreciceFlow::Reset_OldState( bool *StopCalc, double *dt ){
  unsigned long iPoint;
  unsigned short iDim, jDim;

  for (iPoint = 0; iPoint < nPoint; iPoint++) {
if (geometry->node[iPoint]->GetGlobalIndex()==11){
/*cout << "sol node 300 x before reset  " << solver[FLOW_SOL]->node[iPoint]->GetSolution(0)<< "   " << solver[FLOW_SOL]->node[iPoint]->GetSolution(1) <<endl;
cout << "soltn node 300 x before reset  " << solver[FLOW_SOL]->node[iPoint]->GetSolution_time_n(0)<< "   " << solver[FLOW_SOL]->node[iPoint]->GetSolution_time_n(1) <<endl;
cout << "soltn1 node 300 x before reset  " << (solver[FLOW_SOL]->node[iPoint]->GetSolution_time_n1())[0] << "   " <<(solver[FLOW_SOL]->node[iPoint]->GetSolution_time_n1())[1] <<endl;*/
//cout << "coord node 11 before reset  " << geometry->node[iPoint]->GetCoord(0)<< "   " << geometry->node[iPoint]->GetCoord(1) <<endl;
/*cout << "coordn node 300 x before reset  " << geometry->node[iPoint]->GetCoord_n(0)<< "   " << geometry->node[iPoint]->GetCoord_n(1) <<endl;
cout << "coordn1 node 300 x before reset  " << geometry->node[iPoint]->GetCoord_n1(0)<< "   " << geometry->node[iPoint]->GetCoord_n1(1) <<endl;
cout << "coordp1 node 300 x before reset  " << geometry->node[iPoint]->GetCoord_p1()<< "   " << geometry->node[iPoint]->GetCoord_p1() <<endl;
cout << "grid_vel node 300 x before reset  " << geometry->node[iPoint]->GetGridVel()<< "   " << geometry->node[iPoint]->GetGridVel() <<endl;*/
}
    /*--- Reset solutions from the last and current time step ---*/
    solver[FLOW_SOL]->node[iPoint]->SetSolution(solution_Saved[iPoint]);
    if (config->GetUnsteady_Simulation() != NO){
      solver[FLOW_SOL]->node[iPoint]->Set_Solution_time_n(solution_time_n_Saved[iPoint]);
      solver[FLOW_SOL]->node[iPoint]->Set_Solution_time_n1(solution_time_n1_Saved[iPoint]);
    }

    /*--- Reset solutions from the last, current and next time step ---*/    
    if (config->GetUnsteady_Simulation() != NO){
      geometry->node[iPoint]->SetCoord(Coord_n1_Saved[iPoint]);
      geometry->node[iPoint]->SetCoord_n();
      geometry->node[iPoint]->SetCoord_n1();
      geometry->node[iPoint]->SetCoord(Coord_n_Saved[iPoint]);
      geometry->node[iPoint]->SetCoord_n();
      geometry->node[iPoint]->SetCoord_p1(Coord_p1_Saved[iPoint]);
    }
    geometry->node[iPoint]->SetCoord(Coord_Saved[iPoint]);

    /*--- Reload grid velocity and gradient ---*/
    if (config->GetUnsteady_Simulation() != NO){
      geometry->node[iPoint]->SetGridVel(GridVel_Saved[iPoint]);
      for (iDim = 0; iDim < nDim; iDim++) {
        for (jDim = 0; jDim < nDim; jDim++) {
          geometry->node[iPoint]->SetGridVel_Grad(iDim, jDim, GridVel_Grad_Saved[iPoint][iDim][jDim]);
        }
      }
    }
if (geometry->node[iPoint]->GetGlobalIndex()==11){
/*cout << "sol node 300 x after reset  " << solver[FLOW_SOL]->node[iPoint]->GetSolution(0)<< "   " << solver[FLOW_SOL]->node[iPoint]->GetSolution(1) <<endl;
cout << "soltn node 300 x after reset  " << solver[FLOW_SOL]->node[iPoint]->GetSolution_time_n(0)<< "   " << solver[FLOW_SOL]->node[iPoint]->GetSolution_time_n(1) <<endl;
cout << "soltn1 node 300 x after reset  " << (solver[FLOW_SOL]->node[iPoint]->GetSolution_time_n1())[0] << "   " <<(solver[FLOW_SOL]->node[iPoint]->GetSolution_time_n1())[1] <<endl;*/
//cout << "coord node 11 after reset  " << geometry->node[iPoint]->GetCoord(0)<< "   " << geometry->node[iPoint]->GetCoord(1) <<endl;
/*cout << "coordn node 300 x after reset  " << geometry->node[iPoint]->GetCoord_n(0)<< "   " << geometry->node[iPoint]->GetCoord_n(1) <<endl;
cout << "coordn1 node 300 x after reset  " << geometry->node[iPoint]->GetCoord_n1(0)<< "   " << geometry->node[iPoint]->GetCoord_n1(1) <<endl;
cout << "coordp1 node 300 x after reset  " << geometry->node[iPoint]->GetCoord_p1()<< "   " << geometry->node[iPoint]->GetCoord_p1() <<endl;
cout << "grid_vel node 300 x after reset  " << geometry->node[iPoint]->GetGridVel()<< "   " << geometry->node[iPoint]->GetGridVel() <<endl;*/
}
  }

  /*--- Reset the convergence state and the time step of the last simulation ---*/
  *StopCalc = StopCalc_savedState;
  *dt = dt_savedState;

  /*--- Inform preCICE that the reading task has been fulfilled ---*/
  solverInterface.fulfilledAction(coric);

};

CPreciceFEA::CPreciceFEA( int processRank, int processSize, CGeometry**** geometry_container, CSolver***** solver_container,
                            CConfig** config_container, CVolumetricMovement*** grid_movement )
                           : CPrecice(processRank, processSize, geometry_container, solver_container, config_container, grid_movement)
{


  unsigned long iPoint;
  unsigned short iDim;

  /*--- Number of variables ---*/
  nVar = solver[FEA_SOL]->GetnVar();

  dt_savedState = 0.0;
  StopCalc_savedState = false;

  /*--- Initialize structures ---*/
  vertexIDs = NULL;
  forceID = NULL;
  displDeltaID = NULL;
  meshID = NULL;

  
  solution_Saved = new su2double*[nPoint];
  solution_Old_Saved = new su2double*[nPoint];
  solution_time_n_Saved = new su2double*[nPoint];
  solution_vel_Saved = new su2double*[nPoint];
  solution_vel_time_n_Saved = new su2double*[nPoint];
  solution_accel_Saved = new su2double*[nPoint];
  solution_accel_time_n_Saved = new su2double*[nPoint];

  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    solution_Saved[iPoint] = new su2double[nVar];
    solution_Old_Saved[iPoint] = new su2double[nVar];
    solution_time_n_Saved[iPoint] = new su2double[nVar];
    solution_vel_Saved[iPoint] = new su2double[nVar];
    solution_vel_time_n_Saved[iPoint] = new su2double[nVar];
    solution_accel_Saved[iPoint] = new su2double[nVar];
    solution_accel_time_n_Saved[iPoint] = new su2double[nVar];
  }

  /*--- Initialize magnitudes ---*/
  DisplacementDonor = NULL;
  DisplacementDonor_Prev = NULL;

  /*--- Configure preCICE with its configuration file ---*/
  Configure(config->GetpreCICE_ConfigFileName());

};

CPreciceFEA::~CPreciceFEA() {

  unsigned long iPoint;
  unsigned short iDim;

  for (iPoint = 0; iPoint < nPoint; iPoint++) {

    if (solution_Saved[iPoint]         != NULL) delete [] solution_Saved[iPoint];
    if (solution_time_n_Saved[iPoint]  != NULL) delete [] solution_time_n_Saved[iPoint];
    if (solution_vel_Saved[iPoint]  != NULL) delete [] solution_vel_Saved[iPoint];
    if (solution_vel_time_n_Saved[iPoint]  != NULL) delete [] solution_vel_time_n_Saved[iPoint];
    if (solution_accel_Saved[iPoint]  != NULL) delete [] solution_accel_Saved[iPoint];
    if (solution_accel_time_n_Saved[iPoint]  != NULL) delete [] solution_accel_time_n_Saved[iPoint];
    
  }

  if (solution_Saved         != NULL) delete [] solution_Saved;
  if (solution_time_n_Saved  != NULL) delete [] solution_time_n_Saved;
  if (solution_vel_Saved  != NULL) delete [] solution_vel_Saved;
  if (solution_vel_time_n_Saved  != NULL) delete [] solution_vel_time_n_Saved;
  if (solution_accel_Saved  != NULL) delete [] solution_accel_Saved;
  if (solution_accel_time_n_Saved  != NULL) delete [] solution_accel_time_n_Saved;

  if (vertexIDs              != NULL) delete [] vertexIDs;
  if (forceID                != NULL) delete [] forceID;
  if (displDeltaID           != NULL) delete [] displDeltaID;
  if (meshID                 != NULL) delete [] meshID;

};

su2double CPreciceFEA::Initialize() {
  
  unsigned long iSurface, jSurface;
  unsigned long iPoint;
  unsigned short iDim;
  int iVertex;
  su2double precice_dt;

  /*--- Allocate structures to interchange data with preCICE ---*/
  meshID = new int[nWetSurfaces];
  forceID = new int[nWetSurfaces];
  displDeltaID = new int[nWetSurfaces];

  /*--- Recover preCICE mesh IDs ---*/
  for (iSurface = 0; iSurface < nWetSurfaces; iSurface++) {
    meshID[iSurface] = solverInterface.getMeshID("SU2_Mesh_Solid" + to_string(iSurface));
  }

  /*--- Store the number of surfaces belonging to this rank ---*/
  for (iSurface = 0; iSurface < nWetSurfaces; iSurface++) {
    if (config->GetMarker_All_TagBound(config->GetpreCICE_WetSurfaceMarkerName() + to_string(iSurface)) != -1) {
      nWetSurfacesDomain++;
    }
  }

  /*--- Disable the processes in this rank if there are no wet surfaces that belong to it ---*/
  if (nWetSurfacesDomain < 1)  activeProcess = false;

  /*--- When this rank is active (that is, it holds a piece of the wet boundary) ---*/
  if (activeProcess) {

    /*--- Store the wet surface marker IDs in an array ---*/
    valueMarkerWet = new short[nWetSurfacesDomain];
    markerLocalToGlobal = new unsigned long[nWetSurfacesDomain];
    jSurface = 0;
    for (iSurface = 0; iSurface < nWetSurfaces; iSurface++) {
      if (config->GetMarker_All_TagBound(config->GetpreCICE_WetSurfaceMarkerName() + to_string(iSurface)) != -1) {
        valueMarkerWet[jSurface] = config->GetMarker_All_TagBound(config->GetpreCICE_WetSurfaceMarkerName() + to_string(iSurface));
        markerLocalToGlobal[jSurface] = iSurface;
        jSurface++;
      }
    }

    vertexIDs = new int*[nWetSurfacesDomain];
    nVertex = new int[nWetSurfacesDomain];

    for (iSurface = 0; iSurface < nWetSurfacesDomain; iSurface++) {

      nVertex[iSurface] = geometry->nVertex[valueMarkerWet[iSurface]];

      /*--- Allocate memory to store all the node coordinates in the surface ---*/
      su2double **nodeCoord;
      nodeCoord = new su2double*[nVertex[iSurface]];
      for (iVertex = 0; iVertex < nVertex[iSurface]; iVertex++)
        nodeCoord[iVertex] = new su2double[nDim];

      /*--- Loop over the nodes of each surface ---*/
      for (iVertex = 0; iVertex < nVertex[iSurface]; iVertex++) {

        /*--- Store the point index for the vertex nodes ---*/
        iPoint = geometry->vertex[valueMarkerWet[iSurface]][iVertex]->GetNode();

        /*--- Retrieve the coordinates ---*/
        for (iDim = 0; iDim < nDim; iDim++) {
          nodeCoord[iVertex][iDim] = geometry->node[iPoint]->GetCoord(iDim);
          if (verbose) cout << "Process #" << processRank << "/" << processSize-1 << ": Initial coordinates of node (local index, global index, node color): (" << iVertex << ", " << iPoint << ", " << geometry->node[iPoint]->GetColor() << "): " << nodeCoord[iVertex][iDim] << endl; /*--- for debugging purposes ---*/
        }
      }

      /*--- Passive double structure to pass on the vector of coordinates into preCICE ---*/
      passivedouble *coord;
      coord = new passivedouble[nVertex[iSurface]*nDim];

      /*--- Convert the coordinates to a passive double structure that will be communicated to preCICE ---*/
      for (iVertex = 0; iVertex < nVertex[iSurface]; iVertex++) {
        for (iDim = 0; iDim < nDim; iDim++) {
          coord[iVertex*nDim + iDim] = SU2_TYPE::GetValue(nodeCoord[iVertex][iDim]);
        }
      }

      /*--- Variable internal to preCICE ---*/
      vertexIDs[iSurface] = new int[nVertex[iSurface]];

      /*--- Interface to preCICE ---*/
      solverInterface.setMeshVertices(meshID[markerLocalToGlobal[iSurface]], nVertex[iSurface], coord, vertexIDs[iSurface]);
      forceID[markerLocalToGlobal[iSurface]] = solverInterface.getDataID("Forces" + to_string(markerLocalToGlobal[iSurface]), meshID[markerLocalToGlobal[iSurface]]);
      displDeltaID[markerLocalToGlobal[iSurface]] = solverInterface.getDataID("DisplacementDeltas" + to_string(markerLocalToGlobal[iSurface]), meshID[markerLocalToGlobal[iSurface]]);

      /*--- Deallocate the containers of the coordinates ---*/
      if (coord != NULL) delete [] coord;
      for (iVertex = 0; iVertex < nVertex[iSurface]; iVertex++){
        if (nodeCoord[iVertex]!= NULL) delete [] nodeCoord[iVertex];
      }
      if (nodeCoord!= NULL) delete [] nodeCoord;

    }
    for (iSurface = 0; iSurface < nWetSurfaces; iSurface++) {
      for (jSurface = 0; jSurface < nWetSurfacesDomain; jSurface++) {
        if (markerLocalToGlobal[jSurface] != iSurface) {
          solverInterface.setMeshVertices(meshID[iSurface], 0, NULL, NULL);
          forceID[iSurface] = solverInterface.getDataID("Forces" + to_string(iSurface), meshID[iSurface]);
          displDeltaID[iSurface] = solverInterface.getDataID("DisplacementDeltas" + to_string(iSurface), meshID[iSurface]);
        }
      }
    }
  } else {
    /*--- If the surface is not active ---*/
    for (iSurface = 0; iSurface < nWetSurfaces; iSurface++) {
      solverInterface.setMeshVertices(meshID[iSurface], 0, NULL, NULL);
      forceID[iSurface] = solverInterface.getDataID("Forces" + to_string(iSurface), meshID[iSurface]);
      displDeltaID[iSurface] = solverInterface.getDataID("DisplacementDeltas" + to_string(iSurface), meshID[iSurface]);
    }
  }

  /*--- Finally, initiallize preCICE with the information that we have sent to it ---*/
  precice_dt = solverInterface.initialize();

  /*--- readData ---*/
  if (solverInterface.isReadDataAvailable()) {

    /*--- Retrieve the forces from preCICE for the wet surfaces that belong to the active process---*/
    for (iSurface = 0; iSurface < nWetSurfacesDomain; iSurface++) {

      /*--- Allocate memory to store all the forces in the wet surface ---*/
      su2double **forces_su2;
      forces_su2 = new su2double*[nVertex[iSurface]];
      for (iVertex = 0; iVertex < nVertex[iSurface]; iVertex++)
        forces_su2[iVertex] = new su2double[nDim];

      /*--- Allocate memory to retrieve all the forces in the wet surface coming from preCICE ---*/
      passivedouble *forces;
      forces = new passivedouble[nVertex[iSurface]*nDim];

      /*--- Read the force data from preCICE ---*/
      solverInterface.readBlockVectorData(forceID[markerLocalToGlobal[iSurface]], nVertex[iSurface], vertexIDs[iSurface], forces);

      /*--- Store the forces coming from preCICE in the container for SU2 ---*/
      for (iVertex = 0; iVertex < nVertex[iSurface]; iVertex++) {
        for (iDim = 0; iDim < nDim; iDim++) {
          forces_su2[iVertex][iDim] = forces[iVertex*nDim + iDim];
//if (geometry->node[geometry->vertex[valueMarkerWet[iSurface]][iVertex]->GetNode()]->GetGlobalIndex() == 36){
//cout << "ForcesInFea 36  "<<  geometry->node[geometry->vertex[valueMarkerWet[iSurface]][iVertex]->GetNode()]->GetGlobalIndex()   << "   " << forces[iVertex*nDim + iDim] << endl;
//}
	}
      }
	
      /*--- declare Point_Struct ---*/
      unsigned long Point_Struct;
      /*--- Loop over vertices of coupled boundary ---*/
      for (iVertex = 0; iVertex < nVertex[iSurface]; iVertex++) {
        Point_Struct = geometry->vertex[valueMarkerWet[iSurface]][iVertex]->GetNode();
        /*--- Add the forces to the Structural Solver ---*/
        solver[FEA_SOL]->node[Point_Struct]->Set_FlowTraction(forces_su2[iVertex]);
//if (geometry->node[Point_Struct]->GetGlobalIndex()==36){
//cout << "flow tractions node 36 x   " << solver[FEA_SOL]->node[Point_Struct]->Get_FlowTraction(0) <<endl;
//}
      }

      /*--- Deallocate the containers of the force ---*/
      if (forces != NULL)  delete [] forces;
      for (iVertex = 0; iVertex < nVertex[iSurface]; iVertex++){
        if (forces_su2[iVertex]!= NULL) delete [] forces_su2[iVertex];
      }
      if (forces_su2!= NULL) delete [] forces_su2;

    }
  }

  return precice_dt;

};

su2double CPreciceFEA::Advance( su2double computedTimestep ) { 

  unsigned long iSurface;
  unsigned short iDim, jDim;
  int iVertex;
  su2double max_precice_dt;

  if (activeProcess) {

   
    /*--- Compute the displacements for the wet surfaces that belong to the active process ---*/
    for (iSurface = 0; iSurface < nWetSurfacesDomain; iSurface++) {

      /*--- Allocate memory to store all the displacements in the surface ---*/
      su2double **displacements_su2;
      displacements_su2 = new su2double*[nVertex[iSurface]];
      for (iVertex = 0; iVertex < nVertex[iSurface]; iVertex++)
        displacements_su2[iVertex] = new su2double[nDim];

      /*--- Allocate a passive double structure to store the displacements as they are transferred to preCICE ---*/
      passivedouble *displacementDeltas;
      displacementDeltas = new passivedouble[nVertex[iSurface]*nDim];

      
      unsigned long Point_Struct;

      /*--- Loop over vertices of coupled boundary ---*/
      for (iVertex = 0; iVertex < nVertex[iSurface]; iVertex++) {

        /*--- Get the vertex index  ---*/
        Point_Struct = geometry->vertex[valueMarkerWet[iSurface]][iVertex]->GetNode();

        /*--- Retrieve the values of the displacements ---*/
	DisplacementDonor = solver[FEA_SOL]->node[Point_Struct]->GetSolution();
	DisplacementDonor_Prev = solver[FEA_SOL]->node[Point_Struct]->GetSolution_time_n();
/*	if (geometry->node[Point_Struct]->GetGlobalIndex()==36){
cout << "solution node 36   " << solver[FEA_SOL]->node[Point_Struct]->GetSolution(0) <<endl;
cout << "solution time n node 36   " << (solver[FEA_SOL]->node[Point_Struct]->Get_femSolution_time_n())[0] <<endl;
}  */

        /*--- Compute the displacements ---*/
        for (iDim = 0; iDim < nDim; iDim++){
	  //displacements_su2[iVertex][iDim] = DisplacementDonor[iDim];
          displacements_su2[iVertex][iDim] = DisplacementDonor[iDim] - DisplacementDonor_Prev[iDim] ;
	}
      }

      /*--- Store the displacements in a passive double structure, for transfer into preCICE ---*/
      for (iVertex = 0; iVertex < nVertex[iSurface]; iVertex++) {
        for (iDim = 0; iDim < nDim; iDim++) {
	  /*TODO under investigation... -> consistent mapping?!*/
          /*if (geometry->node[Point_Struct]->GetColor() == processRank)*/
            displacementDeltas[iVertex*nDim + iDim] = SU2_TYPE::GetValue(displacements_su2[iVertex][iDim]);
          /*else
            displacementDeltas[iVertex*nDim + iDim] = 0.0;*/
//if (geometry->node[geometry->vertex[valueMarkerWet[iSurface]][iVertex]->GetNode()]->GetGlobalIndex() == 36){
//cout << "DispOutFea 36   "<< geometry->node[geometry->vertex[valueMarkerWet[iSurface]][iVertex]->GetNode()]->GetGlobalIndex()  << "   " << displacementDeltas[iVertex*nDim + iDim] << endl;
//}
        }
      }

      /*--- Send the data to preCICE ---*/
      solverInterface.writeBlockVectorData(displDeltaID[markerLocalToGlobal[iSurface]], nVertex[iSurface], vertexIDs[iSurface], displacementDeltas);

      /*--- Deallocate the containers of the displacements ---*/
      if (displacementDeltas != NULL) delete [] displacementDeltas;
      for (iVertex = 0; iVertex < nVertex[iSurface]; iVertex++){
        if (displacements_su2[iVertex]!= NULL) delete [] displacements_su2[iVertex];
      }
      if (displacements_su2!= NULL) delete [] displacements_su2;

    }

    /*--- Advance the interface in preCICE ---*/
    max_precice_dt = solverInterface.advance( SU2_TYPE::GetValue(computedTimestep) );

    /*--- Retrieve the forces from preCICE for the wet surfaces that belong to the active process---*/
    for (iSurface = 0; iSurface < nWetSurfacesDomain; iSurface++) {

      /*--- Allocate memory to store all the forces in the wet surface ---*/
      su2double **forces_su2;
      forces_su2 = new su2double*[nVertex[iSurface]];
      for (iVertex = 0; iVertex < nVertex[iSurface]; iVertex++)
        forces_su2[iVertex] = new su2double[nDim];

      /*--- Allocate memory to retrieve all the forces in the wet surface coming from preCICE ---*/
      passivedouble *forces;
      forces = new passivedouble[nVertex[iSurface]*nDim];

      /*--- Read the force data from preCICE ---*/
      solverInterface.readBlockVectorData(forceID[markerLocalToGlobal[iSurface]], nVertex[iSurface], vertexIDs[iSurface], forces);

      /*--- Store the forces coming from preCICE in the container for SU2 ---*/
      for (iVertex = 0; iVertex < nVertex[iSurface]; iVertex++) {
        for (iDim = 0; iDim < nDim; iDim++) {
          forces_su2[iVertex][iDim] = forces[iVertex*nDim + iDim];
//if (geometry->node[geometry->vertex[valueMarkerWet[iSurface]][iVertex]->GetNode()]->GetGlobalIndex() == 36){
//cout << "ForcesInFea 36  "<< geometry->node[geometry->vertex[valueMarkerWet[iSurface]][iVertex]->GetNode()]->GetGlobalIndex()  << "   " << forces[iVertex*nDim + iDim] << endl;
//}
	}
      }
	
      /*--- declare Point_Struct ---*/
      unsigned long Point_Struct;
      /*--- Loop over vertices of coupled boundary ---*/
      for (iVertex = 0; iVertex < nVertex[iSurface]; iVertex++) {
        Point_Struct = geometry->vertex[valueMarkerWet[iSurface]][iVertex]->GetNode();
        /*--- Add the forces to the Structural Solver ---*/
        solver[FEA_SOL]->node[Point_Struct]->Set_FlowTraction(forces_su2[iVertex]);
//      if (geometry->node[Point_Struct]->GetGlobalIndex()==36){
//cout << "flow tractions node 36 x   " << solver[FEA_SOL]->node[Point_Struct]->Get_FlowTraction(0) <<endl;
//}  
      }

      /*--- Deallocate the containers of the force ---*/
      if (forces != NULL)  delete [] forces;
      for (iVertex = 0; iVertex < nVertex[iSurface]; iVertex++){
        if (forces_su2[iVertex]!= NULL) delete [] forces_su2[iVertex];
      }
      if (forces_su2!= NULL) delete [] forces_su2;

    }

 return max_precice_dt;
  }
  else {
    /*--- Advance the interface (as nothing is transferred from this process) ---*/
    max_precice_dt = solverInterface.advance( SU2_TYPE::GetValue(computedTimestep) );
    return max_precice_dt;
  }

};


void CPreciceFEA::Set_OldState( bool *StopCalc, double *dt ) { 

  unsigned long iPoint;
  unsigned short iVar, iDim, jDim;
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
if (geometry->node[iPoint]->GetGlobalIndex()==3){
//cout << "sol node 3 before set  " << solver[FEA_SOL]->node[iPoint]->GetSolution(0)<< "   " << solver[FEA_SOL]->node[iPoint]->GetSolution(1) <<endl;
/*cout << "soltn node 36 x before set  " << solver[FEA_SOL]->node[iPoint]->GetSolution_time_n(0)<< "   " << solver[FEA_SOL]->node[iPoint]->GetSolution_time_n(1) <<endl;
cout << "solvel node 36 x before set  " << solver[FEA_SOL]->node[iPoint]->GetSolution_Vel(0)<< "   " << solver[FEA_SOL]->node[iPoint]->GetSolution_Vel(1) <<endl;
cout << "solveltn node 36 x before set  " << solver[FEA_SOL]->node[iPoint]->GetSolution_Vel_time_n(0)<< "   " << solver[FEA_SOL]->node[iPoint]->GetSolution_Vel_time_n(1) <<endl;
cout << "solaccel node 36 x before set  " << solver[FEA_SOL]->node[iPoint]->GetSolution_Accel(0)<< "   " << solver[FEA_SOL]->node[iPoint]->GetSolution_Accel(1) <<endl;
cout << "solacceltn node 36 x before set  " << solver[FEA_SOL]->node[iPoint]->GetSolution_Accel_time_n(0)<< "   " << solver[FEA_SOL]->node[iPoint]->GetSolution_Accel_time_n(1) <<endl;*/
}
    for (iVar = 0; iVar < nVar; iVar++) {
      /*--- Save solutions at last and current time step ---*/
      solution_Saved[iPoint][iVar] = (solver[FEA_SOL]->node[iPoint]->GetSolution())[iVar];
      if (config->GetDynamic_Analysis() == DYNAMIC){
	solution_time_n_Saved[iPoint][iVar] = (solver[FEA_SOL]->node[iPoint]->GetSolution_time_n(iVar));
        solution_vel_Saved[iPoint][iVar] = (solver[FEA_SOL]->node[iPoint]->GetSolution_Vel())[iVar];
        solution_vel_time_n_Saved[iPoint][iVar] = (solver[FEA_SOL]->node[iPoint]->GetSolution_Vel_time_n())[iVar];
        solution_accel_Saved[iPoint][iVar] = (solver[FEA_SOL]->node[iPoint]->GetSolution_Accel())[iVar];
        solution_accel_time_n_Saved[iPoint][iVar] = (solver[FEA_SOL]->node[iPoint]->GetSolution_Accel_time_n())[iVar];
      }    
    }
	if (geometry->node[iPoint]->GetGlobalIndex()==3){
//cout << "sol node 3 after set  " << solver[FEA_SOL]->node[iPoint]->GetSolution(0)<< "   " << solver[FEA_SOL]->node[iPoint]->GetSolution(1) <<endl;
/*cout << "soltn node 36 x after set  " << solver[FEA_SOL]->node[iPoint]->GetSolution_time_n(0)<< "   " << solver[FEA_SOL]->node[iPoint]->GetSolution_time_n(1) <<endl;
cout << "solvel node 36 x after set  " << solver[FEA_SOL]->node[iPoint]->GetSolution_Vel(0)<< "   " << solver[FEA_SOL]->node[iPoint]->GetSolution_Vel(1) <<endl;
cout << "solveltn node 36 x after set  " << solver[FEA_SOL]->node[iPoint]->GetSolution_Vel_time_n(0)<< "   " << solver[FEA_SOL]->node[iPoint]->GetSolution_Vel_time_n(1) <<endl;
cout << "solaccel node 36 x after set  " << solver[FEA_SOL]->node[iPoint]->GetSolution_Accel(0)<< "   " << solver[FEA_SOL]->node[iPoint]->GetSolution_Accel(1) <<endl;
cout << "solacceltn node 36 x after set  " << solver[FEA_SOL]->node[iPoint]->GetSolution_Accel_time_n(0)<< "   " << solver[FEA_SOL]->node[iPoint]->GetSolution_Accel_time_n(1) <<endl;*/  
        }
  }

  /*--- Store the convergence state and the time step of the current simulation ---*/
  StopCalc_savedState = *StopCalc;
  dt_savedState = *dt;

  /*--- Inform preCICE that the writing task has been fulfilled ---*/
  solverInterface.fulfilledAction(cowic);

};

void CPreciceFEA::Reset_OldState( bool *StopCalc, double *dt ) {
  unsigned long iPoint;
  unsigned short iDim, jDim;
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
if (geometry->node[iPoint]->GetGlobalIndex()==3){
//cout << "sol node 3 before reset  " << solver[FEA_SOL]->node[iPoint]->GetSolution(0)<< "   " << solver[FEA_SOL]->node[iPoint]->GetSolution(1) <<endl;
/*cout << "soltn node 36 x before reset  " << solver[FEA_SOL]->node[iPoint]->GetSolution_time_n(0)<< "   " << solver[FEA_SOL]->node[iPoint]->GetSolution_time_n(1) <<endl;
cout << "solvel node 36 x before reset  " << solver[FEA_SOL]->node[iPoint]->GetSolution_Vel(0)<< "   " << solver[FEA_SOL]->node[iPoint]->GetSolution_Vel(1) <<endl;
cout << "solveltn node 36 x before reset  " << solver[FEA_SOL]->node[iPoint]->GetSolution_Vel_time_n(0)<< "   " << solver[FEA_SOL]->node[iPoint]->GetSolution_Vel_time_n(1) <<endl;
cout << "solaccel node 36 x before reset  " << solver[FEA_SOL]->node[iPoint]->GetSolution_Accel(0)<< "   " << solver[FEA_SOL]->node[iPoint]->GetSolution_Accel(1) <<endl;
cout << "solacceltn node 36 x before reset  " << solver[FEA_SOL]->node[iPoint]->GetSolution_Accel_time_n(0)<< "   " << solver[FEA_SOL]->node[iPoint]->GetSolution_Accel_time_n(1) <<endl;*/
}
      /*--- Reset solutions from the last and current time step ---*/
    solver[FEA_SOL]->node[iPoint]->SetSolution(solution_Saved[iPoint]); 
    if (config->GetDynamic_Analysis() == DYNAMIC){    
      solver[FEA_SOL]->node[iPoint]->SetSolution_time_n(solution_time_n_Saved[iPoint]);
      solver[FEA_SOL]->node[iPoint]->SetSolution_Vel(solution_vel_Saved[iPoint]);
      solver[FEA_SOL]->node[iPoint]->SetSolution_Vel_time_n(solution_vel_time_n_Saved[iPoint]);
      solver[FEA_SOL]->node[iPoint]->SetSolution_Accel(solution_accel_Saved[iPoint]);
      solver[FEA_SOL]->node[iPoint]->SetSolution_Accel_time_n(solution_accel_time_n_Saved[iPoint]);
    }
if (geometry->node[iPoint]->GetGlobalIndex()==3){
//cout << "sol node 3 after reset  " << solver[FEA_SOL]->node[iPoint]->GetSolution(0)<< "   " << solver[FEA_SOL]->node[iPoint]->GetSolution(1) << iPoint <<endl;
/*cout << "soltn node 36 x after reset  " << solver[FEA_SOL]->node[iPoint]->GetSolution_time_n(0)<< "   " << solver[FEA_SOL]->node[iPoint]->GetSolution_time_n(1) <<endl;
cout << "solvel node 36 x after reset  " << solver[FEA_SOL]->node[iPoint]->GetSolution_Vel(0)<< "   " << solver[FEA_SOL]->node[iPoint]->GetSolution_Vel(1) <<endl;
cout << "solveltn node 36 x after reset  " << solver[FEA_SOL]->node[iPoint]->GetSolution_Vel_time_n(0)<< "   " << solver[FEA_SOL]->node[iPoint]->GetSolution_Vel_time_n(1) <<endl;
cout << "solaccel node 36 x after reset  " << solver[FEA_SOL]->node[iPoint]->GetSolution_Accel(0)<< "   " << solver[FEA_SOL]->node[iPoint]->GetSolution_Accel(1) <<endl;
cout << "solacceltn node 36 x after reset  " << solver[FEA_SOL]->node[iPoint]->GetSolution_Accel_time_n(0)<< "   " << solver[FEA_SOL]->node[iPoint]->GetSolution_Accel_time_n(1) <<endl; */
 }
  }

  /*--- Reset the convergence state and the time step of the last simulation ---*/
  *StopCalc = StopCalc_savedState;
  *dt = dt_savedState;

  /*--- Inform preCICE that the reading task has been fulfilled ---*/
  solverInterface.fulfilledAction(coric);
};
