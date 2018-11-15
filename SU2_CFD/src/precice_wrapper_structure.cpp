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
solverInterface( "SU2_CFD", processRank, processSize),
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

  forces = NULL;
  displacementDeltas = NULL;

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
  if (forces                 != NULL) delete [] forces;
  if (displacementDeltas     != NULL) delete [] displacementDeltas;

};

void CPreciceFlow::Configure( const string& configurationFilename ) {

  if (verbose) cout << "Process #" << processRank << "/" << processSize-1 << ": Configuring preCICE..." << endl;

  solverInterface.configure( configurationFilename );

  if (verbose) cout << "Process #" << processRank << "/" << processSize-1 << ": ...done configuring preCICE!" << endl;

  //Checking for dimensional consistency of SU2 and preCICE - Exit if not consistent
  if(solverInterface.getDimensions() != nDim){
    cout << "Dimensions of SU2 and preCICE are not equal! Now exiting..." << endl;
    exit(EXIT_FAILURE);
  }

  //Checking for number of wet surfaces - Exit if not cat least one wet surface defined
  if(nWetSurfaces < 1){
    cout << "There must be at least one wet surface! Now exiting..." << endl;
    exit(EXIT_FAILURE);
  }

};

su2double CPreciceFlow::Initialize() {

  unsigned long iSurface, jSurface;
  unsigned long iPoint;
  unsigned short iDim;
  int iVertex;
  su2double precice_dt;

  if (verbose) cout << "Process #" << processRank << "/" << processSize-1 << ": Initializing preCICE..." << endl;

  meshID = new int[nWetSurfaces];
  forceID = new int[nWetSurfaces];
  displDeltaID = new int[nWetSurfaces];

  /*--- Recover preCICE mesh IDs ---*/
  for (iSurface = 0; iSurface < nWetSurfaces; iSurface++) {
    meshID[iSurface] = solverInterface.getMeshID("SU2_Mesh" + to_string(iSurface));
  }

  /*--- Store the number of surfaces belonging to this rank ---*/
  for (iSurface = 0; iSurface < nWetSurfaces; iSurface++) {
    if (config->GetMarker_All_TagBound(config->GetpreCICE_WetSurfaceMarkerName() + to_string(iSurface)) == -1) {
      cout << "Process #" << processRank << "/" << processSize-1 << ": Does not work on " << config->GetpreCICE_WetSurfaceMarkerName() << iSurface << endl;
    } else {
      nWetSurfacesDomain++;
    }
  }

  /*--- Disable the processes in this rank if there are no wet surfaces that belong to it ---*/
  if (nWetSurfacesDomain < 1) {
    cout << "Process #" << processRank << "/" << processSize-1 << ": Does not work on the wet surface at all." << endl;
    activeProcess = false;
  }

  if (activeProcess) {

    //Store the wet surface marker values in an array, which has the size equal to the number of wet surfaces actually being worked on by this process
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
        for (int iDim = 0; iDim < nDim; iDim++) {
          nodeCoord[iVertex][iDim] = geometry->node[iPoint]->GetCoord(iDim);
          if (verbose) {
            cout << "Process #" << processRank << "/" << processSize-1 << ": Initial coordinates of node (local index, global index, node color): (" << iVertex << ", " << iPoint << ", " << geometry->node[iPoint]->GetColor() << "): " << nodeCoord[iVertex][iDim] << endl; /*--- for debugging purposes ---*/
          }
        }
      }

      /*--- Passive double structure to pass on the vector of coordinates into preCICE ---*/
      passivedouble *coord;
      coord = new passivedouble[nVertex[iSurface]*nDim];
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

    }
    for (iSurface = 0; iSurface < nWetSurfaces; iSurface++) {
      bool flag = false;
      for (jSurface = 0; jSurface < nWetSurfacesDomain; jSurface++) {
        if (markerLocalToGlobal[jSurface] == iSurface) {
          flag = true;
        }
      }
      if (!flag) {
        solverInterface.setMeshVertices(meshID[iSurface], 0, NULL, NULL);
        forceID[iSurface] = solverInterface.getDataID("Forces" + to_string(iSurface), meshID[iSurface]);
        displDeltaID[iSurface] = solverInterface.getDataID("DisplacementDeltas" + to_string(iSurface), meshID[iSurface]);
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

  if (verbose) {
    cout << "Process #" << processRank << "/" << processSize-1 << ": There is grid movement (expected: 1): " << config->GetGrid_Movement() << endl; /*--- for debugging purposes ---*/
    cout << "Process #" << processRank << "/" << processSize-1 << ": Kind of grid movement (expected: 13): " << config->GetKind_GridMovement(ZONE_0) << endl;  /*--- for debugging purposes ---*/
  }

  precice_dt = solverInterface.initialize();
  if (verbose) {
    cout << "Process #" << processRank << "/" << processSize-1 << ": ...done initializing preCICE!" << endl;
  }
  return precice_dt;

};

su2double CPreciceFlow::Advance( su2double computedTimestep ) {

  if (activeProcess) {
    if (verbose) {
      cout << "Process #" << processRank << "/" << processSize-1 << ": Advancing preCICE..." << endl;
    }

    //Get physical simulation information
    bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
    bool viscous_flow = ((config->GetKind_Solver() == NAVIER_STOKES) || (config->GetKind_Solver() == RANS));

    //Compute factorForces for redimensionalizing forces ("ND" = Non-Dimensional)
    Velocity_Real = config->GetVelocity_FreeStream();
    Density_Real = config->GetDensity_FreeStream();
    Velocity_ND = config->GetVelocity_FreeStreamND();
    Density_ND = config->GetDensity_FreeStreamND();
    Velocity2_Real = 0.0;  /*--- denotes squared real velocity ---*/
    Velocity2_ND = 0.0;  /*--- denotes squared non-dimensional velocity ---*/
    //Compute squared values
    for (int iDim = 0; iDim < nDim; iDim++){
      Velocity2_Real += Velocity_Real[iDim]*Velocity_Real[iDim];
      Velocity2_ND += Velocity_ND[iDim]*Velocity_ND[iDim];
    }
    //Compute factor for redimensionalizing forces
    double factorForces = Density_Real*Velocity2_Real/(Density_ND*Velocity2_ND);
    if (verbose) {
      cout << "Process #" << processRank << "/" << processSize-1 << ": Factor for (non-/re-)dimensionalization of forces: " << factorForces << endl;  /*--- for debugging purposes ---*/
    }

    for (int i = 0; i < nWetSurfacesDomain; i++) {
      if (verbose) {
        //1. Compute forces
        cout << "Process #" << processRank << "/" << processSize-1 << ": Advancing preCICE: Computing forces for " << config->GetpreCICE_WetSurfaceMarkerName() << markerLocalToGlobal[i] << "..." << endl;
      }
      //Some variables to be used:
      unsigned long nodeVertex[nVertex[i]];
      double normalsVertex[nVertex[i]][nDim];
      double normalsVertex_Unit[nVertex[i]][nDim];
      double Area;
      double Pn = 0.0;  /*--- denotes pressure at a node ---*/
      double Pinf = 0.0;  /*--- denotes environmental (farfield) pressure ---*/
      double** Grad_PrimVar = NULL; /*--- denotes (u.A. velocity) gradients needed for computation of viscous forces ---*/
      double Viscosity = 0.0;
      double Tau[3][3];
      double TauElem[3];
      double forces_su2[nVertex[i]][nDim];  /*--- forces will be stored such, before converting to simple array ---*/

      /*--- Loop over vertices of coupled boundary ---*/
      for (int iVertex = 0; iVertex < nVertex[i]; iVertex++) {
        //Get node number (= index) to vertex (= node)
        nodeVertex[iVertex] = geometry->vertex[valueMarkerWet[i]][iVertex]->GetNode(); /*--- Store all nodes (indices) in a vector ---*/
        // Get normal vector
        for (int iDim = 0; iDim < nDim; iDim++){
          normalsVertex[iVertex][iDim] = (geometry->vertex[valueMarkerWet[i]][iVertex]->GetNormal())[iDim];
        }
        // Unit normals
        Area = 0.0;
        for (int iDim = 0; iDim < nDim; iDim++) {
          Area += normalsVertex[iVertex][iDim]*normalsVertex[iVertex][iDim];
        }
        Area = sqrt(Area);
        for (int iDim = 0; iDim < nDim; iDim++) {
          normalsVertex_Unit[iVertex][iDim] = normalsVertex[iVertex][iDim]/Area;
        }
        // Get the values of pressure and viscosity
        Pn = solver[FLOW_SOL]->node[nodeVertex[iVertex]]->GetPressure();
        Pinf = solver[FLOW_SOL]->GetPressure_Inf();
        if (viscous_flow){
          Grad_PrimVar = solver[FLOW_SOL]->node[nodeVertex[iVertex]]->GetGradient_Primitive();
          Viscosity = solver[FLOW_SOL]->node[nodeVertex[iVertex]]->GetLaminarViscosity();
        }

        // Calculate the forces_su2 in the nodes for the inviscid term --> Units of force (non-dimensional).
        for (int iDim = 0; iDim < nDim; iDim++) {
          forces_su2[iVertex][iDim] = -(Pn-Pinf)*normalsVertex[iVertex][iDim];
        }
        // Calculate the forces_su2 in the nodes for the viscous term
        if (viscous_flow){
          // Divergence of the velocity
          double div_vel = 0.0;
          for (int iDim = 0; iDim < nDim; iDim++){
            div_vel += Grad_PrimVar[iDim+1][iDim];
          }
          if (incompressible){
            div_vel = 0.0;  /*--- incompressible flow is divergence-free ---*/
          }
          for (int iDim = 0; iDim < nDim; iDim++) {
            for (int jDim = 0 ; jDim < nDim; jDim++) {
              // Dirac delta
              double Delta = 0.0;
              if (iDim == jDim){
                Delta = 1.0;
              }
              // Viscous stress
              Tau[iDim][jDim] = Viscosity*(Grad_PrimVar[jDim+1][iDim] + Grad_PrimVar[iDim+1][jDim]) -
              2/3*Viscosity*div_vel*Delta;
              // Add Viscous component in the forces_su2 vector --> Units of force (non-dimensional).
              forces_su2[iVertex][iDim] += Tau[iDim][jDim]*normalsVertex[iVertex][jDim];
            }
          }
        }
        // Rescale forces_su2 to SI units
        for (int iDim = 0; iDim < nDim; iDim++) {
          forces_su2[iVertex][iDim] = forces_su2[iVertex][iDim]*factorForces;
        }
      }
      //convert forces_su2 into forces
      forces = new double[nVertex[i]*nDim];
      for (int iVertex = 0; iVertex < nVertex[i]; iVertex++) {
        for (int iDim = 0; iDim < nDim; iDim++) {
          //Do not write forces for duplicate nodes! -> Check wether the color of the node matches the MPI-rank of this process. Only write forces, if node originally belongs to this process.
          if (geometry->node[nodeVertex[iVertex]]->GetColor() == processRank) {
            forces[iVertex*nDim + iDim] = forces_su2[iVertex][iDim];
          }
          else{
            forces[iVertex*nDim + iDim] = 0;
          }
        }
      }
      if (verbose) {
        cout << "Process #" << processRank << "/" << processSize-1 << ": Advancing preCICE: ...done computing forces for " << config->GetpreCICE_WetSurfaceMarkerName() << markerLocalToGlobal[i] << endl;
      }

      //2. Write forces
      if (verbose) {
        cout << "Process #" << processRank << "/" << processSize-1 << ": Advancing preCICE: Writing forces for " << config->GetpreCICE_WetSurfaceMarkerName() << markerLocalToGlobal[i] << "..." << endl;
      }
      //Load Ramping functionality: Reduce force vector before transferring by a ramping factor, which increases with the number of elapsed time steps; Achtung: ExtIter beginnt bei 0 (ohne Restart) und bei einem Restart (Startlösung) nicht bei 0, sondern bei der Startiterationsnummer
      if (config->GetpreCICE_LoadRamping() && ((config->GetExtIter() - config->GetUnst_RestartIter()) < config->GetpreCICE_LoadRampingDuration())) {
        if (verbose) {
          cout << "Process #" << processRank << "/" << processSize-1 << ": Load ramping factor in preCICE: " << config->GetExtIter() - config->GetUnst_RestartIter() + 1 << "/" << config->GetpreCICE_LoadRampingDuration() << endl;
        }
        *forces = *forces * ((config->GetExtIter() - config->GetUnst_RestartIter()) + 1) / config->GetpreCICE_LoadRampingDuration();
      }
      solverInterface.writeBlockVectorData(forceID[markerLocalToGlobal[i]], nVertex[i], vertexIDs[i], forces);
      if (verbose) {
        cout << "Process #" << processRank << "/" << processSize-1 << ": Advancing preCICE: ...done writing forces for " << config->GetpreCICE_WetSurfaceMarkerName() << markerLocalToGlobal[i] << "." << endl;
      }
      if (forces != NULL){
        delete [] forces;
      }
    }

    //3. Advance solverInterface
    if (verbose) cout << "Process #" << processRank << "/" << processSize-1 << ": Advancing preCICE: Advancing SolverInterface..." << endl;

    su2double max_precice_dt;

    max_precice_dt = solverInterface.advance( SU2_TYPE::GetValue(computedTimestep) );

    if (verbose) cout << "Process #" << processRank << "/" << processSize-1 << ": Advancing preCICE: ...done advancing SolverInterface." << endl;

    // displacements = new double[nVertex*nDim]; //TODO: Delete later
    for (int i = 0; i < nWetSurfacesDomain; i++) {
      //4. Read displacements/displacementDeltas
      if (verbose) {
        cout << "Process #" << processRank << "/" << processSize-1 << ": Advancing preCICE: Reading displacement deltas for " << config->GetpreCICE_WetSurfaceMarkerName() << markerLocalToGlobal[i] << "..." << endl;
      }
      double displacementDeltas_su2[nVertex[i]][nDim]; /*--- displacementDeltas will be stored such, before converting to simple array ---*/
      displacementDeltas = new double[nVertex[i]*nDim];
      solverInterface.readBlockVectorData(displDeltaID[markerLocalToGlobal[i]], nVertex[i], vertexIDs[i], displacementDeltas);
      if (verbose) {
        cout << "Process #" << processRank << "/" << processSize-1 << ": Advancing preCICE: ...done reading displacement deltas for " << config->GetpreCICE_WetSurfaceMarkerName() << markerLocalToGlobal[i] << "." << endl;
      }

      //5. Set displacements/displacementDeltas
      if (verbose) {
        cout << "Process #" << processRank << "/" << processSize-1 << ": Advancing preCICE: Setting displacement deltas for " << config->GetpreCICE_WetSurfaceMarkerName() << markerLocalToGlobal[i] << "..." << endl;
      }
      //convert displacementDeltas into displacementDeltas_su2
      for (int iVertex = 0; iVertex < nVertex[i]; iVertex++) {
        for (int iDim = 0; iDim < nDim; iDim++) {
          displacementDeltas_su2[iVertex][iDim] = displacementDeltas[iVertex*nDim + iDim];
        }
      }
      if (displacementDeltas != NULL) {
        delete [] displacementDeltas;
      }

      //Set change of coordinates (i.e. displacementDeltas)
      for (int iVertex = 0; iVertex < nVertex[i]; iVertex++) {
        geometry->vertex[valueMarkerWet[i]][iVertex]->SetVarCoord(displacementDeltas_su2[iVertex]);
      }
      if (verbose) {
        cout << "Process #" << processRank << "/" << processSize-1 << ": Advancing preCICE: ...done setting displacement deltas for " << config->GetpreCICE_WetSurfaceMarkerName() << markerLocalToGlobal[i] << "." << endl;
      }
    }

    if (verbose) {
      cout << "Process #" << processRank << "/" << processSize-1 << ": ...done advancing preCICE!" << endl;
    }
    return max_precice_dt;
  }
  else{
    if (verbose) cout << "Process #" << processRank << "/" << processSize-1 << ": Advancing preCICE..." << endl;

    //3. Advance solverInterface
    if (verbose) cout << "Process #" << processRank << "/" << processSize-1 << ": Advancing preCICE: Advancing SolverInterface..." << endl;
    su2double max_precice_dt;
    max_precice_dt = solverInterface.advance( SU2_TYPE::GetValue(computedTimestep) );

    if (verbose) {
      cout << "Process #" << processRank << "/" << processSize-1 << ": Advancing preCICE: ...done advancing SolverInterface." << endl;
      cout << "Process #" << processRank << "/" << processSize-1 << ": ...done advancing preCICE!" << endl;
    }
    return max_precice_dt;
  }

};

void CPreciceFlow::Set_OldState( bool *StopCalc, double *dt ) {

  unsigned long iPoint;
  unsigned short iVar;

  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    for (iVar = 0; iVar < nVar; iVar++) {
      //Save solutions at last and current time step
      solution_Saved[iPoint][iVar] = (solver[FLOW_SOL]->node[iPoint]->GetSolution())[iVar];
      solution_time_n_Saved[iPoint][iVar] = (solver[FLOW_SOL]->node[iPoint]->GetSolution_time_n())[iVar];
      solution_time_n1_Saved[iPoint][iVar] = (solver[FLOW_SOL]->node[iPoint]->GetSolution_time_n1())[iVar];
    }
    for (int iDim = 0; iDim < nDim; iDim++) {
      //Save coordinates at last, current and next time step
      Coord_Saved[iPoint][iDim] = (geometry->node[iPoint]->GetCoord())[iDim];
      Coord_n_Saved[iPoint][iDim] = (geometry->node[iPoint]->GetCoord_n())[iDim];
      Coord_n1_Saved[iPoint][iDim] = (geometry->node[iPoint]->GetCoord_n1())[iDim];
      Coord_p1_Saved[iPoint][iDim] = (geometry->node[iPoint]->GetCoord_p1())[iDim];
      //Save grid velocity
      GridVel_Saved[iPoint][iDim] = (geometry->node[iPoint]->GetGridVel())[iDim];
      for (int jDim = 0; jDim < nDim; jDim++) {
        //Save grid velocity gradient
        GridVel_Grad_Saved[iPoint][iDim][jDim] = (geometry->node[iPoint]->GetGridVel_Grad())[iDim][jDim];
      }
    }
  }

  //Save wether simulation should be stopped after the current iteration
  StopCalc_savedState = *StopCalc;
  //Save the time step size
  dt_savedState = *dt;
  //Writing task has been fulfilled successfully
  solverInterface.fulfilledAction(cowic);

};

void CPreciceFlow::Reset_OldState( bool *StopCalc, double *dt )
{
  unsigned long iPoint;
  unsigned short iDim, jDim;

  for (iPoint = 0; iPoint < nPoint; iPoint++) {

    //Reload solutions at last and current time step
    solver[FLOW_SOL]->node[iPoint]->SetSolution(solution_Saved[iPoint]);
    solver[FLOW_SOL]->node[iPoint]->Set_Solution_time_n(solution_time_n_Saved[iPoint]);
    solver[FLOW_SOL]->node[iPoint]->Set_Solution_time_n1(solution_time_n1_Saved[iPoint]);

    //Reload coordinates at last, current and next time step
    geometry->node[iPoint]->SetCoord(Coord_n1_Saved[iPoint]);
    geometry->node[iPoint]->SetCoord_n();
    geometry->node[iPoint]->SetCoord_n1();
    geometry->node[iPoint]->SetCoord(Coord_n_Saved[iPoint]);
    geometry->node[iPoint]->SetCoord_n();
    geometry->node[iPoint]->SetCoord_p1(Coord_p1_Saved[iPoint]);
    geometry->node[iPoint]->SetCoord(Coord_Saved[iPoint]);

    //Reload grid velocity
    geometry->node[iPoint]->SetGridVel(GridVel_Saved[iPoint]);

    //Reload grid velocity gradient
    for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        geometry->node[iPoint]->SetGridVel_Grad(iDim, jDim, GridVel_Grad_Saved[iPoint][iDim][jDim]);
      }
    }
  }

  //Reload wether simulation should be stopped after current iteration
  *StopCalc = StopCalc_savedState;
  //Reload the time step size
  *dt = dt_savedState;
  //Reading task has been fulfilled successfully
  solverInterface.fulfilledAction(coric);

};

void CPreciceFlow::Finalize()
{
  cout << "Process #" << processRank << "/" << processSize-1 << ": Finalizing preCICE..." << endl;
  solverInterface.finalize();
  cout << "Process #" << processRank << "/" << processSize-1 << ": Done finalizing preCICE!" << endl;
};

