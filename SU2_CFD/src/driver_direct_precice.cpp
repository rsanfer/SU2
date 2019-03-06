/*!
 * \file driver_direct_precice.cpp
 * \brief The main subroutines for driving problems using preCICE.
 * \author R. Sanchez
 * \version 6.1.0 "Falcon"
 *
 */

#include "../include/driver_structure.hpp"
#include "../include/definition_structure.hpp"
#include "../include/precice_structure.hpp"

#ifdef VTUNEPROF
#include <ittnotify.h>
#endif

CPreciceDriver::CPreciceDriver(char* confFile, unsigned short val_nZone,
                               unsigned short val_nDim, bool val_periodic,
                               SU2_Comm MPICommunicator) : CSinglezoneDriver(confFile,
                                                                   val_nZone,
                                                                   val_nDim,
                                                                   val_periodic,
                                                                   MPICommunicator) {

  /*--- Initialize the precice object depending on the kind of problem at hand ---*/
  switch (config_container[ZONE_0]->GetpreCICE_Subproblem()){
  case PRECICE_FLOW:
    precice = new CPreciceFlow(rank, size, geometry_container, solver_container, config_container, grid_movement);
    /*--- Load time step ---*/
    dt = new double(config_container[ZONE_0]->GetDelta_UnstTimeND());
    break;
  case PRECICE_FEA:
    precice = new CPreciceFEA(rank, size, geometry_container, solver_container, config_container, grid_movement);
    /*--- Load time step ---*/
    dt = new double(config_container[ZONE_0]->GetDelta_DynTime());
    break;
  default:
    precice = new CPreciceFlow(rank, size, geometry_container, solver_container, config_container, grid_movement);
    /*--- Load time step ---*/
    dt = new double(config_container[ZONE_0]->GetDelta_UnstTimeND());
    break;
  }

  /*--- Load time step ---*/
  

  /*--- Initialize preCICE ---*/
  max_precice_dt = new double(precice->Initialize());

  //DynamicMeshUpdate(0);

}

CPreciceDriver::~CPreciceDriver(void) { }

void CPreciceDriver::StartSolver(){

#ifdef VTUNEPROF
  __itt_resume();
#endif

  /*--- Main external loop of the solver. Within this loop, each iteration ---*/

  if (rank == MASTER_NODE)
    cout << endl <<"-------------------------- Begin preCICE Solver -------------------------" << endl;

  while ((ExtIter < config_container[ZONE_0]->GetnExtIter()) && precice->ActiveCoupling()) {

    /*--- Store the old state for implicit coupling ---*/
    if(precice->ActionRequired(precice->getCowic()))
      precice->Set_OldState(&StopCalc, dt);

    /*--- Set the time step ---*/
    dt = min(max_precice_dt,dt);
    config_container[ZONE_0]->SetDelta_UnstTimeND(*dt);

    /*--- Perform some external iteration preprocessing. ---*/
    PreprocessExtIter(ExtIter);

    /*--- Perform a dynamic mesh update if required. ---*/
    if (!fem_solver) {
      DynamicMeshUpdate(ExtIter);
    }

    /*--- Run a single iteration of the problem. ---*/

    Run();

    /*--- Terminate the simulation if only the Jacobian must be computed. ---*/
    if (config_container[ZONE_0]->GetJacobian_Spatial_Discretization_Only()) break;

    /*--- Advance preCICE ---*/
    *max_precice_dt = precice->Advance(*dt);

    bool output_solution = true;

    /*--- Update the solution for dual time stepping strategy ---*/	
    Update();

    /*--- Monitor the computations after each iteration. ---*/
    Monitor(ExtIter);

    /*--- Test if the preCICE is converged ---*/
    if(precice->ActionRequired(precice->getCoric())){
      /*--- If unconverged, reset to the old state ---*/
      precice->Reset_OldState(&StopCalc, dt);
      output_solution = false;
    }
    else{
    ExtIter++;
    }

    /*--- Output the solution files. ---*/
    if (output_solution) Output(ExtIter);

    /*--- If the convergence criteria has been met, terminate the simulation. ---*/
    if (StopCalc) break;
  }
#ifdef VTUNEPROF
  __itt_pause();
#endif
}

void CPreciceDriver::Finalize(){

  if(precice_usage){
    precice->Finalize();
    if (dt != NULL) {
        delete dt;
    }
    if (max_precice_dt != NULL) {
        delete max_precice_dt;
    }
    if (precice != NULL) {
        delete precice;
    }
  }

}
