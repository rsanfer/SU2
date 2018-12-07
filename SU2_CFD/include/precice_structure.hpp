/*!
 * \file precice_structure.hpp
 * \brief Headers of the structure to couple SU2 via precice.
 *        Based on the precice adaptor file (precice.hpp): https://github.com/precice/su2-adapter.git
 * \author R. Sanchez, based on the work by Alexander Rusch and the preCICE community
 * \version 6.1.0 "Falcon"
 */

#pragma once

#include "SolverInterface.hpp"
#include "SU2_CFD.hpp"
#include <string>
#include <stdlib.h>
using namespace precice;
using namespace std;


/*!
 * \class CPrecice
 * \brief Parent class for preCICE coupling.
 *        There will be a child class for each particular solver (Euler, Navier-Stokes, etc.)
 * \author R. Sanchez, based on the work by Alexander Rusch and the preCICE community
 */

class CPrecice {

protected:

  int processRank, processSize;         /*--- Store the current rank and the total size of the parallelization ---*/
  SolverInterface solverInterface;      /*--- Coupling object ---*/

  CGeometry* geometry;                  /*--- Stores the geometry of the problem ---*/
  CSolver** solver;                     /*--- Stores the current solution of the problem ---*/
  CConfig* config;                      /*--- Stores the configuration of the problem ---*/
  CVolumetricMovement* grid_movement;   /*--- Stores the grid movement container of the current problem ---*/

  unsigned short nDim, nVar;            /*--- Store problem dimensionality ---*/
  unsigned long nPoint;                 /*--- Store number of points in the discretization ---*/

  int *nVertex;              /*--- Store number of vertices per marker in the wet surface ---*/
  short *valueMarkerWet;      /*--- ---*/

  const string& coric;
  const string& cowic;

  bool activeProcess;          /*--- Determines if the process is active in the current surface ---*/
  unsigned long nWetSurfaces,           /*--- Number of wet surfaces in the problem ---*/
                nWetSurfacesDomain;     /*--- Number of wet surfaces in the domain ---*/

  unsigned long* markerLocalToGlobal;            /*--- Mapping of the wet marker from local to global ---*/

  bool verbose;                         /*--- Verbose output to screen ---*/

public:

  /*!
  * \brief Constructor of the class CPrecice.
  * \param[in] processRank - Index of this solver Process.
  * \param[in] processSize - Overall number of solver Processes.
  */
  CPrecice( int processRank, int processSize, CGeometry**** geometry_container, CSolver***** solver_container, CConfig** config_container, CVolumetricMovement*** grid_movement );

  /*!
  * \brief Destructor of the class CPrecice.
  */
  ~CPrecice();

  /*!
  * \brief Check whether preCICE coupling is ongoing.
  * \return State of coupling (True=ongoing, False=completed)
  */
  bool ActiveCoupling();

  /*!
  * \brief Check whether preCICE requires an action (read and/or write checkpoint).
  * \return State of coupling (True=ongoing, False=completed)
  */
  bool ActionRequired(const string& action);

  /*!
  * \brief Determine if an iteration checkpoint needs to be written.
  * \return Coric string.
  */
  const string& getCowic();

  /*!
  * \brief Determine if an iteration checkpoint needs to be read.
  * \return Coric string.
  */
  const string& getCoric();

  /*!
  * \brief Configures preCICE from the given xml file.
  * \param[in] configurationFilename - Name (with path) of the xml configuration file.
  */
  void Configure( const string& configurationFilename );

  /*!
  * \brief Finalizes preCICE.
  */
  void Finalize();

  /*!
  * \brief Initialize preCICE.
  * \return Maximum length of first timestep to be computed by the solver.
  */
  virtual su2double Initialize();

  /*!
  * \brief Advances preCICE after the solver has computed one timestep.
  * \param[in] computedTimestep - Timestep computed by solver.
  * \return Maximum length of next timestep to be computed by solver.
  */
  virtual su2double Advance( su2double computedTimestep );

  /*!
  * \brief Store the current state.
  * \param[in] StopCalc - Determines if the SU2 simulation is completed (convergence or max number of iterations).
  * \param[in] dt - Current time step size.
  */
  virtual void Set_OldState( bool *StopCalc, double *dt );

  /*!
  * \brief Reload the old state.
  * \param[in] StopCalc - Determines if the SU2 simulation is completed (convergence or max number of iterations).
  * \param[in] dt - Current time step size.
  */

  virtual void Reset_OldState( bool *StopCalc, double *dt );

};

/*!
 * \class CPreciceFlow
 * \brief Class for preCICE coupling of flow sub-problems.
 * \author R. Sanchez, based on the work by Alexander Rusch and the preCICE community
 */

class CPreciceFlow : public CPrecice {

protected:


  int **vertexIDs;
  int *forceID;
  int *displDeltaID;
  int *meshID;

  su2double **Coord_Saved,
            **Coord_n_Saved,
            **Coord_n1_Saved,
            **Coord_p1_Saved,
            **GridVel_Saved,
            ***GridVel_Grad_Saved;

  su2double **solution_Saved,
            **solution_time_n_Saved,
            **solution_time_n1_Saved;

  su2double dt_savedState;
  bool StopCalc_savedState;

  /*--- Computation of the redimensionalization ---*/
  su2double* Velocity_Real;
  su2double* Velocity_ND;
  su2double Density_Real,
            Density_ND,
            Velocity2_Real,
            Velocity2_ND;
  su2double factorForces;

public:

  /*!
  * \brief Constructor of the class CPreciceFlow.
  * \param[in] processRank - Index of this solver Process.
  * \param[in] processSize - Overall number of solver Processes.
  */
  CPreciceFlow( int processRank, int processSize, CGeometry**** geometry_container, CSolver***** solver_container, CConfig** config_container, CVolumetricMovement*** grid_movement );

  /*!
  * \brief Destructor of the class CPreciceFlow.
  */
  ~CPreciceFlow();

  /*!
  * \brief Initialize preCICE.
  * \return Maximum length of first timestep to be computed by the solver.
  */
  su2double Initialize();

  /*!
  * \brief Advances preCICE after the solver has computed one timestep.
  * \param[in] computedTimestep - Timestep computed by solver.
  * \return Maximum length of next timestep to be computed by solver.
  */
  su2double Advance( su2double computedTimestep );

  /*!
  * \brief Store the current state.
  * \param[in] StopCalc - Determines if the SU2 simulation is completed (convergence or max number of iterations).
  * \param[in] dt - Current time step size.
  */
  void Set_OldState( bool *StopCalc, double *dt );

  /*!
  * \brief Reload the old state.
  * \param[in] StopCalc - Determines if the SU2 simulation is completed (convergence or max number of iterations).
  * \param[in] dt - Current time step size.
  */

  void Reset_OldState( bool *StopCalc, double *dt );


};

/*!
 * \class CPreciceFEA
 * \brief Class for preCICE coupling of solid FEA sub-problems.
 * \author R. Sanchez
 */

class CPreciceFEA : public CPrecice {

protected:


  int **vertexIDs;
  int *forceID;
  int *displDeltaID;
  int *meshID;

  su2double **Coord_Saved,
            **Coord_n_Saved,
            **Coord_n1_Saved,
            **Coord_p1_Saved,
            **GridVel_Saved,
            ***GridVel_Grad_Saved;

  su2double **solution_Saved,
            **solution_time_n_Saved,
            **solution_time_n1_Saved;

  su2double dt_savedState;
  bool StopCalc_savedState;

  su2double* DisplacementDonor;
  su2double* DisplacementDonor_Prev;
  

public:

  /*!
  * \brief Constructor of the class CPreciceFEA.
  * \param[in] processRank - Index of this solver Process.
  * \param[in] processSize - Overall number of solver Processes.
  */
  CPreciceFEA( int processRank, int processSize, CGeometry**** geometry_container, CSolver***** solver_container, CConfig** config_container, CVolumetricMovement*** grid_movement );

  /*!
  * \brief Destructor of the class CPreciceFEA.
  */
  ~CPreciceFEA();

  /*!
  * \brief Initialize preCICE.
  * \return Maximum length of first timestep to be computed by the solver.
  */
  su2double Initialize();

  /*!
  * \brief Advances preCICE after the solver has computed one timestep.
  * \param[in] computedTimestep - Timestep computed by solver.
  * \return Maximum length of next timestep to be computed by solver.
  */
  su2double Advance( su2double computedTimestep );

  /*!
  * \brief Store the current state.
  * \param[in] StopCalc - Determines if the SU2 simulation is completed (convergence or max number of iterations).
  * \param[in] dt - Current time step size.
  */
  void Set_OldState( bool *StopCalc, double *dt );

  /*!
  * \brief Reload the old state.
  * \param[in] StopCalc - Determines if the SU2 simulation is completed (convergence or max number of iterations).
  * \param[in] dt - Current time step size.
  */

  void Reset_OldState( bool *StopCalc, double *dt );


};

#include "precice_structure.inl"
