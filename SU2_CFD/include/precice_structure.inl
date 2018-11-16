/*!
 * \file precice_structure.hpp
 * \brief Inline routine of the structure to couple SU2 via preCICE.
 *        Based on the precice adaptor file (precice.hpp): https://github.com/precice/su2-adapter.git
 * \author R. Sanchez, based on the work by Alexander Rusch and the preCICE community
 * \version 6.1.0 "Falcon"
 */

#pragma once


inline bool CPrecice::ActiveCoupling() { return solverInterface.isCouplingOngoing(); }

inline bool CPrecice::ActionRequired(const string& action) { return solverInterface.isActionRequired(action); }

inline const string& CPrecice::getCowic() { return cowic; }

inline const string& CPrecice::getCoric() { return coric; }

inline void CPrecice::Configure( const string& configurationFilename ) { }

inline su2double CPrecice::Initialize() { return 0.0; }

inline su2double CPrecice::Advance( su2double computedTimestep ) { return 0.0; }

inline void CPrecice::Set_OldState( bool *StopCalc, double *dt ) { }

inline void CPrecice::Reset_OldState( bool *StopCalc, double *dt ) { }

inline void CPrecice::Finalize() { }
