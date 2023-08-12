/*---------------------------------------------------------------------------*\

Class
    Foam::TabulatedGasPropertiesObject

Description
    TabulatedGasPropertiesObject interpolation class

\*---------------------------------------------------------------------------*/

#ifndef TabulatedGasPropertiesObject_H
#define TabulatedGasPropertiesObject_H

#include "scalar.H"
#include "RectangularMatrix.H"
// #include <vector>

//#include "typeInfo.H"
//#include "autoPtr.H"
//#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       TabulatedGasPropertiesObject class Declaration
\*---------------------------------------------------------------------------*/

class TabulatedGasPropertiesObject
{

// private

  RectangularMatrix<scalar> gasPropertyTable_; // gasPropertyTable itself
  int indexPiMin_;                             // First index of Pi : Pi <= P
  int indexPiMax_;                             // Last index of Pi : Pi <= P
  int indexPiPlus1Min_;                        // First index of P_{i+1} : P_{i+1} >= P
  int indexPiPlus1Max_;                        // Last index of P_{i+1} : P_{i+1} >= P
  scalar f_p_;                                 // Interpolation coefficient for the pression
  int indexTi_;                                // Index of Ti : Ti <= T
  int indexTiPlus1_;                           // Index of T_{i+1} : T_{i+1} >= T
  scalar f_T_;                                 // Interpolation coefficient for T
  scalar M_;                                   // M   (molar mass)
  // scalar cp_;                                  // cpg (gas calorific capacity)
  // scalar gamma_;                               // gamma   (cp/cv)
  scalar h_;                                   // hg  (gas enthalpy)
  scalar mu_;                                  // mu  (viscosity)
  // scalar rho_;                                 // rho (density)
  scalar T_temp_;                              // temporarily stores the last interpolation temperature (to avoid re-interpolating for each variable)
  scalar p_temp_;                              // temporarily stores the last interpolation pressure (to avoid re-interpolating for each variable)
  scalar small_;                               // small number

 public:

  // Constructors
  TabulatedGasPropertiesObject();

  //- Construct from components
  TabulatedGasPropertiesObject(RectangularMatrix<scalar>& gasPropertyTable);


  // Destructor

  virtual ~TabulatedGasPropertiesObject();

  // Access Functions
  int indexPiMin() const {
    return indexPiMin_;
  }
  int indexPiMax() const {
    return indexPiMax_;
  }
  int indexPiPlus1Min() const {
    return indexPiPlus1Min_;
  }
  int indexPiPlus1Max() const {
    return indexPiPlus1Max_;
  }
  int indexTi() const {
    return indexTi_;
  }
  int indexTiPlus1() const {
    return indexTiPlus1_;
  }
  scalar Ti() const {
    return gasPropertyTable_[indexTi_][1];
  }
  scalar TiPlus1() const {
    return gasPropertyTable_[indexTiPlus1_][1];
  }
  scalar M() const {
    return M_;
  }
  // scalar cp() const { return cp_; }
  // scalar gamma() const { return gamma_; }
  scalar h() const {
    return h_;
  }
  scalar mu() const {
    return mu_;
  }
  // scalar rho() const { return rho_; }

  // Member Functions
  // returns the Temperature interpolation factor and updates the indexes of Ti and Ti+1
  scalar indexTi(scalar T_, int i_min_, int i_max_);
  // returns the pressure interpolation factor and updates the indexes of the pressure subtables
  scalar indexPi(scalar p_);
  // returns M, molar mass (in kg/m³)
  scalar M(scalar p_, scalar T_);
  // scalar cp(scalar p_, scalar T_);
  // scalar gamma(scalar p_, scalar T_);
  scalar h(scalar p_, scalar T_);
  scalar mu(scalar p_, scalar T_);
  // scalar rho(scalar p_, scalar T_);


}; // End class TabulatedGasPropertiesObject
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
