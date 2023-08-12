/*

Class
    Foam::LinearInterpolationMaterialPropertiesObject

Description
    LinearInterpolationMaterialPropertiesObject interpolation class

\*---------------------------------------------------------------------------*/

#ifndef LinearInterpolationMaterialPropertiesObject_H
#define LinearInterpolationMaterialPropertiesObject_H

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
                            Bprime class Declaration
\*---------------------------------------------------------------------------*/

class LinearInterpolationMaterialPropertiesObject
{

 protected:

  RectangularMatrix<scalar> materialPropertyTable_; // materialPropertyTable itself
  int indexPiMin_;                                  // First index of Pi : Pi <= P
  int indexPiMax_;                                  // Last index of Pi : Pi <= P
  int indexPiPlus1Min_;                             // First index of P_{i+1} : P_{i+1} >= P
  int indexPiPlus1Max_;                             // Last index of P_{i+1} : P_{i+1} >= P
  scalar f_p_;                                      // Interpolation coefficient for the pression
  int indexTi_;                                     // Index of Ti : Ti <= T
  int indexTiPlus1_;                                // Index of T_{i+1} : T_{i+1} >= T
  scalar f_T_;                                      // Interpolation coefficient for T
  scalar cp_;                                       // calorific capacity
  scalar h_;                                        // enthalpy
  scalar ki_;                                       // conductivity - direction i
  scalar kj_;                                       // conductivity - direction j
  scalar kk_;                                       // conductivity - direction k
  scalar eps_;                                      // emissivity
  scalar alpha_;                                    // absorptivity
  scalar hs_;                                       // sensible enthalpy
  scalar T_temp_;                                   // temporarily stores the last interpolation temperature (to avoid re-interpolating for each variable)
  scalar p_temp_;                                   // temporarily stores the last interpolation pressure (to avoid re-interpolating for each variable)
  scalar small_;                                    // small number

 public:

  // Constructors
  LinearInterpolationMaterialPropertiesObject();
  //- Construct from components
  LinearInterpolationMaterialPropertiesObject(RectangularMatrix<scalar>& materialPropertyTable);

  // Destructor
  virtual ~LinearInterpolationMaterialPropertiesObject();

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
    return materialPropertyTable_[indexTi_][1];
  }
  scalar TiPlus1() const {
    return materialPropertyTable_[indexTiPlus1_][1];
  }
  scalar cp() const {
    return cp_;
  }
  scalar h() const {
    return h_;
  }
  scalar ki() const {
    return ki_;
  }
  scalar kj() const {
    return kj_;
  }
  scalar kk() const {
    return kk_;
  }
  scalar eps() const {
    return eps_;
  }
  scalar alpha() const {
    return alpha_;
  }
  scalar hs() const {
    return hs_;
  }

  // Member Functions
  virtual void update(scalar p_, scalar T_); // update all data

  // return the Temperature interpolation factor and updates the indexes of Ti and Ti+1
  scalar indexTi(scalar T_, int i_min_, int i_max_);

  // return the pressure interpolation factor and updates the indexes of the pressure subtables
  scalar indexPi(scalar p_);
  // return the values of the variables
  scalar M(scalar p_, scalar T_);
  scalar cp(scalar p_, scalar T_);
  scalar h(scalar p_, scalar T_);
  scalar ki(scalar p_, scalar T_);
  scalar kj(scalar p_, scalar T_);
  scalar kk(scalar p_, scalar T_);
  scalar eps(scalar p_, scalar T_);
  scalar alpha(scalar p_, scalar T_);
  scalar hs(scalar p_, scalar T_);
}; // End class LinearInterpolationMaterialPropertiesObject
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
