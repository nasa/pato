/*

Class
    Foam::LinearInterpolationMaterialPropertiesObject

Description
    LinearInterpolationMaterialPropertiesObject interpolation class

\*---------------------------------------------------------------------------*/

#ifndef LinearInterpolation_Factor_MaterialPropertiesObject_H
#define LinearInterpolation_Factor_MaterialPropertiesObject_H

#include "LinearInterpolationMaterialPropertiesObject.H"
#include "List.H"
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

class LinearInterpolation_Factor_MaterialPropertiesObject : public LinearInterpolationMaterialPropertiesObject
{

 private:

  List<scalar> factors_;

 public:

  // Constructors
  LinearInterpolation_Factor_MaterialPropertiesObject();

  //- Construct from Parent
  LinearInterpolation_Factor_MaterialPropertiesObject(const LinearInterpolationMaterialPropertiesObject& obj, List<scalar>& factors);

  // Destructor
  virtual ~LinearInterpolation_Factor_MaterialPropertiesObject();

  // Member Functions
  virtual void update(scalar p_, scalar T_); // update all data

}; // End class LinearInterpolation_Factor_MaterialPropertiesObject
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
