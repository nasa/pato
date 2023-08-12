/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::maximumStressModel

Description

SourceFiles
    maximumStressModel.H
    maximumStressModel.C

\*---------------------------------------------------------------------------*/

#ifndef maximumStressModel_H
#define maximumStressModel_H

#define namespaceFoam
#include "simpleMaterialFailureCriteriaModel.H"
#include "simpleMaterialsModel.H"
#include "IOFunctions_extend.H"
#include "wedgePolyPatch.H"
#include "simpleSolidMechanicsModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class maximumStressModel Declaration
\*---------------------------------------------------------------------------*/

class maximumStressMaterialFailureCriteriaModel
:
public simpleMaterialFailureCriteriaModel
{
 public:

  //- Runtime type information
  TypeName("maximumStress"); // Disable material failure

  // Constructors

  //- Construct from mesh and phase name
  maximumStressMaterialFailureCriteriaModel(const fvMesh& mesh,const word& dictName);

  //- Destructor
  virtual ~maximumStressMaterialFailureCriteriaModel();

  //- Member Functions

  //- Update material failure criteria
  virtual void updateFailureCriteria();

  //- Create the list of boundary types for the failure criteria
  wordList failureCriteriaBC();

  //- Clean the fields
  virtual void cleanFields();

 private:

  //- List of boundary types for the failure criteria
  wordList failureCriteriaBC_;

  // Solid Mechanics Model
  simpleSolidMechanicsModel& solidMechanicsModel_;

  // Failure criteria field
  volScalarField& failureCriteria_;

  // Foam Extend mesh reference
  const Foam_extend_::fvMesh& mesh_extend_;

  // Stress field reference
  const Foam_extend_::volSymmTensorField& sigma_;

  // Failure criteria
  const int SUCCEEDS = 0;
  const int FAILS = 1;

  // Failure modes
  const int NMODES = 7;
  enum mode {
      TENSION_TTT,
      COMPRESSION_TTT,
      TENSION_IP,
      COMPRESSION_IP,
      SHEAR_TTT_TTT,
      SHEAR_IP_IP,
      SHEAR_TTT_IP
  };

  // TTT: throught the thickness
  const double strength_tension_TTT_;
  const double strength_compression_TTT_;
  // IP: in plane
  const double strength_tension_IP_;
  const double strength_compression_IP_;
  // Shear
  const double strength_shear_;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
