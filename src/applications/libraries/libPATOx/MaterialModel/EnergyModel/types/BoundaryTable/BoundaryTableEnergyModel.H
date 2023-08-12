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
    along with OpenFOAM.  If BoundaryTablet, see <http://www.gnu.org/licenses/>.

Class
    Foam::BoundaryTableEnergyModel

Description

SourceFiles
    BoundaryTableEnergyModel.H
    BoundaryTableEnergyModel.C

\*---------------------------------------------------------------------------*/

#ifndef BoundaryTableEnergyModel_H
#define BoundaryTableEnergyModel_H

#include "simpleEnergyModel.H"
#include "uniformFixedValueFvPatchField.H"
#include "boundaryMappingFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class BoundaryTableEnergyModel Declaration
\*---------------------------------------------------------------------------*/

class BoundaryTableEnergyModel
:
public simpleEnergyModel
{

 protected:

  //- Mesh reference
  const fvMesh& mesh_;

  //- Temperature [K]
  volScalarField& T_;

  //- uniformFixedValue boundaryField
  fvPatchScalarField& T_bf_;

 public:

  //- Runtime type information
  TypeName("BoundaryTable");

  // Constructors

  //- Construct from mesh and phase name
  BoundaryTableEnergyModel(const fvMesh& mesh,const word& regionName);

  //- Destructor
  virtual ~BoundaryTableEnergyModel();

  //- Member Functions

  //- Update
  virtual void update();

  //- Find the ID of the patch of type "uniformFixedValueFvPatchField"
  label findTableBoundaryField();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
