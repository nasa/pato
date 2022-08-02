/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "BoundaryTableEnergyModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::BoundaryTableEnergyModel::BoundaryTableEnergyModel
(
    const fvMesh& mesh,
    const word& dictName
)
  :
simpleEnergyModel(mesh, dictName),
mesh_(mesh),
T_(meshLookupOrConstructScalar(mesh, "Ta")),
T_bf_(T_.boundaryFieldRef()[findTableBoundaryField()])
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::BoundaryTableEnergyModel::~BoundaryTableEnergyModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::BoundaryTableEnergyModel::update()
{
  T_.correctBoundaryConditions();

  forAll(T_,cellI) {
    if (T_bf_.size()>0) {
      T_[cellI] = T_bf_[0];
    }
  }

  T_.correctBoundaryConditions();

}

label Foam::BoundaryTableEnergyModel::findTableBoundaryField()
{
  label bf_ = -1;
  forAll(mesh_.boundaryMesh(), patchI) {
    const fvPatchScalarField& T_BF = T_.boundaryFieldRef()[patchI];

    if(isA<uniformFixedValueFvPatchField<scalar> >(T_BF)) {
      bf_ = patchI;
    }
  }
  if (bf_ < 0 ) {
    FatalErrorInFunction << "uniformFixedValue boundary field not found for Ta."
                         << exit(FatalError) ;
  }
  return bf_;
}

// ************************************************************************* //
