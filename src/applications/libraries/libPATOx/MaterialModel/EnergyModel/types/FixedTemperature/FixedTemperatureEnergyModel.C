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
    along with OpenFOAM.  If FixedTemperaturet, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "FixedTemperatureEnergyModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FixedTemperatureEnergyModel::FixedTemperatureEnergyModel
(
    const fvMesh& mesh,
    const word& dictName
)
  :
simpleEnergyModel(mesh, dictName),
intT_(readScalar(simpleEnergyModel::materialDict_.subDict("Energy").lookup("internalField_T"))),
bcT_(simpleEnergyModel::materialDict_.subDict("Energy").lookup("boundaryField_T")),
T_(meshLookupOrConstructScalar(mesh, "Ta"))
{
  if ( bcT_.size() != mesh.boundaryMesh().size()) {
    FatalErrorInFunction << "boundaryField_T must have all the boundaries. boundaryField_T.size()->"<<
                         bcT_.size() << " != mesh.boundaryMesh().size()->" <<  mesh.boundaryMesh().size()
                         << nl << "Mesh boundary names = " <<  mesh.boundaryMesh().names()
                         << exit(FatalError);
  }
  labelList bcID_;
  forAll(bcT_, bcI) {
    label id_ = mesh.boundaryMesh().findPatchID(bcT_[bcI].first());
    if (id_<0) {
      FatalErrorInFunction << bcT_[bcI].first() << " boundary not found in T." << exit(FatalError);
    } else {
      bcID_.append(id_);
    }
  }
  forAll(T_, cellI) {
    T_[cellI]=intT_;
  }
  forAll(bcID_, bcI) {
    label patchI = bcID_[bcI];
    if(isA<fixedValueFvPatchScalarField>( T_.boundaryField()[patchI])) {
      forAll(mesh.boundaryMesh()[patchI], faceI) {
        T_.boundaryFieldRef()[patchI][faceI] = bcT_[bcI].second();
      }
    }
  }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::FixedTemperatureEnergyModel::~FixedTemperatureEnergyModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::FixedTemperatureEnergyModel::update()
{}

// ************************************************************************* //
