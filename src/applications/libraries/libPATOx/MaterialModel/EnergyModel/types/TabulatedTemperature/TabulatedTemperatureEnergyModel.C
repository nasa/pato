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

#include "TabulatedTemperatureEnergyModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::TabulatedTemperatureEnergyModel::TabulatedTemperatureEnergyModel
(
    const fvMesh& mesh,
    const word& dictName
)
  :
simpleEnergyModel(mesh, dictName),
originPatchName_(simpleEnergyModel::materialDict_.subDict("Energy").lookup("originPatchName")),
intTable_fileName_(simpleEnergyModel::materialDict_.subDict("Energy").lookup("internalFieldTable_fileName")),
intTable_data_(readFileData(intTable_fileName_)),
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
  int originPatchID = mesh.boundaryMesh().findPatchID(originPatchName_);
  if (originPatchID<0) {
    FatalErrorInFunction << originPatchName_ << " patch not found." << exit(FatalError);
  }
  List<scalar> distance_originPatch;
  forAll(T_, cellI) {
    scalar distance = -1;
    forAll(mesh.boundaryMesh()[originPatchID], faceI) {
      scalar d = dist(mesh.Cf().boundaryField()[originPatchID][faceI], mesh.C()[cellI]);
      if ( d < distance || distance < 0) {
        distance = d;
      }
    }
    distance_originPatch.append(distance);
  }
  if (intTable_data_.size()!=2) {
    FatalErrorInFunction <<  intTable_fileName_ << " table size is different than 2." << exit(FatalError);
  }
  forAll(T_, cellI) {
    T_[cellI]=linearInterpolation(intTable_data_[0], intTable_data_[1], distance_originPatch[cellI]);
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

Foam::TabulatedTemperatureEnergyModel::~TabulatedTemperatureEnergyModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::TabulatedTemperatureEnergyModel::update()
{
  Info <<  "time = " << this->mesh_.time().value() << " C()[0]=" << this->mesh_.C()[0] << " C()[1]=" << this->mesh_.C()[1]  << " Cf()[0]=" << this->mesh_.Cf()[0] << endl;
}

// ************************************************************************* //
