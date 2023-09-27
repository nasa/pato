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
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "FixedPressureGammaMassModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FixedPressureGammaMassModel::FixedPressureGammaMassModel
(
    const fvMesh& mesh,
    const word& dictName
)
  :
simpleMassModel(mesh, dictName),
intP_(readScalar(simpleMassModel::materialDict_.subDict("Mass").lookup("internalField_P"))),
bcP_(simpleMassModel::materialDict_.subDict("Mass").lookup("boundaryField_P")),
materialChemistryModel_(refModel<simpleMaterialChemistryModel>()),
P_(createVolFieldIfNotFound<scalar>(materialChemistryModel_,"p","yes")),
intGamma_(readScalar(simpleMassModel::materialDict_.subDict("Mass").lookup("internalField_Gamma"))),
bcGamma_(simpleMassModel::materialDict_.subDict("Mass").lookup("boundaryField_Gamma")),
Gamma(createVolField<tensor>("Gamma",dimensionedTensor("0", dimTime, tensor(intGamma_,0,0,0,intGamma_,0,0,0,intGamma_)))),
Gamma_symm(createVolField<symmTensor>("Gamma_symm",symm(Gamma)))
{
  forAll(P_, cellI) {
    P_[cellI]=intP_;
  }
  if ( bcP_.size() != mesh.boundaryMesh().size()) {
    FatalErrorInFunction << "boundaryField_P must have all the boundaries. boundaryField_P.size()->"<<
                         bcP_.size() << " != mesh.boundaryMesh().size()->" <<  mesh.boundaryMesh().size()
                         << nl << "Mesh boundary names = " <<  mesh.boundaryMesh().names()
                         << exit(FatalError);
  }
  labelList bcID_;
  forAll(bcP_, bcI) {
    label id_ = mesh.boundaryMesh().findPatchID(bcP_[bcI].first());
    if (id_<0) {
      FatalErrorInFunction << bcP_[bcI].first() << " boundary not found in T." << exit(FatalError);
    } else {
      bcID_.append(id_);
    }
  }
  forAll(bcID_, bcI) {
    label patchI = bcID_[bcI];
    if(isA<fixedValueFvPatchScalarField>( P_.boundaryField()[patchI])) {
      forAll(mesh.boundaryMesh()[patchI], faceI) {
        P_.boundaryFieldRef()[patchI][faceI] = bcP_[bcI].second();
      }
    }
  }

  labelList bcID2_;

  forAll(bcGamma_, bcI) {
    label id_ = mesh.boundaryMesh().findPatchID(bcP_[bcI].first());
    if (id_<0) {
      FatalErrorInFunction << bcP_[bcI].first() << " boundary not found in T." << exit(FatalError);
    } else {
      bcID2_.append(id_);
    }
  }
  forAll(bcID2_, bcI) {
    label patchI = bcID2_[bcI];
    if(isA<fixedValueFvPatchScalarField>( Gamma.boundaryField()[patchI])) {
      forAll(mesh.boundaryMesh()[patchI], faceI) {
        Gamma.boundaryFieldRef()[patchI][faceI] = bcGamma_[bcI].second();
      }
    }
  }
  Gamma_symm=symm(Gamma);
  modelInitialized();
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::FixedPressureGammaMassModel::~FixedPressureGammaMassModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::FixedPressureGammaMassModel::update()
{}

// ************************************************************************* //
