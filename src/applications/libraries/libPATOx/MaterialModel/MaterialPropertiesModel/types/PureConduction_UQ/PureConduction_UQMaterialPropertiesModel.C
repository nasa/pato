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
    along with OpenFOAM.  If PureConduction_UQ_t, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "PureConduction_UQMaterialPropertiesModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PureConduction_UQMaterialPropertiesModel::PureConduction_UQMaterialPropertiesModel
(
    const fvMesh& mesh,
    const word& dictName
)
  :
simpleMaterialPropertiesModel(mesh, dictName),
materialPropertiesDirectory(fileName(simpleMaterialPropertiesModel::materialDict_.subDict("MaterialProperties").lookup("MaterialPropertiesDirectory")).expand()),
constantPropertiesDictionary
(
    IOobject
    (
        materialPropertiesDirectory+"/constantProperties",
        mesh.time().db().parent(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
),
k_(simpleMaterialPropertiesModel::k_),
cp_(simpleMaterialPropertiesModel::cp_),
rho_s_(simpleMaterialPropertiesModel::rho_s_),
kConst_(readScalar(constantPropertiesDictionary.lookup("k_const"))),
cpConst_(readScalar(constantPropertiesDictionary.lookup("cp_const"))),
rho_sConst_(readScalar(constantPropertiesDictionary.lookup("rho_s_const")))
{
  forAll(cp_, cellI) {
    cp_[cellI] = cpConst_;
    k_[cellI] = tensor(kConst_,0,0,0,kConst_,0,0,0,kConst_);
    rho_s_[cellI] = rho_sConst_;
  }
  forAll(cp_.boundaryField(), patchI) {
    forAll(cp_.boundaryField()[patchI], faceI) {
      cp_.boundaryFieldRef()[patchI][faceI] = cpConst_;
      k_.boundaryFieldRef()[patchI][faceI] = tensor(kConst_,0,0,0,kConst_,0,0,0,kConst_);
      rho_s_.boundaryFieldRef()[patchI][faceI] = rho_sConst_;
    }
  }
  cp_.correctBoundaryConditions();
  k_.correctBoundaryConditions();
  rho_s_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PureConduction_UQMaterialPropertiesModel::~PureConduction_UQMaterialPropertiesModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::PureConduction_UQMaterialPropertiesModel::update()
{}

// ************************************************************************* //
