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
    along with OpenFOAM.  If Porous_const_k_UQt, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "Porous_const_k_UQMaterialPropertiesModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Porous_const_k_UQMaterialPropertiesModel::Porous_const_k_UQMaterialPropertiesModel
(
    const fvMesh& mesh,
    const word& dictName
)
  :
PorousMaterialPropertiesModel(mesh, dictName),
k_const_v_(readScalar(constantPropertiesDictionary.lookup("k_const_v"))),
k_const_c_(readScalar(constantPropertiesDictionary.lookup("k_const_c")))
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Porous_const_k_UQMaterialPropertiesModel::~Porous_const_k_UQMaterialPropertiesModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::Porous_const_k_UQMaterialPropertiesModel::update()
{
  if(simpleMaterialPropertiesModel::debug_) {
    Info << "--- update K --- Foam::porousMaterialPropertiesModel::update()" << endl;
  }
  // update permeability (used in MassModel)
  K_ = Kc_ + (Kv_ - Kc_) * tau_;

  if(simpleMaterialPropertiesModel::debug_) {
    Info << "--- update cp --- Foam::porousMaterialPropertiesModel::update()" << endl;
  }
  // update solid specific heat (used in EnergyModel)
  PorousMaterialPropertiesModel::updateCp();

  if(simpleMaterialPropertiesModel::debug_) {
    Info << "--- update k --- Foam::porousMaterialPropertiesModel::update()" << endl;
  }
  // update solid thermal conductivity (used in EnergyModel)
  updatek();

  if(simpleMaterialPropertiesModel::debug_) {
    Info << "--- update h_bar --- Foam::porousMaterialPropertiesModel::update()" << endl;
  }
  // update averaged solid enthalpy (used in PyrolysisModel)
  PorousMaterialPropertiesModel::updateHbar();

  if(simpleMaterialPropertiesModel::debug_) {
    Info << "--- update emissivity absorptivity --- Foam::porousMaterialPropertiesModel::update()" << endl;
  }
  // update emissivity and absorptivity
  PorousMaterialPropertiesModel::updateEmissivityAbsorptivity();

  if(simpleMaterialPropertiesModel::debug_) {
    Info << "--- update h_c --- Foam::porousMaterialPropertiesModel::update()" << endl;
  }
  // update charred enthalpy
  PorousMaterialPropertiesModel::updateHc();

  // update solid density
  rho_s_ *= 0;
  forAll(solidRho_, phaseI) {
    rho_s_ += solidRho_[phaseI] * solidEps_[phaseI];
  }
}

void Foam::Porous_const_k_UQMaterialPropertiesModel::updatek()
{
  forAll(T_, i) {
    scalar kConst_ =
        k_const_v_ *
        tau_[i] * rho_v_[i] / max(rho_s_[i], rho_c_[i])
        +
        k_const_c_ *
        (1. - tau_[i] * rho_v_[i] / max(rho_s_[i], rho_c_[i]))
        ;

    kijk_[i] = tensor(kiCoef_.value() * kConst_,0,0,0, kjCoef_.value() * kConst_,0,0,0, kkCoef_.value() * kConst_);
  }
  forAll(T_.boundaryField(), patchI) {
    forAll(T_.boundaryField()[patchI], faceI) {
      scalar kConst_ =
          k_const_v_ *
          tau_.boundaryField()[patchI][faceI] * rho_v_.boundaryField()[patchI][faceI] / max(rho_s_.boundaryField()[patchI][faceI], rho_c_.boundaryField()[patchI][faceI])
          +
          k_const_c_ * (1. - tau_.boundaryField()[patchI][faceI] * rho_v_.boundaryField()[patchI][faceI] / max(rho_s_.boundaryField()[patchI][faceI], rho_c_.boundaryField()[patchI][faceI]))
          ;
      kijk_.boundaryFieldRef()[patchI][faceI] =  tensor(kiCoef_.value() * kConst_,0,0,0, kjCoef_.value() * kConst_,0,0,0, kkCoef_.value() * kConst_);
    }
  }
  kijk_.correctBoundaryConditions();



  k_ = tP_.value() & kijk_ & tP_.value().T(); // passage from basis ijk (main directions of the conductivity tensor) to basis xyz (mesh)
  if (multiRegions) {
    k_abl_sym_ = symm(k_);
  }
}

// ************************************************************************* //
