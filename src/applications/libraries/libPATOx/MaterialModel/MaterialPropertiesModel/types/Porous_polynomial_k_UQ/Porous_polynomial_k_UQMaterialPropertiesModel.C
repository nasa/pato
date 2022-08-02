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

#include "Porous_polynomial_k_UQMaterialPropertiesModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Porous_polynomial_k_UQMaterialPropertiesModel::Porous_polynomial_k_UQMaterialPropertiesModel
(
    const fvMesh& mesh,
    const word& dictName
)
  :
PorousMaterialPropertiesModel(mesh, dictName),
kCondVirgin_(simpleMaterialPropertiesModel::materialDict_.subDict("MaterialProperties").found("kCondVirgin")?readScalar(simpleMaterialPropertiesModel::materialDict_.subDict("MaterialProperties").lookup("kCondVirgin")):readScalar(constantPropertiesDictionary.lookup("kCondVirgin"))),
kCondChar_(simpleMaterialPropertiesModel::materialDict_.subDict("MaterialProperties").found("kCondChar")?readScalar(simpleMaterialPropertiesModel::materialDict_.subDict("MaterialProperties").lookup("kCondChar")):readScalar(constantPropertiesDictionary.lookup("kCondChar"))),
kRadVirgin_(simpleMaterialPropertiesModel::materialDict_.subDict("MaterialProperties").found("kRadVirgin")?readScalar(simpleMaterialPropertiesModel::materialDict_.subDict("MaterialProperties").lookup("kRadVirgin")):readScalar(constantPropertiesDictionary.lookup("kRadVirgin"))),
kRadChar_(simpleMaterialPropertiesModel::materialDict_.subDict("MaterialProperties").found("kRadChar")?readScalar(simpleMaterialPropertiesModel::materialDict_.subDict("MaterialProperties").lookup("kRadChar")):readScalar(constantPropertiesDictionary.lookup("kRadChar")))
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Porous_polynomial_k_UQMaterialPropertiesModel::~Porous_polynomial_k_UQMaterialPropertiesModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::Porous_polynomial_k_UQMaterialPropertiesModel::update()
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

void Foam::Porous_polynomial_k_UQMaterialPropertiesModel::updatek()
{
  //porous_polynomial_k_SA_UQMaterialPropertiesObject& virgin = virginObject();
  //porous_polynomial_k_SA_UQMaterialPropertiesObject& charred = charObject();

  forAll(k_, i) {
    scalar k_v_ =  kCondVirgin_ + kRadVirgin_ * pow3(T_[i]);
    scalar k_c_ =  kCondChar_ + kRadChar_ * pow3(T_[i]);
    k_[i].xx() =  kiCoef_.value() *
                  (
                      k_v_ *
                      tau_[i] * rho_v_[i] / max(rho_s_[i], rho_c_[i])
                      +
                      k_c_ *
                      (1. - tau_[i] * rho_v_[i] / max(rho_s_[i], rho_c_[i]))
                  );
    k_[i].yy() =
        kjCoef_.value() *
        (
            k_v_ *
            tau_[i] * rho_v_[i] / max(rho_s_[i], rho_c_[i])
            +
            k_c_ *
            (1. - tau_[i] * rho_v_[i] / max(rho_s_[i], rho_c_[i]))
        );
    k_[i].zz() =
        kkCoef_.value() *
        (
            k_v_ *
            tau_[i] * rho_v_[i] / max(rho_s_[i], rho_c_[i])
            +
            k_c_ *
            (1. - tau_[i] * rho_v_[i] / max(rho_s_[i], rho_c_[i]))
        );
  }

  forAll(mesh_.boundaryMesh(), patchI) {
    forAll(k_.boundaryFieldRef()[patchI], faceI) {
      scalar k_v_ =  kCondVirgin_ + kRadVirgin_ * pow3(T_.boundaryField()[patchI][faceI]);
      scalar k_c_ =  kCondChar_ + kRadChar_ * pow3(T_.boundaryField()[patchI][faceI]);
      k_.boundaryFieldRef()[patchI][faceI].xx() =
          kiCoef_.value() *
          (
              k_v_ *
              tau_.boundaryField()[patchI][faceI] * rho_v_.boundaryField()[patchI][faceI] / max(rho_s_.boundaryField()[patchI][faceI], rho_c_.boundaryField()[patchI][faceI])
              +
              k_c_ *
              (1. - tau_.boundaryField()[patchI][faceI] * rho_v_.boundaryField()[patchI][faceI] / max(rho_s_.boundaryField()[patchI][faceI], rho_c_.boundaryField()[patchI][faceI]))
          );
      k_.boundaryFieldRef()[patchI][faceI].yy() =
          kjCoef_.value() *
          (
              k_v_ *
              tau_.boundaryField()[patchI][faceI] * rho_v_.boundaryField()[patchI][faceI] / max(rho_s_.boundaryField()[patchI][faceI], rho_c_.boundaryField()[patchI][faceI])
              +
              k_c_ *
              (1. - tau_.boundaryField()[patchI][faceI] * rho_v_.boundaryField()[patchI][faceI] / max(rho_s_.boundaryField()[patchI][faceI], rho_c_.boundaryField()[patchI][faceI]))
          );
      k_.boundaryFieldRef()[patchI][faceI].zz() =
          kkCoef_.value() *
          (
              k_v_ *
              tau_.boundaryField()[patchI][faceI] * rho_v_.boundaryField()[patchI][faceI] / max(rho_s_.boundaryField()[patchI][faceI], rho_c_.boundaryField()[patchI][faceI])
              +
              k_c_ *
              (1. - tau_.boundaryField()[patchI][faceI] * rho_v_.boundaryField()[patchI][faceI] / max(rho_s_.boundaryField()[patchI][faceI], rho_c_.boundaryField()[patchI][faceI]))
          );
    }
  }


  k_ = tP_.value() & k_ & tP_.value().T(); // passage from basis ijk (main directions of the conductivity tensor) to basis xyz (mesh)
  k_abl_sym_ = symm(k_);
}


// ************************************************************************* //
