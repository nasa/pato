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
startModelInit_(startModelInit())
{
  Switch found_k=false;
  forAll(listMapFields, mapI) {
    if (listMapFields[mapI].fieldName_=="k") {
      listMapFields.set(mapI,newMapField_Porous_polynomial_k_UQ("k","volTensorField",init_k,update_k));
      found_k=true;
      break;
    }
  }
  if (!found_k) {
    FatalError << "Field k not found in listMapFields" << exit(FatalError);
  }
  modelInitialized();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Porous_polynomial_k_UQMaterialPropertiesModel::~Porous_polynomial_k_UQMaterialPropertiesModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::Porous_polynomial_k_UQMaterialPropertiesModel::init_k(word name)
{

  PorousMaterialPropertiesModel::initVirginCharFiles(name);
  if(name=="k") {
    // Thermal conductivity coefficients: k(T) = kCond + kRad * T^3
    scalarMatProps_.insert("kCondVirgin",createScalarProp("kCondVirgin"));
    scalarMatProps_.insert("kRadVirgin",createScalarProp("kRadVirgin"));
    scalarMatProps_.insert("kCondChar",createScalarProp("kCondChar"));
    scalarMatProps_.insert("kRadChar",createScalarProp("kRadChar"));
  }
}

void Foam::Porous_polynomial_k_UQMaterialPropertiesModel::update_k(word name)
{
  const scalar& kCondVirgin_ = scalarMatProps_["kCondVirgin"];
  const scalar& kRadVirgin_ = scalarMatProps_["kRadVirgin"];
  const scalar& kCondChar_ = scalarMatProps_["kCondChar"];
  const scalar& kRadChar_ = scalarMatProps_["kRadChar"];
  volTensorField& k_ = tensorFields_[name];
  volSymmTensorField& k_abl_sym_ = symmTensorFields_["k_abl_sym"];
  volTensorField& kijk_ = tensorFields_["kijk"];
  const dimensionedScalar& kiCoef_ = dimScalarMatProps_["kiCoef"];
  const dimensionedScalar& kjCoef_ = dimScalarMatProps_["kjCoef"];
  const dimensionedScalar& kkCoef_ = dimScalarMatProps_["kkCoef"];
  const dimensionedTensor& tP_ = dimTensorMatProps_["tP"];

  forAll(k_, i) {
    scalar k_v_ = kCondVirgin_ + kRadVirgin_*pow3(T_[i]);
    scalar k_c_ = kCondChar_ + kRadChar_*pow3(T_[i]);
    k_[i].xx() = kiCoef_.value()
                 *(k_v_*tau_[i]*rho_v_[i]/max(rho_s_[i], rho_c_[i]) + k_c_
                   *(1. - tau_[i]*rho_v_[i]/max(rho_s_[i], rho_c_[i])));
    k_[i].yy() =
        kjCoef_.value()
        *(k_v_*tau_[i]*rho_v_[i]/max(rho_s_[i], rho_c_[i]) + k_c_
          *(1. - tau_[i]*rho_v_[i]/max(rho_s_[i], rho_c_[i])));
    k_[i].zz() =
        kkCoef_.value()
        *(k_v_*tau_[i]*rho_v_[i]/max(rho_s_[i], rho_c_[i]) + k_c_
          *(1. - tau_[i]*rho_v_[i]/max(rho_s_[i], rho_c_[i])));
  }

  forAll(mesh_.boundaryMesh(), patchI) {
    forAll(k_.boundaryFieldRef()[patchI], faceI) {
      scalar k_v_ =
          kCondVirgin_ + kRadVirgin_
          *pow3(T_.boundaryField()[patchI][faceI]);
      scalar k_c_ =
          kCondChar_ + kRadChar_*pow3(T_.boundaryField()[patchI][faceI]);
      k_.boundaryFieldRef()[patchI][faceI].xx() =
          kiCoef_.value()
          *(k_v_*tau_.boundaryField()[patchI][faceI]
            *rho_v_.boundaryField()[patchI][faceI]
            /max
            (
                rho_s_.boundaryField()[patchI][faceI],
                rho_c_.boundaryField()[patchI][faceI]
            )
            + k_c_*(1. - tau_.boundaryField()[patchI][faceI]
                    *rho_v_.boundaryField()[patchI][faceI]
                    /max
                    (
                        rho_s_.boundaryField()[patchI][faceI],
                        rho_c_.boundaryField()[patchI][faceI])
                   ));
      k_.boundaryFieldRef()[patchI][faceI].yy() =
          kjCoef_.value()
          *(k_v_*tau_.boundaryField()[patchI][faceI]
            *rho_v_.boundaryField()[patchI][faceI]
            /max
            (
                rho_s_.boundaryField()[patchI][faceI],
                rho_c_.boundaryField()[patchI][faceI]
            )
            + k_c_ *(1. - tau_.boundaryField()[patchI][faceI]
                     *rho_v_.boundaryField()[patchI][faceI]
                     /max
                     (
                         rho_s_.boundaryField()[patchI][faceI],
                         rho_c_.boundaryField()[patchI][faceI]
                     )
                    ));
      k_.boundaryFieldRef()[patchI][faceI].zz() =
          kkCoef_.value()
          *(k_v_*tau_.boundaryField()[patchI][faceI]
            *rho_v_.boundaryField()[patchI][faceI]
            /max
            (
                rho_s_.boundaryField()[patchI][faceI],
                rho_c_.boundaryField()[patchI][faceI]
            )
            + k_c_*(1. - tau_.boundaryField()[patchI][faceI]
                    *rho_v_.boundaryField()[patchI][faceI]
                    /max
                    (
                        rho_s_.boundaryField()[patchI][faceI],
                        rho_c_.boundaryField()[patchI][faceI]
                    )
                   ));
    }
  }

  /*passage from basis ijk (main directions of the conductivity tensor)
    to basis xyz (mesh)*/
  k_ = tP_.value() & k_ & tP_.value().T();
  k_abl_sym_ = symm(k_);
}


// ************************************************************************* //
