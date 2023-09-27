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
    const word& regionName
)
  :
PorousMaterialPropertiesModel(mesh, regionName),
startModelInit_(startModelInit())
{
  Switch found_k=false;
  forAll(listMapFields, mapI) {
    if (listMapFields[mapI].fieldName_=="k") {
      listMapFields.set(mapI,newMapField_Porous_const_k_UQ("k","volTensorField",init_k,update_k));
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

Foam::Porous_const_k_UQMaterialPropertiesModel::~Porous_const_k_UQMaterialPropertiesModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::Porous_const_k_UQMaterialPropertiesModel::init_k(word name)
{
  PorousMaterialPropertiesModel::initVirginCharFiles(name);
  if(name=="k") {
    //- Constant virgin thermal conductivity (UQ) [W/m/K]
    scalarMatProps_.insert("k_const_v",createScalarProp("k_const_v"));
    //- Constant charred thermal conductivity (UQ) [W/m/K]
    scalarMatProps_.insert("k_const_c",createScalarProp("k_const_c"));
  }
}

void Foam::Porous_const_k_UQMaterialPropertiesModel::update_k(word name)
{
  const scalar& k_const_v_ = scalarMatProps_["k_const_v"];
  const scalar& k_const_c_ = scalarMatProps_["k_const_c"];
  volTensorField& k_ = tensorFields_[name];
  volSymmTensorField& k_abl_sym_ = symmTensorFields_["k_abl_sym"];
  volTensorField& kijk_ = tensorFields_["kijk"];
  const dimensionedScalar& kiCoef_ = dimScalarMatProps_["kiCoef"];
  const dimensionedScalar& kjCoef_ = dimScalarMatProps_["kjCoef"];
  const dimensionedScalar& kkCoef_ = dimScalarMatProps_["kkCoef"];
  const dimensionedTensor& tP_ = dimTensorMatProps_["tP"];

  forAll(T_, i) {
    scalar kConst_ =
        k_const_v_
        *tau_[i]*rho_v_[i]/max(rho_s_[i], rho_c_[i])
        + k_const_c_
        *(1. - tau_[i]*rho_v_[i]/max(rho_s_[i], rho_c_[i]));

    kijk_[i] =
        tensor
        (
            kiCoef_.value()*kConst_,0,0,
            0, kjCoef_.value()*kConst_,0,
            0,0, kkCoef_.value()*kConst_
        );
  }
  forAll(T_.boundaryField(), patchI) {
    forAll(T_.boundaryField()[patchI], faceI) {
      scalar kConst_ =
          k_const_v_
          *tau_.boundaryField()[patchI][faceI]
          *rho_v_.boundaryField()[patchI][faceI]
          /max
          (
              rho_s_.boundaryField()[patchI][faceI],
              rho_c_.boundaryField()[patchI][faceI]
          )
          + k_const_c_
          *(1. - tau_.boundaryField()[patchI][faceI]
            *rho_v_.boundaryField()[patchI][faceI]
            /max
            (
                rho_s_.boundaryField()[patchI][faceI],
                rho_c_.boundaryField()[patchI][faceI])
           );
      kijk_.boundaryFieldRef()[patchI][faceI] =
          tensor
          (
              kiCoef_.value()*kConst_,0,0,
              0, kjCoef_.value()*kConst_,0,
              0,0, kkCoef_.value()*kConst_
          );
    }
  }
  kijk_.correctBoundaryConditions();

  /*passage from basis ijk (main directions of the conductivity tensor)
    to basis xyz (mesh)*/
  k_ = tP_.value() & kijk_ & tP_.value().T();
  if (multiRegions) {
    k_abl_sym_ = symm(k_);
  }
}

// ************************************************************************* //
