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

#include "Porous_constMaterialPropertiesModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Porous_constMaterialPropertiesModel::Porous_constMaterialPropertiesModel
(
    const fvMesh& mesh,
    const word& regionName
)
  :
simpleMaterialPropertiesModel(mesh, regionName),
CommonMaterialPropertiesModel(*this, modelName,"Porous_const")
{
  listMapFields= {
    newMapField_Porous_const("K","volTensorField",init_tensor,update_empty),
    newMapField_Porous_const("Beta","volTensorField",init_tensor,update_empty),
    newMapField_Porous_const("nu","volScalarField",init_scalar,update_empty),
    newMapField_Porous_const("E","volScalarField",init_scalar,update_empty),
    newMapField_Porous_const("alpha","volScalarField",init_scalar,update_empty),
    newMapField_Porous_const("xi","volScalarField",init_scalar,update_empty),
    newMapField_Porous_const("cp","volScalarField",init_scalar,update_empty),
    newMapField_Porous_const("k","volTensorField",init_tensor,update_empty),
    newMapField_Porous_const("pyrolysisFlux","volScalarField",init_scalar,update_empty),
    newMapField_Porous_const("emissivity","volScalarField",init_scalar,update_empty),
    newMapField_Porous_const("absorptivity","volScalarField",init_scalar,update_empty),
    newMapField_Porous_const("h_c","volScalarField",init_scalar,update_empty),
    newMapField_Porous_const("rho_s","volScalarField",init_scalar,update_empty),
  };
  modelInitialized();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Porous_constMaterialPropertiesModel::~Porous_constMaterialPropertiesModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::Porous_constMaterialPropertiesModel::initialize()
{
  CommonMaterialPropertiesModel::initialize();
}

void Foam::Porous_constMaterialPropertiesModel::update()
{}

void Foam::Porous_constMaterialPropertiesModel::init_scalar(word name)
{
  scalarMatProps_.insert(name+"_const",createScalarProp(name+"_const"));
  volScalarField& f = scalarFields_[name];
  scalar& f_const = scalarMatProps_[name+"_const"];
  forAll(f, cellI) {
    f[cellI] = f_const;
  }
  forAll(f.boundaryField(), patchI) {
    forAll(f.boundaryField()[patchI], faceI) {
      f.boundaryFieldRef()[patchI][faceI] = f_const;
    }
  }
  f.correctBoundaryConditions();
}

void Foam::Porous_constMaterialPropertiesModel::init_tensor(word name)
{
  scalarMatProps_.insert(name+"_const",createScalarProp(name+"_const"));
  volTensorField& f = tensorFields_[name];
  scalar& f_const = scalarMatProps_[name+"_const"];
  forAll(f, cellI) {
    f[cellI] = tensor(f_const,0,0,0,f_const,0,0,0,f_const);
  }
  forAll(f.boundaryField(), patchI) {
    forAll(f.boundaryField()[patchI], faceI) {
      f.boundaryFieldRef()[patchI][faceI] = tensor(f_const,0,0,0,f_const,0,0,0,f_const);
    }
  }
  f.correctBoundaryConditions();
}

void Foam::Porous_constMaterialPropertiesModel::init_k(word name)
{
  init_tensor(name);
  dimensionedTensor kI = dimensionedTensor("kI",dimensionSet(1, 1, -3, -1, 0, 0, 0),I_);
  symmTensorFields_.insert("k_abl_sym",createVolField<symmTensor>("k_abl_sym",symm(kI)));
  volSymmTensorField& k_abl_sym_ = symmTensorFields_["k_abl_sym"];
  k_abl_sym_ = symm(tensorFields_[name]);
  k_abl_sym_.correctBoundaryConditions();
}

void Foam::Porous_constMaterialPropertiesModel::update_empty(word name)
{
  FatalError << "Not implemented: Foam::Porous_constMaterialPropertiesModel::update_empty(" << name << ")" << exit(FatalError);
}

// ************************************************************************* //
