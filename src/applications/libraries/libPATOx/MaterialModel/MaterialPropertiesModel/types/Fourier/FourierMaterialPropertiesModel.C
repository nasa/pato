/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2020      PATO
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

#include "FourierMaterialPropertiesModel.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FourierMaterialPropertiesModel::FourierMaterialPropertiesModel
(
    const fvMesh& mesh,
    const word& regionName
)
  :
simpleMaterialPropertiesModel(mesh,regionName),
CommonMaterialPropertiesModel(*this,modelName,"Fourier"),
nCoeffs_(createScalarProp("nCoeffs","yes",5)),
energyModel_(refModel<simpleEnergyModel>()),
T_(energyModel_.refVolField<scalar>("Ta")),
fourierPropDict_(
    IOobject
    (
        materialDirectory_+"/FourierProperties",
        mesh.time().db().parent(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    )
)
{
  listMapFields= {
    newMapField("rho_s","volScalarField",initCoeffs,updateScalarField),
    newMapField("cp","volScalarField",initCoeffs,updateScalarField),
    newMapField("nu","volScalarField",initCoeffs,updateScalarField),
    newMapField("E","volScalarField",initCoeffs,updateScalarField),
    newMapField("alpha","volScalarField",initCoeffs,updateScalarField),
    newMapField("xi","volScalarField",initCoeffs,updateScalarField),
    newMapField("k","volTensorField",init_k_coeffs,update_k),
  };
  // Transfer fourierPropDict_ data to constantPropertiesDictionary_
  wordList fourierPropDictToc_ = fourierPropDict_.toc();
  forAll(fourierPropDictToc_, i) {
    constantPropertiesDictionary_.add(fourierPropDictToc_[i],fourierPropDict_.lookup(fourierPropDictToc_[i]));
  }
  modelInitialized();
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::FourierMaterialPropertiesModel::~FourierMaterialPropertiesModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::FourierMaterialPropertiesModel::update()
{
  CommonMaterialPropertiesModel::update();
}

void Foam::FourierMaterialPropertiesModel::initialize()
{
  CommonMaterialPropertiesModel::initialize();
}

void Foam::FourierMaterialPropertiesModel::initCoeffs(word name)
{
  wordList info_(nCoeffs_);
  forAll(info_, infoI) {
    info_[infoI]=" * T^"+std::to_string(infoI) ;
    if (infoI < nCoeffs_-1) {
      info_[infoI]+=" +";
    }
  }
  coeffs_.insert(name,new scalarList(nCoeffs_));
  forAll(coeffs_[name], coefI) {
    word var_name=name;
    if (name == "rho_s") var_name = "rho";
    if (name == "emissivity") var_name = "e";
    coeffs_[name][coefI] = createScalarProp(var_name+"_sub_n["+std::to_string(coefI)+"]","yes",0.0);
  }
  Info << getTabLevel() << "Fourier coefficients: "<< name << "(T)= ";
  forAll(coeffs_[name], coefI) {
    Info  << coeffs_[name][coefI] << info_[coefI] << " " ;
  }
  Info << endl;
}

void Foam::FourierMaterialPropertiesModel::init_k_coeffs(word name)
{
  dimensionedTensor kI = dimensionedTensor("kI",dimensionSet(1, 1, -3, -1, 0, 0, 0),I_);
  symmTensorFields_.insert("k_abl_sym",createVolField<symmTensor>("k_abl_sym",symm(kI)));
  initCoeffs(name);
}

void Foam::FourierMaterialPropertiesModel::updateScalarField(word name)
{
  fourierScalarFunction(scalarFields_[name],coeffs_[name],T_);
}

void Foam::FourierMaterialPropertiesModel::updateTensorField(word name)
{
  fourierTensorFunction(tensorFields_[name],coeffs_[name],T_);
}

void Foam::FourierMaterialPropertiesModel::update_k(word name)
{
  updateTensorField(name);
  symmTensorFields_["k_abl_sym"] = symm(tensorFields_[name]);
}

void Foam::FourierMaterialPropertiesModel::fourierScalarFunction(volScalarField& field, scalarList& coeffs, volScalarField& T)
{
  forAll(field, fieldI) {
    field[fieldI] = 0;
    forAll(coeffs, coefI) {
      field[fieldI] += coeffs[coefI]*pow(T[fieldI],coefI);
    }
  }
  forAll(mesh_.boundaryMesh(), patchi) {
    if (!isA<emptyPolyPatch>(mesh_.boundaryMesh()[patchi])) {
      forAll(mesh_.boundaryMesh()[patchi], facei) {
        field.boundaryFieldRef()[patchi][facei] =  0;
        forAll(coeffs, coefI) {
          field.boundaryFieldRef()[patchi][facei] +=
              coeffs[coefI]
              *pow(T.boundaryFieldRef()[patchi][facei],coefI);
        }
      }
    }
  }
}

void Foam::FourierMaterialPropertiesModel::fourierTensorFunction(volTensorField& field, scalarList& coeffs, volScalarField& T)
{
  forAll(field, fieldI) {
    field[fieldI] = tensor::zero;
    forAll(coeffs, coefI) {
      field[fieldI].xx() += coeffs[coefI]*pow(T[fieldI],coefI);
    }
    field[fieldI].yy() = field[fieldI].xx();
    field[fieldI].zz() = field[fieldI].xx();
  }
  forAll(mesh_.boundaryMesh(), patchi) {
    if (!isA<emptyPolyPatch>(mesh_.boundaryMesh()[patchi])) {
      forAll(mesh_.boundaryMesh()[patchi], facei) {
        field.boundaryFieldRef()[patchi][facei] = tensor::zero;
        forAll(coeffs, coefI) {
          field.boundaryFieldRef()[patchi][facei].xx() += coeffs[coefI]*pow(T.boundaryFieldRef()[patchi][facei],coefI);
        }
        field.boundaryFieldRef()[patchi][facei].yy() = field.boundaryFieldRef()[patchi][facei].xx();
        field.boundaryFieldRef()[patchi][facei].zz() = field.boundaryFieldRef()[patchi][facei].xx();
      }
    }
  }
}



// ************************************************************************* //
