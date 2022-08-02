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

#include "FourierMaterialPropertiesModel.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FourierMaterialPropertiesModel::FourierMaterialPropertiesModel
(
    const fvMesh& mesh,
    const word& dictName
)
  :
simpleMaterialPropertiesModel(mesh, dictName),
materialPropertiesDirectory(fileName(simpleMaterialPropertiesModel::materialDict_.subDict("MaterialProperties").lookup("MaterialPropertiesDirectory")).expand()),
materialPropertiesDictionary
(
    IOobject
    (
        materialPropertiesDirectory+"/FourierProperties",
        mesh.time().db().parent(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
),
nCoefs_(5),
rhoCoefs_(nCoefs_),
cpCoefs_(nCoefs_),
kCoefs_(nCoefs_),
cp_(simpleMaterialPropertiesModel::cp_),
k_(simpleMaterialPropertiesModel::k_),
rho_s_(simpleMaterialPropertiesModel::rho_s_),
T_(meshLookupOrConstructScalar(mesh, "Ta")),
kijk_
(
    IOobject
    (
        "kijk",
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedTensor("kI",dimensionSet(1, 1, -3, -1, 0, 0, 0),tensor(1, 0, 0, 0, 1, 0, 0, 0, 1))
),
k_abl_sym_
(
    IOobject
    (
        "k_abl_sym",
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    symm(kijk_)

)

{

  Info << "Reading " << materialPropertiesDirectory << endl;

  wordList info_(nCoefs_);
  forAll(info_, infoI) {
    info_[infoI]=" * T^"+std::to_string(infoI) ;
    if (infoI < nCoefs_-1) {
      info_[infoI]+=" +";
    }
  }

  Info << "Fourier coefficients:" << endl;
  Info << "rho(T)= ";
  forAll(rhoCoefs_, coefI) {
    rhoCoefs_[coefI]= materialPropertiesDictionary.lookupOrDefault<scalar>("rho_sub_n["+std::to_string(coefI)+"]", 0.0);
    Info  << rhoCoefs_[coefI] << info_[coefI] << " " ;
  }

  Info << endl;
  Info << "cp(T) = ";
  forAll(cpCoefs_, coefI) {
    cpCoefs_[coefI]= materialPropertiesDictionary.lookupOrDefault<scalar>("cp_sub_n["+std::to_string(coefI)+"]", 0.0);
    Info <<  cpCoefs_[coefI] << info_[coefI] << " " ;
  }

  Info << endl;
  Info << "k(T) = ";
  forAll(kCoefs_, coefI) {
    kCoefs_[coefI]= materialPropertiesDictionary.lookupOrDefault<scalar>("k_sub_n["+std::to_string(coefI)+"]", 0.0);
    Info  << kCoefs_[coefI]  << info_[coefI] << " " ;
  }

  Info << endl;
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::FourierMaterialPropertiesModel::~FourierMaterialPropertiesModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::FourierMaterialPropertiesModel::update()
{
  fourierFunction(rho_s_, rhoCoefs_, T_);
  fourierFunction(cp_, cpCoefs_, T_);
  forAll(k_, fieldI) {
    tensor zeroTens(0, 0, 0, 0, 0, 0, 0, 0, 0);
    k_[fieldI]= zeroTens;
    forAll(kCoefs_, coefI) {
      k_[fieldI].xx()+=kCoefs_[coefI]*pow(T_[fieldI],coefI);
    }

    k_[fieldI].xy()=0;
    k_[fieldI].xz()=0;
    k_[fieldI].yx()=0;
    k_[fieldI].yy()=k_[fieldI].xx();
    k_[fieldI].yz()=0;
    k_[fieldI].zx()=0;
    k_[fieldI].zy()=0;
    k_[fieldI].zz()=k_[fieldI].xx();
  }

  forAll(mesh_.boundaryMesh(), patchI) {
    forAll(mesh_.boundaryMesh()[patchI], faceI) {
      tensor zeroTens(0, 0, 0, 0, 0, 0, 0, 0, 0);
      k_.boundaryFieldRef()[patchI][faceI] = zeroTens;
      forAll(kCoefs_, coefI) {
        k_.boundaryFieldRef()[patchI][faceI].xx() +=kCoefs_[coefI]*pow(T_.boundaryFieldRef()[patchI][faceI],coefI);
      }
      k_.boundaryFieldRef()[patchI][faceI].xy()=0;
      k_.boundaryFieldRef()[patchI][faceI].xz()=0;
      k_.boundaryFieldRef()[patchI][faceI].yx()=0;
      k_.boundaryFieldRef()[patchI][faceI].yy()=k_.boundaryFieldRef()[patchI][faceI].xx();
      k_.boundaryFieldRef()[patchI][faceI].yz()=0;
      k_.boundaryFieldRef()[patchI][faceI].zx()=0;
      k_.boundaryFieldRef()[patchI][faceI].zy()=0;
      k_.boundaryFieldRef()[patchI][faceI].zz()=k_.boundaryFieldRef()[patchI][faceI].xx();
    }
  }

  k_abl_sym_ = symm(k_);
}

void Foam::FourierMaterialPropertiesModel::fourierFunction(volScalarField& field, scalarList& coefs, volScalarField& T)
{
  forAll(field, fieldI) {
    field[fieldI]=0;
    forAll(coefs, coefI) {
      field[fieldI]+=coefs[coefI]*pow(T[fieldI],coefI);
    }
  }
  forAll(mesh_.boundaryMesh(), patchI) {
    forAll(mesh_.boundaryMesh()[patchI], faceI) {
      field.boundaryFieldRef()[patchI][faceI] =  0;
      forAll(coefs, coefI) {
        field.boundaryFieldRef()[patchI][faceI] +=coefs[coefI]*pow(T.boundaryFieldRef()[patchI][faceI],coefI);
      }
    }
  }
}



// ************************************************************************* //
