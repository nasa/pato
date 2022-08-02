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

#include "Fourier_RadiationMaterialPropertiesModel.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Fourier_RadiationMaterialPropertiesModel::Fourier_RadiationMaterialPropertiesModel
(
    const fvMesh& mesh,
    const word& dictName
)
  :
FourierMaterialPropertiesModel(mesh, dictName),
emissivityCoefs_(nCoefs_),
emissivity_(meshLookupOrConstructScalar(mesh, "emissivity",dimensionedScalar("0", dimless, scalar(0.0)))),
Qr_(meshLookupOrConstructScalar(mesh, "Qr",dimensionedScalar("0", dimMass/ pow3(dimTime)/ dimLength, scalar(0.0)))),
sigmaSB(::constant::physicoChemical::sigma),
Tbackground_(meshLookupOrConstructScalar(mesh, "Tbackground", dimensionedScalar("0", dimTemperature, scalar(0.0))))
{
  wordList info_(nCoefs_);
  forAll(info_, infoI) {
    info_[infoI]=" * T^"+std::to_string(infoI) ;
    if (infoI < nCoefs_-1) {
      info_[infoI]+=" +";
    }
  }
  Info << "e(T) = ";
  forAll(emissivityCoefs_, coefI) {
    emissivityCoefs_[coefI]= materialPropertiesDictionary.lookupOrDefault<scalar>("e_sub_n["+std::to_string(coefI)+"]", 0.0);
    Info  << emissivityCoefs_[coefI]  << info_[coefI] << " " ;
  }

  Info << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Fourier_RadiationMaterialPropertiesModel::~Fourier_RadiationMaterialPropertiesModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::Fourier_RadiationMaterialPropertiesModel::update()
{
  FourierMaterialPropertiesModel::update();
  FourierMaterialPropertiesModel::fourierFunction(emissivity_, emissivityCoefs_, T_);

  forAll(mesh_.boundaryMesh(), patchI) {
    forAll(mesh_.boundaryMesh()[patchI], faceI) {
      Qr_.boundaryFieldRef()[patchI][faceI] = - emissivity_.boundaryFieldRef()[patchI][faceI]*
                                              sigmaSB.value()
                                              *(pow4(T_.boundaryFieldRef()[patchI][faceI])
                                                  -pow4(Tbackground_.boundaryFieldRef()[patchI][faceI]));
    }
  }
}

// ************************************************************************* //
