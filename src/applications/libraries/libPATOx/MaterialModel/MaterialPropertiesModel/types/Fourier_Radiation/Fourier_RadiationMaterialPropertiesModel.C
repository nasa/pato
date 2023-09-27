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
    const word& regionName
)
  :
FourierMaterialPropertiesModel(mesh, regionName),
startModelInit_(startModelInit()),
Qr_(createVolField<scalar>("Qr",dimensionedScalar("0", dimMass/ pow3(dimTime)/ dimLength, scalar(0.0)))),
emissivity_(createVolFieldIfNotFound<scalar>(energyModel_,"emissivity",dimensionedScalar("0", dimless, scalar(0.0)))),
sigmaSB(::constant::physicoChemical::sigma),
Tbackground_(createVolFieldIfNotFound<scalar>(energyModel_,"Tbackground", dimensionedScalar("0", dimTemperature, scalar(0.0))))
{
  listMapFields.append(newMapField("emissivity","volScalarField",initCoeffs,updateScalarField));
  modelInitialized();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Fourier_RadiationMaterialPropertiesModel::~Fourier_RadiationMaterialPropertiesModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::Fourier_RadiationMaterialPropertiesModel::update()
{
  FourierMaterialPropertiesModel::update();
  forAll(mesh_.boundaryMesh(), patchi) {
    if (!isA<emptyPolyPatch>(mesh_.boundaryMesh()[patchi])) {
      forAll(mesh_.boundaryMesh()[patchi], facei) {
        Qr_.boundaryFieldRef()[patchi][facei] =
            - emissivity_.boundaryFieldRef()[patchi][facei]
            * sigmaSB.value()
            * (pow4(T_.boundaryFieldRef()[patchi][facei])
               - pow4(Tbackground_.boundaryFieldRef()[patchi][facei]));
      }
    }
  }
}

// ************************************************************************* //
