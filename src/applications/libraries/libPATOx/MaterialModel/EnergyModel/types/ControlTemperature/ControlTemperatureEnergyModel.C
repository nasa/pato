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
    along with OpenFOAM.  If PureConductiont, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "ControlTemperatureEnergyModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ControlTemperatureEnergyModel::ControlTemperatureEnergyModel
(
    const fvMesh& mesh,
    const word& regionName
)
  :
simpleEnergyModel(mesh, regionName),
T_(createVolField<scalar>("Ta")),
rho_s_(createVolField<scalar>("rho_s",dimensionedScalar("0",dimMass/dimVolume,0))),
coeff1
(
    "coeff1",
    dimensionSet(0,0,0,1,0,0,0),
    simpleEnergyModel::materialDict_.subDict("Energy").template lookupOrDefault<int>("coeff1",300)
),
coeff2
(
    "coeff2",
    dimensionSet(0,0,0,1,0,0,0),
    simpleEnergyModel::materialDict_.subDict("Energy").template lookupOrDefault<int>("coeff2",0)
)
{
  modelInitialized();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ControlTemperatureEnergyModel::~ControlTemperatureEnergyModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ControlTemperatureEnergyModel::update()
{
  // solve energy conservation
  const Time& runTime_ = mesh_.time();

  scalar t=runTime_.value();

  T_ = coeff1 + coeff2*t;

}

// ************************************************************************* //
