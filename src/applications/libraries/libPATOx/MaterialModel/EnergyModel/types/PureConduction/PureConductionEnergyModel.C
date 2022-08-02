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

#include "PureConductionEnergyModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PureConductionEnergyModel::PureConductionEnergyModel
(
    const fvMesh& mesh,
    const word& dictName
)
  :
simpleEnergyModel(mesh, dictName),
mesh_(mesh),
dictName_(dictName),
materialPropertiesModel_(meshLookupOrConstructModel<simpleMaterialPropertiesModel>(mesh,dictName,"MaterialProperties")),
T_(meshLookupOrConstructScalar(mesh, "Ta")),
rho_s_(materialPropertiesModel_.rho_s()),
cp_(materialPropertiesModel_.cp()),
k_(materialPropertiesModel_.k())
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PureConductionEnergyModel::~PureConductionEnergyModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::PureConductionEnergyModel::update()
{
  // solve energy conservation
  if(!this->dynamicMesh_) {
    solve
    (
        fvm::ddt(cp_*rho_s_, T_) -
        fvm::laplacian(k_, T_)
    );
  } else {
    solve
    (
        cp_*rho_s_ * (fvm::ddt(T_) - fvc::div(mesh_.phi(),T_))
        - fvm::laplacian(k_, T_)
    );
  }
}

// ************************************************************************* //
