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
T_(createVolField<scalar>("Ta")),
rho_s_(createVolField<scalar>("rho_s",dimensionedScalar("0",dimMass/dimVolume,0))),
k_(createVolField<tensor>("k",dimensionedTensor("0",dimensionSet(1, 1, -3, -1, 0, 0, 0),tensor(1,0,0,0,1,0,0,0,1)))),
cp_(createVolField<scalar>("cp",dimensionedScalar("0", dimensionSet(0,2,-2,-1,0,0,0), scalar(0))))
{
  modelInitialized();
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
