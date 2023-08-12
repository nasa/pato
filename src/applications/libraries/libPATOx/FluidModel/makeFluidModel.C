/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Description
    Creates modelName instances

\*---------------------------------------------------------------------------*/

#include "basicFluidModel.H"
#include "noFluidModel.H"
#include "chtMultiRegionFoam.H"
#include "simpleFoam.H"
#include "pimpleFoam.H"
#include "fireFoam.H"
#include "reactingFoam.H"
#include "rhoCentralFoamFluidModel.H"
#include "rhoCentralFoamUserThermo.H"
#include "pureThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(noFluidModel, 0);
addToRunTimeSelectionTable(basicFluidModel, noFluidModel, Time);

defineTypeNameAndDebug(chtMultiRegionFoam, 0);
addToRunTimeSelectionTable(basicFluidModel, chtMultiRegionFoam, Time);

defineTypeNameAndDebug(simpleFoam, 0);
addToRunTimeSelectionTable(basicFluidModel, simpleFoam, Time);

defineTypeNameAndDebug(pimpleFoam, 0);
addToRunTimeSelectionTable(basicFluidModel, pimpleFoam, Time);

defineTypeNameAndDebug(fireFoam, 0);
addToRunTimeSelectionTable(basicFluidModel, fireFoam, Time);

defineTypeNameAndDebug(reactingFoam, 0);
addToRunTimeSelectionTable(basicFluidModel, reactingFoam, Time);

defineTypeNameAndDebug(rhoCentralFoamFluidModel, 0);
addToRunTimeSelectionTable(basicFluidModel, rhoCentralFoamFluidModel, Time);

defineTypeNameAndDebug(rhoCentralFoamUserThermo, 0);
addToRunTimeSelectionTable(basicFluidModel, rhoCentralFoamUserThermo, Time);

defineTypeNameAndDebug(pureThermo, 0);
addToRunTimeSelectionTable(basicFluidModel, pureThermo, Time);
}

// ************************************************************************* //
