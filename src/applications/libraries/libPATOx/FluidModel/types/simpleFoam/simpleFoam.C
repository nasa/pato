/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "simpleFoam.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::simpleFoam::simpleFoam
(
    Time& runTime
)
  :
basicFluidModel(runTime)
{
  // ADDED IN MAIN
#include "simpleFoam/createFluidMeshes.H"
#include "simpleFoam/createFluidFields.H"
#include "simpleFoam/initContinuityErrs.H"
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::simpleFoam::~simpleFoam()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::simpleFoam::updateBefore()
{
#include "simpleFoam/readSIMPLEControls.H"
}

void Foam::simpleFoam::updateAfter()
{
  forAll(fluidRegions, i) {
#include "simpleFoam/storeOldFluidFields.H"
    Info<< "\nSolving for fluid region "
        << fluidRegions[i].name() << endl;
#include "simpleFoam/setRegionFluidFields.H"
#include "simpleFoam/readFluidMultiRegionSIMPLEControls.H"
#include "simpleFoam/solveFluid.H"
  }
}

// ************************************************************************* //

