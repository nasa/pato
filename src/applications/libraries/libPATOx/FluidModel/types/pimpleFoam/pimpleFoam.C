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

#include "pimpleFoam.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pimpleFoam::pimpleFoam
(
    Time& runTime
)
  :
basicFluidModel(runTime)
{
#include "pimpleFoam/createFluidMeshes.H"
#include "pimpleFoam/createFluidFields.H"
#include "pimpleFoam/initContinuityErrs.H"
  // From pimpleFoam
#include "pimpleFoam/createTimeControls.H"
#include "pimpleFoam/multiRegionCourantNo.H"
#include "pimpleFoam/setInitialDeltaT.H"
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pimpleFoam::~pimpleFoam()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pimpleFoam::updateBefore()
{
#include "pimpleFoam/readPIMPLEControls.H"
#include "pimpleFoam/readDyMControls.H"
#include "pimpleFoam/multiRegionCourantNo.H"
#include "pimpleFoam/setFluidDeltaT.H"
}

void Foam::pimpleFoam::updateAfter()
{
  // if (nOuterCorr != 1)
  // {
  //     forAll(fluidRegions, i)
  //     {
  //             #include "pimpleFoam/storeOldFluidFields.H"
  //     }
  // }

  // --- PIMPLE Loop
  for (int oCorr = 0; oCorr < nOuterCorr; oCorr++) {
    bool finalIter = oCorr == nOuterCorr - 1;

    forAll(fluidRegions, i) {
      Info<< "\nSolving for fluid region "
          << fluidRegions[i].name() << endl;
#include "pimpleFoam/setRegionFluidFields.H"
#include "pimpleFoam/readFluidMultiRegionPIMPLEControls.H"
#include "pimpleFoam/solveFluid.H"
    }
  }
}

// ************************************************************************* //
