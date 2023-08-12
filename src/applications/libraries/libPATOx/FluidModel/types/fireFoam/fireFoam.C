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

#include "fireFoam.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fireFoam::fireFoam
(
    Time& runTime
)
  :
basicFluidModel(runTime)
{
  // ADDED IN MAIN
#include "fireFoam/createFluidMeshes.H"
#include "fireFoam/createFluidFields.H"
#include "fireFoam/initContinuityErrs.H"

  //  from fireFoam
#include "fireFoam/createTimeControls.H"
#include "fireFoam/readTimeControls.H"
#include "fireFoam/compressibleMultiRegionCourantNo.H"
#include "fireFoam/setInitialMultiRegionDeltaT.H"
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fireFoam::~fireFoam()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fireFoam::updateBefore()
{
#include "fireFoam/readTimeControls.H"
#include "fireFoam/readPIMPLEControls.H"
#include "fireFoam/compressibleMultiRegionCourantNo.H"
#include "fireFoam/setFluidDeltaT.H"
}

void Foam::fireFoam::updateAfter()
{
  if (nOuterCorr != 1) {
    forAll(fluidRegions, i) {
#include "fireFoam/storeOldFluidFields.H"
    }
  }

  // --- PIMPLE loop
  for (int oCorr=0; oCorr<nOuterCorr; oCorr++) {
    bool finalIter = oCorr == nOuterCorr-1;

    forAll(fluidRegions, i) {
      Info<< "\nSolving for fluid region "
          << fluidRegions[i].name() << endl;
#include "fireFoam/setRegionFluidFields.H"
#include "fireFoam/readFluidMultiRegionPIMPLEControls.H"
#include "fireFoam/solveFluid.H"
    }
  }
}


// ************************************************************************* //
