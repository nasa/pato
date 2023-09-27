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

#include "chtMultiRegionFoam.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::chtMultiRegionFoam::chtMultiRegionFoam
(
    Time& runTime
)
  :
basicFluidModel(runTime)
{
  // ADDED IN MAIN
#include "chtMultiRegionFoam/createFluidMeshes.H"
#include "chtMultiRegionFoam/createFluidFields.H"
#include "chtMultiRegionFoam/initContinuityErrs.H"

  //  from chtMultiRegionFoam
#include "chtMultiRegionFoam/createTimeControls.H"
#include "chtMultiRegionFoam/readTimeControls.H"
#include "chtMultiRegionFoam/compressibleMultiRegionCourantNo.H"
#include "chtMultiRegionFoam/setInitialMultiRegionDeltaT.H"
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::chtMultiRegionFoam::~chtMultiRegionFoam()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::chtMultiRegionFoam::updateBefore()
{
#include "chtMultiRegionFoam/readTimeControls.H"
#include "chtMultiRegionFoam/readPIMPLEControls.H"
#include "chtMultiRegionFoam/compressibleMultiRegionCourantNo.H"
#include "chtMultiRegionFoam/setFluidDeltaT.H"
}

void Foam::chtMultiRegionFoam::updateAfter()
{
  if (nOuterCorr != 1) {
    forAll(fluidRegions, i) {
#include "chtMultiRegionFoam/storeOldFluidFields.H"
    }
  }

  // --- PIMPLE loop
  for (int oCorr=0; oCorr<nOuterCorr; oCorr++) {
    bool finalIter = oCorr == nOuterCorr-1;

    forAll(fluidRegions, i) {
      Info<< "\nSolving for fluid region "
          << fluidRegions[i].name() << endl;
#include "chtMultiRegionFoam/setRegionFluidFields.H"
#include "chtMultiRegionFoam/readFluidMultiRegionPIMPLEControls.H"
#include "chtMultiRegionFoam/solveFluid.H"
    }
  }
}

// ************************************************************************* //
