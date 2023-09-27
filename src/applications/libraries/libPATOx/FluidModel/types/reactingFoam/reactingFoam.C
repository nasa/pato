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

#include "reactingFoam.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reactingFoam::reactingFoam
(
    Time& runTime
)
  :
basicFluidModel(runTime)
{
  // ADDED IN MAIN
#include "reactingFoam/createFluidMeshes.H"
#include "reactingFoam/createFluidFields.H"
#include "reactingFoam/initContinuityErrs.H"
#include "reactingFoam/createTimeControls.H"
#include "reactingFoam/readTimeControls.H"
#include "reactingFoam/compressibleMultiRegionCourantNo.H"
#include "reactingFoam/setInitialMultiRegionDeltaT.H"
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::reactingFoam::~reactingFoam()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::reactingFoam::updateBefore()
{
#include "reactingFoam/readTimeControls.H"
#include "reactingFoam/readPIMPLEControls.H"
#include "reactingFoam/compressibleMultiRegionCourantNo.H"
#include "reactingFoam/setFluidDeltaT.H"
}

void Foam::reactingFoam::updateAfter()
{
//    forAll(fluidRegions, i){
//        dynamicFvMesh& mesh = (dynamicFvMesh&)(fluidRegions[i]);
//        mesh.update();
//        const Time& runTime = mesh.time();
//        #include "volContinuity.H" // checks the topology of the mesh after the motion
//    }

  if (nOuterCorr != 1) {
    forAll(fluidRegions, i) {
#include "reactingFoam/storeOldFluidFields.H"
    }
  }

  // --- PIMPLE loop
  for (int oCorr=0; oCorr<nOuterCorr; oCorr++) {
    bool finalIter = oCorr == nOuterCorr-1;

    forAll(fluidRegions, i) {
      Info<< "\nSolving for fluid region "
          << fluidRegions[i].name() << endl;
#include "reactingFoam/setRegionFluidFields.H"
#include "reactingFoam/readFluidMultiRegionPIMPLEControls.H"
#include "reactingFoam/solveFluid.H"
    }

  }

}

// ************************************************************************* //
