/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
  \\    /   O peration     |
  \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

  Utility
  PATOx

  Description
  Run the material and fluid models, update the different model types.
  Stiching coupling strategy based on time scale separation.



  \*---------------------------------------------------------------------------*/

#include "fvCFD.H" // General OpenFOAM library                                                                                                                                                                                                               
#include "PATOx.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
int main(int argc, char *argv[])
{
#include "postProcess.H" // post process initialization 
#include "setRootCase.H" // Create args from argc and argv
#include "createTime.H" // Create runTime object                                                                                                                                                                                                            
  autoPtr<basicFluidModel> fluidPATOxModel(basicFluidModel::New(runTime));
  simpleMaterialsModel matPATOxModel(runTime);

  scalar deltat = 0;
  scalar maxStitchIter = 0;

  Info << "Start main loop" << endl << endl;
  while(runTime.run()) {

    while ( ((fluidPATOxModel().maxDeltaTw() > fluidPATOxModel().fluidStitchTolerance()) && (maxStitchIter < 1000)) || (runTime.value() < fluidPATOxModel().delaySolid()) ) {
      if(runTime.run()) {
        fluidPATOxModel().updateBefore();
        runTime++;
        Info << "runTime = " << runTime.timeName() << " s  Time step = " << runTime.deltaT().value() << " s"
             << "   maxDeltaTw = " << fluidPATOxModel().maxDeltaTw() << " K/s " << endl;
        fluidPATOxModel().updateAfter();
        maxStitchIter = maxStitchIter + 1;
        Info << " ------------------------------------------------------------------ " << endl;
        Info << endl;
      } else {
        return 0;
      }
    }

    maxStitchIter = 0;
    fluidPATOxModel().maxDeltaTw() = fluidPATOxModel().fluidStitchTolerance() * 2;

    while ( (matPATOxModel.minDeltaTw() < fluidPATOxModel().solidStitchTolerance()) || (deltat < fluidPATOxModel().forceFluidUpdate()) ) {
      if(runTime.run()) {
        matPATOxModel.update();
        deltat = deltat + runTime.deltaT().value();
        runTime++;
        runTime.write();
        Info << "runTime = " << runTime.timeName() << " s  Time step = " << runTime.deltaT().value() << " s";
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
        Info << " ------------------------------------------------------------------ " << endl;
        Info << endl;
      } else {
        return 0;
      }
    }

    deltat = 0;

  }
  return 0;
}


// ************************************************************************* //

//  LocalWords:  dt endl
