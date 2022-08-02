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

  Info << "Start main loop" << endl << endl;
  while(runTime.run()) { // main loop
    if (runTime.value() >= fluidPATOxModel().delayFluid()) {
    fluidPATOxModel().updateBefore(); // update fluid model before runTime++
    }
    runTime++;
    if (runTime.value() >= fluidPATOxModel().delayFluid()) {
    Info << "runTime = " << runTime.timeName() << " s  Time step = " << runTime.deltaT().value() << " s" << endl;
    fluidPATOxModel().updateAfter(); // update fluid model after runTime++ 
    }
    if (runTime.value() >= fluidPATOxModel().delaySolid()) {
      matPATOxModel.update(); // update material model when runtime > delay solid time  
    }
    runTime.write();
    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;
    Info << " ------------------------------------------------------------------ " << endl;
    Info << endl;
  }
  return 0;
}


// ************************************************************************* //

//  LocalWords:  dt endl
