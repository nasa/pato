/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Global
    setMultiRegionDeltaT

Description
    Reset the timestep to maintain a constant maximum courant and
    diffusion Numbers. Reduction of time-step is immediate, but
    increase is damped to avoid unstable oscillations.

\*---------------------------------------------------------------------------*/

/**** Initialization ****/

scalar maxCo = runTime.controlDict().lookupOrDefault<scalar>("maxCo", 1.0);
scalar maxDeltaT = runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);
scalar minDeltaT = runTime.controlDict().lookupOrDefault<scalar>("minDeltaT", SMALL);
scalar adjustStartTime = runTime.controlDict().lookupOrDefault<scalar>("adjustStartTime", 0);

// ***** Adjust time step ***** //
//CoNum=0.0;

if (adjustTimeStep && (runTime.value() > adjustStartTime))
{
  if (CoNum == -great) {
    CoNum = small;
  }

  if (DiNum == -great) {
    DiNum = small;
  }

  scalar TFactorFluid = maxCo/(CoNum + SMALL);
  scalar TFactorSolid = maxDi/(DiNum + SMALL);

//  scalar TFactorFilm = maxCo/(surfaceFilm.CourantNumber() + SMALL);

  scalar dt0 = runTime.deltaTValue();

  scalar deltaTphys =
      max
      (
          min
          (
              dt0*min(min(TFactorFluid,TFactorSolid), 1.2),
              maxDeltaT
          ),
          minDeltaT
      );

  scalar deltaTFact = min(min(TFactorFluid, 1.0 + 0.1*TFactorFluid), 1.2);

  runTime.setDeltaT
  (
      min
      (
          deltaTFact*deltaTphys,
          maxDeltaT
      )
  );

//  scalar maxDeltaTFact = maxCo / (CoNum + SMALL);
//  scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.2 * maxDeltaTFact), 1.05);
//  scalar deltaTphys =
//      max
//      (
//          min
//          (
//              deltaTFact * runTime.deltaT().value(),
//              maxDeltaT
//          ),
//          minDeltaT
//      );
//
//  runTime.setDeltaT(deltaTphys);
}

// ************************************************************************* //
