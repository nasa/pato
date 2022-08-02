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

Switch adjustTimeStep
(
    runTime.controlDict().lookup("adjustTimeStep")
);

scalar maxCo
(
    readScalar(runTime.controlDict().lookup("maxCo"))
);

scalar maxDeltaT = GREAT;

if (runTime.controlDict().found("maxDeltaT"))
{
  maxDeltaT = readScalar(runTime.controlDict().lookup("maxDeltaT"));
}

scalar minDeltaT = SMALL;

if (runTime.controlDict().found("minDeltaT"))
{
  minDeltaT = readScalar(runTime.controlDict().lookup("minDeltaT"));
}

scalar adjustStartTime
(
    readScalar(runTime.controlDict().lookup("adjustStartTime"))
);

// ***** Adjust time step ***** //

CoNum=0.0;
if (adjustTimeStep && (runTime.value() > adjustStartTime))
{
  scalar maxDeltaTFact = maxCo / (CoNum + SMALL);
  scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.2 * maxDeltaTFact), 1.05);
  scalar deltaTphys =
      max
      (
          min
          (
              deltaTFact * runTime.deltaT().value(),
              maxDeltaT
          ),
          minDeltaT
      );

  runTime.setDeltaT(deltaTphys);
}

// ************************************************************************* //
