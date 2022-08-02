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

#include "makeFiniteRateChemistrySolverTypes.H"

#include "thermoPhysicsTypes.H"
#include "psiReactionThermo.H"
#include "rhoReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// Chemistry solvers based on sensibleEnthalpy
makeFiniteRateChemistrySolverTypes(psiReactionThermo, constGasHThermoPhysics);
makeFiniteRateChemistrySolverTypes(psiReactionThermo, gasHThermoPhysics);
makeFiniteRateChemistrySolverTypes
(
    psiReactionThermo,
    constIncompressibleGasHThermoPhysics
);
makeFiniteRateChemistrySolverTypes
(
    psiReactionThermo,
    incompressibleGasHThermoPhysics
);
makeFiniteRateChemistrySolverTypes(psiReactionThermo, icoPoly8HThermoPhysics);
makeFiniteRateChemistrySolverTypes(psiReactionThermo, constFluidHThermoPhysics);
makeFiniteRateChemistrySolverTypes
(
    psiReactionThermo,
    constAdiabaticFluidHThermoPhysics
);
makeFiniteRateChemistrySolverTypes(psiReactionThermo, constHThermoPhysics);

makeFiniteRateChemistrySolverTypes(rhoReactionThermo, constGasHThermoPhysics);
makeFiniteRateChemistrySolverTypes(rhoReactionThermo, gasHThermoPhysics);
makeFiniteRateChemistrySolverTypes
(
    rhoReactionThermo,
    constIncompressibleGasHThermoPhysics
);
makeFiniteRateChemistrySolverTypes
(
    rhoReactionThermo,
    incompressibleGasHThermoPhysics
);
makeFiniteRateChemistrySolverTypes(rhoReactionThermo, icoPoly8HThermoPhysics);
makeFiniteRateChemistrySolverTypes(rhoReactionThermo, constFluidHThermoPhysics);
makeFiniteRateChemistrySolverTypes
(
    rhoReactionThermo,
    constAdiabaticFluidHThermoPhysics
);
makeFiniteRateChemistrySolverTypes(rhoReactionThermo, constHThermoPhysics);

// Chemistry solvers based on sensibleInternalEnergy
makeFiniteRateChemistrySolverTypes(psiReactionThermo, constGasEThermoPhysics);
makeFiniteRateChemistrySolverTypes(psiReactionThermo, gasEThermoPhysics);
makeFiniteRateChemistrySolverTypes
(
    psiReactionThermo,
    constIncompressibleGasEThermoPhysics
);
makeFiniteRateChemistrySolverTypes
(
    psiReactionThermo,
    incompressibleGasEThermoPhysics
);
makeFiniteRateChemistrySolverTypes(psiReactionThermo, icoPoly8EThermoPhysics);
makeFiniteRateChemistrySolverTypes(psiReactionThermo, constFluidEThermoPhysics);
makeFiniteRateChemistrySolverTypes
(
    psiReactionThermo,
    constAdiabaticFluidEThermoPhysics
);
makeFiniteRateChemistrySolverTypes(psiReactionThermo, constEThermoPhysics);

makeFiniteRateChemistrySolverTypes(rhoReactionThermo, constGasEThermoPhysics);
makeFiniteRateChemistrySolverTypes(rhoReactionThermo, gasEThermoPhysics);
makeFiniteRateChemistrySolverTypes
(
    rhoReactionThermo,
    constIncompressibleGasEThermoPhysics
);
makeFiniteRateChemistrySolverTypes
(
    rhoReactionThermo,
    incompressibleGasEThermoPhysics
);
makeFiniteRateChemistrySolverTypes(rhoReactionThermo, icoPoly8EThermoPhysics);
makeFiniteRateChemistrySolverTypes(rhoReactionThermo, constFluidEThermoPhysics);
makeFiniteRateChemistrySolverTypes
(
    rhoReactionThermo,
    constAdiabaticFluidEThermoPhysics
);
makeFiniteRateChemistrySolverTypes(rhoReactionThermo, constEThermoPhysics);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
