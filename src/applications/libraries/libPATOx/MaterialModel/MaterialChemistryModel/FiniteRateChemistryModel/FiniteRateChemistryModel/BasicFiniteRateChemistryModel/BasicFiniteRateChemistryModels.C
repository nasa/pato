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

InClass
    Foam::psiChemistryModel

Description
    Creates chemistry model instances templated on the type of thermodynamics

\*---------------------------------------------------------------------------*/

#include "makeFiniteRateChemistryModel.H"

#include "psiReactionThermo.H"
#include "rhoReactionThermo.H"

#include "StandardFiniteRateChemistryModel.H"
#include "TDACFiniteRateChemistryModel.H"
#include "thermoPhysicsTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// Make base types
makeChemistryModel(psiReactionThermo);
makeChemistryModel(rhoReactionThermo);

// Chemistry moldels based on sensibleEnthalpy
makeChemistryModelType
(
    StandardFiniteRateChemistryModel,
    psiReactionThermo,
    constGasHThermoPhysics
);

makeChemistryModelType
(
    StandardFiniteRateChemistryModel,
    psiReactionThermo,
    gasHThermoPhysics
);

makeChemistryModelType
(
    StandardFiniteRateChemistryModel,
    psiReactionThermo,
    constIncompressibleGasHThermoPhysics
);

makeChemistryModelType
(
    StandardFiniteRateChemistryModel,
    psiReactionThermo,
    incompressibleGasHThermoPhysics
);

makeChemistryModelType
(
    StandardFiniteRateChemistryModel,
    psiReactionThermo,
    icoPoly8HThermoPhysics
);

makeChemistryModelType
(
    StandardFiniteRateChemistryModel,
    psiReactionThermo,
    constFluidHThermoPhysics
);

makeChemistryModelType
(
    StandardFiniteRateChemistryModel,
    psiReactionThermo,
    constAdiabaticFluidHThermoPhysics
);

makeChemistryModelType
(
    StandardFiniteRateChemistryModel,
    psiReactionThermo,
    constHThermoPhysics
);


makeChemistryModelType
(
    StandardFiniteRateChemistryModel,
    rhoReactionThermo,
    constGasHThermoPhysics
);

makeChemistryModelType
(
    StandardFiniteRateChemistryModel,
    rhoReactionThermo,
    gasHThermoPhysics
);

makeChemistryModelType
(
    StandardFiniteRateChemistryModel,
    rhoReactionThermo,
    constIncompressibleGasHThermoPhysics
);

makeChemistryModelType
(
    StandardFiniteRateChemistryModel,
    rhoReactionThermo,
    incompressibleGasHThermoPhysics
);

makeChemistryModelType
(
    StandardFiniteRateChemistryModel,
    rhoReactionThermo,
    icoPoly8HThermoPhysics
);

makeChemistryModelType
(
    StandardFiniteRateChemistryModel,
    rhoReactionThermo,
    constFluidHThermoPhysics
);

makeChemistryModelType
(
    StandardFiniteRateChemistryModel,
    rhoReactionThermo,
    constAdiabaticFluidHThermoPhysics
);

makeChemistryModelType
(
    StandardFiniteRateChemistryModel,
    rhoReactionThermo,
    constHThermoPhysics
);


makeChemistryModelType
(
    TDACFiniteRateChemistryModel,
    psiReactionThermo,
    constGasHThermoPhysics
);

makeChemistryModelType
(
    TDACFiniteRateChemistryModel,
    psiReactionThermo,
    gasHThermoPhysics
);

makeChemistryModelType
(
    TDACFiniteRateChemistryModel,
    psiReactionThermo,
    constIncompressibleGasHThermoPhysics
);

makeChemistryModelType
(
    TDACFiniteRateChemistryModel,
    psiReactionThermo,
    incompressibleGasHThermoPhysics
);

makeChemistryModelType
(
    TDACFiniteRateChemistryModel,
    psiReactionThermo,
    icoPoly8HThermoPhysics
);

makeChemistryModelType
(
    TDACFiniteRateChemistryModel,
    psiReactionThermo,
    constFluidHThermoPhysics
);

makeChemistryModelType
(
    TDACFiniteRateChemistryModel,
    psiReactionThermo,
    constAdiabaticFluidHThermoPhysics
);

makeChemistryModelType
(
    TDACFiniteRateChemistryModel,
    psiReactionThermo,
    constHThermoPhysics
);


makeChemistryModelType
(
    TDACFiniteRateChemistryModel,
    rhoReactionThermo,
    constGasHThermoPhysics
);

makeChemistryModelType
(
    TDACFiniteRateChemistryModel,
    rhoReactionThermo,
    gasHThermoPhysics
);

makeChemistryModelType
(
    TDACFiniteRateChemistryModel,
    rhoReactionThermo,
    constIncompressibleGasHThermoPhysics
);

makeChemistryModelType
(
    TDACFiniteRateChemistryModel,
    rhoReactionThermo,
    incompressibleGasHThermoPhysics
);

makeChemistryModelType
(
    TDACFiniteRateChemistryModel,
    rhoReactionThermo,
    icoPoly8HThermoPhysics
);

makeChemistryModelType
(
    TDACFiniteRateChemistryModel,
    rhoReactionThermo,
    constFluidHThermoPhysics
);

makeChemistryModelType
(
    TDACFiniteRateChemistryModel,
    rhoReactionThermo,
    constAdiabaticFluidHThermoPhysics
);

makeChemistryModelType
(
    TDACFiniteRateChemistryModel,
    rhoReactionThermo,
    constHThermoPhysics
);


// Chemistry moldels based on sensibleInternalEnergy
makeChemistryModelType
(
    StandardFiniteRateChemistryModel,
    psiReactionThermo,
    constGasEThermoPhysics
);

makeChemistryModelType
(
    StandardFiniteRateChemistryModel,
    psiReactionThermo,
    gasEThermoPhysics
);

makeChemistryModelType
(
    StandardFiniteRateChemistryModel,
    psiReactionThermo,
    constIncompressibleGasEThermoPhysics
);

makeChemistryModelType
(
    StandardFiniteRateChemistryModel,
    psiReactionThermo,
    incompressibleGasEThermoPhysics
);

makeChemistryModelType
(
    StandardFiniteRateChemistryModel,
    psiReactionThermo,
    icoPoly8EThermoPhysics
);

makeChemistryModelType
(
    StandardFiniteRateChemistryModel,
    psiReactionThermo,
    constFluidEThermoPhysics
);

makeChemistryModelType
(
    StandardFiniteRateChemistryModel,
    psiReactionThermo,
    constAdiabaticFluidEThermoPhysics
);

makeChemistryModelType
(
    StandardFiniteRateChemistryModel,
    psiReactionThermo,
    constEThermoPhysics
);



makeChemistryModelType
(
    StandardFiniteRateChemistryModel,
    rhoReactionThermo,
    constGasEThermoPhysics
);

makeChemistryModelType
(
    StandardFiniteRateChemistryModel,
    rhoReactionThermo,
    gasEThermoPhysics
);

makeChemistryModelType
(
    StandardFiniteRateChemistryModel,
    rhoReactionThermo,
    constIncompressibleGasEThermoPhysics
);

makeChemistryModelType
(
    StandardFiniteRateChemistryModel,
    rhoReactionThermo,
    incompressibleGasEThermoPhysics
);

makeChemistryModelType
(
    StandardFiniteRateChemistryModel,
    rhoReactionThermo,
    icoPoly8EThermoPhysics
);

makeChemistryModelType
(
    StandardFiniteRateChemistryModel,
    rhoReactionThermo,
    constFluidEThermoPhysics
);

makeChemistryModelType
(
    StandardFiniteRateChemistryModel,
    rhoReactionThermo,
    constAdiabaticFluidEThermoPhysics
);

makeChemistryModelType
(
    StandardFiniteRateChemistryModel,
    rhoReactionThermo,
    constEThermoPhysics
);


makeChemistryModelType
(
    TDACFiniteRateChemistryModel,
    psiReactionThermo,
    constGasEThermoPhysics
);

makeChemistryModelType
(
    TDACFiniteRateChemistryModel,
    psiReactionThermo,
    gasEThermoPhysics
);

makeChemistryModelType
(
    TDACFiniteRateChemistryModel,
    psiReactionThermo,
    constIncompressibleGasEThermoPhysics
);

makeChemistryModelType
(
    TDACFiniteRateChemistryModel,
    psiReactionThermo,
    incompressibleGasEThermoPhysics
);

makeChemistryModelType
(
    TDACFiniteRateChemistryModel,
    psiReactionThermo,
    icoPoly8EThermoPhysics
);

makeChemistryModelType
(
    TDACFiniteRateChemistryModel,
    psiReactionThermo,
    constFluidEThermoPhysics
);

makeChemistryModelType
(
    TDACFiniteRateChemistryModel,
    psiReactionThermo,
    constAdiabaticFluidEThermoPhysics
);

makeChemistryModelType
(
    TDACFiniteRateChemistryModel,
    psiReactionThermo,
    constEThermoPhysics
);


makeChemistryModelType
(
    TDACFiniteRateChemistryModel,
    rhoReactionThermo,
    constGasEThermoPhysics
);

makeChemistryModelType
(
    TDACFiniteRateChemistryModel,
    rhoReactionThermo,
    gasEThermoPhysics
);

makeChemistryModelType
(
    TDACFiniteRateChemistryModel,
    rhoReactionThermo,
    constIncompressibleGasEThermoPhysics
);

makeChemistryModelType
(
    TDACFiniteRateChemistryModel,
    rhoReactionThermo,
    incompressibleGasEThermoPhysics
);

makeChemistryModelType
(
    TDACFiniteRateChemistryModel,
    rhoReactionThermo,
    icoPoly8EThermoPhysics
);

makeChemistryModelType
(
    TDACFiniteRateChemistryModel,
    rhoReactionThermo,
    constFluidEThermoPhysics
);

makeChemistryModelType
(
    TDACFiniteRateChemistryModel,
    rhoReactionThermo,
    constAdiabaticFluidEThermoPhysics
);

makeChemistryModelType
(
    TDACFiniteRateChemistryModel,
    rhoReactionThermo,
    constEThermoPhysics
);
}

// ************************************************************************* //
