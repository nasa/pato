/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Description
    Creates MassModel instances

\*---------------------------------------------------------------------------*/

#include "simpleMassModel.H"
#include "noMassModel.H"
#include "DarcyLaw_HeterogeneousMassModel.H"
#include "DarcyLawMassModel.H"
#include "DarcyForchheimerLawMassModel.H"
#include "DarcyLaw2TMassModel.H"
#include "DarcyForchheimerLaw2TMassModel.H"
#include "FixedPressureGammaMassModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(noMassModel, 0);
addToRunTimeSelectionTable(simpleMassModel, noMassModel, fvMesh);

defineTypeNameAndDebug(DarcyLaw_HeterogeneousMassModel, 0);
addToRunTimeSelectionTable(simpleMassModel, DarcyLaw_HeterogeneousMassModel, fvMesh);

defineTypeNameAndDebug(DarcyLawMassModel, 0);
addToRunTimeSelectionTable(simpleMassModel, DarcyLawMassModel, fvMesh);

defineTypeNameAndDebug(DarcyLaw2TMassModel, 0);
addToRunTimeSelectionTable(simpleMassModel, DarcyLaw2TMassModel, fvMesh);

defineTypeNameAndDebug(FixedPressureGammaMassModel, 0);
addToRunTimeSelectionTable(simpleMassModel, FixedPressureGammaMassModel, fvMesh);

defineTypeNameAndDebug(DarcyForchheimerLawMassModel, 0);
addToRunTimeSelectionTable(simpleMassModel, DarcyForchheimerLawMassModel, fvMesh);

defineTypeNameAndDebug(DarcyForchheimerLaw2TMassModel, 0);
addToRunTimeSelectionTable(simpleMassModel, DarcyForchheimerLaw2TMassModel, fvMesh);
}

// ************************************************************************* //
