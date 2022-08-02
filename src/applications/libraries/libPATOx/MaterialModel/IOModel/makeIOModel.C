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
    Creates IOModel instances

\*---------------------------------------------------------------------------*/

#include "simpleIOModel.H"
#include "noIOModel.H"
#include "massIOModel.H"
#include "WriteControlIOModel.H"
#include "InverseProblemIOModel.H"
#include "LinearInterpolationIOModel.H"
#include "ProfileIOModel.H"
#include "PurePyrolysisIOModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(noIOModel, 0);
addToRunTimeSelectionTable(simpleIOModel, noIOModel, fvMesh);

defineTypeNameAndDebug(InverseProblemIOModel, 0);
addToRunTimeSelectionTable(simpleIOModel, InverseProblemIOModel, fvMesh);

defineTypeNameAndDebug(LinearInterpolationIOModel, 0);
addToRunTimeSelectionTable(simpleIOModel, LinearInterpolationIOModel, fvMesh);

defineTypeNameAndDebug(ProfileIOModel, 0);
addToRunTimeSelectionTable(simpleIOModel, ProfileIOModel, fvMesh);

defineTypeNameAndDebug(PurePyrolysisIOModel, 0);
addToRunTimeSelectionTable(simpleIOModel, PurePyrolysisIOModel, fvMesh);

defineTypeNameAndDebug(WriteControlIOModel, 0);
addToRunTimeSelectionTable(simpleIOModel, WriteControlIOModel, fvMesh);
defineTypeNameAndDebug(massIOModel, 0);
addToRunTimeSelectionTable(simpleIOModel, massIOModel, fvMesh);
}

// ************************************************************************* //
