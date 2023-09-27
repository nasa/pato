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
    Creates PyrolysisModel instances

\*---------------------------------------------------------------------------*/

#include "simplePyrolysisModel.H"
#include "noPyrolysisModel.H"
#include "FIATPyrolysisModel.H"
#include "LinearArrheniusPyrolysisModel.H"
#include "virginPyrolysisModel.H"
#include "charPyrolysisModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(noPyrolysisModel, 0);
addToRunTimeSelectionTable(simplePyrolysisModel, noPyrolysisModel, fvMesh);

defineTypeNameAndDebug(FIATPyrolysisModel, 0);
addToRunTimeSelectionTable(simplePyrolysisModel, FIATPyrolysisModel, fvMesh);

defineTypeNameAndDebug(LinearArrheniusPyrolysisModel, 0);
addToRunTimeSelectionTable(simplePyrolysisModel, LinearArrheniusPyrolysisModel, fvMesh);

defineTypeNameAndDebug(virginPyrolysisModel, 0);
addToRunTimeSelectionTable(simplePyrolysisModel, virginPyrolysisModel, fvMesh);

defineTypeNameAndDebug(charPyrolysisModel, 0);
addToRunTimeSelectionTable(simplePyrolysisModel, charPyrolysisModel, fvMesh);
}

// ************************************************************************* //
