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
    Creates GasPropertiesModel instances

\*---------------------------------------------------------------------------*/

#include "simpleGasPropertiesModel.H"
#include "noGasPropertiesModel.H"
#include "EquilibriumGasPropertiesModel.H"
#include "FiniteRateGasPropertiesModel.H"
#include "TabulatedGasPropertiesModel.H"
#include "Tabulated2TGasPropertiesModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(noGasPropertiesModel, 0);
addToRunTimeSelectionTable(simpleGasPropertiesModel, noGasPropertiesModel, fvMesh);

defineTypeNameAndDebug(EquilibriumGasPropertiesModel, 0);
addToRunTimeSelectionTable(simpleGasPropertiesModel, EquilibriumGasPropertiesModel, fvMesh);

defineTypeNameAndDebug(FiniteRateGasPropertiesModel, 0);
addToRunTimeSelectionTable(simpleGasPropertiesModel, FiniteRateGasPropertiesModel, fvMesh);

defineTypeNameAndDebug(TabulatedGasPropertiesModel, 0);
addToRunTimeSelectionTable(simpleGasPropertiesModel, TabulatedGasPropertiesModel, fvMesh);

defineTypeNameAndDebug(Tabulated2TGasPropertiesModel, 0);
addToRunTimeSelectionTable(simpleGasPropertiesModel, Tabulated2TGasPropertiesModel, fvMesh);

}

// ************************************************************************* //
