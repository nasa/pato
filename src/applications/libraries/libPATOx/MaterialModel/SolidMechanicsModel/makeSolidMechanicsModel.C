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
    Creates SolidMechanicsModel instances

\*---------------------------------------------------------------------------*/

#include "simpleSolidMechanicsModel.H"
#include "noSolidMechanicsModel.H"
#if defined(FOAM_EXTEND)
#include "elasticSolidFoamSolidMechanicsModel.H"
#include "elasticOrthoSolidFoamSolidMechanicsModel.H"
#endif
#include "DisplacementSolidMechanicsModel.H"
#include "ElasticThermalSolidMechanicsModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(noSolidMechanicsModel, 0);
addToRunTimeSelectionTable
(
    simpleSolidMechanicsModel,
    noSolidMechanicsModel,
    fvMesh
);

#if defined(FOAM_EXTEND)
defineTypeNameAndDebug(elasticSolidFoamSolidMechanicsModel, 0);
addToRunTimeSelectionTable
(
    simpleSolidMechanicsModel,
    elasticSolidFoamSolidMechanicsModel,
    fvMesh
);

defineTypeNameAndDebug(elasticOrthoSolidFoamSolidMechanicsModel, 0);
addToRunTimeSelectionTable
(
    simpleSolidMechanicsModel,
    elasticOrthoSolidFoamSolidMechanicsModel,
    fvMesh
);
#endif

defineTypeNameAndDebug(DisplacementSolidMechanicsModel, 0);
addToRunTimeSelectionTable
(
    simpleSolidMechanicsModel,
    DisplacementSolidMechanicsModel,
    fvMesh
);

defineTypeNameAndDebug(ElasticThermalSolidMechanicsModel, 0);
addToRunTimeSelectionTable
(
    simpleSolidMechanicsModel,
    ElasticThermalSolidMechanicsModel,
    fvMesh
);
}

// ************************************************************************* //
