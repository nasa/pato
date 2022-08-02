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
    Creates BoundaryMappingModel instances

\*---------------------------------------------------------------------------*/

#include "simpleBoundaryMappingModel.H"
#include "twoDaxiFluxMapBoundaryMappingModel.H"
#include "twoDaxiPressureMapBoundaryMappingModel.H"
#include "twoDaxiFluxMapTimeBoundaryMappingModel.H"
#include "twoDaxiPressureMapTimeBoundaryMappingModel.H"
#include "twoDaxiBoundaryMappingModel.H"
#include "threeDtecplotBoundaryMappingModel.H"
#include "constantBoundaryMappingModel.H"
#include "gaussianBoundaryMappingModel.H"
#include "gaussianMixtureBoundaryMappingModel.H"
#include "polynomialBoundaryMappingModel.H"
#include "sigmoidBoundaryMappingModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(twoDaxiFluxMapBoundaryMappingModel, 0);
addToRunTimeSelectionTable(simpleBoundaryMappingModel, twoDaxiFluxMapBoundaryMappingModel, fvMesh);

defineTypeNameAndDebug(twoDaxiPressureMapBoundaryMappingModel, 0);
addToRunTimeSelectionTable(simpleBoundaryMappingModel, twoDaxiPressureMapBoundaryMappingModel, fvMesh);

defineTypeNameAndDebug(twoDaxiFluxMapTimeBoundaryMappingModel, 0);
addToRunTimeSelectionTable(simpleBoundaryMappingModel, twoDaxiFluxMapTimeBoundaryMappingModel, fvMesh);

defineTypeNameAndDebug(twoDaxiPressureMapTimeBoundaryMappingModel, 0);
addToRunTimeSelectionTable(simpleBoundaryMappingModel, twoDaxiPressureMapTimeBoundaryMappingModel, fvMesh);

defineTypeNameAndDebug(twoDaxiBoundaryMappingModel, 0);
addToRunTimeSelectionTable(simpleBoundaryMappingModel, twoDaxiBoundaryMappingModel, fvMesh);

defineTypeNameAndDebug(threeDtecplotBoundaryMappingModel, 0);
addToRunTimeSelectionTable(simpleBoundaryMappingModel, threeDtecplotBoundaryMappingModel, fvMesh);

defineTypeNameAndDebug(constantBoundaryMappingModel, 0);
addToRunTimeSelectionTable(simpleBoundaryMappingModel, constantBoundaryMappingModel, fvMesh);

defineTypeNameAndDebug(gaussianBoundaryMappingModel, 0);
addToRunTimeSelectionTable(simpleBoundaryMappingModel, gaussianBoundaryMappingModel, fvMesh);

defineTypeNameAndDebug(gaussianMixtureBoundaryMappingModel, 0);
addToRunTimeSelectionTable(simpleBoundaryMappingModel, gaussianMixtureBoundaryMappingModel, fvMesh);

defineTypeNameAndDebug(polynomialBoundaryMappingModel, 0);
addToRunTimeSelectionTable(simpleBoundaryMappingModel, polynomialBoundaryMappingModel, fvMesh);

defineTypeNameAndDebug(sigmoidBoundaryMappingModel, 0);
addToRunTimeSelectionTable(simpleBoundaryMappingModel, sigmoidBoundaryMappingModel, fvMesh);
}




// ************************************************************************* //
