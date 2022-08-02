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
    Creates MaterialPropertiesModel instances

\*---------------------------------------------------------------------------*/

#include "simpleMaterialPropertiesModel.H"
#include "noMaterialPropertiesModel.H"
#include "Porous_const_k_UQMaterialPropertiesModel.H"
#include "Porous_polynomial_k_UQMaterialPropertiesModel.H"
#include "PorousMaterialPropertiesModel.H"
#include "Porous_factorMaterialPropertiesModel.H"
#include "GradedPorousMaterialPropertiesModel.H"
#include "PureConductionMaterialPropertiesModel.H"
#include "PureConduction_UQMaterialPropertiesModel.H"
#include "FourierMaterialPropertiesModel.H"
#include "Fourier_RadiationMaterialPropertiesModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(noMaterialPropertiesModel, 0);
addToRunTimeSelectionTable(simpleMaterialPropertiesModel, noMaterialPropertiesModel, fvMesh);

defineTypeNameAndDebug(Porous_const_k_UQMaterialPropertiesModel, 0);
addToRunTimeSelectionTable(simpleMaterialPropertiesModel, Porous_const_k_UQMaterialPropertiesModel, fvMesh);

defineTypeNameAndDebug(Porous_polynomial_k_UQMaterialPropertiesModel, 0);
addToRunTimeSelectionTable(simpleMaterialPropertiesModel, Porous_polynomial_k_UQMaterialPropertiesModel, fvMesh);

defineTypeNameAndDebug(PorousMaterialPropertiesModel, 0);
addToRunTimeSelectionTable(simpleMaterialPropertiesModel, PorousMaterialPropertiesModel, fvMesh);

defineTypeNameAndDebug(PureConductionMaterialPropertiesModel, 0);
addToRunTimeSelectionTable(simpleMaterialPropertiesModel, PureConductionMaterialPropertiesModel, fvMesh);

defineTypeNameAndDebug(PureConduction_UQMaterialPropertiesModel, 0);
addToRunTimeSelectionTable(simpleMaterialPropertiesModel, PureConduction_UQMaterialPropertiesModel, fvMesh);

defineTypeNameAndDebug(FourierMaterialPropertiesModel, 0);
addToRunTimeSelectionTable(simpleMaterialPropertiesModel, FourierMaterialPropertiesModel, fvMesh);

defineTypeNameAndDebug(Fourier_RadiationMaterialPropertiesModel, 0);
addToRunTimeSelectionTable(simpleMaterialPropertiesModel, Fourier_RadiationMaterialPropertiesModel, fvMesh);

defineTypeNameAndDebug(GradedPorousMaterialPropertiesModel, 0);
addToRunTimeSelectionTable(simpleMaterialPropertiesModel, GradedPorousMaterialPropertiesModel, fvMesh);

defineTypeNameAndDebug(Porous_factorMaterialPropertiesModel, 0);
addToRunTimeSelectionTable(simpleMaterialPropertiesModel, Porous_factorMaterialPropertiesModel, fvMesh);
}

// ************************************************************************* //
