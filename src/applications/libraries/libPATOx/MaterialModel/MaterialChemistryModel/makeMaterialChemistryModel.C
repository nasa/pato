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
    Creates MaterialChemistryModel instances

\*---------------------------------------------------------------------------*/

#include "simpleMaterialChemistryModel.H"
#include "noMaterialChemistryModel.H"
#include "ConstantEquilibriumMaterialChemistryModel.H"
#include "ConstantFiniteRateMaterialChemistryModel.H"
#include "EquilibriumElementMaterialChemistryModel.H"
#include "OnlyFiniteRateMaterialChemistryModel.H"
#include "SpeciesConservationMaterialChemistryModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(noMaterialChemistryModel, 0);
addToRunTimeSelectionTable(simpleMaterialChemistryModel, noMaterialChemistryModel, fvMesh);

defineTypeNameAndDebug(ConstantEquilibriumMaterialChemistryModel, 0);
addToRunTimeSelectionTable(simpleMaterialChemistryModel, ConstantEquilibriumMaterialChemistryModel, fvMesh);

defineTypeNameAndDebug(ConstantFiniteRateMaterialChemistryModel, 0);
addToRunTimeSelectionTable(simpleMaterialChemistryModel, ConstantFiniteRateMaterialChemistryModel, fvMesh);

defineTypeNameAndDebug(EquilibriumElementMaterialChemistryModel, 0);
addToRunTimeSelectionTable(simpleMaterialChemistryModel, EquilibriumElementMaterialChemistryModel, fvMesh);

defineTypeNameAndDebug(OnlyFiniteRateMaterialChemistryModel, 0);
addToRunTimeSelectionTable(simpleMaterialChemistryModel, OnlyFiniteRateMaterialChemistryModel, fvMesh);

defineTypeNameAndDebug(SpeciesConservationMaterialChemistryModel, 0);
addToRunTimeSelectionTable(simpleMaterialChemistryModel, SpeciesConservationMaterialChemistryModel, fvMesh);

}

// ************************************************************************* //
