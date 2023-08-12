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
    Creates EnergyModel instances

\*---------------------------------------------------------------------------*/

#include "simpleEnergyModel.H"
#include "noEnergyModel.H"
#include "BoundaryTableEnergyModel.H"
#include "CorrectBCEnergyModel.H"
#include "FixedTemperatureEnergyModel.H"
#include "TabulatedTemperatureEnergyModel.H"
#include "PureConductionEnergyModel.H"
#include "PyrolysisEnergyModel.H"
#include "ForchheimerPyrolysisEnergyModel.H"
#include "Pyrolysis2TEnergyModel.H"
#include "Darcy2TEnergyModel.H"
#include "ForchheimerPyrolysis2TEnergyModel.H"
#include "Forchheimer2TEnergyModel.H"
#include "Pyrolysis_Heterogeneous_SpeciesDiffusionEnergyModel.H"
#include "addToRunTimeSelectionTable.H"
#include "ControlTemperatureEnergyModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(noEnergyModel, 0);
addToRunTimeSelectionTable(simpleEnergyModel, noEnergyModel, fvMesh);

defineTypeNameAndDebug(BoundaryTableEnergyModel, 0);
addToRunTimeSelectionTable(simpleEnergyModel, BoundaryTableEnergyModel, fvMesh);

defineTypeNameAndDebug(CorrectBCEnergyModel, 0);
addToRunTimeSelectionTable(simpleEnergyModel, CorrectBCEnergyModel, fvMesh);

defineTypeNameAndDebug(FixedTemperatureEnergyModel, 0);
addToRunTimeSelectionTable(simpleEnergyModel, FixedTemperatureEnergyModel, fvMesh);

defineTypeNameAndDebug(PureConductionEnergyModel, 0);
addToRunTimeSelectionTable(simpleEnergyModel, PureConductionEnergyModel, fvMesh);

defineTypeNameAndDebug(PyrolysisEnergyModel, 0);
addToRunTimeSelectionTable(simpleEnergyModel, PyrolysisEnergyModel, fvMesh);

defineTypeNameAndDebug(Pyrolysis2TEnergyModel, 0);
addToRunTimeSelectionTable(simpleEnergyModel, Pyrolysis2TEnergyModel, fvMesh);

defineTypeNameAndDebug(Pyrolysis_Heterogeneous_SpeciesDiffusionEnergyModel, 0);
addToRunTimeSelectionTable(simpleEnergyModel, Pyrolysis_Heterogeneous_SpeciesDiffusionEnergyModel, fvMesh);

defineTypeNameAndDebug(TabulatedTemperatureEnergyModel, 0);
addToRunTimeSelectionTable(simpleEnergyModel, TabulatedTemperatureEnergyModel, fvMesh);

defineTypeNameAndDebug(ForchheimerPyrolysisEnergyModel, 0);
addToRunTimeSelectionTable(simpleEnergyModel, ForchheimerPyrolysisEnergyModel, fvMesh);

defineTypeNameAndDebug(ForchheimerPyrolysis2TEnergyModel, 0);
addToRunTimeSelectionTable(simpleEnergyModel, ForchheimerPyrolysis2TEnergyModel, fvMesh);

defineTypeNameAndDebug(Darcy2TEnergyModel, 0);
addToRunTimeSelectionTable(simpleEnergyModel, Darcy2TEnergyModel, fvMesh);

defineTypeNameAndDebug(Forchheimer2TEnergyModel, 0);
addToRunTimeSelectionTable(simpleEnergyModel, Forchheimer2TEnergyModel, fvMesh);

defineTypeNameAndDebug(ControlTemperatureEnergyModel, 0);
addToRunTimeSelectionTable(simpleEnergyModel, ControlTemperatureEnergyModel, fvMesh);
}

// ************************************************************************* //
