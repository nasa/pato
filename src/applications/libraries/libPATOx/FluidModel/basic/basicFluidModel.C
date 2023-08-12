/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "basicFluidModel.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
defineTypeNameAndDebug(basicFluidModel, 0);
defineRunTimeSelectionTable(basicFluidModel, Time);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicFluidModel::basicFluidModel
(
    Time& runTime_
)
  :
rp(runTime_),
runTime(runTime_),
delaySolid_(runTime_.controlDict().lookupOrDefault<scalar>("delaySolid",0)),
delayFluid_(runTime_.controlDict().lookupOrDefault<scalar>("delayFluid",0)),
solidStitchTolerance_(runTime_.controlDict().lookupOrDefault<scalar>("solidStitchTolerance",0)),
fluidStitchTolerance_(runTime_.controlDict().lookupOrDefault<scalar>("fluidStitchTolerance",0)),
forceFluidUpdate_(runTime_.controlDict().lookupOrDefault<scalar>("forceFluidUpdate",0)),
maxDeltaTw_(fluidStitchTolerance_*2)
{
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


Foam::autoPtr<Foam::basicFluidModel> Foam::basicFluidModel::New
(
    Time& runTime_
)
{
  IOdictionary setCaseDict_
  (
      IOobject
      (
          "setCase",
          runTime_.constant(),
          runTime_.db(),
          IOobject::READ_IF_PRESENT,
          IOobject::NO_WRITE,
          false
      )
  );

  word typeName(setCaseDict_.lookupOrDefault<word>("fluidType","noFluidModel"));
  Info << "Selecting fluid type \"" << typeName << "\"" << endl << endl;

  typename basicFluidModel::TimeConstructorTable::iterator cstrIter =
      basicFluidModel::TimeConstructorTablePtr_->find(typeName);

  if (cstrIter == basicFluidModel::TimeConstructorTablePtr_->end()) {

    wordList list_available =  basicFluidModel::TimeConstructorTablePtr_->sortedToc();

    FatalErrorInFunction
        << "Unknown \"" << basicFluidModel::typeName << "\" fluid type in constant/setCase "
        << typeName << nl << nl
        << "Valid fluid types are:"
        << list_available
        << exit(FatalError);
  }

  return autoPtr<basicFluidModel>(cstrIter()(runTime_));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::basicFluidModel::~basicFluidModel()
{}


// * * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * * * //


// ************************************************************************* //

